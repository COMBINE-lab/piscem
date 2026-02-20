use std::ffi::{OsStr, OsString};
use std::io;
use std::os::raw::{c_char, c_int};
use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand};
use serde_json::json;
use tracing::{Level, error, info, warn};

mod piscem_commands;
use piscem_commands::*;

use piscem_rs::cli::build as rs_build;
use piscem_rs::cli::poison as rs_poison;
use piscem_rs::cli::map_bulk;
use piscem_rs::cli::map_scrna;
use piscem_rs::cli::map_scatac;
use piscem_rs::cli::build::BuildArgs;
use piscem_rs::cli::poison::BuildPoisonArgs;
use piscem_rs::cli::map_scrna::MapScrnaArgs;
use piscem_rs::cli::map_bulk::MapBulkArgs;
use piscem_rs::cli::map_scatac::MapScatacArgs;

#[link(name = "cfcore_static", kind = "static", modifiers = "+whole-archive")]
unsafe extern "C" {
    pub fn cf_build(args: c_int, argsv: *const *const c_char) -> c_int;
}

/// Indexing and mapping to compacted colored de Bruijn graphs
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
struct Cli {
    /// be quiet (no effect yet for cDBG building phase of indexing).
    #[arg(short, long)]
    quiet: bool,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Index a reference sequence
    #[command(arg_required_else_help = true)]
    Build(BuildOpts),

    /// map reads for single-cell processing
    #[command(arg_required_else_help = true)]
    MapSC(MapSCOpts),

    /// map reads for bulk processing
    #[command(arg_required_else_help = true)]
    MapBulk(MapBulkOpts),

    /// map reads for scAtac processing
    #[command(arg_required_else_help = true)]
    MapSCAtac(MapSCAtacOpts),
}

// from: https://stackoverflow.com/questions/74322541/how-to-append-to-pathbuf
fn append_to_path(p: impl Into<OsString>, s: impl AsRef<OsStr>) -> PathBuf {
    let mut p = p.into();
    p.push(s);
    p.into()
}

fn main() -> Result<(), anyhow::Error> {
    let cli_args = Cli::parse();

    let quiet = cli_args.quiet;
    if quiet {
        tracing_subscriber::fmt()
            .with_max_level(Level::WARN)
            .with_writer(io::stderr)
            .init();
    } else {
        tracing_subscriber::fmt()
            .with_max_level(Level::INFO)
            .with_writer(io::stderr)
            .init();
    }

    let ncpus = num_cpus::get();

    match cli_args.command {
        Commands::Build(BuildOpts {
            ref_seqs,
            ref_lists,
            ref_dirs,
            klen,
            mlen,
            threads,
            output,
            keep_intermediate_dbg,
            work_dir,
            overwrite,
            no_ec_table,
            polya_clip_length,
            decoy_paths,
            seed,
        }) => {
            use std::ffi::CString;
            use std::os::raw::c_char;

            info!("starting piscem build");
            if threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    threads
                );
            }
            if threads > ncpus {
                bail!(
                    "the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    threads,
                    ncpus
                );
            }
            if mlen >= klen {
                bail!(
                    "minimizer length ({}) must be < k-mer length ({})",
                    mlen,
                    klen
                );
            }

            let idxthreads = if !threads.is_power_of_two() {
                let idxthreads = 1_usize.max(threads.next_power_of_two() / 2);
                warn!(
                    r#"The number of threads used for the indexing step must be a power of 2.
Using {} for cDBG construction and {} for indexing.
NOTE: This is a temporary restriction and should be lifted in a future version of piscem."#,
                    threads, idxthreads
                );
                idxthreads
            } else {
                threads
            };

            // if the decoy sequences are provided, ensure they are valid paths
            if let Some(ref decoys) = decoy_paths {
                for d in decoys {
                    match d.try_exists() {
                        Ok(true) => {}
                        Ok(false) => {
                            bail!(
                                "Path for decoy file {} seems not to point to a valid file",
                                d.display()
                            );
                        }
                        Err(e) => {
                            bail!(
                                "Error {} when checking the existence of decoy file {}",
                                e,
                                d.display()
                            );
                        }
                    }
                }
            }

            // Ensure the parent directory of the output prefix exists.
            if let Some(parent) = output.parent() {
                if !parent.as_os_str().is_empty() && !parent.exists() {
                    std::fs::create_dir_all(parent)?;
                    info!(
                        "created output directory {} (did not previously exist)",
                        parent.display()
                    );
                }
            }

            let mut args: Vec<CString> = vec![];

            let cf_out = PathBuf::from(output.as_path().to_string_lossy().into_owned() + "_cfish");
            let cf_base_path = cf_out.as_path();
            let seg_file = append_to_path(cf_base_path, ".cf_seg");
            let seq_file = append_to_path(cf_base_path, ".cf_seq");
            let struct_file = append_to_path(cf_base_path, ".json");
            let build_ret;

            if overwrite {
                if struct_file.exists() {
                    std::fs::remove_file(struct_file.clone())?;
                }
                if seg_file.exists() {
                    std::fs::remove_file(seg_file.clone())?;
                }
                if seq_file.exists() {
                    std::fs::remove_file(seq_file.clone())?;
                }
            }

            if struct_file.exists() && (!seq_file.exists() || !seg_file.exists()) {
                warn!(
                    "The prefix you have chosen for output already corresponds to an existing cDBG structure file {:?}.",
                    struct_file.display()
                );
                warn!(
                    "However, the corresponding seq and seg files do not exist. Please either delete this structure file, choose another output prefix, or use the --overwrite flag."
                );
                bail!("Cannot write over existing index without the --overwrite flag.");
            }

            args.push(CString::new("cdbg_builder").unwrap());

            let mut has_input = false;

            if let Some(seqs) = ref_seqs {
                if !seqs.is_empty() {
                    let out_stem =
                        PathBuf::from(output.as_path().to_string_lossy().into_owned() + ".sigs");
                    let configs = prepare_fasta::RecordParseConfig {
                        input: seqs.clone(),
                        output_stem: out_stem,
                        polya_clip_length,
                    };
                    info!("Computing and recording reference signatures...");
                    prepare_fasta::parse_records(configs)?;
                    info!("done.");
                    args.push(CString::new("--seq").unwrap());
                    let reflist = seqs.join(",");
                    args.push(CString::new(reflist.as_str()).unwrap());
                    has_input = true;
                }
            }

            if let Some(lists) = ref_lists {
                if !lists.is_empty() {
                    args.push(CString::new("--list").unwrap());
                    let reflist = lists.join(",");
                    args.push(CString::new(reflist.as_str()).unwrap());
                    has_input = true;
                }
            }

            if let Some(dirs) = ref_dirs {
                if !dirs.is_empty() {
                    args.push(CString::new("--dir").unwrap());
                    let reflist = dirs.join(",");
                    args.push(CString::new(reflist.as_str()).unwrap());
                    has_input = true;
                }
            }

            assert!(
                has_input,
                "Input (via --ref-seqs, --ref-lists, or --ref-dirs) must be provided."
            );

            args.push(CString::new("-k").unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new("--track-short-seqs").unwrap());
            args.push(CString::new("--poly-N-stretch").unwrap());
            args.push(CString::new("--collate-output-in-mem").unwrap());
            // check if the provided work directory exists.
            // If not, then try and create it.
            match work_dir.try_exists() {
                Ok(true) => {
                    info!(
                        "will use {} as the work directory for temporary files.",
                        work_dir.display()
                    );
                }
                Ok(false) => {
                    // try to create it
                    match std::fs::create_dir_all(&work_dir) {
                        Ok(_) => {}
                        Err(e) => {
                            error!(
                                "when attempting to create working directory {}, encountered error {:#?}",
                                &work_dir.display(),
                                e
                            );
                            bail!(
                                "Failed to create working directory {} for index construction : {:#?}",
                                &work_dir.display(),
                                e
                            );
                        }
                    }
                }
                Err(e) => {
                    error!("when checking existence of working directory {:#?}", e);
                    bail!(
                        "Failed to create working directory for index construction : {:#?}",
                        e
                    );
                }
            }

            // check if the provided output path is more than just a prefix
            // if so, check if the specified directory exists and create it
            // if it doesn't.
            if let Some(parent_path) = cf_out.parent() {
                if !parent_path.exists() {
                    std::fs::create_dir_all(parent_path)?;
                    info!(
                        "directory {} did not already exist; creating it.",
                        parent_path.display()
                    );
                }
            }

            args.push(CString::new("-o").unwrap());
            args.push(CString::new(cf_out.as_path().to_string_lossy().into_owned()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(threads.to_string()).unwrap());
            // output format
            args.push(CString::new("-f").unwrap());
            args.push(CString::new("3").unwrap());
            // work directory
            args.push(CString::new("-w").unwrap());
            args.push(CString::new(work_dir.as_path().to_string_lossy().into_owned()).unwrap());

            info!("args = {:?}", args);
            {
                let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                let args_len: c_int = args.len() as c_int;
                build_ret = unsafe { cf_build(args_len, arg_ptrs.as_ptr()) };
            }

            if build_ret != 0 {
                bail!(
                    "cDBG constructor returned exit code {}; failure.",
                    build_ret
                );
            }

            // Build piscem-rs index
            rs_build::run(BuildArgs {
                input: cf_out.clone(),
                output: output.clone(),
                klen,
                mlen,
                threads: idxthreads,
                build_ec_table: !no_ec_table,
                canonical: true,
                seed,
                single_mphf: false,
            })?;

            // Build poison table if decoys were provided
            if let Some(decoy_pathbufs) = decoy_paths {
                rs_poison::run(BuildPoisonArgs {
                    index: output.clone(),
                    decoys: decoy_pathbufs,
                    threads,
                })?;
            }

            if !keep_intermediate_dbg {
                info!("removing intermediate cdBG files produced by cuttlefish.");

                match std::fs::remove_file(seg_file.clone()) {
                    Ok(_) => {
                        info!("removed segment file {}", seg_file.display());
                    }
                    Err(e) => {
                        warn!(
                            "cannot remove {}, encountered error {:?}!",
                            seg_file.display(),
                            e
                        );
                    }
                };

                match std::fs::remove_file(seq_file.clone()) {
                    Ok(_) => {
                        info!("removed tiling file {}", seq_file.display());
                    }
                    Err(e) => {
                        warn!(
                            "cannot remove {}, encountered error {:?}!",
                            seq_file.display(),
                            e
                        );
                    }
                };
                // for now, let the json file stick around. It's
                // generally very small and may contain useful information
                // about the references being indexed.
            }

            let piscem_rs_ver = piscem_rs::VERSION;
            let cuttlefish_ver = env!("cuttlefish-ver");
            let piscem_ver = clap::crate_version!();

            let version_json = json!({
                "piscem-rs": piscem_rs_ver,
                "cuttlefish": cuttlefish_ver,
                "piscem": piscem_ver
            });

            let mut fname = output.as_path().to_string_lossy().into_owned();
            fname.push_str("_ver.json");

            let fname = PathBuf::from(fname);
            let ver_file = std::fs::File::create(fname)?;
            serde_json::to_writer_pretty(ver_file, &version_json)?;

            info!("piscem build finished.");
        }

        Commands::MapSC(sc_opts) => {
            if sc_opts.threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    sc_opts.threads
                );
            }
            if sc_opts.threads > ncpus {
                bail!(
                    "the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    sc_opts.threads,
                    ncpus
                );
            }

            if sc_opts.struct_constraints {
                warn!("--struct-constraints is not supported by piscem-rs and will be ignored");
            }

            let args = MapScrnaArgs {
                index: PathBuf::from(&sc_opts.index),
                read1: sc_opts.read1.iter().map(PathBuf::from).collect(),
                read2: sc_opts.read2.iter().map(PathBuf::from).collect(),
                geometry: sc_opts.geometry,
                output: sc_opts.output.clone(),
                threads: sc_opts.threads,
                skipping_strategy: sc_opts.skipping_strategy,
                no_poison: sc_opts.no_poison,
                ignore_ambig_hits: sc_opts.ignore_ambig_hits,
                max_ec_card: sc_opts.max_ec_card,
                max_hit_occ: sc_opts.max_hit_occ as usize,
                max_hit_occ_recover: sc_opts.max_hit_occ_recover as usize,
                max_read_occ: sc_opts.max_read_occ as usize,
                with_position: sc_opts.with_position,
                quiet,
            };
            map_scrna::run(args)?;
        }

        Commands::MapSCAtac(scatac_opts) => {
            if scatac_opts.threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    scatac_opts.threads
                );
            }
            if scatac_opts.threads > ncpus {
                bail!(
                    "the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    scatac_opts.threads,
                    ncpus
                );
            }

            for (flag, name) in [
                (scatac_opts.sam_format, "--sam-format"),
                (scatac_opts.bed_format, "--bed-format"),
                (scatac_opts.use_chr, "--use-chr"),
                (scatac_opts.check_kmer_orphan, "--check-kmer-orphan"),
                (scatac_opts.struct_constraints, "--struct-constraints"),
            ] {
                if flag {
                    warn!("{} is not supported by piscem-rs and will be ignored", name);
                }
            }

            let barcode = scatac_opts.barcode.unwrap_or_default();
            let args = MapScatacArgs {
                index: PathBuf::from(&scatac_opts.index),
                reads: scatac_opts.reads.unwrap_or_default().iter().map(PathBuf::from).collect(),
                read1: scatac_opts.read1.unwrap_or_default().iter().map(PathBuf::from).collect(),
                read2: scatac_opts.read2.unwrap_or_default().iter().map(PathBuf::from).collect(),
                barcode: barcode.iter().map(PathBuf::from).collect(),
                output: scatac_opts.output.clone(),
                threads: scatac_opts.threads,
                bc_len: scatac_opts.bclen as usize,
                no_tn5_shift: scatac_opts.no_tn5_shift,
                no_poison: scatac_opts.no_poison,
                check_ambig_hits: !scatac_opts.ignore_ambig_hits,
                max_ec_card: scatac_opts.max_ec_card,
                max_hit_occ: scatac_opts.max_hit_occ as usize,
                max_hit_occ_recover: scatac_opts.max_hit_occ_recover as usize,
                max_read_occ: scatac_opts.max_read_occ as usize,
                end_cache_capacity: scatac_opts.end_cache_capacity,
                bin_size: scatac_opts.bin_size as u64,
                bin_overlap: scatac_opts.bin_overlap as u64,
                thr: scatac_opts.thr,
                min_overlap: 30,
                skipping_strategy: None,
                quiet,
            };
            map_scatac::run(args)?;
        }

        Commands::MapBulk(bulk_opts) => {
            if bulk_opts.threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    bulk_opts.threads
                );
            }
            if bulk_opts.threads > ncpus {
                bail!(
                    "the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    bulk_opts.threads,
                    ncpus
                );
            }

            if bulk_opts.struct_constraints {
                warn!("--struct-constraints is not supported by piscem-rs and will be ignored");
            }

            let args = MapBulkArgs {
                index: PathBuf::from(&bulk_opts.index),
                reads: bulk_opts.reads.unwrap_or_default().iter().map(PathBuf::from).collect(),
                read1: bulk_opts.read1.unwrap_or_default().iter().map(PathBuf::from).collect(),
                read2: bulk_opts.read2.unwrap_or_default().iter().map(PathBuf::from).collect(),
                output: bulk_opts.output.clone(),
                threads: bulk_opts.threads,
                skipping_strategy: bulk_opts.skipping_strategy,
                no_poison: bulk_opts.no_poison,
                ignore_ambig_hits: bulk_opts.ignore_ambig_hits,
                max_ec_card: bulk_opts.max_ec_card,
                max_hit_occ: bulk_opts.max_hit_occ as usize,
                max_hit_occ_recover: bulk_opts.max_hit_occ_recover as usize,
                max_read_occ: bulk_opts.max_read_occ as usize,
                quiet,
            };
            map_bulk::run(args)?;
        }
    }
    Ok(())
}
