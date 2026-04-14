use mimalloc::MiMalloc;
use std::io;
use std::path::PathBuf;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use anyhow::{Context, Result, bail};
use cf1_rs::{CfInput, cf_build as cf1_build};
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
            dict,
        }) => {
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

            // Build CfInput from CLI args, using native variants where possible.
            let cf_input = if let Some(ref seqs) = ref_seqs {
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
                    CfInput::Files(seqs.iter().map(PathBuf::from).collect())
                } else {
                    bail!("--ref-seqs provided but empty");
                }
            } else if let Some(ref lists) = ref_lists {
                if lists.len() == 1 {
                    CfInput::ListFile(PathBuf::from(&lists[0]))
                } else {
                    // Multiple list files — resolve all to individual paths.
                    let mut files = Vec::new();
                    for list_path in lists {
                        let contents = std::fs::read_to_string(list_path)
                            .with_context(|| format!("reading list file {}", list_path))?;
                        for line in contents.lines() {
                            let line = line.trim();
                            if !line.is_empty() {
                                files.push(PathBuf::from(line));
                            }
                        }
                    }
                    CfInput::Files(files)
                }
            } else if let Some(ref dirs) = ref_dirs {
                if dirs.len() == 1 {
                    CfInput::Directory(PathBuf::from(&dirs[0]))
                } else {
                    // Multiple directories — resolve all to individual file paths.
                    let mut files = Vec::new();
                    for dir_path in dirs {
                        for entry in std::fs::read_dir(dir_path)
                            .with_context(|| format!("reading directory {}", dir_path))?
                        {
                            let path = entry?.path();
                            if path.is_file() {
                                files.push(path);
                            }
                        }
                    }
                    CfInput::Files(files)
                }
            } else {
                bail!("Input (via --ref-seqs, --ref-lists, or --ref-dirs) must be provided.");
            };

            let cf_out = PathBuf::from(output.as_path().to_string_lossy().into_owned() + "_cfish");
            let cf_out_str = cf_out.as_path().to_string_lossy();
            let seg_file = PathBuf::from(format!("{}.cf_seg", cf_out_str));
            let seq_file = PathBuf::from(format!("{}.cf_seq", cf_out_str));
            let struct_file = PathBuf::from(format!("{}.json", cf_out_str));

            if overwrite {
                for f in [&struct_file, &seg_file, &seq_file] {
                    if f.exists() {
                        std::fs::remove_file(f)?;
                    }
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

            // Ensure work directory exists.
            match work_dir.try_exists() {
                Ok(true) => {
                    info!(
                        "will use {} as the work directory for temporary files.",
                        work_dir.display()
                    );
                }
                Ok(false) => {
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

            // Ensure the output parent directory exists.
            if let Some(parent_path) = cf_out.parent() {
                if !parent_path.as_os_str().is_empty() && !parent_path.exists() {
                    std::fs::create_dir_all(parent_path)?;
                    info!(
                        "directory {} did not already exist; creating it.",
                        parent_path.display()
                    );
                }
            }

            info!("starting cDBG construction with cf1-rs");
            let cf_result = cf1_build()
                .input(cf_input)
                .output_prefix(cf_out.clone())
                .k(klen)
                .threads(threads)
                .work_dir(work_dir.clone())
                .call()?;

            info!(
                "cDBG construction complete: {} unitigs, {} vertices",
                cf_result.unitig_count, cf_result.vertex_count
            );

            // Build piscem-rs index from cuttlefish output
            rs_build::run(BuildArgs {
                input: cf_out.clone(),
                output: output.clone(),
                klen,
                mlen,
                threads,
                build_ec_table: !no_ec_table,
                canonical: true,
                seed,
                single_mphf: false,
                dict,
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
                info!("removing intermediate cdBG files produced by cf1-rs.");
                for (path, label) in [
                    (&cf_result.seg_file, "segment"),
                    (&cf_result.seq_file, "tiling"),
                ] {
                    match std::fs::remove_file(path) {
                        Ok(_) => info!("removed {} file {}", label, path.display()),
                        Err(e) => warn!("cannot remove {}, error: {:?}", path.display(), e),
                    }
                }
                // Keep the json file — it's small and may contain useful info.
            }

            let piscem_rs_ver = piscem_rs::VERSION;
            let cf1_rs_ver = env!("CF1_RS_VERSION");
            let piscem_ver = clap::crate_version!();

            let version_json = json!({
                "piscem-rs": piscem_rs_ver,
                "cf1-rs": cf1_rs_ver,
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

            let args = MapScrnaArgs {
                index: PathBuf::from(&sc_opts.index),
                read1: sc_opts.read1.iter().map(PathBuf::from).collect(),
                read2: sc_opts.read2.iter().map(PathBuf::from).collect(),
                geometry: sc_opts.geometry,
                output: sc_opts.output.clone(),
                threads: sc_opts.threads,
                skipping_strategy: sc_opts.skipping_strategy,
                no_poison: sc_opts.no_poison,
                struct_constraints: sc_opts.struct_constraints,
                ignore_ambig_hits: sc_opts.ignore_ambig_hits,
                max_ec_card: sc_opts.max_ec_card,
                max_hit_occ: sc_opts.max_hit_occ as usize,
                max_hit_occ_recover: sc_opts.max_hit_occ_recover as usize,
                max_read_occ: sc_opts.max_read_occ as usize,
                with_position: sc_opts.with_position,
                quiet,
                dict: sc_opts.dict,
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
                dict: scatac_opts.dict,
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

            let args = MapBulkArgs {
                index: PathBuf::from(&bulk_opts.index),
                reads: bulk_opts.reads.unwrap_or_default().iter().map(PathBuf::from).collect(),
                read1: bulk_opts.read1.unwrap_or_default().iter().map(PathBuf::from).collect(),
                read2: bulk_opts.read2.unwrap_or_default().iter().map(PathBuf::from).collect(),
                output: bulk_opts.output.clone(),
                threads: bulk_opts.threads,
                skipping_strategy: bulk_opts.skipping_strategy,
                no_poison: bulk_opts.no_poison,
                struct_constraints: bulk_opts.struct_constraints,
                ignore_ambig_hits: bulk_opts.ignore_ambig_hits,
                max_ec_card: bulk_opts.max_ec_card,
                max_hit_occ: bulk_opts.max_hit_occ as usize,
                max_hit_occ_recover: bulk_opts.max_hit_occ_recover as usize,
                max_read_occ: bulk_opts.max_read_occ as usize,
                quiet,
                dict: bulk_opts.dict,
            };
            map_bulk::run(args)?;
        }
    }
    Ok(())
}
