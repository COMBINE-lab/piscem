use std::ffi::CString;
use std::io;
use std::os::raw::{c_char, c_int};
use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{bail, Result};
use clap::{ArgGroup, Args, Parser, Subcommand};
use num_cpus;
use tracing::{error, info, warn, Level};

#[link(name = "pesc_static", kind = "static")]
extern "C" {
    pub fn run_pesc_sc(args: c_int, argsv: *const *const c_char) -> c_int;
    pub fn run_pesc_bulk(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "build_static", kind = "static")]
extern "C" {
    pub fn run_build(args: c_int, argsv: *const *const c_char) -> c_int;
    pub fn run_build_poison_table(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "cfcore_static", kind = "static", modifiers = "+whole-archive")]
extern "C" {
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

pub trait AsArgv {
    fn as_argv(&self) -> Vec<CString>;
}

#[derive(Args, Clone, Debug)]
pub struct MapSCOpts {
    /// input index prefix
    #[arg(short, long)]
    index: String,

    /// geometry of barcode, umi and read
    #[arg(short, long)]
    geometry: String,

    /// path to list of read 1 files
    #[arg(short = '1', long, value_delimiter = ',', required = true)]
    read1: Vec<String>,

    /// path to list of read 2 files
    #[arg(short = '2', long, value_delimiter = ',', required = true)]
    read2: Vec<String>,

    /// number of threads to use
    #[arg(short, long)]
    threads: usize,

    /// path to output directory
    #[arg(short, long)]
    output: String,

    /// do not consider poison k-mers, even if the underlying index contains them.
    /// In this case, the mapping results will be identical to those obtained as if
    /// no poison table was added to the index.
    #[arg(long)]
    no_poison: bool,

    /// apply structural constraints when performing mapping.
    #[arg(short = 'c', long)]
    struct_constraints: bool,

    // skipping strategy to use for k-mer collection
    #[arg(long, default_value = "permissive", value_parser = clap::builder::PossibleValuesParser::new(["permissive", "strict"]))]
    skipping_strategy: String,

    /// skip checking of the equivalence classes of k-mers that were too
    /// ambiguous to be otherwise considered (passing this flag can speed up
    /// mapping slightly, but may reduce specificity).
    #[arg(long)]
    ignore_ambig_hits: bool,

    /// determines the maximum cardinality equivalence class
    /// (number of (txp, orientation status) pairs) to examine if performing check-ambig-hits.
    #[arg(
        long,
        short,
        requires = "check_ambig_hits",
        default_value_t = 4096,
        conflicts_with = "ignore_ambig_hits",
        help_heading = "Advanced options"
    )]
    max_ec_card: u32,

    /// in the first pass, consider only k-mers having <= --max-hit-occ hits.
    #[arg(long, default_value_t = 256, help_heading = "Advanced options")]
    max_hit_occ: u32,

    /// if all k-mers have > --max-hit-occ hits, then make a second pass and consider k-mers
    /// having <= --max-hit-occ-recover hits.
    #[arg(long, default_value_t = 1024, help_heading = "Advanced options")]
    max_hit_occ_recover: u32,

    /// reads with more than this number of mappings will not have
    /// their mappings reported.
    #[arg(long, default_value_t = 2500, help_heading = "Advanced options")]
    max_read_occ: u32,
}

#[derive(Args, Clone, Debug)]
#[command(group(
        ArgGroup::new("read_source")
        .required(true)
        .args(["read1", "reads"])
))]
pub struct MapBulkOpts {
    /// input index prefix
    #[arg(short, long)]
    index: String,

    /// path to list of read 1 files
    #[arg(short = '1', long, value_delimiter = ',', requires = "read2")]
    read1: Option<Vec<String>>,

    /// path to list of read 2 files
    #[arg(short = '2', long, value_delimiter = ',', requires = "read1")]
    read2: Option<Vec<String>>,

    /// path to list of read unpaired read files
    #[arg(short = 'r', long, value_delimiter = ',', conflicts_with_all = ["read1", "read2"])]
    reads: Option<Vec<String>>,

    /// number of threads to use
    #[arg(short, long)]
    threads: usize,

    /// path to output directory
    #[arg(short, long)]
    output: String,

    /// do not consider poison k-mers, even if the underlying index contains them.
    /// In this case, the mapping results will be identical to those obtained as if
    /// no poison table was added to the index.
    #[arg(long)]
    no_poison: bool,

    /// apply structural constraints when performing mapping.
    #[arg(short = 'c', long)]
    struct_constraints: bool,

    // skipping strategy to use for k-mer collection
    #[arg(long, default_value = "permissive", value_parser = clap::builder::PossibleValuesParser::new(["permissive", "strict"]))]
    skipping_strategy: String,

    /// skip checking of the equivalence classes of k-mers that were too
    /// ambiguous to be otherwise considered (passing this flag can speed up
    /// mapping slightly, but may reduce specificity).
    #[arg(long)]
    ignore_ambig_hits: bool,

    /// determines the maximum cardinality equivalence class
    /// (number of (txp, orientation status) pairs) to examine if performing check-ambig-hits.
    #[arg(
        long,
        short,
        requires = "check_ambig_hits",
        default_value_t = 4096,
        conflicts_with = "ignore_ambig_hits",
        help_heading = "Advanced options"
    )]
    max_ec_card: u32,

    /// in the first pass, consider only k-mers having <= --max-hit-occ hits.
    #[arg(long, default_value_t = 256, help_heading = "Advanced options")]
    max_hit_occ: u32,

    /// if all k-mers have > --max-hit-occ hits, then make a second pass and consider k-mers
    /// having <= --max-hit-occ-recover hits.
    #[arg(long, default_value_t = 1024, help_heading = "Advanced options")]
    max_hit_occ_recover: u32,

    /// reads with more than this number of mappings will not have
    /// their mappings reported.
    #[arg(long, default_value_t = 2500, help_heading = "Advanced options")]
    max_read_occ: u32,
}

impl AsArgv for MapSCOpts {
    fn as_argv(&self) -> Vec<CString> {
        let r1_string = self.read1.join(",");
        let r2_string = self.read2.join(",");

        let mut args: Vec<CString> = vec![
            CString::new("sc_ref_mapper").unwrap(),
            CString::new("-i").unwrap(),
            CString::new(self.index.clone()).unwrap(),
            CString::new("-g").unwrap(),
            CString::new(self.geometry.clone()).unwrap(),
            CString::new("-1").unwrap(),
            CString::new(r1_string.as_str()).unwrap(),
            CString::new("-2").unwrap(),
            CString::new(r2_string.as_str()).unwrap(),
            CString::new("-t").unwrap(),
            CString::new(self.threads.to_string()).unwrap(),
            CString::new("-o").unwrap(),
            CString::new(self.output.as_str()).unwrap(),
        ];

        if self.ignore_ambig_hits {
            args.push(CString::new("--ignore-ambig-hits").unwrap());
        } else {
            args.push(CString::new("--max-ec-card").unwrap());
            args.push(CString::new(self.max_ec_card.to_string()).unwrap());
        }

        if self.no_poison {
            args.push(CString::new("--no-poison").unwrap());
        }

        args.push(CString::new("--skipping-strategy").unwrap());
        args.push(CString::new(self.skipping_strategy.to_string()).unwrap());

        if self.struct_constraints {
            args.push(CString::new("--struct-constraints").unwrap());
        }

        args.push(CString::new("--max-hit-occ").unwrap());
        args.push(CString::new(self.max_hit_occ.to_string()).unwrap());

        args.push(CString::new("--max-hit-occ-recover").unwrap());
        args.push(CString::new(self.max_hit_occ_recover.to_string()).unwrap());

        args.push(CString::new("--max-read-occ").unwrap());
        args.push(CString::new(self.max_read_occ.to_string()).unwrap());

        args
    }
}

impl AsArgv for MapBulkOpts {
    fn as_argv(&self) -> Vec<CString> {

        let mut args: Vec<CString> = vec![
            CString::new("bulk_ref_mapper").unwrap(),
            CString::new("-i").unwrap(),
            CString::new(self.index.clone()).unwrap(),
            CString::new("-t").unwrap(),
            CString::new(self.threads.to_string()).unwrap(),
            CString::new("-o").unwrap(),
            CString::new(self.output.as_str()).unwrap(),
        ];

        if let Some(ref unpaired_reads) = &self.reads {
            let r_string = unpaired_reads.clone().join(",");
            args.push(CString::new("-r").unwrap());
            args.push(CString::new(r_string.as_str()).unwrap());
        } else if let (Some(ref r1), Some(ref r2)) = (&self.read1, &self.read2) {
            let r1_string = r1.clone().join(",");
            let r2_string = r2.clone().join(",");
            args.push(CString::new("-1").unwrap());
            args.push(CString::new(r1_string.as_str()).unwrap());
            args.push(CString::new("-2").unwrap());
            args.push(CString::new(r2_string.as_str()).unwrap());
        }

        if self.ignore_ambig_hits {
            args.push(CString::new("--ignore-ambig-hits").unwrap());
        } else {
            args.push(CString::new("--max-ec-card").unwrap());
            args.push(CString::new(self.max_ec_card.to_string()).unwrap());
            //idx_suffixes.push("ectab".into());
        }

        if self.no_poison {
            args.push(CString::new("--no-poison").unwrap());
        }

        args.push(CString::new("--skipping-strategy").unwrap());
        args.push(CString::new(self.skipping_strategy.to_string()).unwrap());

        if self.struct_constraints {
            args.push(CString::new("--struct-constraints").unwrap());
        }

        args.push(CString::new("--max-hit-occ").unwrap());
        args.push(CString::new(self.max_hit_occ.to_string()).unwrap());

        args.push(CString::new("--max-hit-occ-recover").unwrap());
        args.push(CString::new(self.max_hit_occ_recover.to_string()).unwrap());

        args.push(CString::new("--max-read-occ").unwrap());
        args.push(CString::new(self.max_read_occ.to_string()).unwrap());

        //if self.quiet {
        //args.push(CString::new("--quiet").unwrap());
        //}
        args
    }
}


#[derive(Debug, Subcommand)]
enum Commands {
    /// Index a reference sequence
    #[command(arg_required_else_help = true)]
    #[command(group(
            ArgGroup::new("ref-input")
            .required(true)
            .args(&["ref_seqs", "ref_lists", "ref_dirs"]),
            ))]
    Build {
        /// ',' separated list of reference FASTA files
        #[arg(short = 's', long, value_delimiter = ',', required = true)]
        ref_seqs: Option<Vec<String>>,

        /// ',' separated list of files (each listing input FASTA files)
        #[arg(short = 'l', long, value_delimiter = ',', required = true)]
        ref_lists: Option<Vec<String>>,

        /// ',' separated list of directories (all FASTA files in each directory will be indexed,
        /// but not recursively).
        #[arg(short = 'd', long, value_delimiter = ',', required = true)]
        ref_dirs: Option<Vec<String>>,

        /// length of k-mer to use
        #[arg(short, long)]
        klen: usize,

        /// length of minimizer to use
        #[arg(short, long)]
        mlen: usize,

        /// number of threads to use
        #[arg(short, long)]
        threads: usize,

        /// output file stem
        #[arg(short, long)]
        output: PathBuf,

        /// retain the reduced format GFA files produced by cuttlefish that
        /// describe the reference cDBG (the default is to remove these).
        #[arg(long)]
        keep_intermediate_dbg: bool,

        /// working directory where temporary files should be placed.
        #[arg(short = 'w', long, default_value_os_t = PathBuf::from("."))]
        work_dir: PathBuf,

        /// overwite an existing index if the output path is the same.
        #[arg(long)]
        overwrite: bool,

        /// skip the construction of the equivalence class lookup table
        /// when building the index (not recommended).
        #[arg(long)]
        no_ec_table: bool,

        /// path to (optional) decoy sequence used to insert poison
        /// k-mer information into the index.
        #[arg(long)]
        decoy_paths: Option<Vec<PathBuf>>,
    },

    /// map reads for single-cell processing
    #[command(arg_required_else_help = true)]
    MapSC(MapSCOpts),

    /// map reads for bulk processing
    #[command(arg_required_else_help = true)]
    MapBulk(MapBulkOpts),
}

fn main() -> Result<(), anyhow::Error> {
    let cli_args = Cli::parse();
    //env_logger::Builder::from_env(Env::default().default_filter_or("warn")).init();

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
        Commands::Build {
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
            decoy_paths,
        } => {
            info!("starting piscem build");
            if !(threads > 0) {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    threads
                );
            }
            if threads > ncpus {
                bail!("the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    threads, ncpus);
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

            let mut args: Vec<CString> = vec![];

            let cf_out = PathBuf::from(output.as_path().to_string_lossy().into_owned() + "_cfish");
            let cf_base_path = cf_out.as_path();
            let seg_file = cf_base_path.with_extension("cf_seg");
            let seq_file = cf_base_path.with_extension("cf_seq");
            let struct_file = cf_base_path.with_extension("json");
            let mut build_ret;

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

            if struct_file.exists() {
                if !seq_file.exists() || !seg_file.exists() {
                    warn!("The prefix you have chosen for output already corresponds to an existing cDBG structure file {:?}.", struct_file.display());
                    warn!("However, the corresponding seq and seg files do not exist. Please either delete this structure file, choose another output prefix, or use the --overwrite flag.");
                    bail!("Cannot write over existing index without the --overwrite flag.");
                }
            }

            args.push(CString::new("cdbg_builder").unwrap());

            // We can treat the different input options independently
            // here because the argument parser should have enforced
            // their exclusivity.
            let mut has_input = false;

            if let Some(seqs) = ref_seqs {
                if !seqs.is_empty() {
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
                            error!("when attempting to create working directory {}, encountered error {:#?}", &work_dir.display(), e);
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

            args.clear();
            args.push(CString::new("ref_index_builder").unwrap());

            args.push(CString::new("-i").unwrap());
            args.push(CString::new(cf_out.as_path().to_string_lossy().into_owned()).unwrap());
            args.push(CString::new("-k").unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new("-m").unwrap());
            args.push(CString::new(mlen.to_string()).unwrap()); // minimizer length

            args.push(CString::new("--canonical-parsing").unwrap());
            if !no_ec_table {
                args.push(CString::new("--build-ec-table").unwrap());
            }
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(output.as_path().to_string_lossy().into_owned()).unwrap());

            args.push(CString::new("-d").unwrap());
            args.push(CString::new(work_dir.as_path().to_string_lossy().into_owned()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(threads.to_string()).unwrap());

            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            {
                println!("{:?}", args);
                let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                let args_len: c_int = args.len() as c_int;
                build_ret = unsafe { run_build(args_len, arg_ptrs.as_ptr()) };
            }

            if build_ret != 0 {
                bail!("indexer returned exit code {}; failure.", build_ret);
            }

            // now, build the poison table if there are decoys
            if let Some(decoy_pathbufs) = decoy_paths {
                args.clear();
                args.push(CString::new("poison_table_builder").unwrap());

                // index is the one we just built
                args.push(CString::new("-i").unwrap());
                args.push(CString::new(output.as_path().to_string_lossy().into_owned()).unwrap());

                args.push(CString::new("-t").unwrap());
                args.push(CString::new(threads.to_string()).unwrap());

                if overwrite {
                    args.push(CString::new("--overwrite").unwrap());
                }

                if quiet {
                    args.push(CString::new("--quiet").unwrap());
                }

                let path_args = decoy_pathbufs
                    .into_iter()
                    .map(|x| x.to_string_lossy().into_owned())
                    .collect::<Vec<String>>()
                    .join(",");
                args.push(CString::new("-d").unwrap());
                args.push(CString::new(path_args).unwrap());

                {
                    println!("{:?}", args);
                    let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                    let args_len: c_int = args.len() as c_int;
                    build_ret = unsafe { run_build_poison_table(args_len, arg_ptrs.as_ptr()) };
                }
                if build_ret != 0 {
                    bail!(
                        "building poison table returned exit code {}; failure.",
                        build_ret
                    );
                }
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

            info!("piscem build finished");
        }

        Commands::MapSC(sc_opts) => {
            if !(sc_opts.threads > 0) {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    sc_opts.threads
                );
            }
            if sc_opts.threads > ncpus {
                bail!("the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    sc_opts.threads, ncpus);
            }

            let mut idx_suffixes: Vec<String> =
                vec!["sshash".into(), "ctab".into(), "refinfo".into()];

            if !sc_opts.ignore_ambig_hits {
                idx_suffixes.push("ectab".into());
            }

            let idx_path = PathBuf::from_str(&sc_opts.index)?;
            for s in idx_suffixes {
                let req_file = idx_path.with_extension(s);
                if !req_file.exists() {
                    bail!("To load the index with the specified prefix {}, piscem expects the file {} to exist, but it does not!", &sc_opts.index, req_file.display());
                }
            }

            let mut args = sc_opts.as_argv();
            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            info!("cmd: {:?}", args);
            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            let map_ret = unsafe { run_pesc_sc(args_len, arg_ptrs.as_ptr()) };
            if map_ret != 0 {
                bail!("mapper returned exit code {}; failure", map_ret);
            }
        }

        Commands::MapBulk(bulk_opts) => {
            if !(bulk_opts.threads > 0) {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    bulk_opts.threads
                );
            }
            if bulk_opts.threads > ncpus {
                bail!("the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    bulk_opts.threads, ncpus);
            }

            let mut idx_suffixes: Vec<String> =
                vec!["sshash".into(), "ctab".into(), "refinfo".into()];

            if !bulk_opts.ignore_ambig_hits {
                idx_suffixes.push("ectab".into());
            }

            let idx_path = PathBuf::from_str(&bulk_opts.index)?;
            for s in idx_suffixes {
                let req_file = idx_path.with_extension(s);
                if !req_file.exists() {
                    bail!("To load the index with the specified prefix {}, piscem expects the file {} to exist, but it does not!", &bulk_opts.index, req_file.display());
                }
            }
 
            let mut args = bulk_opts.as_argv();

            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            let map_ret = unsafe { run_pesc_bulk(args_len, arg_ptrs.as_ptr()) };
            if map_ret != 0 {
                bail!("mapper returned exit code {}; failure", map_ret);
            }
        }
    }
    Ok(())
}
