use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::path::PathBuf;

use anyhow::{bail, Result};
use clap::{ArgGroup, Parser, Subcommand};
use env_logger::Env;
use log::info;

#[link(name = "pesc_static", kind = "static")]
extern "C" {
    pub fn run_pesc_sc(args: c_int, argsv: *const *const c_char) -> c_int;
    pub fn run_pesc_bulk(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "build_static", kind = "static")]
extern "C" {
    pub fn run_build(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "cfcore_static", kind = "static", modifiers = "+whole-archive")]
extern "C" {
    pub fn cf_build(args: c_int, argsv: *const *const c_char) -> c_int;
}

/// Indexing and mapping to compacted colored de Bruijn graphs
#[derive(Debug, Parser)]
#[clap(version)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Index a reference sequence
    #[clap(arg_required_else_help = true)]
    #[clap(group(
            ArgGroup::new("ref-input")
            .required(true)
            .args(&["ref-seqs", "ref-lists", "ref-dirs"]),
            ))]
    Build {
        /// reference FASTA location
        #[clap(
            short = 's',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        ref_seqs: Option<Vec<String>>,

        #[clap(
            short = 'l',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        ref_lists: Option<Vec<String>>,

        #[clap(
            short = 'd',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        ref_dirs: Option<Vec<String>>,
        /// length of k-mer to use
        #[clap(short, long, value_parser)]
        klen: usize,

        /// length of minimizer to use
        #[clap(short, long, value_parser)]
        mlen: usize,

        /// number of threads to use
        #[clap(short, long, value_parser)]
        threads: usize,

        /// output file stem
        #[clap(short, long, value_parser)]
        output: PathBuf,

        /// be quiet during the indexing phase (no effect yet for cDBG building).
        #[clap(short, action)]
        quiet: bool,
    },

    /// map sc reads
    #[clap(arg_required_else_help = true)]
    MapSC {
        /// input index prefix
        #[clap(short, long, value_parser)]
        index: String,

        /// geometry of barcode, umi and read
        #[clap(short, long, value_parser)]
        geometry: String,

        /// path to list of read 1 files
        #[clap(
            short = '1',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        read1: Vec<String>,

        /// path to list of read 1 files
        #[clap(
            short = '2',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        read2: Vec<String>,

        /// number of threads to use
        #[clap(short, long, value_parser)]
        threads: usize,

        /// path to output directory
        #[clap(short, long, value_parser)]
        output: String,

        /// be quiet during mapping
        #[clap(short, action)]
        quiet: bool,
    },

    /// map bulk reads
    #[clap(arg_required_else_help = true)]
    MapBulk {
        /// input index prefix
        #[clap(short, long, value_parser)]
        index: String,

        /// path to list of read 1 files
        #[clap(
            short = '1',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        read1: Vec<String>,

        /// path to list of read 1 files
        #[clap(
            short = '2',
            long,
            value_parser,
            value_delimiter = ',',
            required = true
        )]
        read2: Vec<String>,

        /// number of threads to use
        #[clap(short, long, value_parser)]
        threads: usize,

        /// path to output directory
        #[clap(short, long, value_parser)]
        output: String,

        /// be quiet during mapping
        #[clap(short, action)]
        quiet: bool,
    },
}

fn main() -> Result<(), anyhow::Error> {
    let cli_args = Cli::parse();
    env_logger::Builder::from_env(Env::default().default_filter_or("warn")).init();

    match cli_args.command {
        Commands::Build {
            ref_seqs,
            ref_lists,
            ref_dirs,
            klen,
            mlen,
            threads,
            output,
            quiet,
        } => {
            info!("starting piscem build");
            assert!(
                mlen < klen,
                "minimizer length ({}) >= k-mer len ({})",
                mlen,
                klen
            );

            let mut args: Vec<CString> = vec![];

            let cf_out = PathBuf::from(output.as_path().to_string_lossy().into_owned() + "_cfish");
            let mut build_ret;

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
            args.push(CString::new(cf_out.as_path().to_string_lossy().into_owned()).unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new(mlen.to_string()).unwrap()); // minimizer length
            args.push(CString::new("--canonical-parsing").unwrap());
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(output.as_path().to_string_lossy().into_owned()).unwrap());
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

            info!("piscem build finished");
        }

        Commands::MapSC {
            index,
            geometry,
            read1,
            read2,
            threads,
            output,
            quiet,
        } => {
            let r1_string = read1.join(",");
            let r2_string = read2.join(",");
            let mut args: Vec<CString> = vec![
                CString::new("sc_ref_mapper").unwrap(),
                CString::new("-i").unwrap(),
                CString::new(index).unwrap(),
                CString::new("-g").unwrap(),
                CString::new(geometry).unwrap(),
                CString::new("-1").unwrap(),
                CString::new(r1_string.as_str()).unwrap(),
                CString::new("-2").unwrap(),
                CString::new(r2_string.as_str()).unwrap(),
                CString::new("-t").unwrap(),
                CString::new(threads.to_string()).unwrap(),
                CString::new("-o").unwrap(),
                CString::new(output.as_str()).unwrap(),
            ];

            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            let map_ret = unsafe { run_pesc_sc(args_len, arg_ptrs.as_ptr()) };
            if map_ret != 0 {
                bail!("mapper returned exit code {}; failure", map_ret);
            }
        }

        Commands::MapBulk {
            index,
            read1,
            read2,
            threads,
            output,
            quiet,
        } => {
            let r1_string = read1.join(",");
            let r2_string = read2.join(",");

            let mut args: Vec<CString> = vec![
                CString::new("bulk_ref_mapper").unwrap(),
                CString::new("-i").unwrap(),
                CString::new(index).unwrap(),
                CString::new("-1").unwrap(),
                CString::new(r1_string.as_str()).unwrap(),
                CString::new("-2").unwrap(),
                CString::new(r2_string.as_str()).unwrap(),
                CString::new("-t").unwrap(),
                CString::new(threads.to_string()).unwrap(),
                CString::new("-o").unwrap(),
                CString::new(output.as_str()).unwrap(),
            ];

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
