use std::ffi::CString;
use std::os::raw::{c_char, c_int};

use anyhow::{bail, Result};
use clap::{Parser, Subcommand};

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
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Index a reference sequence
    #[clap(arg_required_else_help = true)]
    Build {
        /// reference FASTA location
        #[clap(short, long, value_parser)]
        reference: String,

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
        output: String,

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
        #[clap(short = '1', long, value_parser)]
        read1: Vec<String>,

        /// path to list of read 1 files
        #[clap(short = '2', long, value_parser)]
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
        #[clap(short = '1', long, value_parser)]
        read1: Vec<String>,

        /// path to list of read 1 files
        #[clap(short = '2', long, value_parser)]
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

    match cli_args.command {
        Commands::Build {
            reference,
            klen,
            mlen,
            threads,
            output,
            quiet,
        } => {
            assert!(
                mlen < klen,
                "minimizer length ({}) >= k-mer len ({})",
                mlen,
                klen
            );

            let mut args: Vec<CString> = vec![];

            let cf_out = output.clone() + "_cfish";
            let mut build_ret;

            args.push(CString::new("cdbg_builder").unwrap());
            args.push(CString::new("--seq").unwrap());
            args.push(CString::new(reference.as_str()).unwrap());
            args.push(CString::new("-k").unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(cf_out.as_str()).unwrap());
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
            args.push(CString::new(cf_out.as_str()).unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new(mlen.to_string()).unwrap()); // minimizer length
            args.push(CString::new("--canonical-parsing").unwrap());
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(output.as_str()).unwrap());
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
        },

        Commands::MapSC {
            index,
            geometry,
            read1,
            read2,
            threads,
            output,
            quiet,
        } => {
            let mut args: Vec<CString> = vec![];

            args.push(CString::new("sc_ref_mapper").unwrap());
            args.push(CString::new("-i").unwrap());
            args.push(CString::new(index).unwrap());
            args.push(CString::new("-g").unwrap());
            args.push(CString::new(geometry).unwrap());

            args.push(CString::new("-1").unwrap());
            let r1_string = read1.join(",");
            args.push(CString::new(r1_string.as_str()).unwrap());

            args.push(CString::new("-2").unwrap());
            let r2_string = read2.join(",");
            args.push(CString::new(r2_string.as_str()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(threads.to_string()).unwrap());

            args.push(CString::new("-o").unwrap());
            args.push(CString::new(output.as_str()).unwrap());
            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            let map_ret = unsafe { run_pesc_sc(args_len, arg_ptrs.as_ptr()) };
            if map_ret != 0 {
                bail!("mapper returned exit code {}; failure", map_ret);
            }
        },

        Commands::MapBulk {
            index,
            read1,
            read2,
            threads,
            output,
            quiet,
        } => {
            let mut args: Vec<CString> = vec![];

            args.push(CString::new("bulk_ref_mapper").unwrap());
            args.push(CString::new("-i").unwrap());
            args.push(CString::new(index).unwrap());

            args.push(CString::new("-1").unwrap());
            let r1_string = read1.join(",");
            args.push(CString::new(r1_string.as_str()).unwrap());

            args.push(CString::new("-2").unwrap());
            let r2_string = read2.join(",");
            args.push(CString::new(r2_string.as_str()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(threads.to_string()).unwrap());

            args.push(CString::new("-o").unwrap());
            args.push(CString::new(output.as_str()).unwrap());
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
