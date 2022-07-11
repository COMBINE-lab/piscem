use std::ffi::CString;
use std::os::raw::{c_char, c_int};

use clap::{arg, command, Command}; 
use anyhow::{bail, Result};

#[link(name = "pesc_static", kind = "static")]
extern "C" {
    pub fn run_pesc(args: c_int, argsv: *const *const c_char) -> c_int;
    pub fn run_build(args: c_int, argsv: *const *const c_char) -> c_int;
}

fn print_help() {
    eprintln!("piscem version 0.0.1: ");
    eprintln!("valid commands are {{build, map}}, for more information ");
    eprintln!("try build -h or map -h.");
}


/*
Options:
  -h,--help                   Print this help message and exit
  -i,--index TEXT REQUIRED    input index prefix
  -1,--read1 TEXT ... REQUIRED
                              path to list of read 1 files
  -2,--read2 TEXT ... REQUIRED
                              path to list of read 2 files
  -o,--output TEXT REQUIRED   path to output directory
  -g,--geometry TEXT REQUIRED geometry of barcode, umi and read
  -t,--threads UINT [16]      An integer that specifies the number of threads to use
 */

fn main() -> Result<(), anyhow::Error> {

    let matches = command!()
        .propagate_version(true)
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("build")
            .about("build the piscem index")
            .args(&[
                  arg!(-r --ref <REF> "reference FASTA location"),
                  arg!(-k --klen <K> "length of k-mer to use"),
                  arg!(-t --threads <THREADS> "number of threads to use"),
                  arg!(-o --output <OUTPUT> "output file stem")
            ])
        )
        .subcommand(
            Command::new("map")
            .about("map reads")
            .args(&[
                  arg!(-i --index <INDEX> "input index prefix"),
                  arg!(-'1' --read1 <READ1> "path to list of read 1 files").require_value_delimiter(true), 
                  arg!(-'2' --read2 <READ2> "path to list of read 2 files").require_value_delimiter(true), 
                  arg!(-o --output <OUTPUT> "path to output directory"),
                  arg!(-g --geometry <GEO> "geometry of barcode, umi and read"),
                  arg!(-t --threads <THREADS> "an interger specifying the number of threads to use")
            ])
        ).get_matches();

    if std::env::args().len() < 2 {
        print_help();
        bail!("program exited abnormally.");
    } else {
        // treat the 1st argument as the "command"
        let cmd = std::env::args().nth(1).expect("no command given");
        // from
        // https://stackoverflow.com/questions/69437925/problem-with-calling-c-function-that-receive-command-line-arguments-from-rust
        let mut args = std::env::args()
            .map(|arg| CString::new(arg).unwrap())
            .collect::<Vec<CString>>();
        _ = args.remove(1);
        let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
        let args_len: c_int = args.len() as c_int;

        match cmd.as_str() {
            "map" => {
                unsafe { run_pesc(args_len, arg_ptrs.as_ptr()) };
            }
            "build" => {
                unsafe { run_build(args_len, arg_ptrs.as_ptr()) };
            }
            "--version" => {
                eprintln!("piscem version 0.0.1");
            }
            "-h" | "--help" => {
                print_help();
            }
            c => {
                eprintln!(
                    "{} is not a valid command; valid commands are {{build, map}}.",
                    c
                );
                bail!("program exited abnormally.");
            }
        }
    }
    Ok(())
}
