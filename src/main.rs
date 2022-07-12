use std::ffi::CString;
use std::os::raw::{c_char, c_int};

use anyhow::{bail, Result};
use clap::{arg, command, value_parser, Command};

#[link(name = "pesc_static", kind = "static")]
extern "C" {
    pub fn run_pesc(args: c_int, argsv: *const *const c_char) -> c_int;
}


#[link(name = "build_static", kind = "static")]
extern "C" {
    pub fn run_build(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "cfcore_static", kind = "static", modifiers = "+whole-archive")]
extern "C" {
    pub fn cf_build(args: c_int, argsv: *const *const c_char) -> c_int;
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
                    arg!(-r --reference <REFERENCE> "reference FASTA location"),
                    arg!(-k --klen <K> "length of k-mer to use").value_parser(value_parser!(usize)),
                    arg!(-m --mlen <M> "length of the minimizer use")
                        .value_parser(value_parser!(usize)),
                    arg!(-t --threads <THREADS> "number of threads to use")
                        .value_parser(value_parser!(usize)),
                    arg!(-o --output <OUTPUT> "output file stem"),
                ]),
        )
        .subcommand(
            Command::new("map").about("map reads").args(&[
                arg!(-i --index <INDEX> "input index prefix"),
                arg!(-'1' --read1 <READ1> "path to list of read 1 files")
                    .require_value_delimiter(true),
                arg!(-'2' --read2 <READ2> "path to list of read 2 files")
                    .require_value_delimiter(true),
                arg!(-o --output <OUTPUT> "path to output directory"),
                arg!(-g --geometry <GEO> "geometry of barcode, umi and read"),
                arg!(-t --threads <THREADS> "an interger specifying the number of threads to use")
                    .value_parser(value_parser!(u32)),
            ]),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("build", sub_matches)) => {
            let mut args: Vec<CString> = vec![];
            let k: usize = *sub_matches.get_one("klen").expect("K should be an integer");
            let m: usize = *sub_matches.get_one("mlen").expect("M should be an integer");
            let r: String = sub_matches
                .get_one::<String>("reference")
                .expect("REFERENCE missing")
                .to_string();
            let t: usize = *sub_matches
                .get_one("threads")
                .expect("THREADS should be an integer");
            let o: String = sub_matches
                .get_one::<String>("output")
                .expect("OUTPUT missing")
                .to_string();

            assert!(m < k, "minimizer length ({}) >= k-mer len ({})", m, k);

            let cf_out = o.clone() + "_cfish";
            let mut build_ret = 0;
            /*
            args.push(CString::new("cdbg_builder").unwrap());
            args.push(CString::new("--seq").unwrap());
            args.push(CString::new(r.as_str()).unwrap());
            args.push(CString::new("-k").unwrap());
            args.push(CString::new(k.to_string()).unwrap());
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(cf_out.as_str()).unwrap());
            args.push(CString::new("-t").unwrap());
            args.push(CString::new(t.to_string()).unwrap());
            // format
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
            */
            args.push(CString::new("ref_index_builder").unwrap());
            args.push(CString::new(cf_out.as_str()).unwrap());
            args.push(CString::new(k.to_string()).unwrap());
            args.push(CString::new(m.to_string()).unwrap()); // minimizer length
            args.push(CString::new("--canonical-parsing").unwrap());
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(o.as_str()).unwrap());
            {
                println!("{:?}", args);
                let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                let args_len: c_int = args.len() as c_int;
                build_ret = unsafe { run_build(args_len, arg_ptrs.as_ptr()) };
            }

            if build_ret != 0 {
                bail!("indexer returned exit code {}; failure.", build_ret);
            }
        }
        Some(("map", sub_matches)) => {
            let mut args: Vec<CString> = vec![];
            let r1: String = sub_matches
                .get_many("read1")
                .expect("read2 empty")
                .cloned()
                .collect::<Vec<String>>()
                .join(",");
            let r2: String = sub_matches
                .get_many("read2")
                .expect("read1 empty")
                .cloned()
                .collect::<Vec<String>>()
                .join(",");
            let i: String = sub_matches
                .get_one::<String>("index")
                .expect("INDEX missing")
                .to_string();
            let t: u32 = *sub_matches
                .get_one("threads")
                .expect("THREADS should be an integer");
            let g: String = sub_matches
                .get_one::<String>("geometry")
                .expect("OUTPUT missing")
                .to_string();
            let o: String = sub_matches
                .get_one::<String>("output")
                .expect("OUTPUT missing")
                .to_string();

            args.push(CString::new("ref_mapper").unwrap());
            args.push(CString::new("-i").unwrap());
            args.push(CString::new(i.as_str()).unwrap());
            args.push(CString::new("-g").unwrap());
            args.push(CString::new(g.to_string()).unwrap());

            args.push(CString::new("-1").unwrap());
            args.push(CString::new(r1.as_str()).unwrap());

            args.push(CString::new("-2").unwrap());
            args.push(CString::new(r2.as_str()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(t.to_string()).unwrap());

            args.push(CString::new("-o").unwrap());
            args.push(CString::new(o.as_str()).unwrap());

            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            unsafe { run_pesc(args_len, arg_ptrs.as_ptr()) };
        }
        Some((cmd, &_)) => {
            bail!("Invalid command {}; program exited abnormally.", cmd);
        }
        None => {
            bail!("missing command; program exited abnormally.");
        }
    }

    Ok(())
}
