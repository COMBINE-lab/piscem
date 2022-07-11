use std::ffi::CString;
use std::os::raw::{c_char, c_int};

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

fn main() -> Result<(), anyhow::Error> {
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
