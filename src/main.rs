use std::os::raw::{c_char, c_int};
use std::ffi::CString;

#[link(name = "pesc_static", kind = "static")]
extern "C" {
    pub fn run_pesc(args: c_int, argsv: *const *const c_char) -> c_int;
}

fn main() {
    println!("Hello, world, from Rust!");
    // from
    // https://stackoverflow.com/questions/69437925/problem-with-calling-c-function-that-receive-command-line-arguments-from-rust
    let args = std::env::args()
        .map(|arg| CString::new(arg).unwrap())
        .collect::<Vec<CString>>();
    let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
    let args_len: c_int = args.len() as c_int;
    unsafe { run_pesc(args_len, arg_ptrs.as_ptr()) };
}
