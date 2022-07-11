use cmake::Config;
use std::env;

fn main() {
    let custom_cc = env::var("CC");
    let custom_cxx = env::var("CXX");

    let mut cfg = Box::new(Config::new("piscem-cpp"));
    match custom_cc {
        Ok(cc_var) => {
            (*cfg).define("CMAKE_C_COMPILER", cc_var);
        }
        Err(_) => {}
    }

    match custom_cxx {
        Ok(cxx_var) => {
            (*cfg).define("CMAKE_CXX_COMPILER", cxx_var);
        }
        Err(_) => {}
    }

    let dst = (*cfg).build();

    println!(
        "cargo:rustc-link-search=native={}",
        dst.join("lib").display()
    );
    println!("cargo:rustc-link-lib=static=pesc_static");
    println!("cargo:rustc-link-lib=static=build_static");
    println!("cargo:rustc-link-lib=static=sshash_static");
    println!("cargo:rustc-link-lib=static=z");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
