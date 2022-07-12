use cmake::Config;
use std::env;

fn main() {
    let custom_cc = env::var("CC");
    let custom_cxx = env::var("CXX");

    let mut cfg_piscem_cpp = Box::new(Config::new("piscem-cpp"));
    let mut cfg_cf = Box::new(Config::new("cuttlefish"));

    (*cfg_cf).define("INSTANCE_COUNT", "32");
    match custom_cc {
        Ok(cc_var) => {
            (*cfg_piscem_cpp).define("CMAKE_C_COMPILER", cc_var.clone());
            (*cfg_cf).define("CMAKE_C_COMPILER", cc_var);
        }
        Err(_) => {}
    }

    match custom_cxx {
        Ok(cxx_var) => {
            (*cfg_piscem_cpp).define("CMAKE_CXX_COMPILER", cxx_var.clone());
            (*cfg_cf).define("CMAKE_CXX_COMPILER", cxx_var);
        }
        Err(_) => {}
    }

    let dst_piscem_cpp = (*cfg_piscem_cpp).build();
    let dst_cf = (*cfg_cf).build();

    println!(
        "cargo:rustc-link-search=native={}",
        dst_cf.join("lib").display()
    );
    println!(
        "cargo:rustc-link-search=native={}",
        dst_piscem_cpp.join("lib").display()
    );
    

    println!("cargo:rustc-link-lib=static=kmc_core");
    println!("cargo:rustc-link-lib=static=pesc_static");
    println!("cargo:rustc-link-lib=static=build_static");
    println!("cargo:rustc-link-lib=static=sshash_static");
    println!("cargo:rustc-link-lib=static=z");
    println!("cargo:rustc-link-lib=static=bz2");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
