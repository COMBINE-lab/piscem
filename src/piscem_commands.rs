use anyhow::{bail, Result};
use clap::{ArgGroup, Args};
use std::ffi::CString;
use std::path::{Path, PathBuf};
use std::str::FromStr;

trait DefaultMappingParams {
    const MAX_EC_CARD: u32;
    const MAX_HIT_OCC: u32;
    const MAX_HIT_OCC_RECOVER: u32;
    const MAX_READ_OCC: u32;
    const SKIPPING_STRATEGY: &'static str;
}

struct DefaultParams;

impl DefaultMappingParams for DefaultParams {
    const MAX_EC_CARD: u32 = 4096;
    const MAX_HIT_OCC: u32 = 256;
    const MAX_HIT_OCC_RECOVER: u32 = 1024;
    const MAX_READ_OCC: u32 = 2500;
    const SKIPPING_STRATEGY: &'static str = "permissive";
}

/// Trait to produce a proper set of command-line arguments
/// from a populated struct.  There is a single method,
/// `as_argv`, which produces a Vec<CString> that can be parsed
/// and passed to a C function as the `char** argv` parameter.
pub trait AsArgv {
    fn as_argv(&self) -> Result<Vec<CString>>;
}

#[derive(Args, Clone, Debug)]
#[command(arg_required_else_help = true)]
#[command(group(
    ArgGroup::new("ref-input")
    .required(true)
    .args(&["ref_seqs", "ref_lists", "ref_dirs"]),
))]
pub(crate) struct BuildOpts {
    /// ',' separated list of reference FASTA files
    #[arg(short = 's', long, value_delimiter = ',', required = true)]
    pub ref_seqs: Option<Vec<String>>,

    /// ',' separated list of files (each listing input FASTA files)
    #[arg(short = 'l', long, value_delimiter = ',', required = true)]
    pub ref_lists: Option<Vec<String>>,

    /// ',' separated list of directories (all FASTA files in each directory will be indexed,
    /// but not recursively).
    #[arg(short = 'd', long, value_delimiter = ',', required = true)]
    pub ref_dirs: Option<Vec<String>>,

    /// length of k-mer to use
    #[arg(short, long)]
    pub klen: usize,

    /// length of minimizer to use
    #[arg(short, long)]
    pub mlen: usize,

    /// number of threads to use
    #[arg(short, long)]
    pub threads: usize,

    /// output file stem
    #[arg(short, long)]
    pub output: PathBuf,

    /// retain the reduced format GFA files produced by cuttlefish that
    /// describe the reference cDBG (the default is to remove these).
    #[arg(long)]
    pub keep_intermediate_dbg: bool,

    /// working directory where temporary files should be placed.
    #[arg(short = 'w', long, default_value_os_t = PathBuf::from("."))]
    pub work_dir: PathBuf,

    /// overwite an existing index if the output path is the same.
    #[arg(long)]
    pub overwrite: bool,

    /// skip the construction of the equivalence class lookup table
    /// when building the index (not recommended).
    #[arg(long)]
    pub no_ec_table: bool,

    /// path to (optional) ',' sparated list of decoy sequences used to insert poison
    /// k-mer information into the index.
    #[arg(long, value_delimiter = ',')]
    pub decoy_paths: Option<Vec<PathBuf>>,
}

#[derive(Args, Clone, Debug)]
pub(crate) struct MapSCOpts {
    /// input index prefix
    #[arg(short, long)]
    pub index: String,

    /// geometry of barcode, umi and read
    #[arg(short, long)]
    pub geometry: String,

    /// path to list of read 1 files
    #[arg(short = '1', long, value_delimiter = ',', required = true)]
    pub read1: Vec<String>,

    /// path to list of read 2 files
    #[arg(short = '2', long, value_delimiter = ',', required = true)]
    pub read2: Vec<String>,

    /// number of threads to use
    #[arg(short, long, default_value_t = 16)]
    pub threads: usize,

    /// path to output directory
    #[arg(short, long)]
    pub output: String,

    /// do not consider poison k-mers, even if the underlying index contains them.
    /// In this case, the mapping results will be identical to those obtained as if
    /// no poison table was added to the index.
    #[arg(long)]
    pub no_poison: bool,

    /// apply structural constraints when performing mapping.
    #[arg(short = 'c', long)]
    pub struct_constraints: bool,

    /// the skipping strategy to use for k-mer collection
    #[arg(long, default_value = &DefaultParams::SKIPPING_STRATEGY, value_parser = clap::builder::PossibleValuesParser::new(["permissive", "strict"]))]
    pub skipping_strategy: String,

    /// skip checking of the equivalence classes of k-mers that were too
    /// ambiguous to be otherwise considered (passing this flag can speed up
    /// mapping slightly, but may reduce specificity).
    #[arg(long)]
    pub ignore_ambig_hits: bool,

    /// determines the maximum cardinality equivalence class
    /// (number of (txp, orientation status) pairs) to examine (cannot be used with
    /// --ignore-ambig-hits).
    #[arg(
        long,
        short,
        default_value_t = DefaultParams::MAX_EC_CARD,
        conflicts_with = "ignore_ambig_hits",
        help_heading = "Advanced options"
    )]
    pub max_ec_card: u32,

    /// in the first pass, consider only k-mers having <= --max-hit-occ hits.
    #[arg(long, default_value_t = DefaultParams::MAX_HIT_OCC, help_heading = "Advanced options")]
    pub max_hit_occ: u32,

    /// if all k-mers have > --max-hit-occ hits, then make a second pass and consider k-mers
    /// having <= --max-hit-occ-recover hits.
    #[arg(long, default_value_t = DefaultParams::MAX_HIT_OCC_RECOVER, help_heading = "Advanced options")]
    pub max_hit_occ_recover: u32,

    /// reads with more than this number of mappings will not have
    /// their mappings reported.
    #[arg(long, default_value_t = DefaultParams::MAX_READ_OCC, help_heading = "Advanced options")]
    pub max_read_occ: u32,
}

#[derive(Args, Clone, Debug)]
#[command(group(
        ArgGroup::new("read_source")
        .required(true)
        .args(["read1", "reads"])
))]
pub(crate) struct MapBulkOpts {
    /// input index prefix
    #[arg(short, long)]
    pub index: String,

    /// path to list of read 1 files
    #[arg(short = '1', long, value_delimiter = ',', requires = "read2")]
    pub read1: Option<Vec<String>>,

    /// path to list of read 2 files
    #[arg(short = '2', long, value_delimiter = ',', requires = "read1")]
    pub read2: Option<Vec<String>>,

    /// path to list of read unpaired read files
    #[arg(short = 'r', long, value_delimiter = ',', conflicts_with_all = ["read1", "read2"])]
    pub reads: Option<Vec<String>>,

    /// number of threads to use
    #[arg(short, long, default_value_t = 16)]
    pub threads: usize,

    /// path to output directory
    #[arg(short, long)]
    pub output: String,

    /// do not consider poison k-mers, even if the underlying index contains them.
    /// In this case, the mapping results will be identical to those obtained as if
    /// no poison table was added to the index.
    #[arg(long)]
    pub no_poison: bool,

    /// apply structural constraints when performing mapping.
    #[arg(short = 'c', long)]
    pub struct_constraints: bool,

    /// skipping strategy to use for k-mer collection
    #[arg(long, default_value = &DefaultParams::SKIPPING_STRATEGY, value_parser = clap::builder::PossibleValuesParser::new(["permissive", "strict"]))]
    pub skipping_strategy: String,

    /// skip checking of the equivalence classes of k-mers that were too
    /// ambiguous to be otherwise considered (passing this flag can speed up
    /// mapping slightly, but may reduce specificity).
    #[arg(long)]
    pub ignore_ambig_hits: bool,

    /// determines the maximum cardinality equivalence class
    /// (number of (txp, orientation status) pairs) to examine (cannot be used with
    /// --ignore-ambig-hits).
    #[arg(
        long,
        short,
        requires = "check_ambig_hits",
        default_value_t = DefaultParams::MAX_EC_CARD,
        conflicts_with = "ignore_ambig_hits",
        help_heading = "Advanced options"
    )]
    pub max_ec_card: u32,

    /// in the first pass, consider only k-mers having <= --max-hit-occ hits.
    #[arg(long, default_value_t = DefaultParams::MAX_HIT_OCC, help_heading = "Advanced options")]
    pub max_hit_occ: u32,

    /// if all k-mers have > --max-hit-occ hits, then make a second pass and consider k-mers
    /// having <= --max-hit-occ-recover hits.
    #[arg(long, default_value_t = DefaultParams::MAX_HIT_OCC_RECOVER, help_heading = "Advanced options")]
    pub max_hit_occ_recover: u32,

    /// reads with more than this number of mappings will not have
    /// their mappings reported.
    #[arg(long, default_value_t = DefaultParams::MAX_READ_OCC, help_heading = "Advanced options")]
    pub max_read_occ: u32,
}

impl AsArgv for MapSCOpts {
    fn as_argv(&self) -> Result<Vec<CString>> {
        // first check if the relevant index files exist
        let mut idx_suffixes: Vec<String> = vec!["sshash".into(), "ctab".into(), "refinfo".into()];

        if !self.ignore_ambig_hits {
            idx_suffixes.push("ectab".into());
        }

        {
            let idx_path = get_index_path(&self.index)?;
            for s in idx_suffixes {
                let req_file = idx_path.with_extension(s);
                if !req_file.exists() {
                    bail!("To load the index with the specified prefix {}, piscem expects the file {} to exist, but it does not!", &self.index, req_file.display());
                }
            }
        }

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

        Ok(args)
    }
}

fn get_index_path(base: &str) -> Result<PathBuf> {
    if Path::new(base).exists() {
        bail!(
            concat!("The path {} was provided as the base path for the index, but this corresponds ",
                    "to a specific existing file. The provided path should be the file stem (e.g. without the extension)."),
            base);
    }

    if let Some(_ext) = Path::new(base).extension() {
        Ok(PathBuf::from_str(&format!("{}.dummy", base))?)
    } else {
        Ok(PathBuf::from_str(base)?)
    }
}

impl AsArgv for MapBulkOpts {
    fn as_argv(&self) -> Result<Vec<CString>> {
        let mut idx_suffixes: Vec<String> = vec!["sshash".into(), "ctab".into(), "refinfo".into()];

        if !self.ignore_ambig_hits {
            idx_suffixes.push("ectab".into());
        }

        {
            let idx_path = get_index_path(&self.index)?;
            for s in idx_suffixes {
                let req_file = idx_path.with_extension(s);
                if !req_file.exists() {
                    bail!("To load the index with the specified prefix {}, piscem expects the file {} to exist, but it does not!", &self.index, req_file.display());
                }
            }
        }

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

        Ok(args)
    }
}
