use anyhow::{Result, anyhow, bail};
use clap::{ArgGroup, Args};
use std::path::PathBuf;

trait DefaultMappingParams {
    const MAX_EC_CARD: u32;
    const MAX_HIT_OCC: u32;
    const MAX_HIT_OCC_RECOVER: u32;
    const MAX_READ_OCC: u32;
    const SKIPPING_STRATEGY: &'static str;
    const THRESHOLD: f32;
    const BIN_SIZE: u32;
    const BIN_OVERLAP: u32;
    const BCLEN: u16;
    const END_CACHE_CAPACITY: usize;
}

struct DefaultParams;

impl DefaultMappingParams for DefaultParams {
    const MAX_EC_CARD: u32 = 4096;
    const MAX_HIT_OCC: u32 = 256;
    const MAX_HIT_OCC_RECOVER: u32 = 1024;
    const MAX_READ_OCC: u32 = 2500;
    const SKIPPING_STRATEGY: &'static str = "permissive";
    const THRESHOLD: f32 = 0.7;
    const BIN_SIZE: u32 = 1000;
    const BIN_OVERLAP: u32 = 300;
    const BCLEN: u16 = 16;
    const END_CACHE_CAPACITY: usize = 5_000_000;
}

fn klen_is_good(s: &str) -> Result<usize> {
    let k: usize = s
        .parse()
        .map_err(|_| anyhow!("`{s}` can't be parsed as a number"))?;
    if k > 31 {
        bail!("klen = {k} must be <= 31");
    } else if (k & 1) == 0 {
        bail!("klen = {k} must be odd");
    } else {
        Ok(k)
    }
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
    #[arg(short = 's', long, help_heading = "Input", value_delimiter = ',')]
    pub ref_seqs: Option<Vec<String>>,

    /// ',' separated list of files (each listing input FASTA files)
    #[arg(short = 'l', long, help_heading = "Input", value_delimiter = ',')]
    pub ref_lists: Option<Vec<String>>,

    /// ',' separated list of directories (all FASTA files in each directory will be indexed,
    /// but not recursively).
    #[arg(short = 'd', long, help_heading = "Input", value_delimiter = ',')]
    pub ref_dirs: Option<Vec<String>>,

    /// length of k-mer to use, must be <= 31 and odd
    #[arg(short, long, help_heading = "Index Construction Parameters", default_value_t = 31, value_parser = klen_is_good)]
    pub klen: usize,

    /// length of minimizer to use; must be < `klen`
    #[arg(
        short,
        long,
        help_heading = "Index Construction Parameters",
        default_value_t = 19
    )]
    pub mlen: usize,

    /// number of threads to use
    #[arg(short, long, help_heading = "Index Construction Parameters")]
    pub threads: usize,

    /// output file stem
    #[arg(short, long)]
    pub output: PathBuf,

    /// retain the reduced format GFA files produced by cuttlefish that
    /// describe the reference cDBG (the default is to remove these).
    #[arg(long, help_heading = "Indexing Details")]
    pub keep_intermediate_dbg: bool,

    /// working directory where temporary files should be placed.
    #[arg(short = 'w', long, help_heading = "Indexing Details", default_value_os_t = PathBuf::from("./workdir.noindex"))]
    pub work_dir: PathBuf,

    /// overwite an existing index if the output path is the same.
    #[arg(long, help_heading = "Indexing Details")]
    pub overwrite: bool,

    /// skip the construction of the equivalence class lookup table
    /// when building the index (not recommended).
    #[arg(long, help_heading = "Index Construction Parameters")]
    pub no_ec_table: bool,

    /// If provided (default is not to clip polyA), then reference sequences
    /// ending with polyA tails of length greater than or equal to this value
    /// will be clipped.
    #[arg(long, help_heading = "Index Construction Parameters")]
    pub polya_clip_length: Option<usize>,

    /// path to (optional) ',' sparated list of decoy sequences used to insert poison
    /// k-mer information into the index.
    #[arg(long, value_delimiter = ',')]
    pub decoy_paths: Option<Vec<PathBuf>>,

    /// index construction seed (seed value passed to SSHash index construction; useful if empty
    /// buckets occur).
    #[arg(
        long = "seed",
        help_heading = "Index Construction Parameters",
        default_value_t = 1
    )]
    pub seed: u64,
}

#[derive(Args, Clone, Debug)]
pub(crate) struct MapSCOpts {
    /// input index prefix
    #[arg(short, long, help_heading = "Input")]
    pub index: String,

    /// geometry of barcode, umi and read
    #[arg(short, long)]
    pub geometry: String,

    /// path to a ',' separated list of read 1 files
    #[arg(
        short = '1',
        long,
        help_heading = "Input",
        value_delimiter = ',',
        required = true
    )]
    pub read1: Vec<String>,

    /// path to a ',' separated list of read 2 files
    #[arg(
        short = '2',
        long,
        help_heading = "Input",
        value_delimiter = ',',
        required = true
    )]
    pub read2: Vec<String>,

    /// number of threads to use
    #[arg(short, long, default_value_t = 16)]
    pub threads: usize,

    /// path to output directory
    #[arg(short, long)]
    pub output: PathBuf,

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

    /// includes the positions of each mapped read in the resulting RAD file.
    #[arg(long)]
    pub with_position: bool,

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
    #[arg(short, long, help_heading = "Input")]
    pub index: String,

    /// path to a comma-separated list of read 1 files
    #[arg(
        short = '1',
        long,
        help_heading = "Input",
        value_delimiter = ',',
        requires = "read2"
    )]
    pub read1: Option<Vec<String>>,

    /// path to a ',' separated list of read 2 files
    #[arg(
        short = '2',
        long,
        help_heading = "Input",
        value_delimiter = ',',
        requires = "read1"
    )]
    pub read2: Option<Vec<String>>,

    /// path to a ',' separated list of read unpaired read files
    #[arg(short = 'r', long, help_heading = "Input", value_delimiter = ',', conflicts_with_all = ["read1", "read2"])]
    pub reads: Option<Vec<String>>,

    /// number of threads to use
    #[arg(short, long, default_value_t = 16)]
    pub threads: usize,

    /// path to output directory
    #[arg(short, long)]
    pub output: PathBuf,

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
pub(crate) struct MapSCAtacOpts {
    /// input index prefix
    #[arg(short, long, help_heading = "Input")]
    pub index: String,

    /// path to a ',' separated list of read 1 files
    #[arg(
        short = '1',
        long,
        help_heading = "Input",
        value_delimiter = ',',
        requires = "read2"
    )]
    pub read1: Option<Vec<String>>,

    /// path to a ',' separated list of read 2 files
    #[arg(
        short = '2',
        long,
        help_heading = "Input",
        value_delimiter = ',',
        requires = "read1"
    )]
    pub read2: Option<Vec<String>>,

    #[arg(short = 'r', long, help_heading = "Input", value_delimiter = ',', conflicts_with_all = ["read1", "read2"])]
    pub reads: Option<Vec<String>>,

    /// path to a ',' separated list of barcode files
    #[arg(
        short = 'b',
        long,
        help_heading = "Input",
        value_delimiter = ','
    )]
    pub barcode: Option<Vec<String>>,

    /// number of threads to use
    #[arg(short, long, default_value_t = 16)]
    pub threads: usize,

    /// path to output directory
    #[arg(short, long)]
    pub output: PathBuf,

    /// skip checking of the equivalence classes of k-mers that were too
    /// ambiguous to be otherwise considered (passing this flag can speed up
    /// mapping slightly, but may reduce specificity).
    #[arg(long)]
    pub ignore_ambig_hits: bool,

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

    /// output mappings in sam format
    #[arg(long)]
    pub sam_format: bool,

    /// output mappings in bed format
    #[arg(long)]
    pub bed_format: bool,

    /// use chromosomes as color
    #[arg(long)]
    pub use_chr: bool,

    /// threshold to be considered for pseudoalignment, default set to 0.7
    #[arg(long, default_value_t = DefaultParams::THRESHOLD)]
    pub thr: f32,

    /// size of virtual color, default set to 1000
    #[arg(long, default_value_t = DefaultParams::BIN_SIZE)]
    pub bin_size: u32,

    /// size for bin overlap, default set to 300
    #[arg(long, default_value_t = DefaultParams::BIN_OVERLAP)]
    pub bin_overlap: u32,

    /// do not apply Tn5 shift to mapped positions
    #[arg(long)]
    pub no_tn5_shift: bool,

    /// Check if any mapping kmer exist for a mate which is not mapped,
    /// but there exists mapping for the other read.
    /// If set to true and a mapping kmer exists, then the pair would not be mapped (default false)
    #[arg(long)]
    pub check_kmer_orphan: bool,

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

    /// the length of the barcode sequence
    #[arg(long, default_value_t = DefaultParams::BCLEN, help_heading = "Advanced options")]
    pub bclen: u16,

    /// the capacity of the cache used to provide fast lookup for k-mers at the ends of unitigs
    #[arg(long, default_value_t = DefaultParams::END_CACHE_CAPACITY, help_heading = "Advanced options")]
    pub end_cache_capacity: usize,
}
