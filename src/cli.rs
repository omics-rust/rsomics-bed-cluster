use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_cluster::cluster;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// Cluster overlapping BED intervals — bedtools cluster equivalent.
///
/// Reads a coordinate-sorted BED file and appends a 1-based cluster-ID column
/// to each record. Two records belong to the same cluster when they overlap
/// (or, with -d, when their gap is ≤ d bases).
///
/// Input must be sorted by chromosome then start (and by strand when using -s).
#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-cluster", disable_help_flag = true)]
pub struct Cli {
    /// Input BED file (default: stdin).
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Output file (default: stdout).
    #[arg(long = "out", short = 'o')]
    pub output: Option<PathBuf>,

    /// Maximum gap (bp) allowed between intervals in the same cluster.
    /// 0 = require actual overlap (default).
    #[arg(short = 'd', long = "dist", default_value = "0")]
    pub dist: i64,

    /// Require same strand for clustering (input must be sorted by strand too).
    #[arg(short = 's', long = "strand")]
    pub strand: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };
        match self.input {
            Some(ref p) => {
                let f = File::open(p).map_err(RsomicsError::Io)?;
                cluster(f, out, self.dist, self.strand)
            }
            None => {
                let stdin = io::stdin();
                cluster(stdin.lock(), out, self.dist, self.strand)
            }
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Cluster overlapping BED intervals (bedtools cluster equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdin"),
                description: "Input BED file (default: stdin)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "out",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output file (default: stdout)",
                why_default: None,
            },
            FlagSpec {
                short: Some('d'),
                long: "dist",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("i64"),
                required: false,
                default: Some("0"),
                description: "Maximum gap (bp) between intervals in the same cluster; 0 = require actual overlap",
                why_default: None,
            },
            FlagSpec {
                short: Some('s'),
                long: "strand",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Require same strand for clustering",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Cluster a sorted BED file",
            command: "rsomics-bed-cluster -i sorted.bed",
        },
        Example {
            description: "Cluster with 50 bp gap tolerance",
            command: "rsomics-bed-cluster -d 50 -i sorted.bed",
        },
        Example {
            description: "Strand-specific clustering",
            command: "rsomics-bed-cluster -s -i sorted.bed",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
