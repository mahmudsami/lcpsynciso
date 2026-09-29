//! lcpsynciso command line: dispatch to the `predict`, `cluster` and `find-isoforms`
//! subcommands.

use lcpsynciso::{cluster, isoform, predict};

fn print_usage() {
    eprintln!(
        "lcpsynciso — cluster long reads by shared LCP-syncmer seeds, then detect isoforms

USAGE:
    lcpsynciso <command> [options]

COMMANDS:
    predict          fused reads-to-isoforms pipeline (file read once, all in RAM)
    cluster          greedy clustering of reads by shared, level-weighted seeds
    find-isoforms    split each cluster into isoforms by internal structure

Run `lcpsynciso <command> --help` for command-specific options."
    );
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let cmd = match argv.first() {
        Some(c) => c.as_str(),
        None => {
            print_usage();
            std::process::exit(2);
        }
    };
    let rest = &argv[1..];
    match cmd {
        "predict" => predict::run(predict::parse_args(rest)),
        "cluster" => cluster::run(cluster::parse_args(rest)),
        "find-isoforms" | "isoforms" | "isofinder" => {
            isoform::run(isoform::parse_args(rest))
        }
        "-h" | "--help" => print_usage(),
        other => {
            eprintln!("unknown command: {other}\n");
            print_usage();
            std::process::exit(2);
        }
    }
}
