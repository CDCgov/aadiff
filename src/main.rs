#![feature(let_chains)]
#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

use clap::Parser;
use either::Either;
use std::{
    fs::OpenOptions,
    io::{BufReader, BufWriter, Write, stdin, stdout},
    path::PathBuf,
};
use zoe::prelude::*;

// TODO: we want to change the delimiter

/// Holds arguments for the APD CLI.
#[derive(Debug, Parser)]
#[command(about = "Tool for calculating amino acid difference tables")]
pub struct APDArgs {
    #[arg(short = 'i', long)]
    /// Optional input fasta
    input_fasta: Option<PathBuf>,

    #[arg(short = 'o', long)]
    /// Optional output delimited
    output_tsv: Option<PathBuf>,

    #[arg(short = 't', long)]
    /// Turns on simultaneous multi-threading
    multi_threaded: bool,

    #[arg(short = 'e', long)]
    /// Use unix line-endings instead of Windows ones
    unix_line_endings: bool,
}

/// CLI entry point
fn main() {
    let args = APDArgs::parse();
    let line_ending = if args.unix_line_endings { "" } else { "\r" };

    let mut reader = if let Some(ref file_path) = args.input_fasta {
        FastaReader::new(BufReader::new(Either::Left(
            OpenOptions::new().read(true).open(file_path).expect("File opening error"),
        )))
    } else {
        FastaReader::new(BufReader::new(Either::Right(stdin())))
    };

    let mut writer = if let Some(ref file_path) = args.output_tsv {
        BufWriter::new(Either::Left(
            OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .open(file_path)
                .expect("File write error"),
        ))
    } else {
        BufWriter::new(Either::Right(stdout()))
    };

    let Some(Ok(dna_reference)) = reader.next() else {
        eprintln!("No first record available!");
        std::process::exit(1);
    };
    // TODO: recode to IUPAC if not in IUPAC.
    let reference = {
        let mut r = dna_reference.translate();
        r.sequence.replace_all_bytes(b'~', b'X');
        r
    };

    let other_sequences = reader
        .map(|record|
            // TODO: recode to IUPAC
            // TODO: don't translate, instead defer until later
            record.map(|r| {let mut aa = r.translate(); aa.sequence.replace_all_bytes(b'~',b'X'); aa}))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_die("Could not process other data.");

    let mut buffer = format!(",{}", reference.name);
    for query_header in other_sequences.iter().map(|f| f.name.as_str()) {
        buffer.push(',');
        buffer.push_str(query_header);
    }
    writeln!(&mut writer, "{buffer}{line_ending}").unwrap_or_fail();

    for (i, ref_aa) in reference.sequence.into_iter().enumerate() {
        let mut differences_found = false;
        buffer.clear();

        // TODO: grab codon in the future instead and consider not resolvable ambiguations
        // TODO: check all sequences have same length earlier
        // TODO: make sure you use uppercase for both
        // TODO: do we annotate deletions?
        for query_aa in other_sequences.iter().map(|f| f.sequence[i]) {
            if ref_aa != query_aa && query_aa.is_ascii_alphabetic() {
                buffer.push_str(r#",""#);
                buffer.push(query_aa as char);
                buffer.push('"');
                differences_found = true;
            } else {
                buffer.push(',');
            }
        }

        if differences_found {
            writeln!(&mut writer, "{p},{aa}{buffer}{line_ending}", p = i + 1, aa = ref_aa as char).unwrap_or_fail();
        }
    }

    writer.flush().unwrap_or_fail();
}
