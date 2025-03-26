#![feature(let_chains)]
#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

use clap::Parser;
use either::Either;
use std::{
    collections::HashMap,
    fs::OpenOptions,
    io::{BufReader, BufWriter, Write, stdin, stdout},
    ops::Range,
    path::PathBuf,
    sync::LazyLock,
};
use zoe::{data::fasta::FastaNT, prelude::*};

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

    // This appeared to be available in the legacy program.
    #[arg(short = 'r', long)]
    /// Restrict to non-ambiguous alignable regions, pairwise.
    restrict_to_pairwise_alignable: bool,

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
    let ref_range = get_valid_range(&reference.sequence, args.restrict_to_pairwise_alignable);

    let other_sequences = reader
        .map(|record|
            // TODO: recode to IUPAC
            // TODO: don't translate, instead defer until later
            record.map(|r| {
                let FastaNT { name, sequence } = r.into();
                let mut residues: AminoAcids = sequence.to_aa_iter().collect();
                residues.replace_all_bytes(b'~',b'X');
                let valid_range = get_valid_range(&residues, args.restrict_to_pairwise_alignable);

                ValidSeq {
                    name, residues,
                    codons: sequence,
                    valid_range
                }
              }))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_die("Could not process other data.");

    let mut buffer = format!(",{}", reference.name);
    for query_header in other_sequences.iter().map(|f| f.name.as_str()) {
        buffer.push(',');
        buffer.push_str(query_header);
    }
    writeln!(&mut writer, "{buffer}{line_ending}").unwrap_or_fail();

    for (i, &ref_aa) in reference.sequence[ref_range].iter().enumerate() {
        let mut differences_found = false;
        buffer.clear();

        // TODO: make sure you use uppercase for both
        for ValidSeq {
            name: _,
            residues,
            codons,
            valid_range,
        } in other_sequences.iter()
        {
            let aa = residues[i];
            let codon = [codons[i * 3], codons[i * 3 + 1], codons[i * 3 + 2]];

            if valid_range.contains(&i) && ref_aa != aa {
                buffer.push_str(r#",""#);

                if aa == b'-' {
                    buffer.push_str("del");
                } else if aa == b'X'
                    && let Some(degen_aa) = GC3.get(&codon)
                {
                    // We currently support degeneracy up to 3 distinct as beyond that it is kind of useless.
                    buffer.push_str(degen_aa);
                } else {
                    buffer.push(aa as char);
                }
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

struct ValidSeq {
    name:        String,
    residues:    AminoAcids,
    codons:      Nucleotides,
    valid_range: std::ops::Range<usize>,
}

fn get_valid_range(aa: &AminoAcids, restrict: bool) -> Range<usize> {
    if restrict {
        let (Some(s), Some(e)) = (
            aa.iter().position(|&aa| aa != b'X' && aa != b'-'),
            aa.iter().rposition(|&aa| aa != b'X' && aa != b'-'),
        ) else {
            eprintln!("Sequence doesn't contain valid data for comparison.");
            std::process::exit(1);
        };

        s..e + 1
    } else {
        0..aa.len()
    }
}

macro_rules! fill_gc3 {
    ($( $key: expr => $val: expr ),*) => {{
        let mut map = HashMap::new();
        $( map.insert(*$key, $val); )*
        map
   }}
}

static GC3: LazyLock<HashMap<[u8; 3], &'static str>> = LazyLock::new(|| {
    fill_gc3!(
        b"AAB" => "K/N",   b"AAD" => "K/N",   b"AAH" => "K/N",   b"AAK" => "K/N",   b"AAM" => "K/N",
        b"AAN" => "K/N",   b"AAS" => "K/N",   b"AAV" => "K/N",   b"AAW" => "K/N",   b"ABA" => "I/R/T",
        b"ABC" => "I/S/T", b"ABG" => "M/R/T", b"ABT" => "I/S/T", b"ABY" => "I/S/T", b"ADA" => "I/K/R",
        b"ADC" => "I/N/S", b"ADG" => "K/M/R", b"ADT" => "I/N/S", b"ADY" => "I/N/S", b"AGB" => "R/S",
        b"AGD" => "R/S",   b"AGH" => "R/S",   b"AGK" => "R/S",   b"AGM" => "R/S",   b"AGN" => "R/S",
        b"AGS" => "R/S",   b"AGV" => "R/S",   b"AGW" => "R/S",   b"AHA" => "I/K/T", b"AHC" => "I/N/T",
        b"AHG" => "K/M/T", b"AHT" => "I/N/T", b"AHY" => "I/N/T", b"AKA" => "I/R",   b"AKC" => "I/S",
        b"AKG" => "M/R",   b"AKH" => "I/R/S", b"AKM" => "I/R/S", b"AKR" => "I/M/R", b"AKT" => "I/S",
        b"AKW" => "I/R/S", b"AKY" => "I/S",   b"AMA" => "K/T",   b"AMB" => "K/N/T", b"AMC" => "N/T",
        b"AMD" => "K/N/T", b"AMG" => "K/T",   b"AMH" => "K/N/T", b"AMK" => "K/N/T", b"AMM" => "K/N/T",
        b"AMN" => "K/N/T", b"AMR" => "K/T",   b"AMS" => "K/N/T", b"AMT" => "N/T",   b"AMV" => "K/N/T",
        b"AMW" => "K/N/T", b"AMY" => "N/T",   b"ARA" => "K/R",   b"ARC" => "N/S",   b"ARG" => "K/R",
        b"ARR" => "K/R",   b"ART" => "N/S",   b"ARY" => "N/S",   b"ASA" => "R/T",   b"ASB" => "R/S/T",
        b"ASC" => "S/T",   b"ASD" => "R/S/T", b"ASG" => "R/T",   b"ASH" => "R/S/T", b"ASK" => "R/S/T",
        b"ASM" => "R/S/T", b"ASN" => "R/S/T", b"ASR" => "R/T",   b"ASS" => "R/S/T", b"AST" => "S/T",
        b"ASV" => "R/S/T", b"ASW" => "R/S/T", b"ASY" => "S/T",   b"ATB" => "I/M",   b"ATD" => "I/M",
        b"ATK" => "I/M",   b"ATN" => "I/M",   b"ATR" => "I/M",   b"ATS" => "I/M",   b"ATV" => "I/M",
        b"AVA" => "K/R/T", b"AVC" => "N/S/T", b"AVG" => "K/R/T", b"AVR" => "K/R/T", b"AVT" => "N/S/T",
        b"AVY" => "N/S/T", b"AWA" => "I/K",   b"AWC" => "I/N",   b"AWG" => "K/M",   b"AWH" => "I/K/N",
        b"AWM" => "I/K/N", b"AWR" => "I/K/M", b"AWT" => "I/N",   b"AWW" => "I/K/N", b"AWY" => "I/N",
        b"AYA" => "I/T",   b"AYB" => "I/M/T", b"AYC" => "I/T",   b"AYD" => "I/M/T", b"AYG" => "M/T",
        b"AYH" => "I/T",   b"AYK" => "I/M/T", b"AYM" => "I/T",   b"AYN" => "I/M/T", b"AYR" => "I/M/T",
        b"AYS" => "I/M/T", b"AYT" => "I/T",   b"AYV" => "I/M/T", b"AYW" => "I/T",   b"AYY" => "I/T",
        b"BAA" => "E/Q/*", b"BAC" => "D/H/Y", b"BAG" => "E/Q/*", b"BAR" => "E/Q/*", b"BAT" => "D/H/Y",
        b"BAY" => "D/H/Y", b"BCA" => "A/P/S", b"BCB" => "A/P/S", b"BCC" => "A/P/S", b"BCD" => "A/P/S",
        b"BCG" => "A/P/S", b"BCH" => "A/P/S", b"BCK" => "A/P/S", b"BCM" => "A/P/S", b"BCN" => "A/P/S",
        b"BCR" => "A/P/S", b"BCS" => "A/P/S", b"BCT" => "A/P/S", b"BCV" => "A/P/S", b"BCW" => "A/P/S",
        b"BCY" => "A/P/S", b"BGA" => "G/R/*", b"BGC" => "C/G/R", b"BGG" => "G/R/W", b"BGT" => "C/G/R",
        b"BGY" => "C/G/R", b"BTA" => "L/V",   b"BTB" => "F/L/V", b"BTC" => "F/L/V", b"BTD" => "F/L/V",
        b"BTG" => "L/V",   b"BTH" => "F/L/V", b"BTK" => "F/L/V", b"BTM" => "F/L/V", b"BTN" => "F/L/V",
        b"BTR" => "L/V",   b"BTS" => "F/L/V", b"BTT" => "F/L/V", b"BTV" => "F/L/V", b"BTW" => "F/L/V",
        b"BTY" => "F/L/V", b"CAB" => "H/Q",   b"CAD" => "H/Q",   b"CAH" => "H/Q",   b"CAK" => "H/Q",
        b"CAM" => "H/Q",   b"CAN" => "H/Q",   b"CAS" => "H/Q",   b"CAV" => "H/Q",   b"CAW" => "H/Q",
        b"CBA" => "L/P/R", b"CBB" => "L/P/R", b"CBC" => "L/P/R", b"CBD" => "L/P/R", b"CBG" => "L/P/R",
        b"CBH" => "L/P/R", b"CBK" => "L/P/R", b"CBM" => "L/P/R", b"CBN" => "L/P/R", b"CBR" => "L/P/R",
        b"CBS" => "L/P/R", b"CBT" => "L/P/R", b"CBV" => "L/P/R", b"CBW" => "L/P/R", b"CBY" => "L/P/R",
        b"CDA" => "L/Q/R", b"CDC" => "H/L/R", b"CDG" => "L/Q/R", b"CDR" => "L/Q/R", b"CDT" => "H/L/R",
        b"CDY" => "H/L/R", b"CHA" => "L/P/Q", b"CHC" => "H/L/P", b"CHG" => "L/P/Q", b"CHR" => "L/P/Q",
        b"CHT" => "H/L/P", b"CHY" => "H/L/P", b"CKA" => "L/R",   b"CKB" => "L/R",   b"CKC" => "L/R",
        b"CKD" => "L/R",   b"CKG" => "L/R",   b"CKH" => "L/R",   b"CKK" => "L/R",   b"CKM" => "L/R",
        b"CKN" => "L/R",   b"CKR" => "L/R",   b"CKS" => "L/R",   b"CKT" => "L/R",   b"CKV" => "L/R",
        b"CKW" => "L/R",   b"CKY" => "L/R",   b"CMA" => "P/Q",   b"CMB" => "H/P/Q", b"CMC" => "H/P",
        b"CMD" => "H/P/Q", b"CMG" => "P/Q",   b"CMH" => "H/P/Q", b"CMK" => "H/P/Q", b"CMM" => "H/P/Q",
        b"CMN" => "H/P/Q", b"CMR" => "P/Q",   b"CMS" => "H/P/Q", b"CMT" => "H/P",   b"CMV" => "H/P/Q",
        b"CMW" => "H/P/Q", b"CMY" => "H/P",   b"CRA" => "Q/R",   b"CRB" => "H/Q/R", b"CRC" => "H/R",
        b"CRD" => "H/Q/R", b"CRG" => "Q/R",   b"CRH" => "H/Q/R", b"CRK" => "H/Q/R", b"CRM" => "H/Q/R",
        b"CRN" => "H/Q/R", b"CRR" => "Q/R",   b"CRS" => "H/Q/R", b"CRT" => "H/R",   b"CRV" => "H/Q/R",
        b"CRW" => "H/Q/R", b"CRY" => "H/R",   b"CSA" => "P/R",   b"CSB" => "P/R",   b"CSC" => "P/R",
        b"CSD" => "P/R",   b"CSG" => "P/R",   b"CSH" => "P/R",   b"CSK" => "P/R",   b"CSM" => "P/R",
        b"CSN" => "P/R",   b"CSR" => "P/R",   b"CSS" => "P/R",   b"CST" => "P/R",   b"CSV" => "P/R",
        b"CSW" => "P/R",   b"CSY" => "P/R",   b"CVA" => "P/Q/R", b"CVC" => "H/P/R", b"CVG" => "P/Q/R",
        b"CVR" => "P/Q/R", b"CVT" => "H/P/R", b"CVY" => "H/P/R", b"CWA" => "L/Q",   b"CWB" => "H/L/Q",
        b"CWC" => "H/L",   b"CWD" => "H/L/Q", b"CWG" => "L/Q",   b"CWH" => "H/L/Q", b"CWK" => "H/L/Q",
        b"CWM" => "H/L/Q", b"CWN" => "H/L/Q", b"CWR" => "L/Q",   b"CWS" => "H/L/Q", b"CWT" => "H/L",
        b"CWV" => "H/L/Q", b"CWW" => "H/L/Q", b"CWY" => "H/L",   b"CYA" => "L/P",   b"CYB" => "L/P",
        b"CYC" => "L/P",   b"CYD" => "L/P",   b"CYG" => "L/P",   b"CYH" => "L/P",   b"CYK" => "L/P",
        b"CYM" => "L/P",   b"CYN" => "L/P",   b"CYR" => "L/P",   b"CYS" => "L/P",   b"CYT" => "L/P",
        b"CYV" => "L/P",   b"CYW" => "L/P",   b"CYY" => "L/P",   b"DAA" => "E/K/*", b"DAC" => "D/N/Y",
        b"DAG" => "E/K/*", b"DAR" => "E/K/*", b"DAT" => "D/N/Y", b"DAY" => "D/N/Y", b"DCA" => "A/S/T",
        b"DCB" => "A/S/T", b"DCC" => "A/S/T", b"DCD" => "A/S/T", b"DCG" => "A/S/T", b"DCH" => "A/S/T",
        b"DCK" => "A/S/T", b"DCM" => "A/S/T", b"DCN" => "A/S/T", b"DCR" => "A/S/T", b"DCS" => "A/S/T",
        b"DCT" => "A/S/T", b"DCV" => "A/S/T", b"DCW" => "A/S/T", b"DCY" => "A/S/T", b"DGA" => "G/R/*",
        b"DGC" => "C/G/S", b"DGG" => "G/R/W", b"DGT" => "C/G/S", b"DGY" => "C/G/S", b"DTA" => "I/L/V",
        b"DTC" => "F/I/V", b"DTG" => "L/M/V", b"DTT" => "F/I/V", b"DTY" => "F/I/V", b"GAB" => "D/E",
        b"GAD" => "D/E",   b"GAH" => "D/E",   b"GAK" => "D/E",   b"GAM" => "D/E",   b"GAN" => "D/E",
        b"GAS" => "D/E",   b"GAV" => "D/E",   b"GAW" => "D/E",   b"GBA" => "A/G/V", b"GBB" => "A/G/V",
        b"GBC" => "A/G/V", b"GBD" => "A/G/V", b"GBG" => "A/G/V", b"GBH" => "A/G/V", b"GBK" => "A/G/V",
        b"GBM" => "A/G/V", b"GBN" => "A/G/V", b"GBR" => "A/G/V", b"GBS" => "A/G/V", b"GBT" => "A/G/V",
        b"GBV" => "A/G/V", b"GBW" => "A/G/V", b"GBY" => "A/G/V", b"GDA" => "E/G/V", b"GDC" => "D/G/V",
        b"GDG" => "E/G/V", b"GDR" => "E/G/V", b"GDT" => "D/G/V", b"GDY" => "D/G/V", b"GHA" => "A/E/V",
        b"GHC" => "A/D/V", b"GHG" => "A/E/V", b"GHR" => "A/E/V", b"GHT" => "A/D/V", b"GHY" => "A/D/V",
        b"GKA" => "G/V",   b"GKB" => "G/V",   b"GKC" => "G/V",   b"GKD" => "G/V",   b"GKG" => "G/V",
        b"GKH" => "G/V",   b"GKK" => "G/V",   b"GKM" => "G/V",   b"GKN" => "G/V",   b"GKR" => "G/V",
        b"GKS" => "G/V",   b"GKT" => "G/V",   b"GKV" => "G/V",   b"GKW" => "G/V",   b"GKY" => "G/V",
        b"GMA" => "A/E",   b"GMB" => "A/D/E", b"GMC" => "A/D",   b"GMD" => "A/D/E", b"GMG" => "A/E",
        b"GMH" => "A/D/E", b"GMK" => "A/D/E", b"GMM" => "A/D/E", b"GMN" => "A/D/E", b"GMR" => "A/E",
        b"GMS" => "A/D/E", b"GMT" => "A/D",   b"GMV" => "A/D/E", b"GMW" => "A/D/E", b"GMY" => "A/D",
        b"GRA" => "E/G",   b"GRB" => "D/E/G", b"GRC" => "D/G",   b"GRD" => "D/E/G", b"GRG" => "E/G",
        b"GRH" => "D/E/G", b"GRK" => "D/E/G", b"GRM" => "D/E/G", b"GRN" => "D/E/G", b"GRR" => "E/G",
        b"GRS" => "D/E/G", b"GRT" => "D/G",   b"GRV" => "D/E/G", b"GRW" => "D/E/G", b"GRY" => "D/G",
        b"GSA" => "A/G",   b"GSB" => "A/G",   b"GSC" => "A/G",   b"GSD" => "A/G",   b"GSG" => "A/G",
        b"GSH" => "A/G",   b"GSK" => "A/G",   b"GSM" => "A/G",   b"GSN" => "A/G",   b"GSR" => "A/G",
        b"GSS" => "A/G",   b"GST" => "A/G",   b"GSV" => "A/G",   b"GSW" => "A/G",   b"GSY" => "A/G",
        b"GVA" => "A/E/G", b"GVC" => "A/D/G", b"GVG" => "A/E/G", b"GVR" => "A/E/G", b"GVT" => "A/D/G",
        b"GVY" => "A/D/G", b"GWA" => "E/V",   b"GWB" => "D/E/V", b"GWC" => "D/V",   b"GWD" => "D/E/V",
        b"GWG" => "E/V",   b"GWH" => "D/E/V", b"GWK" => "D/E/V", b"GWM" => "D/E/V", b"GWN" => "D/E/V",
        b"GWR" => "E/V",   b"GWS" => "D/E/V", b"GWT" => "D/V",   b"GWV" => "D/E/V", b"GWW" => "D/E/V",
        b"GWY" => "D/V",   b"GYA" => "A/V",   b"GYB" => "A/V",   b"GYC" => "A/V",   b"GYD" => "A/V",
        b"GYG" => "A/V",   b"GYH" => "A/V",   b"GYK" => "A/V",   b"GYM" => "A/V",   b"GYN" => "A/V",
        b"GYR" => "A/V",   b"GYS" => "A/V",   b"GYT" => "A/V",   b"GYV" => "A/V",   b"GYW" => "A/V",
        b"GYY" => "A/V",   b"HAA" => "K/Q/*", b"HAC" => "H/N/Y", b"HAG" => "K/Q/*", b"HAR" => "K/Q/*",
        b"HAT" => "H/N/Y", b"HAY" => "H/N/Y", b"HCA" => "P/S/T", b"HCB" => "P/S/T", b"HCC" => "P/S/T",
        b"HCD" => "P/S/T", b"HCG" => "P/S/T", b"HCH" => "P/S/T", b"HCK" => "P/S/T", b"HCM" => "P/S/T",
        b"HCN" => "P/S/T", b"HCR" => "P/S/T", b"HCS" => "P/S/T", b"HCT" => "P/S/T", b"HCV" => "P/S/T",
        b"HCW" => "P/S/T", b"HCY" => "P/S/T", b"HGA" => "R/*",   b"HGC" => "C/R/S", b"HGG" => "R/W",
        b"HGR" => "R/W/*", b"HGT" => "C/R/S", b"HGY" => "C/R/S", b"HTA" => "I/L",   b"HTC" => "F/I/L",
        b"HTG" => "L/M",   b"HTH" => "F/I/L", b"HTM" => "F/I/L", b"HTR" => "I/L/M", b"HTT" => "F/I/L",
        b"HTW" => "F/I/L", b"HTY" => "F/I/L", b"KAA" => "E/*",   b"KAC" => "D/Y",   b"KAG" => "E/*",
        b"KAR" => "E/*",   b"KAT" => "D/Y",   b"KAY" => "D/Y",   b"KCA" => "A/S",   b"KCB" => "A/S",
        b"KCC" => "A/S",   b"KCD" => "A/S",   b"KCG" => "A/S",   b"KCH" => "A/S",   b"KCK" => "A/S",
        b"KCM" => "A/S",   b"KCN" => "A/S",   b"KCR" => "A/S",   b"KCS" => "A/S",   b"KCT" => "A/S",
        b"KCV" => "A/S",   b"KCW" => "A/S",   b"KCY" => "A/S",   b"KGA" => "G/*",   b"KGB" => "C/G/W",
        b"KGC" => "C/G",   b"KGG" => "G/W",   b"KGH" => "C/G/*", b"KGK" => "C/G/W", b"KGM" => "C/G/*",
        b"KGR" => "G/W/*", b"KGS" => "C/G/W", b"KGT" => "C/G",   b"KGW" => "C/G/*", b"KGY" => "C/G",
        b"KRA" => "E/G/*", b"KTA" => "L/V",   b"KTB" => "F/L/V", b"KTC" => "F/V",   b"KTD" => "F/L/V",
        b"KTG" => "L/V",   b"KTH" => "F/L/V", b"KTK" => "F/L/V", b"KTM" => "F/L/V", b"KTN" => "F/L/V",
        b"KTR" => "L/V",   b"KTS" => "F/L/V", b"KTT" => "F/V",   b"KTV" => "F/L/V", b"KTW" => "F/L/V",
        b"KTY" => "F/V",   b"MAA" => "K/Q",   b"MAC" => "H/N",   b"MAG" => "K/Q",   b"MAR" => "K/Q",
        b"MAT" => "H/N",   b"MAY" => "H/N",   b"MCA" => "P/T",   b"MCB" => "P/T",   b"MCC" => "P/T",
        b"MCD" => "P/T",   b"MCG" => "P/T",   b"MCH" => "P/T",   b"MCK" => "P/T",   b"MCM" => "P/T",
        b"MCN" => "P/T",   b"MCR" => "P/T",   b"MCS" => "P/T",   b"MCT" => "P/T",   b"MCV" => "P/T",
        b"MCW" => "P/T",   b"MCY" => "P/T",   b"MGB" => "R/S",   b"MGC" => "R/S",   b"MGD" => "R/S",
        b"MGH" => "R/S",   b"MGK" => "R/S",   b"MGM" => "R/S",   b"MGN" => "R/S",   b"MGS" => "R/S",
        b"MGT" => "R/S",   b"MGV" => "R/S",   b"MGW" => "R/S",   b"MGY" => "R/S",   b"MKA" => "I/L/R",
        b"MKG" => "L/M/R", b"MRA" => "K/Q/R", b"MRG" => "K/Q/R", b"MRR" => "K/Q/R", b"MSA" => "P/R/T",
        b"MSG" => "P/R/T", b"MSR" => "P/R/T", b"MTA" => "I/L",   b"MTB" => "I/L/M", b"MTC" => "I/L",
        b"MTD" => "I/L/M", b"MTG" => "L/M",   b"MTH" => "I/L",   b"MTK" => "I/L/M", b"MTM" => "I/L",
        b"MTN" => "I/L/M", b"MTR" => "I/L/M", b"MTS" => "I/L/M", b"MTT" => "I/L",   b"MTV" => "I/L/M",
        b"MTW" => "I/L",   b"MTY" => "I/L",   b"NGA" => "G/R/*", b"NGG" => "G/R/W", b"NTA" => "I/L/V",
        b"NTG" => "L/M/V", b"RAA" => "E/K",   b"RAC" => "D/N",   b"RAG" => "E/K",   b"RAR" => "E/K",
        b"RAT" => "D/N",   b"RAY" => "D/N",   b"RCA" => "A/T",   b"RCB" => "A/T",   b"RCC" => "A/T",
        b"RCD" => "A/T",   b"RCG" => "A/T",   b"RCH" => "A/T",   b"RCK" => "A/T",   b"RCM" => "A/T",
        b"RCN" => "A/T",   b"RCR" => "A/T",   b"RCS" => "A/T",   b"RCT" => "A/T",   b"RCV" => "A/T",
        b"RCW" => "A/T",   b"RCY" => "A/T",   b"RGA" => "G/R",   b"RGB" => "G/R/S", b"RGC" => "G/S",
        b"RGD" => "G/R/S", b"RGG" => "G/R",   b"RGH" => "G/R/S", b"RGK" => "G/R/S", b"RGM" => "G/R/S",
        b"RGN" => "G/R/S", b"RGR" => "G/R",   b"RGS" => "G/R/S", b"RGT" => "G/S",   b"RGV" => "G/R/S",
        b"RGW" => "G/R/S", b"RGY" => "G/S",   b"RTA" => "I/V",   b"RTB" => "I/M/V", b"RTC" => "I/V",
        b"RTD" => "I/M/V", b"RTG" => "M/V",   b"RTH" => "I/V",   b"RTK" => "I/M/V", b"RTM" => "I/V",
        b"RTN" => "I/M/V", b"RTR" => "I/M/V", b"RTS" => "I/M/V", b"RTT" => "I/V",   b"RTV" => "I/M/V",
        b"RTW" => "I/V",   b"RTY" => "I/V",   b"SAA" => "E/Q",   b"SAC" => "D/H",   b"SAG" => "E/Q",
        b"SAR" => "E/Q",   b"SAT" => "D/H",   b"SAY" => "D/H",   b"SCA" => "A/P",   b"SCB" => "A/P",
        b"SCC" => "A/P",   b"SCD" => "A/P",   b"SCG" => "A/P",   b"SCH" => "A/P",   b"SCK" => "A/P",
        b"SCM" => "A/P",   b"SCN" => "A/P",   b"SCR" => "A/P",   b"SCS" => "A/P",   b"SCT" => "A/P",
        b"SCV" => "A/P",   b"SCW" => "A/P",   b"SCY" => "A/P",   b"SGA" => "G/R",   b"SGB" => "G/R",
        b"SGC" => "G/R",   b"SGD" => "G/R",   b"SGG" => "G/R",   b"SGH" => "G/R",   b"SGK" => "G/R",
        b"SGM" => "G/R",   b"SGN" => "G/R",   b"SGR" => "G/R",   b"SGS" => "G/R",   b"SGT" => "G/R",
        b"SGV" => "G/R",   b"SGW" => "G/R",   b"SGY" => "G/R",   b"STA" => "L/V",   b"STB" => "L/V",
        b"STC" => "L/V",   b"STD" => "L/V",   b"STG" => "L/V",   b"STH" => "L/V",   b"STK" => "L/V",
        b"STM" => "L/V",   b"STN" => "L/V",   b"STR" => "L/V",   b"STS" => "L/V",   b"STT" => "L/V",
        b"STV" => "L/V",   b"STW" => "L/V",   b"STY" => "L/V",   b"TAB" => "Y/*",   b"TAD" => "Y/*",
        b"TAH" => "Y/*",   b"TAK" => "Y/*",   b"TAM" => "Y/*",   b"TAN" => "Y/*",   b"TAS" => "Y/*",
        b"TAV" => "Y/*",   b"TAW" => "Y/*",   b"TBA" => "L/S/*", b"TBC" => "C/F/S", b"TBG" => "L/S/W",
        b"TBT" => "C/F/S", b"TBY" => "C/F/S", b"TDA" => "L/*",   b"TDC" => "C/F/Y", b"TDG" => "L/W/*",
        b"TDR" => "L/W/*", b"TDT" => "C/F/Y", b"TDY" => "C/F/Y", b"TGB" => "C/W",   b"TGD" => "C/W/*",
        b"TGH" => "C/*",   b"TGK" => "C/W",   b"TGM" => "C/*",   b"TGN" => "C/W/*", b"TGR" => "W/*",
        b"TGS" => "C/W",   b"TGV" => "C/W/*", b"TGW" => "C/*",   b"THA" => "L/S/*", b"THC" => "F/S/Y",
        b"THG" => "L/S/*", b"THR" => "L/S/*", b"THT" => "F/S/Y", b"THY" => "F/S/Y", b"TKA" => "L/*",
        b"TKC" => "C/F",   b"TKG" => "L/W",   b"TKR" => "L/W/*", b"TKT" => "C/F",   b"TKY" => "C/F",
        b"TMA" => "S/*",   b"TMB" => "S/Y/*", b"TMC" => "S/Y",   b"TMD" => "S/Y/*", b"TMG" => "S/*",
        b"TMH" => "S/Y/*", b"TMK" => "S/Y/*", b"TMM" => "S/Y/*", b"TMN" => "S/Y/*", b"TMR" => "S/*",
        b"TMS" => "S/Y/*", b"TMT" => "S/Y",   b"TMV" => "S/Y/*", b"TMW" => "S/Y/*", b"TMY" => "S/Y",
        b"TNA" => "L/S/*", b"TRC" => "C/Y",   b"TRG" => "W/*",   b"TRH" => "C/Y/*", b"TRM" => "C/Y/*",
        b"TRR" => "W/*",   b"TRT" => "C/Y",   b"TRW" => "C/Y/*", b"TRY" => "C/Y",   b"TSA" => "S/*",
        b"TSB" => "C/S/W", b"TSC" => "C/S",   b"TSG" => "S/W",   b"TSH" => "C/S/*", b"TSK" => "C/S/W",
        b"TSM" => "C/S/*", b"TSR" => "S/W/*", b"TSS" => "C/S/W", b"TST" => "C/S",   b"TSW" => "C/S/*",
        b"TSY" => "C/S",   b"TTB" => "F/L",   b"TTD" => "F/L",   b"TTH" => "F/L",   b"TTK" => "F/L",
        b"TTM" => "F/L",   b"TTN" => "F/L",   b"TTS" => "F/L",   b"TTV" => "F/L",   b"TTW" => "F/L",
        b"TVA" => "S/*",   b"TVC" => "C/S/Y", b"TVG" => "S/W/*", b"TVR" => "S/W/*", b"TVT" => "C/S/Y",
        b"TVY" => "C/S/Y", b"TWA" => "L/*",   b"TWC" => "F/Y",   b"TWG" => "L/*",   b"TWR" => "L/*",
        b"TWT" => "F/Y",   b"TWY" => "F/Y",   b"TYA" => "L/S",   b"TYB" => "F/L/S", b"TYC" => "F/S",
        b"TYD" => "F/L/S", b"TYG" => "L/S",   b"TYH" => "F/L/S", b"TYK" => "F/L/S", b"TYM" => "F/L/S",
        b"TYN" => "F/L/S", b"TYR" => "L/S",   b"TYS" => "F/L/S", b"TYT" => "F/S",   b"TYV" => "F/L/S",
        b"TYW" => "F/L/S", b"TYY" => "F/S",   b"VAA" => "E/K/Q", b"VAC" => "D/H/N", b"VAG" => "E/K/Q",
        b"VAR" => "E/K/Q", b"VAT" => "D/H/N", b"VAY" => "D/H/N", b"VCA" => "A/P/T", b"VCB" => "A/P/T",
        b"VCC" => "A/P/T", b"VCD" => "A/P/T", b"VCG" => "A/P/T", b"VCH" => "A/P/T", b"VCK" => "A/P/T",
        b"VCM" => "A/P/T", b"VCN" => "A/P/T", b"VCR" => "A/P/T", b"VCS" => "A/P/T", b"VCT" => "A/P/T",
        b"VCV" => "A/P/T", b"VCW" => "A/P/T", b"VCY" => "A/P/T", b"VGA" => "G/R",   b"VGB" => "G/R/S",
        b"VGC" => "G/R/S", b"VGD" => "G/R/S", b"VGG" => "G/R",   b"VGH" => "G/R/S", b"VGK" => "G/R/S",
        b"VGM" => "G/R/S", b"VGN" => "G/R/S", b"VGR" => "G/R",   b"VGS" => "G/R/S", b"VGT" => "G/R/S",
        b"VGV" => "G/R/S", b"VGW" => "G/R/S", b"VGY" => "G/R/S", b"VTA" => "I/L/V", b"VTC" => "I/L/V",
        b"VTG" => "L/M/V", b"VTH" => "I/L/V", b"VTM" => "I/L/V", b"VTT" => "I/L/V", b"VTW" => "I/L/V",
        b"VTY" => "I/L/V", b"WAA" => "K/*",   b"WAC" => "N/Y",   b"WAG" => "K/*",   b"WAR" => "K/*",
        b"WAT" => "N/Y",   b"WAY" => "N/Y",   b"WCA" => "S/T",   b"WCB" => "S/T",   b"WCC" => "S/T",
        b"WCD" => "S/T",   b"WCG" => "S/T",   b"WCH" => "S/T",   b"WCK" => "S/T",   b"WCM" => "S/T",
        b"WCN" => "S/T",   b"WCR" => "S/T",   b"WCS" => "S/T",   b"WCT" => "S/T",   b"WCV" => "S/T",
        b"WCW" => "S/T",   b"WCY" => "S/T",   b"WGA" => "R/*",   b"WGC" => "C/S",   b"WGG" => "R/W",
        b"WGR" => "R/W/*", b"WGT" => "C/S",   b"WGY" => "C/S",   b"WRA" => "K/R/*", b"WSC" => "C/S/T",
        b"WST" => "C/S/T", b"WSY" => "C/S/T", b"WTA" => "I/L",   b"WTC" => "F/I",   b"WTG" => "L/M",
        b"WTH" => "F/I/L", b"WTM" => "F/I/L", b"WTR" => "I/L/M", b"WTT" => "F/I",   b"WTW" => "F/I/L",
        b"WTY" => "F/I",   b"YAA" => "Q/*",   b"YAC" => "H/Y",   b"YAG" => "Q/*",   b"YAR" => "Q/*",
        b"YAT" => "H/Y",   b"YAY" => "H/Y",   b"YCA" => "P/S",   b"YCB" => "P/S",   b"YCC" => "P/S",
        b"YCD" => "P/S",   b"YCG" => "P/S",   b"YCH" => "P/S",   b"YCK" => "P/S",   b"YCM" => "P/S",
        b"YCN" => "P/S",   b"YCR" => "P/S",   b"YCS" => "P/S",   b"YCT" => "P/S",   b"YCV" => "P/S",
        b"YCW" => "P/S",   b"YCY" => "P/S",   b"YGA" => "R/*",   b"YGB" => "C/R/W", b"YGC" => "C/R",
        b"YGG" => "R/W",   b"YGH" => "C/R/*", b"YGK" => "C/R/W", b"YGM" => "C/R/*", b"YGR" => "R/W/*",
        b"YGS" => "C/R/W", b"YGT" => "C/R",   b"YGW" => "C/R/*", b"YGY" => "C/R",   b"YKA" => "L/R/*",
        b"YKG" => "L/R/W", b"YRA" => "Q/R/*", b"YTB" => "F/L",   b"YTC" => "F/L",   b"YTD" => "F/L",
        b"YTH" => "F/L",   b"YTK" => "F/L",   b"YTM" => "F/L",   b"YTN" => "F/L",   b"YTS" => "F/L",
        b"YTT" => "F/L",   b"YTV" => "F/L",   b"YTW" => "F/L",   b"YTY" => "F/L",   b"YWA" => "L/Q/*",
        b"YWG" => "L/Q/*", b"YWR" => "L/Q/*", b"YYA" => "L/P/S", b"YYG" => "L/P/S", b"YYR" => "L/P/S"
    )
});
