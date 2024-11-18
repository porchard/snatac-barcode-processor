use std::fs::File;
use std::collections::{HashSet,HashMap};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{Read,BufReader,BufWriter};
use bio::io::fastq;
use itertools::izip;
use log::info;
use crate::trie::Trie;

fn likelihood_of_errors (uncorrected: &[u8], corrected: &[u8], phred: &[u8]) -> f64 {

    assert_eq!(uncorrected.len(), corrected.len());
    
    let mut l: f64 = 1.0;
    
    for (u, c, p) in izip!(uncorrected, corrected, phred) {
        if u != c {
            let power_base: f64 = 10.0;
            let q: f64 = (*p as f64) - 33.0;
            let byte_l: f64 = power_base.powf(-1.0 * q / 10.0);
            l *= byte_l;
        }
    }

    l
}

/// Correct a non-whitelisted barcode.
/// 
/// Given the uncorrected barcode, it's phred score, a vector of similar whitelisted barcodes (e.g., 
/// whitelisted barcodes w/in Hamming distance two of the uncorrected barcode), and a vector of counts 
/// for the similar whitelisted barcodes (representing how often each of those similar barcodes are 
/// observed in the library; these act as a sort of "prior"), attempts to correct the uncorrected barcode
/// to one of the similar whitelisted barcodes.
fn correct_barcode<'a> (uncorrected: &[u8], uncorrected_phred: &[u8], similar: &Vec<&'a [u8]>, similar_counts: &Vec<&usize>) -> Option<&'a [u8]> {

    if similar.is_empty() {
        return None;
    } else if similar.len() == 1 {
        return Some(similar[0]);
    } else {
        let likelihood_of_errors: Vec<f64> = similar.iter().map(|&s| likelihood_of_errors(uncorrected, s, uncorrected_phred)).collect();
        let likelihood: Vec<f64> = likelihood_of_errors.iter().zip(similar_counts.iter()).map(|(&i, &&j)| i*(j as f64)).collect();
        let norm_factor: f64 = likelihood.iter().sum();
        let norm_likelihoods: Vec<f64> = likelihood.iter().map(|i| i / norm_factor).collect();

        for (i, &correction) in similar.iter().enumerate() {
            if norm_likelihoods[i] >= 0.975 {
                return Some(correction);
            }
        }
    }

    None

}


pub fn correct_barcodes_in_fastq (input_fastq_filename: &str, whitelist_filename: &str, counts_filename: &str, output_fastq_filename: &str, max_edit_distance: usize) {

    // read the whitelist
    let mut whitelist_file = File::open(whitelist_filename).unwrap();
    let mut whitelist: String = String::new();
    whitelist_file.read_to_string(&mut whitelist).unwrap();
    let whitelist: HashSet<&[u8]> = whitelist.split("\n").map(|s| s.trim_end().as_bytes()).collect();

    let mut whitelist_trie = Trie::new();
    for &whitelisted_barcode in whitelist.iter() {
        whitelist_trie.add_word(whitelisted_barcode);
    }

    // read the counts
    let mut counts: HashMap<&[u8], usize> = HashMap::new();
    let mut counts_file = File::open(counts_filename).unwrap();
    let mut counts_string = String::new();
    counts_file.read_to_string(&mut counts_string).unwrap();
    counts_string = counts_string.trim().to_string();
    for i in counts_string.split("\n") {
        let barcode_and_count: Vec<&str> = i.split("\t").collect();
        let barcode = barcode_and_count[0].as_bytes();
        let count = barcode_and_count[1].parse::<usize>().unwrap();
        let e = counts.entry(barcode).or_insert(0);
        *e += count;
    }
    // add pseudocount
    for &whitelisted_barcode in whitelist.iter() {
        if counts.contains_key(whitelisted_barcode) {
            *(counts.get_mut(&whitelisted_barcode).unwrap()) += 1;
        } else {
            counts.insert(whitelisted_barcode, 1);
        }
    }

    let fastq_in = BufReader::new(GzDecoder::new(File::open(input_fastq_filename).unwrap()));
    let fastq_reader = fastq::Reader::from_bufread(fastq_in);

    let fastq_out = BufWriter::new(GzEncoder::new(File::create(output_fastq_filename).unwrap(), Compression::default()));
    let mut fastq_writer = fastq::Writer::from_bufwriter(fastq_out);

    let mut matched_whitelist_before_correction: usize = 0;
    let mut matched_whitelist_after_correction: usize = 0;
    let mut total: usize = 0;

    for result in fastq_reader.records() {
        total += 1;

        let record = result.unwrap();

        if whitelist.contains(&record.seq()) {
            matched_whitelist_before_correction += 1;
            matched_whitelist_after_correction += 1;
            let new_description = format!("CR:Z:{}\tCB:Z:{}\tCY:Z:{}", String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(record.qual().to_vec()).unwrap());

            fastq_writer.write(record.id(), Some(&new_description), record.seq(), record.qual()).unwrap();
        } else {
            let corrections = whitelist_trie.get_words_within_hamming_distance(record.seq(), max_edit_distance);
            let corrections: Vec<&[u8]> = corrections.iter().map(|(s, _c)| s.as_bytes()).collect();
            let corrections_counts: Vec<&usize> = corrections.iter().map(|&s| counts.get(s).unwrap_or(&0)).collect();
            let corrected = correct_barcode(record.seq(), record.qual(), &corrections, &corrections_counts);

            let new_description = match corrected {
                Some(x) => {
                    matched_whitelist_after_correction += 1;
                    format!("CR:Z:{}\tCB:Z:{}\tCY:Z:{}", String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(x.to_vec()).unwrap(), String::from_utf8(record.qual().to_vec()).unwrap())
                },
                None => {
                    format!("CR:Z:{}\tCY:Z:{}", String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(record.qual().to_vec()).unwrap())
                },
            };
            
            fastq_writer.write(record.id(), Some(&new_description), record.seq(), record.qual()).unwrap();
        }
        
        if total % 1000000 == 0 {
            info!("Processed {total} records so far; {matched_whitelist_before_correction} matched whitelist before correction, {matched_whitelist_after_correction} matched whitelist after correction");
        }
    }

    fastq_writer.flush().unwrap();

}