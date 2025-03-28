
use clap::Parser;

use gc_counter::BedData;
use needletail::parse_fastx_file;

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;



#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the fasta file (should be gzipped!)
    #[clap(short, long)]
    file: String,
    /// the output file
    #[clap(short, long)]
    outfile: String,
    /// how many threads to use to analyze this (default 1)
    #[clap( short, long)]
    num_threads: Option<usize>,
    /// genome sizes - a tab separated text file with genome sizes
    #[clap( short, long)]
    genome_sizes: String,
    /// genome bin length
    #[clap( short, long)]
    bin_width: usize,
    
}


fn main() {

    let opts: Opts = Opts::parse();


    let mut expr_file = parse_fastx_file(&opts.file).expect("valid path/file");

    let genome_sizes = read_tsv_to_vec( &opts.genome_sizes ).unwrap();
    let mut data = BedData::new( &genome_sizes, opts.bin_width );
    let delimiter = b"|";

    while let Some(e_record) = expr_file.next() {
        let seqrec = e_record.expect("invalid record");
        let id = std::str::from_utf8( 
                seqrec.id().split(|&x| x == delimiter[0]).collect::<Vec<&[u8]>>()[0]
            ).unwrap().to_string();
        let seq = seqrec.seq();
        // Scan the sequence for "GC" pairs
        for j in 0..seq.len().saturating_sub(1) {
            if seq[j] == b'G' && seq[j + 1] == b'C' {
                data.add( &id, j, 1.0 );
                //println!("Found GC at position {} in sequence {}", j, id);
            }
        }
    }

    data.write_bigwig( &opts.outfile );


}

fn read_tsv_to_vec<P: AsRef<Path>>(filename: P) -> Result<Vec<(String, usize)>, String> {
    let file = File::open(&filename).map_err(|e| format!("Failed to open file {:?}: {}",filename.as_ref(), e))?;
    let reader = io::BufReader::new(file);
    
    let mut result = Vec::new();
    
    for (index, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading line {}: {}", index + 1, e))?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            return Err(format!("Invalid format on line {}: expected 2 columns, found {}", index + 1, parts.len()));
        }
        let num = parts[1].parse::<usize>().map_err(|e| format!("Failed to parse usize on line {}: {}", index + 1, e))?;
        result.push((parts[0].to_string(), num));
    }
    Ok(result)
}