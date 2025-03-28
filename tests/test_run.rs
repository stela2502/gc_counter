use std::process::Command;
use std::fs;
use std::path::Path;

#[test]
fn test_gc_counter_command() {
    // Define the paths for input and output files
    let input_fasta = "tests/data/chrM.fa.gz";
    let output_bw = "tests/data/out.bw";
    let genome_sizes = "tests/data/chrm.size.tsv";
    
    // Run the command using std::process::Command
    let status = Command::new("./target/debug/gc_counter")
        .arg("-f")
        .arg(input_fasta)
        .arg("-o")
        .arg(output_bw)
        .arg("--genome-sizes")
        .arg(genome_sizes)
        .arg("-b")
        .arg("50")
        .status()
        .expect("Failed to execute gc_counter");

    // Check if the command was successful (status code 0)
    assert!(status.success(), "Command failed to execute successfully");

    // Check if the output file exists
    let output_path = Path::new(output_bw);
    assert!(output_path.exists(), "Output file was not created: {}", output_bw);

    // Optional: Clean up the output file after the test (if necessary)
    fs::remove_file(output_bw).expect("Failed to clean up output file");
}
