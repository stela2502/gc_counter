# gc_counter

This is an example implementation of a simple bioinformatics workflow:

Count the GC contents of a fasta file and write the results into a bigwig file counting the CG content of the sequence in valiable bin sizes. The GC content is measured in percent.

This is a tiny project that should introduce the usage of Rust in Bioinformatics to experienced R and Python programmers.

And this repo is the example implementation, not the tutorial on the basics of Rust progamming...

Took me ~ 2h to implement. And a little more to smooth the edges.

The comparison to R took me actually longer and I am sure that one is not efficient.

The resuts are comparable, just that my R does not compact concurrent same values into a bigger area and it actually lacks the true end of the sequence.

# The timings:


## Rust

```
time target/release/gc_counter -f tests/data/chrM.fa.gz -o tests/data/out.bw -b 50 -g tests/data/chrm.size.tsv 
finished in 0 h 0 min 0 sec 2 milli sec

real    0m0,005s
user    0m0,000s
sys 0m0,005s
```


## R

After I had installed all libraries....

```
time Rscript tests/data/Same_in_R.R -i tests/data/chrM.fa.gz -o tests/data/out.R.bw -b 50

BigWig file created at: tests/data/out.R.bw 

real    0m8,869s
user    0m8,259s
sys 0m0,612s
```


But that was unfair from the start...

# Implementation

## lib.rs

The primary goal of this project is to store numeric data corresponding to positions in a genomic sequence. The length of the string (genome) is known, and we want to store values using a windowed approach. This means we need a storage structure for each of these values.

A typical genome contains around 3 gigabases (Gb) of information. With a window size of 50 base pairs (bp), this results in a vector of **60 million** `f32` entries, which occupies approximately **229 MB** of memory. This is not an excessive amount of memory for the task at hand, and this approach is far more efficient than any alternative methods. As such, we need to implement a class that initializes the array based on the genome sizes and provides methods to convert from chromosome and position to index, and vice versa.

The sequences we are working with are typically very long strings. Therefore, it is feasible to store genome sizes in a `HashMap`.

The implementation I followed was inspired by a [previous project of mine](https://github.com/stela2502/bam_tide).

To generate a BigWig file from this class, I am utilizing the Rust `bigtools` package. It is straightforward to export the data using this library once you understand how it is meant to be used. The steps are as follows:

1. Create a class that implements the `Iterator` trait. This involves implementing the `next()` function, which returns the next value to be added to the file and `None` when finished.
2. Set up a Tokio runtime (though I'm not entirely sure what it does—just follow the steps 😉).
3. Use the `BedParserStreamingIterator::wrap_infallible_iter` function from the `bigdata` library.
4. Export the data to the BigWig file.

I honestly wouldn’t have figured this out on my own—[I asked the developer of the `bigtools` library, and they responded the next day](https://github.com/jackh726/bigtools/discussions/74). I hope you would also follow this approach!

## main.rs

The remaining task is handling user input and performing the calculations for the program. For user interaction, I used the Rust `clap` package. Any AI chatbot could suggest a similar solution. 

The FASTA file interaction is handled by the `needletail` library. I also borrowed this approach from one of my old Rust projects 😉. The calculations themselves were written by ChatGPT. It's not particularly complicated—I'm just lazy.

All other data handling has already been implemented in the library.
