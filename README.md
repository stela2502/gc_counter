# gc_counter

This is an example implementation of a simple bioinformatics workflow:

Count the GC contents of a fasta file and write the results into a bigwig file counting the CG content of the sequence in valiable bin sizes.

This is a tiny project that should introduce the usage of Rust in Bioinformatics to experienced R and Python programmers.

And this repo is the example implementation, not the tutorial on the basics of Rust progamming...

Took me ~ 2h to implement. And a little more to break the edges.

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
...


But that was unfair from the start...


