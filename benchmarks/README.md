# Datasets

We used the datasets by Pibiri (2022). The genomes were preprocessed using [BCALM](https://github.com/GATB/bcalm) by Chikhi et al. and [UST](https://github.com/medvedevgroup/UST) by Rahman and Medvedev to remove duplicate k-mers.
The procedure can be reproduced using the script
```
./download_and_preprocess_datasets.sh
```
adapt the path to bcalm and ust on your system accordingly.
The genomes preprecessed by [GGCAT](https://github.com/algbio/ggcat) used by Pibiri and Patro (2026) are available at https://zenodo.org/records/17582116 

# Benchmarks

Run the benchmarks with script
```
./benchmark_rshash.sh
```
It measures the building time and memory consumption using time, query times are measured internally with std::chrono.
The correctness can be checked with test/hashtable that puts all k-mers of a input file in a std::unordered_set<uint64_t> and counts all k-mers in a input query file present in the input sequences.
