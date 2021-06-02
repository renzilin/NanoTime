# NanoTime
A mini tool for nanopore reads (fastq files) extraction based on the start time of sequencing

**The tool is still under development on signals (fast5) files**

## Package
```text
python packages:
    - tqdm
    - multiprocessing
    - Biopython
    - gzip
    - dateutil
```

## Usage
```text
usage: nanotime.py [-h] --fq string --pool_dir string [--duration INT]
                   [--ncpus INT]

Nanotime v0.1. Nanopore reads extraction based on the start time of
sequencing; The process on fast5(signal) file is still under development

optional arguments:
  -h, --help         show this help message and exit
  --fq string        the path of compressed fastq file (e.g. *.fq.gz or
                     *.fastq.gz)
  --pool_dir string  The directory of pool where fastq files were stored (e.g.
                     ./test/20210601-pool1)
  --duration INT     The time duration (hour) of sequencing (e.g. 1)
  --ncpus INT        The number of cores used in scan process for parallel
                     acceleration (e.g. 3)
```

## Example
```sh
git clone https://github.com/renzilin/NanoTime

cd NanoTime/test

# The output of barcode01.fastq should contain 2 reads.
python ../nanotime.py --fq 20210601-pool1/barcode01.fq.gz --pool_dir 20210601-pool1 --duration 1 --ncpus 2

# The output of barcode02 should be empty.
python ../nanotime.py --fq 20210601-pool1/barcode02.fq.gz --pool_dir 20210601-pool1 --duration 1 --ncpus 2

```



