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
usage: nanotime_v2.py [-h] -q FQZ -s SES -l LEN -o OUT

Nanotime v2. Nanopore reads extraction based on the start time of sequencing; The process on fast5(signal) file is still
under development

optional arguments:
  -h, --help         show this help message and exit
  -q FQZ, --fqz FQZ  the path of compressed fastq file (e.g. *.fastq.gz)
  -s SES, --ses SES  the path to sequencing_summary.txt
  -l LEN, --len LEN  the length of run of sequencing (e.g. 1)
  -o OUT, --out OUT  the path of output

```

## Example
```sh
git clone https://github.com/renzilin/NanoTime

python nanotime_v2.py -q <input>/barcode01.fq.gz -s <where-is>/sequencing_summary.txt -l <length> -o <output>/barcode01.fq.gz

```



