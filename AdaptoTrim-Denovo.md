# Trimming adapters from 2bRAD reads and *de novo* reference assembly

Notes on the previous steps for trimming to 36bps and quality trimming can be found [here](https://github.com/ksil91/2016_Notebook/blob/master/TrimQual2brad20160217.md). Most analyses were run on my department's single node, 24 processor high performance computer. For brevity I don't include the .sh scripts to submit the job to the HPC, just the commands that were run.  
Raw sequencing files are found at http://owl.fish.washington.edu/nightingales/O_lurida/2bRAD_Dec2015/  
 ### Required scripts and programs  
 For these steps you need:  
 -  [sequence_processing scripts](https://github.com/Eli-Meyer/sequence_processing) from Eli Meyer's github
 -  [2bRAD_utilities version 2](https://github.com/Eli-Meyer/2brad_utilities/tree/v2.0) from Eli Meyer's github
 -  [raxmlHPC]( http://sco.h-its.org/exelixis/web/software/raxml/cluster.html)
 ### Trimming off adapter sequences and filtering non-informative reads  

The Meyer's pipeline includes a script for fitering out reads that have adapter sequences in them, but it involves a program cross.match that is not simple to download and feels arbitrary in how it filters for adapters as it involves a scoring system and threshold. I decided to use the program [cutadapt]( https://cutadapt.readthedocs.org/en/stable/) as it is commonly used for adapter trimming/filtering, has lots of documentation, and is entirely open-source.  
Previously trimmed and quality filtered reads (.q0) are found [here](owl.fish.washington.edu/wetgenes) in the **2brad_201512_qual** folder.
```sh
for file in *_q0.fastq;
    do cutadapt -a AGATCG -m 34 --too-short-output ${file/_q0.fastq/toosh.fastq} -o ${file/_q0.fastq/ca.fastq} $file;
done
```
This searches for the start of the 3' adapter sequence AGATCG in an error-tolerant way and then trims the adapter sequence and everything following it. It retains reads that are at least 34bp long in a ```ca.fastq``` file 
and puts the rest in a file with the ```toosh.fastq``` suffix.  
The summary from cutadapt for each file can be seen [here](http://owl.fish.washington.edu/wetgenes/2brad_201512_ca/ca_2brad_201512.out). Need to write something to grep the number of reads passing filter for each file and put in a .csv format. The output fastq files are found in the **2brad_201512_ca** folder in the wetgenes/ directory.  
```
$ head -n 
This is cutadapt 1.9.1 with Python 2.7.6
Command line parameters: -a AGATCG -m 34 --too-short-output HC1-10toosh.fastq -o HC1-10ca.fastq HC1-10_q0.fastq
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 17.18 s (14 us/read; 4.28 M reads/minute).
=== Summary ===
Total reads processed:               1,225,170
Reads with adapters:                     5,583 (0.5%)
Reads that were too short:               5,583 (0.5%)
Reads written (passing filters):     1,219,587 (99.5%)
Total basepairs processed:    44,106,120 bp
Total written (filtered):     43,905,132 bp (99.5%)
=== Adapter 1 ===
Sequence: AGATCG; Type: regular 3'; Length: 6; Trimmed: 5583 times.
No. of allowed errors:
0-6 bp: 0
Bases preceding removed adapters:
  A: 23.0%
  C: 41.6%
  G: 17.6%
  T: 8.2%
  none/other: 9.7%
Overview of removed sequences
length	count	expect	max.err	error counts
3	5	19143.3	0	5
4	8	4785.8	0	8
5	1015	1196.5	0	1015
6	561	299.1	0	561
8	85	299.1	0	85
9	27	299.1	0	27
10	111	299.1	0	111
11	95	299.1	0	95
12	244	299.1	0	244
21	148	299.1	0	148
22	341	299.1	0	341
29	1714	299.1	0	1714
30	96	299.1	0	96
31	46	299.1	0	46
32	38	299.1	0	38
33	208	299.1	0	208
34	34	299.1	0	34
35	268	299.1	0	268
36	539	299.1	0	539
```
### Preparing reference  
While a genome is in the works for *O. lurida*, it's not ready yet so I'm following the Meyer's pipeline to assemble a reference of AlfI sites. First I made a combined fastq file from 28 samples across all three populations that had at least 1 million reads passing all filtering steps.  
```sh
head -n 2000000 HC1-11ca.fastq >> combined1215.fastq
head -n 4000000 HC1-15ca.fastq >> combined1215.fastq
head -n 2000000 HC2-12ca.fastq >> combined1215.fastq
head -n 4000000 HC2-16ca.fastq >> combined1215.fastq
head -n 2000000 HC3-10ca.fastq >> combined1215.fastq
head -n 4000000 HC3-11ca.fastq >> combined1215.fastq
head -n 2000000 HC4-1Aca.fastq >> combined1215.fastq
head -n 4000000 HC4-2ca.fastq >> combined1215.fastq
head -n 2000000 HC5-12ca.fastq >> combined1215.fastq
head -n 4000000 HC5-14ca.fastq >> combined1215.fastq
head -n 2000000 NF1-12ca.fastq >> combined1215.fastq
head -n 4000000 NF1-13ca.fastq >> combined1215.fastq
head -n 2000000 NF2-15ca.fastq >> combined1215.fastq
head -n 4000000 NF2-2ca.fastq >> combined1215.fastq
head -n 2000000 NF3-3ca.fastq >> combined1215.fastq
head -n 4000000 NF3-5ca.fastq >> combined1215.fastq
head -n 2000000 NF4-19ca.fastq >> combined1215.fastq
head -n 4000000 NF4-2ca.fastq >> combined1215.fastq
head -n 2000000 NF5-12ca.fastq >> combined1215.fastq
head -n 4000000 NF5-4ca.fastq >> combined1215.fastq
head -n 2000000 SS1-8ca.fastq >> combined1215.fastq
head -n 4000000 SS1-12ca.fastq >> combined1215.fastq
head -n 2000000 SS2-12Aca.fastq >> combined1215.fastq
head -n 4000000 SS2-16ca.fastq >> combined1215.fastq
head -n 2000000 SS3-10ca.fastq >> combined1215.fastq
head -n 4000000 SS3-11ca.fastq >> combined1215.fastq
head -n 4000000 SS4-8ca.fastq >> combined1215.fastq
head -n 4000000 SS5-13Bca.fastq >> combined1215.fastq
```
The had to convert the `combined1215.fastq` file to `.fasta`.
```sh
FastqToFasta.pl combined1215.fastq combined1215.fasta combined1215.qua
```
The the call to build the reference `ref1215.fasta`.
```sh
BuildRef.pl combined1215.fasta combined1215.qual ref1215.fasta
```


