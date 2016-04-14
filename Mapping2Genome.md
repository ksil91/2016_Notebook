# Mapping Subset of 2bRAD Samples to Genome

There were difficulties with the *de novo* assembly of the 2bRAD reads using the Meyer lab scripts, but now we have a draft *O. lurida* genome and so don't need to *de novo* assemble the 2bRAD loci. I am mostly following the [Matz lab's protocol for genotyping with GATK](https://github.com/z0on/2bRAD_GATK) as it includes steps for using recalibrating variant quality scores using technical replicates. We have 2 technical replicates within each library and 3 individuals that were sequenced across all 4 libraries. I had to adapt their scripts in a few places as they were written around job submission to their university's cluster. Most analyses were run on my department's single node, 24 processor high performance computer. For brevity I don't include the .sh scripts to submit the job to the HPC, just the commands that were run. All of the input and output files can be found at: 
http://owl.fish.washington.edu/wetgenes/index.php?dir=subset_genotype%2F


First I selected 20 individuals from each population that had at least 1 million sequencing reads, as well as a few technical repeats for each population.
|South Sound | Hood Canal | Fidalgo Bay |
|------------| --------| --------------|
|SS3-14      | HC2-17A | NF2-15 |
|SS3-15     | HC2-17B | NF5-19 |
|SS3-16   |HC3-1 | NF5-14 |
|SS3-3 | HC3-10 | NF2-12 |
|SS3-20 | HC3-11 | NF3-17 |
|SS2-14 | HC3-7 | NF1-8 |
|SS5-18 | HC1-4 | NF4-12 |
|SS5-7 | HC4-3 | NF2-6A |
|SS1-8 | HC3-4 | NF2-6B |
|SS2-15 | HC4-6 | NF2-6C |
|SS4-8 | HC5-12 | NF1-14A |
|SS3-2 | HC2-1 | NF1-14B |
|SS3-6 | HC4-9 | NF5-3 |
|SS2-12A | HC4-14 | NF2-9 |
|SS2-12B | HC4-4A | NF3-11 |
|SS2-4A | HC4-4B | NF1-19 |
|SS2-4B | HC1-17 | NF4-13 |
|SS3-21 | HC2-8 | NF3-9 |
|SS5-17 | HC5-10 | NF4-9 |
|SS4-13 | HC3-17 | NF5-11 |
|SS3-7 | HC2-14 | NF4-11 |
|SS2-13 | HC2-18 | NF3-5 |
|   |   | NF3-16 |


 ### Required scripts and programs  
 For these steps you need:  
 -  [concatFasta.pl](https://github.com/z0on/2bRAD_GATK/blob/master/concatFasta.pl) from Mikhail Matz's github
 -  bowtie2
 -  Picard Tools
 -  samtools
 ### Preparing reference  
For the reference I'm using scaffolds from the genome assembly by BGI that are over 10KB.
http://de.iplantcollaborative.org/dl/d/5E084D53-E706-420E-AC7D-8620F6F0A535/OlyBGI-scaffold-10k.fa

Concatenating scaffolds into small number of "pseudo-chromosomes":
```sh
concatFasta.pl fasta=OlyBGI-scaffold-10k.fasta
```
Normalize fasta records to ensure all lines are of the same length, using Picard
```sh
java -Xmx1g -jar $PICARD NormalizeFasta INPUT=OlyBGI-scaffold-10k_cc.fasta OUTPUT=OlyBGI-scaffold-10k_ccn.fasta
```
Create genome index
```sh
export GENOME_FASTA=OlyBGI-scaffold-10k_ccn.fasta 
export GENOME_DICT=OlyBGI-scaffold-10k_ccn.dict  
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA
java -jar $PICARD CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT
```
I then wrote a script to automate the creation of job submission shell scripts for mapping reads from each sample to the reference.
```sh
export REF=/home/ksilliman/CommonG/OlyBGI-scaffold-10k_ccn.fasta 

for ffile in *.fastq; do 
    echo "#!/bin/bash"  >> ${ffile/ca.fastq/.sh}
    echo "#PBS -N ${ffile/ca.fastq/m}" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -S /bin/bash" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -l nodes=1:ppn=1"  >> ${ffile/ca.fastq/.sh}
    echo "#PBS -o ${ffile/ca.fastq/m.out}"  >>  ${ffile/ca.fastq/.sh}
    echo "#PBS -e ${ffile/ca.fastq/m.err}" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -l walltime=120:00:00" >> ${ffile/ca.fastq/.sh}
    echo "#PBS -V" >> ${ffile/ca.fastq/.sh}
    echo "cd ~/CommonG/2bRAD_Dec2015/best30_ca" >> ${ffile/ca.fastq/.sh}
    echo "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $REF -U $ffile -S ${ffile/ca.fastq/.bt2.sam}" >> ${ffile/ca.fastq/.sh}
done
```
The line to run bowtie2 is:
```sh
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $REF -U $ffile -S ${ffile/ca.fastq/.bt2.sam}
```
The *.err files list the output of the mapping. It seems like most samples only had ~30% of the reads map to the genome. I suspect that is because we are missing some of the genome with these scaffolds.

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