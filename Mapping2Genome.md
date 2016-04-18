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
 ### Mapping Reads 
I then wrote a script to automate the creation of job submission shell scripts for mapping reads from each sample to the reference. This creates .bt2.sam files for each sample.
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

Making .bam files from .sam files:
```sh
for ffile in *.bt2.sam; do 
    echo "#!/bin/bash"  >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -N ${ffile/bt2.sam/s2b}" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -S /bin/bash" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -l mem=5gb,nodes=1:ppn=1"  >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -o ${ffile/bt2.sam/s2b.out}"  >>  ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -e ${ffile/bt2.sam/s2b.err}" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -l walltime=50:00:00" >> ${ffile/bt2.sam/s2b.sh}
    echo "#PBS -V" >> ${ffile/bt2.sam/s2b.sh}
    echo "cd ~/CommonG/2bRAD_Dec2015/best30_ca" >> ${ffile/bt2.sam/s2b.sh}
    echo "samtools view -bS $ffile > ${ffile/.bt2.sam/_unsorted.bam}" >> ${ffile/bt2.sam/s2b.sh}
    echo "samtools sort ${ffile/.bt2.sam/_unsorted.bam} -f ${ffile/.bt2.sam/_sorted.bam}" >> ${ffile/bt2.sam/s2b.sh}
    echo "java -Xmx5g -jar $PICARD AddOrReplaceReadGroups INPUT=${ffile/.bt2.sam/_sorted.bam} OUTPUT=${ffile/bt2.sam/bam} RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${ffile/bt2.sam/m}" >> ${ffile/bt2.sam/s2b.sh} 
    echo "samtools index ${ffile/bt2.sam/bam}" >> ${ffile/bt2.sam/s2b.sh}

done
```

Tried VQSR:
cat round2.vcf | perl -pe 's/\.m//g' | perl -pe 's/^chrom/chr/' >round2.names.vcf

~/CommonG/z0on/replicatesMatch.pl vcf=round2.names.vcf replicates=testreps.tab > vqsr.vcf

7995 total SNPs
1443 pass hets and match filters
948 show non-reference alleles
738 have alterantive alleles in at least 2 replicate pair(s)
738 have matching heterozygotes in at least 0 replicate pair(s)
317 polymorphic
738 written

VCFtools - 0.1.15
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf vqsr.vcf
	--TsTv-summary

After filtering, kept 67 out of 67 Individuals
Outputting Ts/Tv summary
Ts/Tv ratio: 1.463
After filtering, kept 738 out of a possible 738 Sites
Run Time = 0.00 seconds

vqsr matz:
~/CommonG/z0on/replicatesMatch.pl vcf=round2.names.vcf replicates=testreps.tab polyonly=1 > vqsr_matz.vcf
7995 total SNPs
1443 pass hets and match filters
948 show non-reference alleles
738 have alterantive alleles in at least 2 replicate pair(s)
738 have matching heterozygotes in at least 0 replicate pair(s)
317 polymorphic
317 written

MQ quantiles:
1	8.52
5	14.45
10	17.00
15	19.71
20	20.75
30	24.94
40	28.27
50	31.23
60	35.50
70	38.27
80	40.15
90	41.76

FS quantiles:
1	16.189
5	6.129
10	3.762
15	2.799
20	1.977
30	1.185
40	0.736
50	0.534
60	0.000
70	0.000
80	0.000
90	0.000

DP quantiles:
1	218	16527
5	340	16239
10	551	15544
15	724	15175
20	946	14175
30	1221	12608
40	1513	11155
50	1848	9800
60	2276	8662
70	2916	7140
80	3540	6037
90	4016	5301

QD quantiles:
1	22.87	
5	21.73	
10	19.86	
15	18.62	
20	17.27	
30	13.83	
40	11.28	
50	9.31	
60	7.26	
70	5.30	
80	3.75	
90	2.07	

JOINT quantiles:
1	0.00030
5	0.00225
10	0.00500
15	0.00675
20	0.00840
30	0.01500
40	0.02520
50	0.03840
60	0.06750
70	0.10800
80	0.16000
90	0.33600

------------------------
20.61%	at qual <1 (19.61% gain)
58.16%	at qual <5 (53.16% gain)
67.35%	at qual <10 (57.35% gain)
69.86%	at qual <15 (54.86% gain)
71.32%	at qual <20 (51.32% gain)
75.75%	at qual <30 (45.75% gain)

-nofs
43.04%	at qual <1 (42.04% gain)
66.62%	at qual <5 (61.62% gain)
71.14%	at qual <10 (61.14% gain)
74.11%	at qual <15 (59.11% gain)
75.85%	at qual <20 (55.85% gain)
79.71%	at qual <30 (49.71% gain)

-nomq
------------------------
7.93%	at qual <1 (6.93% gain)
43.30%	at qual <5 (38.30% gain)
57.37%	at qual <10 (47.37% gain)
58.21%	at qual <15 (43.21% gain)
60.53%	at qual <20 (40.53% gain)
69.22%	at qual <30 (39.22% gain)
------------------------
-nodp
------------------------
5.70%	at qual <1 (4.70% gain)
21.50%	at qual <5 (16.50% gain)
31.09%	at qual <10 (21.09% gain)
38.77%	at qual <15 (23.77% gain)
43.49%	at qual <20 (23.49% gain)
56.27%	at qual <30 (26.27% gain)
------------------------
-noqd
------------------------
18.15%	at qual <1 (17.15% gain)
45.03%	at qual <5 (40.03% gain)
58.30%	at qual <10 (48.30% gain)
59.10%	at qual <15 (44.10% gain)
63.61%	at qual <20 (43.61% gain)
68.49%	at qual <30 (38.49% gain)
------------------------


retabvcf.pl vcf=gatk_after_vqsr_matz.vcf tab=/home/ksilliman/CommonG/OlyBGI-scaffold-10k_cc.tab > retab.vcf

# discarding loci with too many heterozygotes, which are likely lumped paralogs
# (by default, fraction of heterozygotes should not exceed maxhet=0.75)
# this step can also filter for the fraction of missing genotypes (default maxmiss=0.5)
hetfilter.pl vcf=retab.vcf > hetfilt_def.vcf
7995 total loci
4351 dropped because fraction of missing genotypes exceeded 0.5
134 dropped because fraction of heterozygotes exceeded 0.75
3510 written

thinning to 1 SNP per tag
thinner.pl vcf=hetfilt.vcf > thin_def.vcf
3510 total loci
257 loci skipped because they were closer than 40
1765 loci selected

# applying filter and selecting polymorphic biallelic loci genotyped in 80% or more individuals
# for non parametric (GATK-based) recalibration: replace --minQ 15 in the following line
# with the quantile of the highest "gain" as reported by recalibrateSNPs_gatk.pl
vcftools --vcf thin_def.vcf --minQ 5 --max-missing 0.8  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filt0
After filtering, kept 67 out of 67 Individuals
After filtering, kept 755 out of a possible 1765 Sites

# genotypic match between pairs of replicates 
# (the most important one is the last one, Heterozygote Discovery Rate)	
repMatchStats.pl vcf=filt0.recode.vcf replicates=testreps.tab 
pair	gtyped	match	[ 00	01	11 ]	HetMatch	HomoHetMismatch	HetNoCall	HetsDiscoveryRate
HC2-17A:HC2-17B	725	660(91.0%)	 [57%	18%	25% ]	119	53	2	0.81	
HC4-4A:HC4-4B	722	672(93.1%)	 [54%	21%	25% ]	143	38	4	0.87	
SS2-12A:SS2-12B	749	533(71.2%)	 [60%	18%	23% ]	95	184	3	0.50	
SS2-4A:SS2-4B	744	687(92.3%)	 [52%	25%	23% ]	174	55		0.86	
NF2-6A:NF2-6B	728	664(91.2%)	 [56%	18%	25% ]	121	43	4	0.84	
NF2-6A:NF2-6C	728	569(78.2%)	 [57%	17%	25% ]	98	56	11	0.75	
NF2-6C:NF2-6B	631	569(90.2%)	 [57%	17%	26% ]	97	53	13	0.75	
NF1-14A:NF1-14B	724	669(92.4%)	 [57%	19%	24% ]	126	44	1	0.85	

------------------------
hets called homos depth: 
lower 25%	11
median		41
upper 75%	116
# creating final filtered file without clones (must list them in the file clones2remove):
# (for non-parametric recalibration, replace --remove-filtered-all with --minQ [quantile of the highest gain] )
vcftools --vcf thin_def.vcf --remove testreps_del.tab --minQ 5 --max-missing 0.7  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out testfinal.7
After filtering, kept 59 out of 67 Individuals
Outputting VCF file...
After filtering, kept 882 out of a possible 1765 Sites

for Tajima's D
vcftools --vcf hetfilt_def.vcf --remove testreps_del.tab --minQ 5 --max-missing 0.7  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out testfinalUnthin.7
After filtering, kept 1855 out of a possible 3510 Sites












