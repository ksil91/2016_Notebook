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
 -  [2bRAD_GATK](https://github.com/z0on/2bRAD_GATK/) from Mikhail Matz's github
 -  bowtie2
 -  Picard Tools
 -  samtools
 -  vcftools
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
Key lines:
```sh
samtools view -bS NF3-16.bt2.sam > NF3-16_unsorted.bam
samtools sort NF3-16_unsorted.bam -f NF3-16_sorted.bam
java -Xmx5g -jar /usr/local/picard/dist/picard.jar AddOrReplaceReadGroups INPUT=NF3-16_sorted.bam OUTPUT=NF3-16.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=NF3-16.m
samtools index NF3-16.bam
```
### Running GATK
 
 First realign around indels. For each sample.bam file:
 ```sh
 export GENOME_REF=/home/ksilliman/CommonG/OlyBGI-scaffold-10k_ccn.fasta
java -Xmx5g -jar /usr/local/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /home/ksilliman/CommonG/OlyBGI-scaffold-10k_ccn.fasta -I sample.bam -o sample.intervals

java -Xmx5g -jar /usr/local/GenomeAnalysisTK.jar -T IndelRealigner -R /home/ksilliman/CommonG/OlyBGI-scaffold-10k_ccn.fasta -targetIntervals sample.intervals -I sample.bam -o sample.real.bam -LOD 0.4
 ```
 Run UnifiedGenotyper for a first pass in order to get base quality recalibration:
  ```sh
java -jar $GATK -T UnifiedGenotyper -R $GENOME_REF -nt 4 -nct 1  \
-I HC2-17A.real.bam \
-I HC2-18.real.bam \
-I HC3-1.real.bam \
-I HC3-10.real.bam \
-I HC3-11.real.bam \
-I HC3-7.real.bam \
-I HC2-17B.real.bam \
-I HC1-4.real.bam \
-I HC4-3.real.bam \
-I HC3-4.real.bam \
-I HC4-6.real.bam \
-I HC5-12.real.bam \
-I HC2-1.real.bam \
-I HC4-9.real.bam \
-I HC4-14.real.bam \
-I HC4-4A.real.bam \
-I HC4-4B.real.bam \
-I HC1-17.real.bam \
-I HC2-8.real.bam \
-I HC5-10.real.bam \
-I HC3-17.real.bam \
-I HC2-14.real.bam \
-I SS3-14.real.bam \
-I SS3-15.real.bam \
-I SS3-16.real.bam \
-I SS3-3.real.bam \
-I SS3-20.real.bam \
-I SS2-14.real.bam \
-I SS5-18.real.bam \
-I SS5-7.real.bam \
-I SS1-8.real.bam \
-I SS2-15.real.bam \
-I SS4-8.real.bam \
-I SS3-2.real.bam \
-I SS3-6.real.bam \
-I SS2-12A.real.bam \
-I SS2-12B.real.bam \
-I SS2-4A.real.bam \
-I SS2-4B.real.bam \
-I SS3-21.real.bam \
-I SS5-17.real.bam \
-I SS4-13.real.bam \
-I SS3-7.real.bam \
-I SS2-13.real.bam \
-I NF2-15.real.bam \
-I NF5-19.real.bam \
-I NF5-14.real.bam \
-I NF2-12.real.bam \
-I NF3-17.real.bam \
-I NF1-8.real.bam \
-I NF4-12.real.bam \
-I NF2-6A.real.bam \
-I NF2-6B.real.bam \
-I NF2-6C.real.bam \
-I NF1-14A.real.bam \
-I NF1-14B.real.bam \
-I NF5-3.real.bam \
-I NF2-9.real.bam \
-I NF3-11.real.bam \
-I NF1-19.real.bam \
-I NF4-13.real.bam \
-I NF3-9.real.bam \
-I NF4-9.real.bam \
-I NF5-11.real.bam \
-I NF4-11.real.bam \
-I NF3-5.real.bam \
-I NF3-16.real.bam \
-o round1.vcf

python ~/CommonG/z0on/GetHighQualVcfs.py -i round1.vcf --percentile 75 -o .
 ```
Base quality score recalibration (BQSR), for each sample:
```sh
java -Xmx20g -jar /usr/local/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/ksilliman/CommonG/OlyBGI-scaffold-10k_ccn.fasta -knownSites sample.m_HQ.vcf -I sample.real.bam -o sample.recalibration_report.grp
java -Xmx10g -jar /usr/local/GenomeAnalysisTK.jar -T PrintReads -R /home/ksilliman/CommonG/OlyBGI-scaffold-10k_ccn.fasta  -I sample.bam -BQSR sample.recalibration_report.grp -o sample.recal.bam
```
2nd iteration of UnifiedGenotyper on quality-recalibrated files:
```sh
java -jar $GATK -T UnifiedGenotyper -R $GENOME_REF -nt 6 -nct 1  --genotype_likelihoods_model SNP \
-I HC2-17A.recal.bam \
-I HC2-18.recal.bam \
-I HC3-1.recal.bam \
-I HC3-10.recal.bam \
-I HC3-11.recal.bam \
-I HC3-7.recal.bam \
-I HC2-17B.recal.bam \
-I HC1-4.recal.bam \
-I HC4-3.recal.bam \
-I HC3-4.recal.bam \
-I HC4-6.recal.bam \
-I HC5-12.recal.bam \
-I HC2-1.recal.bam \
-I HC4-9.recal.bam \
-I HC4-14.recal.bam \
-I HC4-4A.recal.bam \
-I HC4-4B.recal.bam \
-I HC1-17.recal.bam \
-I HC2-8.recal.bam \
-I HC5-10.recal.bam \
-I HC3-17.recal.bam \
-I HC2-14.recal.bam \
-I SS3-14.recal.bam \
-I SS3-15.recal.bam \
-I SS3-16.recal.bam \
-I SS3-3.recal.bam \
-I SS3-20.recal.bam \
-I SS2-14.recal.bam \
-I SS5-18.recal.bam \
-I SS5-7.recal.bam \
-I SS1-8.recal.bam \
-I SS2-15.recal.bam \
-I SS4-8.recal.bam \
-I SS3-2.recal.bam \
-I SS3-6.recal.bam \
-I SS2-12A.recal.bam \
-I SS2-12B.recal.bam \
-I SS2-4A.recal.bam \
-I SS2-4B.recal.bam \
-I SS3-21.recal.bam \
-I SS5-17.recal.bam \
-I SS4-13.recal.bam \
-I SS3-7.recal.bam \
-I SS2-13.recal.bam \
-I NF2-15.recal.bam \
-I NF5-19.recal.bam \
-I NF5-14.recal.bam \
-I NF2-12.recal.bam \
-I NF3-17.recal.bam \
-I NF1-8.recal.bam \
-I NF4-12.recal.bam \
-I NF2-6A.recal.bam \
-I NF2-6B.recal.bam \
-I NF2-6C.recal.bam \
-I NF1-14A.recal.bam \
-I NF1-14B.recal.bam \
-I NF5-3.recal.bam \
-I NF2-9.recal.bam \
-I NF3-11.recal.bam \
-I NF1-19.recal.bam \
-I NF4-13.recal.bam \
-I NF3-9.recal.bam \
-I NF4-9.recal.bam \
-I NF5-11.recal.bam \
-I NF4-11.recal.bam \
-I NF3-5.recal.bam \
-I NF3-16.recal.bam \
-o round2.vcf

#rename sample names to get rid of .m
cat round2.vcf | perl -pe 's/\.m//g' | perl -pe 's/^chrom/chr/' >round2.names.vcf
```
### Using Replicates to recalibrate genotype calls
For these next steps I made a tab-delimited file listing the replicate pairs in the test set called testreps.tab. I tried to run GATK's variant quality score recalibration (VQSR), but it would not work- likely because there aren't enough variants in this test set. Here's the code I tried:
```sh
~/CommonG/z0on/replicatesMatch.pl vcf=round2.names.vcf replicates=testreps.tab > vqsr.vcf
```
```
7995 total SNPs
1443 pass hets and match filters
948 show non-reference alleles
738 have alterantive alleles in at least 2 replicate pair(s)
738 have matching heterozygotes in at least 0 replicate pair(s)
317 polymorphic
738 written
```
```sh
vcftools --vcf vqsr.vcf --TsTv-summary
```
```
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
```
```sh
java -jar $GATK -T VariantRecalibrator -R $GENOME_REF -input round2.names.vcf -nt 6 \
-resource:repmatch,known=true,training=true,truth=true,prior=10  vqsr.vcf \
-an QD -an MQ -an FS -mode SNP --maxGaussians 4 \
--target_titv 1.463 -tranche 85.0 -tranche 90.0 -tranche 95.0 -tranche 99.0 -tranche 100 \
-recalFile round2.recal -tranchesFile recalibrate_SNP.tranches --rscript_file recalibrate_SNP_plots.R 
```
Here's where it would throw an error:
```
##### ERROR MESSAGE: NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --minNumBadVariants 5000, for example).
##### ERROR ------------------------------------------------------------------------------------------

```
So instead I used Matz's non-parametric quantile-based recalibration scripts.
```sh
~/CommonG/z0on/replicatesMatch.pl vcf=round2.names.vcf replicates=testreps.tab polyonly=1 > vqsr_matz.vcf
```
```
7995 total SNPs
1443 pass hets and match filters
948 show non-reference alleles
738 have alterantive alleles in at least 2 replicate pair(s)
738 have matching heterozygotes in at least 0 replicate pair(s)
317 polymorphic
317 written
```
To choose a filtering setting for later, you want the setting with the greatest "gain". After trying the different settings, I chose the one below and a filter of 5:
```sh
recalibrateSNPs_gatk.pl vcf=round2.names.vcf true=vqsr_matz.vcf -nofs >gatk_after_vqsr_matz.vcf
```
```
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
1	0.00100
5	0.00750
10	0.01500
15	0.02250
20	0.02700
30	0.04200
40	0.06400
50	0.09000
60	0.12000
70	0.16800
80	0.25200
90	0.40000

------------------------
43.04%	at qual <1 (42.04% gain)
66.62%	at qual <5 (61.62% gain)
71.14%	at qual <10 (61.14% gain)
74.11%	at qual <15 (59.11% gain)
75.85%	at qual <20 (55.85% gain)
79.71%	at qual <30 (49.71% gain)
------------------------
```
```sh
retabvcf.pl vcf=gatk_after_vqsr_matz.vcf tab=/home/ksilliman/CommonG/OlyBGI-scaffold-10k_cc.tab > retab.vcf
# discarding loci with too many heterozygotes, which are likely lumped paralogs
# (by default, fraction of heterozygotes should not exceed maxhet=0.75)
# this step can also filter for the fraction of missing genotypes (default maxmiss=0.5)
hetfilter.pl vcf=retab.vcf > hetfilt_def.vcf
```
```
7995 total loci
4351 dropped because fraction of missing genotypes exceeded 0.5
134 dropped because fraction of heterozygotes exceeded 0.75
3510 written
```
Thinning to 1 SNP per tag, leaves SNP with highest minor allele frequency:
```sh
thinner.pl vcf=hetfilt_def.vcf > thin_def.vcf
```
```
3510 total loci
257 loci skipped because they were closer than 40
1765 loci selected
```

Applying filter and selecting polymorphic biallelic loci genotyped in 70% or more individuals. --minQ is where you set filter established earlier.
```sh
vcftools --vcf thin_def.vcf --minQ 5 --max-missing 0.7  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out filt0
```
```
After filtering, kept 67 out of 67 Individuals
After filtering, kept 755 out of a possible 1765 Sites
```

Checking genotypic match between pairs of replicates
```sh
repMatchStats.pl vcf=filt0.recode.vcf replicates=testreps.tab 
```
```
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
```
SS2-12A and SS2-12B had lo match so excluding these for now from further analyses. I created a file called testreps_del.tab listing which replicate to delete, based on which replicate had the most mapped reads.
```sh
vcftools --vcf thin_def.vcf --remove testreps_del.tab --minQ 5 --max-missing 0.7  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out testfinal.7
```
```
After filtering, kept 59 out of 67 Individuals
Outputting VCF file...
After filtering, kept 882 out of a possible 1765 Sites
```
for Tajima's D calculation:
```sh
vcftools --vcf hetfilt_def.vcf --remove testreps_del.tab --minQ 5 --max-missing 0.7  --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out testfinalUnthin.7
```
```
After filtering, kept 1855 out of a possible 3510 Sites
```












