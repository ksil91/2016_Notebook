## Trimming and quality filtering Dec2015 2bRAD sequence data  
I am roughly following [Eli Meyer's 2brad geneotyping guide](http://people.oregonstate.edu/~meyere/2bRAD_analysis2.0.html) with some scripts also taken [from Mikhail Matz's 2bRAD github](https://github.com/z0on/2bRAD_GATK). Most analyses were run on my department's single node, 24 processor high performance computer. For brevity I don't include the .sh scripts to submit the job to the HPC, just the commands that were run.  
Raw sequencing files are found at http://owl.fish.washington.edu/nightingales/O_lurida/2bRAD_Dec15/
### Required scripts and programs  
For these steps, you need: 
- Bioperl
- [sequence_processing scripts](https://github.com/Eli-Meyer/sequence_processing) from Eli Meyer's github
- 2bRAD_trim_launch.pl and trim2bRAD.pl from z0on github  
### Trimming adapter sequences  
As the fragments produced by AlfI are only 36bp in length and we sequenced our data to 50bp, the extra 14bp will be adapter sequences. I chose to use Matz's script for this as it checks for the expected adapter sequence on the end of the read (a quality test in itself) before trimming. Also the way that the Meyer protocol checks for adapter sequences in not intuitive for me as it requires parameters for how much a sequence should match the adapter to be discarded.  
```sh
for file in *.fastq.gz;
    do gunzip $file;
done

./2bRAD_trim_launch.pl fastq sampleID=1 site=.{12}GCA.{6}TGC.{12} > trims
```
This produces a file called trims with a call to trim2bRAD.pl for each sample. This can easily be added to a shell script for job submission.
Ex:
```sh
./trim2bRAD.pl SS2-4A_TGACCA-GATCTCT_L001_R1_001.fastq ".{12}GCA.{6}TGC.{12}" "AGATCGGAA" >SS2-4A.tr0
```
After running this, I changed some of the samples names that were incorrect on the original barcode sheet. Trimmed sequences with the correct names are in: http://owl.fish.washington.edu/wetgenes/2brad_201512_trim
### Quality filtering
```sh
for file in *.tr0; 
	do QualFilterFastq.pl $file 25 18 ${ffile/.tr0/_q0.fastq};
done
mv *q0.fastq qual_25_18/ 
```
This is using the Meyer lab's quality filtering script to remove reads that have more than 18 bases with quality scores less than 25. I then moved the results into the folder qual_25_18 in wetgenes. Will check the number of reads retained per sample and then rerun at different cutoffs.


