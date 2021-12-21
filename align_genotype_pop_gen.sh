#alignment, genotyping, and pop gen stats workflow for Woodruff et al. 2021 C. inopinata population genomics paper, "Alignment of genetic differentiation across trophic levels in a fig community"

#if you have any questions, please contact Gavin Woodruff at gcwoodruff@ou.edu

#put working directory here
wkdir="/projects/phillipslab/gavincw/pop_gen_reproduce_june_2021"

#demultiplexed, re-named fastq files are in folder $wkdir/05_reads_to_align/

#get number of reads per sample

cd $wkdir/05_reads_to_align/
for i in *; do echo $(cat $i |wc -l)/4|bc; done

#7193109
#7193109
#3460357
#3460357
#717508
#717508
#19731572
#19731572
#10900527
#10900527
#22514736
#22514736
#26836136
#26836136
#25936266
#25936266
#19438200
#19438200
#3129605
#3129605
#1297288
#1297288
#22758380
#22758380
#660719
#660719
#22720765
#22720765
#15504982
#15504982
#5342526
#5342526
#9187534
#9187534
#19813439
#19813439
#11079441
#11079441
#36402219
#36402219
#9327704
#9327704
#1448396
#1448396
#37482621
#37482621
#8146378
#8146378
	#these values are in a supplemental table

mkdir $wkdir/gsnap_genome_db

gmap_build --dir=$wkdir/gsnap_genome_db -d inopinata_genome_unmasked $wkdir/inopinata_genome/inopinata_genome.fa 


#align reads to reference

mkdir $wkdir/06_gsnap/

gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/A03.1.fq $wkdir/05_reads_to_align/A03.2.fq  > $wkdir/06_gsnap/A03.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/A05.1.fq $wkdir/05_reads_to_align/A05.2.fq  > $wkdir/06_gsnap/A05.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/B03.1.fq $wkdir/05_reads_to_align/B03.2.fq  > $wkdir/06_gsnap/B03.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/dec_2016_D01.1.fq $wkdir/05_reads_to_align/dec_2016_D01.2.fq  > $wkdir/06_gsnap/dec_2016_D01.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/dec_2016_D02.1.fq $wkdir/05_reads_to_align/dec_2016_D02.2.fq  > $wkdir/06_gsnap/dec_2016_D02.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/dec_2016_D03.1.fq $wkdir/05_reads_to_align/dec_2016_D03.2.fq  > $wkdir/06_gsnap/dec_2016_D03.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/dec_2016_D04.1.fq $wkdir/05_reads_to_align/dec_2016_D04.2.fq  > $wkdir/06_gsnap/dec_2016_D04.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/dec_2016_D05.1.fq $wkdir/05_reads_to_align/dec_2016_D05.2.fq  > $wkdir/06_gsnap/dec_2016_D05.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/dec_2016_D10.1.fq $wkdir/05_reads_to_align/dec_2016_D10.2.fq  > $wkdir/06_gsnap/dec_2016_D10.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_A11.1.fq $wkdir/05_reads_to_align/may_2017_A11.2.fq  > $wkdir/06_gsnap/may_2017_A11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_A12.1.fq $wkdir/05_reads_to_align/may_2017_A12.2.fq  > $wkdir/06_gsnap/may_2017_A12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_B11.1.fq $wkdir/05_reads_to_align/may_2017_B11.2.fq  > $wkdir/06_gsnap/may_2017_B11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_B12.1.fq $wkdir/05_reads_to_align/may_2017_B12.2.fq  > $wkdir/06_gsnap/may_2017_B12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_C11.1.fq $wkdir/05_reads_to_align/may_2017_C11.2.fq  > $wkdir/06_gsnap/may_2017_C11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_C12.1.fq $wkdir/05_reads_to_align/may_2017_C12.2.fq  > $wkdir/06_gsnap/may_2017_C12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_D11.1.fq $wkdir/05_reads_to_align/may_2017_D11.2.fq  > $wkdir/06_gsnap/may_2017_D11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_D12.1.fq $wkdir/05_reads_to_align/may_2017_D12.2.fq  > $wkdir/06_gsnap/may_2017_D12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_E11.1.fq $wkdir/05_reads_to_align/may_2017_E11.2.fq  > $wkdir/06_gsnap/may_2017_E11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_E12.1.fq $wkdir/05_reads_to_align/may_2017_E12.2.fq  > $wkdir/06_gsnap/may_2017_E12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_F11.1.fq $wkdir/05_reads_to_align/may_2017_F11.2.fq  > $wkdir/06_gsnap/may_2017_F11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_G12.1.fq $wkdir/05_reads_to_align/may_2017_G12.2.fq  > $wkdir/06_gsnap/may_2017_G12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_H10.1.fq $wkdir/05_reads_to_align/may_2017_H10.2.fq  > $wkdir/06_gsnap/may_2017_H10.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_H11.1.fq $wkdir/05_reads_to_align/may_2017_H11.2.fq  > $wkdir/06_gsnap/may_2017_H11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/05_reads_to_align/may_2017_H12.1.fq $wkdir/05_reads_to_align/may_2017_H12.2.fq  > $wkdir/06_gsnap/may_2017_H12.sam


#extract unique aligning reads

	#gsnap provides flags for that.
	#"XO:Z:CU" for concordant unique
	#"XO:Z:HU" , half mapping unique
	#"XO:Z:PI" , paired unique inversion
	#"XO:Z:PS" , paired unique scramble
	#"XO:Z:PL" , paired unique long
	#"XO:Z:UU" , unpaired unique
		#putting these values in a file for grep , file /projects/phillipslab/gavincw/pop_gen_reproduce_june_2021/additional_files/gsnap_sam_flags (for GCW: /projects/phillipslab/gavincw/pop_gen_4-2019/10_gsnap_pipeline_prev_reads/gsnap_sam_flags)



mkdir $wkdir/07_grep_uniq_alignments
cd $wkdir/07_grep_uniq_alignments
mkdir uniq_alignments
mkdir headers
mkdir cat

#get the unique reads

cd $wkdir/06_gsnap

for i in *; do LC_ALL=C fgrep -w -f $wkdir/additional_files/gsnap_sam_flags $i > $wkdir/07_grep_uniq_alignments/uniq_alignments/$i; done &

#get the headers

cd $wkdir/06_gsnap

for i in *; do LC_ALL=C fgrep "@" $i > $wkdir/07_grep_uniq_alignments/headers/$i; done &

#concatenate headers and unique reads

cd $wkdir/07_grep_uniq_alignments/headers/

for i in *; do cat $i $wkdir/07_grep_uniq_alignments/uniq_alignments/$i > $wkdir/07_grep_uniq_alignments/cat/$i; done &

#get rid of extra stuff

cd $wkdir/07_grep_uniq_alignments/

rm -r uniq_alignments
rm -r headers


#sort alignments
 #for u of oregon's hpc: module load racs-eb SAMtools/1.7-intel-2018a

mkdir $wkdir/08_samtools_sort/

cd $wkdir/07_grep_uniq_alignments/cat/

samtools view -bS A03.sam  | samtools sort > $wkdir/08_samtools_sort/A03.bam
samtools view -bS A05.sam  | samtools sort > $wkdir/08_samtools_sort/A05.bam
samtools view -bS B03.sam  | samtools sort > $wkdir/08_samtools_sort/B03.bam
samtools view -bS dec_2016_D01.sam  | samtools sort > $wkdir/08_samtools_sort/dec_2016_D01.bam
samtools view -bS dec_2016_D02.sam  | samtools sort > $wkdir/08_samtools_sort/dec_2016_D02.bam
samtools view -bS dec_2016_D03.sam  | samtools sort > $wkdir/08_samtools_sort/dec_2016_D03.bam
samtools view -bS dec_2016_D04.sam  | samtools sort > $wkdir/08_samtools_sort/dec_2016_D04.bam
samtools view -bS dec_2016_D05.sam  | samtools sort > $wkdir/08_samtools_sort/dec_2016_D05.bam
samtools view -bS dec_2016_D10.sam  | samtools sort > $wkdir/08_samtools_sort/dec_2016_D10.bam
samtools view -bS may_2017_A11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_A11.bam
samtools view -bS may_2017_A12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_A12.bam
samtools view -bS may_2017_B11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_B11.bam
samtools view -bS may_2017_B12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_B12.bam
samtools view -bS may_2017_C11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_C11.bam
samtools view -bS may_2017_C12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_C12.bam
samtools view -bS may_2017_D11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_D11.bam
samtools view -bS may_2017_D12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_D12.bam
samtools view -bS may_2017_E11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_E11.bam
samtools view -bS may_2017_E12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_E12.bam
samtools view -bS may_2017_F11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_F11.bam
samtools view -bS may_2017_G12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_G12.bam
samtools view -bS may_2017_H10.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_H10.bam
samtools view -bS may_2017_H11.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_H11.bam
samtools view -bS may_2017_H12.sam  | samtools sort > $wkdir/08_samtools_sort/may_2017_H12.bam

#get number of mapped reads per sample

cd $wkdir/08_samtools_sort/
mkdir $wkdir/09_samtools_flagstat/


for i in *; do samtools flagstat $i > $wkdir/09_samtools_flagstat/$i; done &
cd $wkdir/09_samtools_flagstat/

for i in *; do head -1 $i | sed 's/ + .*//g'; done

#13282662
#6259856
#1319264
#33774381
#17775212
#32627339
#36269045
#38122436
#32080204
#5091126
#2259416
#38861026
#1004678
#36905099
#28321158
#9585812
#16374130
#36224016
#19233394
#44181645
#16677504
#2630648
#41035051
#13975100

	#the above values are in a supplemental table

#call genotypes with bcftools mpileup (v 1.9)

#ploidy file in $wkdir/additional_files/ploidy_file.txt
#sample files in $wkdir/additional_files/sample_files/

mkdir $wkdir/10_bcftools_mpileup_call
cd $wkdir/08_samtools_sort/

bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa A03.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/A03.bam -o $wkdir/10_bcftools_mpileup_call/A03.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa A05.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/A05.bam -o $wkdir/10_bcftools_mpileup_call/A05.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa B03.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/B03.bam -o $wkdir/10_bcftools_mpileup_call/B03.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa dec_2016_D01.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/dec_2016_D01.bam -o $wkdir/10_bcftools_mpileup_call/dec_2016_D01.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa dec_2016_D02.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/dec_2016_D02.bam -o $wkdir/10_bcftools_mpileup_call/dec_2016_D02.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa dec_2016_D03.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/dec_2016_D03.bam -o $wkdir/10_bcftools_mpileup_call/dec_2016_D03.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa dec_2016_D04.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/dec_2016_D04.bam -o $wkdir/10_bcftools_mpileup_call/dec_2016_D04.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa dec_2016_D05.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/dec_2016_D05.bam -o $wkdir/10_bcftools_mpileup_call/dec_2016_D05.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa dec_2016_D10.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/dec_2016_D10.bam -o $wkdir/10_bcftools_mpileup_call/dec_2016_D10.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_A11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_A11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_A11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_A12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_A12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_A12.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_B11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_B11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_B11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_B12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_B12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_B12.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_C11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_C11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_C11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_C12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_C12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_C12.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_D11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_D11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_D11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_D12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_D12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_D12.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_E11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_E11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_E11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_E12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_E12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_E12.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_F11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_F11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_F11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_G12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_G12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_G12.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_H10.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_H10.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_H10.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_H11.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_H11.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_H11.bcf
bcftools mpileup -Ou -f $wkdir/inopinata_genome/inopinata_genome.fa may_2017_H12.bam | bcftools call -c -Ob --ploidy-file $wkdir/additional_files/ploidy_file.txt --samples-file $wkdir/additional_files/sample_files/may_2017_H12.bam -o $wkdir/10_bcftools_mpileup_call/may_2017_H12.bcf

#get only sites with >14x coverage

mkdir $wkdir/11_bcftools_view_depth_15


cd $wkdir/10_bcftools_mpileup_call/

for i in *; do bcftools view -i 'DP>=15' $i > $wkdir/11_bcftools_view_depth_15/${i%.bcf}.vcf; done &

#convert back to bcf and index


mkdir $wkdir/12_bcftools_view_bcf/

cd $wkdir/11_bcftools_view_depth_15

for i in *; do bcftools view $i -O b > $wkdir/12_bcftools_view_bcf/${i%.vcf}.bcf; done &

cd $wkdir/12_bcftools_view_bcf/

for i in *; do bcftools index $i; done


#merge all bcf

mkdir $wkdir/13_bcftools_merge/

bcftools merge --info-rules DP:join,MQ0F:join,AF1:join,AC1:join,DP4:join,MQ:join,FQ:join -m snps A03.bcf A05.bcf B03.bcf dec_2016_D01.bcf dec_2016_D02.bcf dec_2016_D03.bcf dec_2016_D04.bcf dec_2016_D05.bcf dec_2016_D10.bcf may_2017_A11.bcf may_2017_A12.bcf may_2017_B11.bcf may_2017_B12.bcf may_2017_C11.bcf may_2017_C12.bcf may_2017_D11.bcf may_2017_D12.bcf may_2017_E11.bcf may_2017_E12.bcf may_2017_F11.bcf may_2017_G12.bcf may_2017_H10.bcf may_2017_H11.bcf may_2017_H12.bcf -o $wkdir/13_bcftools_merge/inopinata_24.vcf


#split autosomes and X for genotyping filter

mkdir $wkdir/14_split_autosomes_X/
mkdir $wkdir/14_split_autosomes_X/autosomes
mkdir $wkdir/14_split_autosomes_X/X

#autosomes

cd $wkdir/13_bcftools_merge/

#put aside header
grep "#" inopinata_24.vcf > $wkdir/14_split_autosomes_X/autosomes/00_header

#get autosomes
grep -v "#" inopinata_24.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 != "X"' > $wkdir/14_split_autosomes_X/autosomes/01_autosomes_no_header

#combine and sort vcf
cd $wkdir/14_split_autosomes_X/autosomes/

cat 00_header 01_autosomes_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_24_autosomes.vcf


#X chromosome

cd $wkdir/13_bcftools_merge/

#put aside header
grep "#" inopinata_24.vcf > $wkdir/14_split_autosomes_X/X/00_header


#get X
grep -v "#" inopinata_24.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 == "X"' > $wkdir/14_split_autosomes_X/X/01_X_no_header 

#combine and sort vcf
cd $wkdir/14_split_autosomes_X/X/

cat 00_header 01_X_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_24_X.vcf 



#just get sites with ~ >80% of samples having genotype calls

mkdir $wkdir/15_bcftools_view_gt_count/

cd $wkdir/14_split_autosomes_X/autosomes/

	#5/24 = 0.2083333
bcftools view -i 'COUNT(GT="mis")<5' inopinata_24_autosomes.vcf > $wkdir/15_bcftools_view_gt_count/inopinata_24_autosomes.vcf &

cd $wkdir/14_split_autosomes_X/X/

	#there are four animals with unknown sex, and only using females for X stats, so the denominator is smaller for the X; 4/18 = 0.22
bcftools view -i 'COUNT(GT="mis")<4' inopinata_24_X.vcf > $wkdir/15_bcftools_view_gt_count/inopinata_24_X.vcf &



#combine autosomes and X

mkdir $wkdir/16_cat/

cd $wkdir/15_bcftools_view_gt_count/
#set aside header
grep "#" inopinata_24_autosomes.vcf > $wkdir/16_cat/00_header &
#autosomes
grep -v "#" inopinata_24_autosomes.vcf >  $wkdir/16_cat/01_inopinata_24_autosomes_no_header &
#x
grep -v "#" inopinata_24_X.vcf > $wkdir/16_cat/02_inopinata_24_X_no_header &

cd $wkdir/16_cat/
#combine vcf
cat 00_header 01_inopinata_24_autosomes_no_header 02_inopinata_24_X_no_header > 03_combined_unsorted.vcf &
#sort
cat 03_combined_unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_24_sorted.vcf &


#next, invariant sites and variant sites to get biallelic snps

#get just the invariant sites....

mkdir $wkdir/17_invariant_sites/

cd $wkdir/16_cat/

grep "#" inopinata_24_sorted.vcf > $wkdir/17_invariant_sites/header &

grep -v "#" inopinata_24_sorted.vcf > $wkdir/17_invariant_sites/inopinata_24_no_header &

cd $wkdir/17_invariant_sites/

awk 'BEGIN {FS="\t"} {OFS="\t"} $5 == "."' inopinata_24_no_header > inopinata_24_no_alt &

cat header inopinata_24_no_alt > inopinata_24_invariant_sites.vcf

#get biallelic snps with at least one copy per allele

mkdir $wkdir/18_bcftools_view_biallelic_snps_maf/

cd $wkdir/16_cat/

bcftools view -m2 -M2 -v snps --min-ac 2:minor inopinata_24_sorted.vcf > $wkdir/18_bcftools_view_biallelic_snps_maf/inopinata_24.vcf &
	#this file is also called inopinata_24_biallelic_snps.vcf and used for pca, clustering



#combine

mkdir $wkdir/19_cat_biallelic_invariant

cd $wkdir/18_bcftools_view_biallelic_snps_maf/

grep "#" inopinata_24.vcf  > $wkdir/19_cat_biallelic_invariant/00_header

grep -v "#" inopinata_24.vcf  >  $wkdir/19_cat_biallelic_invariant/01_inopinata_24_biallelic_snps_no_header

cd $wkdir/17_invariant_sites/

grep -v "#" inopinata_24_invariant_sites.vcf > $wkdir/19_cat_biallelic_invariant/02_inopinata_24_invariant_sites_no_header &

cd $wkdir/19_cat_biallelic_invariant/

cat 00_header 01_inopinata_24_biallelic_snps_no_header 02_inopinata_24_invariant_sites_no_header > 03_combined_unsorted.vcf &
#sort
cat 03_combined_unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_24_sorted.vcf &





#ok, split autosomes, male x and female x for pop gen stats. 

mkdir $wkdir/20_split_autosomes_X_sex

cd $wkdir/20_split_autosomes_X_sex

mkdir 00_autosomes


cd $wkdir/19_cat_biallelic_invariant/

#set aside header
grep "#" inopinata_24_sorted.vcf > $wkdir/20_split_autosomes_X_sex/00_autosomes/00_header

#autosomes
grep -v "#" inopinata_24_sorted.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 != "X"' > $wkdir/20_split_autosomes_X_sex/00_autosomes/01_autosomes_no_header &

cd $wkdir/20_split_autosomes_X_sex/00_autosomes/
#sort
cat 00_header 01_autosomes_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_24_autosomes.vcf &



#separate males and females 

mkdir $wkdir/20_split_autosomes_X_sex/01_split_by_sex

cd $wkdir/19_cat_biallelic_invariant/

#males
bcftools view -s B03.bam,may_2017_H12.bam inopinata_24_sorted.vcf > $wkdir/20_split_autosomes_X_sex/01_split_by_sex/inopinata_males.vcf &

#females
bcftools view -s A03.bam,dec_2016_D01.bam,dec_2016_D02.bam,dec_2016_D03.bam,dec_2016_D04.bam,dec_2016_D05.bam,may_2017_A11.bam,may_2017_A12.bam,may_2017_B11.bam,may_2017_C11.bam,may_2017_C12.bam,may_2017_D11.bam,may_2017_D12.bam,may_2017_E11.bam,may_2017_E12.bam,may_2017_F11.bam,may_2017_H10.bam,may_2017_H11.bam inopinata_24_sorted.vcf > $wkdir/20_split_autosomes_X_sex/01_split_by_sex/inopinata_females.vcf &


#extract X for males and females
mkdir $wkdir/20_split_autosomes_X_sex/02_split_X

cd $wkdir/20_split_autosomes_X_sex/01_split_by_sex/

for i in *; do grep "#" $i > $wkdir/20_split_autosomes_X_sex/02_split_X/$i.header; done &

for i in *; do grep -v "#" $i | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 == "X"' > $wkdir/20_split_autosomes_X_sex/02_split_X/${i%.vcf}.X.vcf_no_header; done &

cd $wkdir/20_split_autosomes_X_sex/02_split_X/

cat inopinata_females.vcf.header inopinata_females.X.vcf_no_header  | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_females_X.vcf &

cat inopinata_males.vcf.header inopinata_males.X.vcf_no_header  | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_males_X.vcf &




#ok, cool, popgenWindows (retrieved 2019-05-28 ; https://github.com/simonhmartin/genomics_general) for estimating pi in 10 kb genomic windows -- autosomes

mkdir $wkdir/21_popgenWindows/
mkdir $wkdir/21_popgenWindows/00_parseVCF/
mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/


#using Biopython v 1.70

python2.7 parseVCF.py -i  $wkdir/20_split_autosomes_X_sex/00_autosomes/inopinata_24_autosomes.vcf | gzip > $wkdir/21_popgenWindows/00_parseVCF/inopinata_24_autosomes.geno.gz &

python2.7 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/21_popgenWindows/00_parseVCF/inopinata_24_autosomes.geno.gz -o  $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_autosomes_stats.csv.gz -f phased  &

cd $wkdir/21_popgenWindows/01_popgenWindows_py/

gunzip inopinata_autosomes_stats.csv.gz


#now, inopinata X -- only using inopinata females for X pop gen stats 

python2.7 parseVCF.py -i $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf | gzip > $wkdir/21_popgenWindows/00_parseVCF/inopinata_females_X.geno.gz &

#now estimate stats

python2.7 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g /$wkdir/21_popgenWindows/00_parseVCF/inopinata_females_X.geno.gz -o  $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_females_X.csv.gz -f phased  &


#combine the files

cd $wkdir/21_popgenWindows/01_popgenWindows_py/

gunzip inopinata_females_X.csv.gz

#remove header
sed '1d' inopinata_females_X.csv > inopinata_females_X_no_header

#combine autosomes and X pi estimates
cat inopinata_autosomes_stats.csv inopinata_females_X_no_header > inopinata_autosomes_females_X_pi.csv

#add species id's and C. elegans data... see elegans.sh for elegans data

awk 'BEGIN {FS=","} {OFS=","} {print $0,"C. inopinata"}' inopinata_autosomes_females_X_pi.csv > inopinata_autosomes_females_X_pi.tmp.csv

#repair header
sed '0,/C. inopinata/{s/C. inopinata/species/}' inopinata_autosomes_females_X_pi.tmp.csv > inopinata_autosomes_females_X_pi.csv 

#comibine with C. elegans data
	#see elegans.sh for elegans data
cat inopinata_autosomes_females_X_pi.csv $wkdir/elegans/11_popgenwindows_py/elegans_stats.csv > ino_elg_pi.csv


#inopinata FST

mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst

cp /projects/phillipslab/gavincw/pop_gen_Oct_2019/16_re-do_M_F_genotyping/11_popgenWindows/01_popgenWindows_py/inopinata_fst/pops_file_islands $wkdir/additional_files/inopinata_pops_file_islands

cp /projects/phillipslab/gavincw/pop_gen_Oct_2019/16_re-do_M_F_genotyping/11_popgenWindows/01_popgenWindows_py/inopinata_fst/pops_file_islands_females $wkdir/additional_files/inopinata_pops_file_islands_females

#make sure popgenwindows.py is in path
#cd /projects/phillipslab/gavincw/download/genomics_general/

	#fst autosomes

python2.7 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 --analysis popPairDist -g $wkdir/21_popgenWindows/00_parseVCF/inopinata_24_autosomes.geno.gz -o $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/inopinata_autosomes_fst_pops_file.csv.gz -f phased -p ishigaki -p iriomote -p yonaguni --popsFile $wkdir/additional_files/inopinata_pops_file_islands &

	#fst female X

python2.7 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 --analysis popPairDist -g $wkdir/21_popgenWindows/00_parseVCF/inopinata_females_X.geno.gz -o $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/inopinata_females_X_fst.csv.gz -f phased -p ishigaki -p iriomote -p yonaguni --popsFile $wkdir/additional_files/inopinata_pops_file_islands_females &


cd $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst

#prepare FST file for stats and figures
gunzip inopinata_autosomes_fst_pops_file.csv.gz
gunzip inopinata_females_X_fst.csv.gz

mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/


# $1,$2,
#$9 = Fst_ishigaki_iriomote
#$10 = Fst_ishigaki_yonaguni
#$11 = Fst_iriomote_yonaguni

#split by population pair to label by population pair

sed '1d' inopinata_autosomes_fst_pops_file.csv | awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3,$9,"ishi-irio"}' > $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/00_ishi-irio_1

sed '1d' inopinata_autosomes_fst_pops_file.csv | awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3,$10,"ishi-yona"}' > $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/01_ishi-yona_1

sed '1d' inopinata_autosomes_fst_pops_file.csv | awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3,$11,"irio-yona"}' > $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/02_irio-yona_1

#recombine labeled data

cd $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/


cat 00_ishi-irio_1 01_ishi-yona_1 02_irio-yona_1 > 04_cat_fst

#remove "nan" ie missing data
grep -v "nan" 04_cat_fst > 05_remove_nan_fst

#extract fst values
awk 'BEGIN {FS=","} {OFS=","} {print $4}' 05_remove_nan_fst > 06_fst_values

##if interested in replacing negative values with zero
##awk '{print($0<0?0:$0)}' 06_fst_values > 07_replace_negative

#get genome coordinates
awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3}' 05_remove_nan_fst > 08_first_three_col

#extract pop pair label
awk 'BEGIN {FS=","} {OFS=","} {print $5}' 05_remove_nan_fst > 09_last_col

#put together coordinates, fst, pop pair label
paste 08_first_three_col 06_fst_values 09_last_col > 10_fst_tmp

sed -i -e 's/\t/,/g' 10_fst_tmp

#add header

echo -e "Chr,BP_start,BP_end,FST,pop_pair" | cat - 10_fst_tmp > autosomes_fst_combined.csv


#X chromosome FST



cd  $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst
mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/X


sed '1d' inopinata_females_X_fst.csv | awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3,$9,"ishi-irio"}' > $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/X/00_ishi-irio_1

sed '1d' inopinata_females_X_fst.csv | awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3,$10,"ishi-yona"}' > $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/X/01_ishi-yona_1


sed '1d' inopinata_females_X_fst.csv | awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3,$11,"irio-yona"}' > $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/X/02_irio-yona_1


#recombine labeled data
cd $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/X

cat 00_ishi-irio_1 01_ishi-yona_1 02_irio-yona_1 > 04_cat_fst

#remove "nan" ie missing data
grep -v "nan" 04_cat_fst > 05_remove_nan_fst


awk 'BEGIN {FS=","} {OFS=","} {print $4}' 05_remove_nan_fst > 06_fst_values

##if interested in replacing negative values with zero
##awk '{print($0<0?0:$0)}' 06_fst_values > 07_replace_negative

#get genome coordinates
awk 'BEGIN {FS=","} {OFS=","} {print $1,$2,$3}' 05_remove_nan_fst > 08_first_three_col

#extract pop pair label
awk 'BEGIN {FS=","} {OFS=","} {print $5}' 05_remove_nan_fst > 09_last_col

#put together coordinates, fst, pop pair label
paste 08_first_three_col 06_fst_values 09_last_col > 10_fst_tmp

sed -i -e 's/\t/,/g' 10_fst_tmp

#combine X and autosomes

cd $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/

cat autosomes_fst_combined.csv $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/X/10_fst_tmp > all_fst_combined.csv



#FIS with stacks

#cp /projects/phillipslab/gavincw/pop_gen_11-2018/12_samtools_mpileup/05_pca/08_population_map_okinawa $wkdir/additional_files/inopinata_autosomes_pop_map_stacks.txt

mkdir $wkdir/22_stacks_populations/
mkdir $wkdir/22_stacks_populations/autosomes/


populations -V $wkdir/20_split_autosomes_X_sex/00_autosomes/inopinata_24_autosomes.vcf -O $wkdir/22_stacks_populations/autosomes/ -M $wkdir/additional_files/inopinata_autosomes_pop_map_stacks.txt  --sigma 3333 --genepop --structure --phylip &


#prep FIS file for stats and figures
mkdir $wkdir/22_stacks_populations/autosomes/site_fis

cd $wkdir/22_stacks_populations/autosomes/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$3,$17}' inopinata_24_autosomes.p.sumstats.tsv > $wkdir/22_stacks_populations/autosomes/site_fis/01_fis_cols

cd $wkdir/22_stacks_populations/autosomes/site_fis


#remove first two lines

sed -i '1d' 01_fis_cols

sed -i '1d' 01_fis_cols

#if interested in replacing negative values with 0
#awk '$4 < 0 { $4=0 }1' 01_fis_cols > 02_fis_cols_replace_negative

#repair bed file (ie, replace spaces with tabs)
sed -i -e 's/ /\t/g' 01_fis_cols


#bedtools to get windowed Fis!



bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.50bp.windows -b 01_fis_cols > 03_fis_50bp_windows &


#remove missing data
awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 03_fis_50bp_windows > 04_fis_50bp_windows_no_missing_data

#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.10kb.windows -b 04_fis_50bp_windows_no_missing_data > 05_fis_10kb_windows &


#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 05_fis_10kb_windows > 06_fis_10kb_windows_no_missing_data


#X chromosome



mkdir $wkdir/22_stacks_populations/X/


populations -V $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf -O $wkdir/22_stacks_populations/X/ -M $wkdir/additional_files/inopinata_X_females_pop_map_stacks.txt  --sigma 3333 --genepop --structure --phylip &


#prep FIS file for stats and figures
mkdir $wkdir/22_stacks_populations/X/site_fis

cd $wkdir/22_stacks_populations/X/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$3,$17}' inopinata_females_X.p.sumstats.tsv > $wkdir/22_stacks_populations/X/site_fis/01_fis_cols

cd $wkdir/22_stacks_populations/X/site_fis


#remove first two lines

sed -i '1d' 01_fis_cols

sed -i '1d' 01_fis_cols

#if interested in replacing negative values with 0
#awk '$4 < 0 { $4=0 }1' 01_fis_cols > 02_fis_cols_replace_negative

#repair bed file (ie, replace spaces with tabs)
sed -i -e 's/ /\t/g' 01_fis_cols


#bedtools to get windowed Fis!


bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.50bp.windows.X -b 01_fis_cols > 03_fis_50bp_windows &

#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 03_fis_50bp_windows > 04_fis_50bp_windows_no_missing_data


#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.10kb.windows.X -b 04_fis_50bp_windows_no_missing_data > 05_fis_10kb_windows &

#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 05_fis_10kb_windows > 06_fis_10kb_windows_no_missing_data




#combine autosomes and X chromosome

mkdir $wkdir/22_stacks_populations/combine_fis/

cat $wkdir/22_stacks_populations/autosomes/site_fis/06_fis_10kb_windows_no_missing_data $wkdir/22_stacks_populations/X/site_fis/06_fis_10kb_windows_no_missing_data > $wkdir/22_stacks_populations/combine_fis/08_fis_aut_and_x_no_missing_data


awk 'BEGIN {OFS="\t"} {print $0, "inopinata"}' $wkdir/22_stacks_populations/combine_fis/08_fis_aut_and_x_no_missing_data > $wkdir/22_stacks_populations/combine_fis/10_fis_aut_and_x_no_missing_data_inopinata

#combine elegans and inopinata Fis data
cat  $wkdir/22_stacks_populations/combine_fis/10_fis_aut_and_x_no_missing_data_inopinata $wkdir/elegans/12_stacks/site_fis/08_elegans_no_header >  $wkdir/22_stacks_populations/combine_fis/10_b_elegans_and_inopinata_fis_no_header


echo -e "Chr\tBP\tFis\tspecies" | cat -  $wkdir/22_stacks_populations/combine_fis/10_b_elegans_and_inopinata_fis_no_header >  $wkdir/22_stacks_populations/combine_fis/fis_elegans_and_inopinata.tsv



#get number of loci , length of loci


mkdir $wkdir/23_bedtools_number_of_loci

mkdir $wkdir/23_bedtools_number_of_loci/autosomes/

mkdir $wkdir/23_bedtools_number_of_loci/X/

#autosomes
#remove header

grep -v "#" $wkdir/20_split_autosomes_X_sex/00_autosomes/inopinata_24_autosomes.vcf > $wkdir/23_bedtools_number_of_loci/autosomes/00_inopinata_24_sorted_no_header &



#make bed
cd $wkdir/23_bedtools_number_of_loci/autosomes

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$2,$4,$5}' 00_inopinata_24_sorted_no_header > 01_inopinata_24_sorted_bed &


#bedtools merge to get contiguous loci

bedtools merge -i 01_inopinata_24_sorted_bed > 02_inopinata_24_sorted_merge_bed &




#X
#remove header

grep -v "#" $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf > $wkdir/23_bedtools_number_of_loci/X/00_inopinata_24_sorted_no_header &



#make bed
cd $wkdir/23_bedtools_number_of_loci/X

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$2,$4,$5}' 00_inopinata_24_sorted_no_header > 01_inopinata_24_sorted_bed &


#bedtools merge to get contiguous loci

bedtools merge -i 01_inopinata_24_sorted_bed > 02_inopinata_24_sorted_merge_bed &

#combine

cat $wkdir/23_bedtools_number_of_loci/autosomes/02_inopinata_24_sorted_merge_bed $wkdir/23_bedtools_number_of_loci/X/02_inopinata_24_sorted_merge_bed > $wkdir/23_bedtools_number_of_loci/02_inopinata_24_sorted_merge_bed

#get number of loci

cd $wkdir/23_bedtools_number_of_loci/

wc -l 02_inopinata_24_sorted_merge_bed

# 77737


#get locus lengths

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,$3-$2+1}' 02_inopinata_24_sorted_merge_bed > 03_inopinata_loci_lengths

echo -e "Chr\tbp_start\tbp_end\tlocus_length" | cat -  03_inopinata_loci_lengths >  inopinata_loci_lengths.tsv


#number of sites, number of snps


#number of sites

grep -v "#" $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf | wc -l 

	#4974872 , 4,974,872 sites

#number of biallelic snps

wc -l $wkdir/19_cat_biallelic_invariant/01_inopinata_24_biallelic_snps_no_header

	#226527 , 226,527 snps



#prep file for splitstree
#get 20,000 random snps for splitstree... 

mkdir $wkdir/25_splitstree


shuf -n 20000 $wkdir/19_cat_biallelic_invariant/01_inopinata_24_biallelic_snps_no_header > $wkdir/25_splitstree/inopinata_24_biallelic_snps_20k_random_snps_no_header

cat $wkdir/19_cat_biallelic_invariant/00_header $wkdir/25_splitstree/inopinata_24_biallelic_snps_20k_random_snps_no_header > $wkdir/25_splitstree/inopinata_24_biallelic_snps_20k_random_snps.vcf

#module load easybuild  GCC/6.3.0-2.27  OpenMPI/2.0.2 Python/3.7.0

cd $wkdir/25_splitstree/

#vcf2phylip.py by Edgardo M. Ortiz
#available at https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py

#vcf to nexus
python $wkdir/additional_files/vcf2phylip.py --input inopinata_24_biallelic_snps_20k_random_snps.vcf -n -f &

#For splitstree figure, file "inopinata_24_biallelic_snps_20k_random_snps.min4.nexus" was loaded in SplitsTree4 (version 4.16.1, built 19 May 2020). Had to edit .nexus file-- 19507 characters (not 20000) were incorporated into the nexus. Replace "NCHAR=20000" with "NCHAR=19507" 




