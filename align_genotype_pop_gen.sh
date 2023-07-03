#alignment, genotyping, and pop gen stats workflow for Woodruff et al. 2021 
#C. inopinata population genomics paper, "Patterns of genomic diversity in a fig-associated 
#close relative of C. elegans." Formerly, "Alignment of genetic differentiation across 
#trophic levels in a fig community"

#if you have any questions, please contact Gavin Woodruff at gcwoodruff@ou.edu

#put working directory here
wkdir="/scratch/gcwoodruff/pop_gen_4-2022/"

#demultiplexed, re-named fastq files are in folder $wkdir/00_reads/

#get number of reads per sample

cd $wkdir/00_reads/
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

#get C. inopinata genome files from WormBase ParaSite

mkdir $wkdir/inopinata_genome
cd $wkdir/inopinata_genome

#unmasked fasta

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_inopinata/PRJDB5687/caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa.gz

#annotation
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_inopinata/PRJDB5687/caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations.gff3.gz

#protein
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_inopinata/PRJDB5687/caenorhabditis_inopinata.PRJDB5687.WBPS16.protein.fa.gz

#full length transcript
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_inopinata/PRJDB5687/caenorhabditis_inopinata.PRJDB5687.WBPS16.mRNA_transcripts.fa.gz

#cds transcript
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_inopinata/PRJDB5687/caenorhabditis_inopinata.PRJDB5687.WBPS16.CDS_transcripts.fa.gz

gunzip *



#make gsnap genome db

mkdir $wkdir/gsnap_genome_db


gmap_build --dir=$wkdir/gsnap_genome_db -d inopinata_genome_unmasked $wkdir/inopinata_genome/caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa


#align reads to reference


mkdir $wkdir/06_gsnap/

gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/A03.1.fq $wkdir/00_reads/A03.2.fq  > $wkdir/06_gsnap/A03.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/A05.1.fq $wkdir/00_reads/A05.2.fq  > $wkdir/06_gsnap/A05.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/B03.1.fq $wkdir/00_reads/B03.2.fq  > $wkdir/06_gsnap/B03.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/dec_2016_D01.1.fq $wkdir/00_reads/dec_2016_D01.2.fq  > $wkdir/06_gsnap/dec_2016_D01.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/dec_2016_D02.1.fq $wkdir/00_reads/dec_2016_D02.2.fq  > $wkdir/06_gsnap/dec_2016_D02.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/dec_2016_D03.1.fq $wkdir/00_reads/dec_2016_D03.2.fq  > $wkdir/06_gsnap/dec_2016_D03.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/dec_2016_D04.1.fq $wkdir/00_reads/dec_2016_D04.2.fq  > $wkdir/06_gsnap/dec_2016_D04.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/dec_2016_D05.1.fq $wkdir/00_reads/dec_2016_D05.2.fq  > $wkdir/06_gsnap/dec_2016_D05.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/dec_2016_D10.1.fq $wkdir/00_reads/dec_2016_D10.2.fq  > $wkdir/06_gsnap/dec_2016_D10.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_A11.1.fq $wkdir/00_reads/may_2017_A11.2.fq  > $wkdir/06_gsnap/may_2017_A11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_A12.1.fq $wkdir/00_reads/may_2017_A12.2.fq  > $wkdir/06_gsnap/may_2017_A12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_B11.1.fq $wkdir/00_reads/may_2017_B11.2.fq  > $wkdir/06_gsnap/may_2017_B11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_B12.1.fq $wkdir/00_reads/may_2017_B12.2.fq  > $wkdir/06_gsnap/may_2017_B12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_C11.1.fq $wkdir/00_reads/may_2017_C11.2.fq  > $wkdir/06_gsnap/may_2017_C11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_C12.1.fq $wkdir/00_reads/may_2017_C12.2.fq  > $wkdir/06_gsnap/may_2017_C12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_D11.1.fq $wkdir/00_reads/may_2017_D11.2.fq  > $wkdir/06_gsnap/may_2017_D11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_D12.1.fq $wkdir/00_reads/may_2017_D12.2.fq  > $wkdir/06_gsnap/may_2017_D12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_E11.1.fq $wkdir/00_reads/may_2017_E11.2.fq  > $wkdir/06_gsnap/may_2017_E11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_E12.1.fq $wkdir/00_reads/may_2017_E12.2.fq  > $wkdir/06_gsnap/may_2017_E12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_F11.1.fq $wkdir/00_reads/may_2017_F11.2.fq  > $wkdir/06_gsnap/may_2017_F11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_G12.1.fq $wkdir/00_reads/may_2017_G12.2.fq  > $wkdir/06_gsnap/may_2017_G12.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_H10.1.fq $wkdir/00_reads/may_2017_H10.2.fq  > $wkdir/06_gsnap/may_2017_H10.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_H11.1.fq $wkdir/00_reads/may_2017_H11.2.fq  > $wkdir/06_gsnap/may_2017_H11.sam
gsnap -d inopinata_genome_unmasked -D $wkdir/gsnap_genome_db/ --trim-mismatch-score=0 --trim-indel-score=0 --format=sam $wkdir/00_reads/may_2017_H12.1.fq $wkdir/00_reads/may_2017_H12.2.fq  > $wkdir/06_gsnap/may_2017_H12.sam

	#the above commands were run in parallel via slurm scripts


#get unique alignments

	#gsnap provides flags for that.
	#"XO:Z:CU" for concordant unique
	#"XO:Z:HU" , half mapping unique
	#"XO:Z:PI" , paired unique inversion
	#"XO:Z:PS" , paired unique scramble
	#"XO:Z:PL" , paired unique long
	#"XO:Z:UU" , unpaired unique
		#putting these values in a file for grep , file $wkdir/additional_files/gsnap_sam_flags 

mkdir $wkdir/07_grep_uniq_alignments
cd $wkdir/07_grep_uniq_alignments


#get the unique reads

cd $wkdir/06_gsnap

for i in *; do LC_ALL=C fgrep -w -f $wkdir/additional_files/gsnap_sam_flags $i > $wkdir/07_grep_uniq_alignments/uniq_alignments/$i; done

#get the headers

cd $wkdir/06_gsnap

for i in *; do LC_ALL=C fgrep "@" $i > $wkdir/07_grep_uniq_alignments/headers/$i; done

#concatenate headers and unique reads

cd $wkdir/07_grep_uniq_alignments/headers/

for i in *; do cat $i $wkdir/07_grep_uniq_alignments/uniq_alignments/$i > $wkdir/07_grep_uniq_alignments/cat/$i; done

#get rid of unneeded files

cd $wkdir/07_grep_uniq_alignments/

rm -r uniq_alignments
rm -r headers


#need samtools next ; OSCER module below
	#SAMtools/1.11-GCC-10.2.0
	#module load SAMtools/1.11-GCC-10.2.0

#sort sam files



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
	#the above commands were run in parallel via slurm scripts

#unmapped reads to get fraction reads aligned

mkdir $wkdir/NNN_unmapped_reads

cd $wkdir/06_gsnap

for i in *; do samtools view -f 4 $i -c > $wkdir/NNN_unmapped_reads/$i; done

cd $wkdir/NNN_unmapped_reads/

cat *

#1576949
#1129884
#148440
#5334168
#4699662
#7468529
#6429788
#5475575
#6443918
#1409635
#405200
#5424368
#357865
#5135714
#3457112
#1427160
#2494577
#4396817
#3483095
#6695332
#2409330
#338267
#8289285
#3141715
	#these were put in the supplemental table and used to determine the fraction of reads aligned

#call genotypes

#ploidy file in $wkdir/additional_files/ploidy_file.txt
#sample files in $wkdir/additional_files/sample_files/

#need bcftools
#module load BCFtools/1.11-GCC-10.2.0
mkdir $wkdir/10_bcftools_mpileup_call

cd $wkdir/08_samtools_sort/


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

	#the above commands were run in parallel via slurm scripts


#get only sites with >14x coverage

mkdir $wkdir/11_bcftools_view_depth_15


cd $wkdir/10_bcftools_mpileup_call/
bcftools view -i 'DP>=15' A03.bcf > $wkdir/11_bcftools_view_depth_15/A03.vcf
bcftools view -i 'DP>=15' A05.bcf > $wkdir/11_bcftools_view_depth_15/A05.vcf
bcftools view -i 'DP>=15' B03.bcf > $wkdir/11_bcftools_view_depth_15/B03.vcf
bcftools view -i 'DP>=15' dec_2016_D01.bcf > $wkdir/11_bcftools_view_depth_15/dec_2016_D01.vcf
bcftools view -i 'DP>=15' dec_2016_D02.bcf > $wkdir/11_bcftools_view_depth_15/dec_2016_D02.vcf
bcftools view -i 'DP>=15' dec_2016_D03.bcf > $wkdir/11_bcftools_view_depth_15/dec_2016_D03.vcf
bcftools view -i 'DP>=15' dec_2016_D04.bcf > $wkdir/11_bcftools_view_depth_15/dec_2016_D04.vcf
bcftools view -i 'DP>=15' dec_2016_D05.bcf > $wkdir/11_bcftools_view_depth_15/dec_2016_D05.vcf
bcftools view -i 'DP>=15' dec_2016_D10.bcf > $wkdir/11_bcftools_view_depth_15/dec_2016_D10.vcf
bcftools view -i 'DP>=15' may_2017_A11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_A11.vcf
bcftools view -i 'DP>=15' may_2017_A12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_A12.vcf
bcftools view -i 'DP>=15' may_2017_B11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_B11.vcf
bcftools view -i 'DP>=15' may_2017_B12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_B12.vcf
bcftools view -i 'DP>=15' may_2017_C11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_C11.vcf
bcftools view -i 'DP>=15' may_2017_C12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_C12.vcf
bcftools view -i 'DP>=15' may_2017_D11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_D11.vcf
bcftools view -i 'DP>=15' may_2017_D12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_D12.vcf
bcftools view -i 'DP>=15' may_2017_E11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_E11.vcf
bcftools view -i 'DP>=15' may_2017_E12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_E12.vcf
bcftools view -i 'DP>=15' may_2017_F11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_F11.vcf
bcftools view -i 'DP>=15' may_2017_G12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_G12.vcf
bcftools view -i 'DP>=15' may_2017_H10.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_H10.vcf
bcftools view -i 'DP>=15' may_2017_H11.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_H11.vcf
bcftools view -i 'DP>=15' may_2017_H12.bcf > $wkdir/11_bcftools_view_depth_15/may_2017_H12.vcf



#convert back to bcf and index

mkdir $wkdir/12_bcftools_view_bcf/

cd $wkdir/11_bcftools_view_depth_15


cd $wkdir/11_bcftools_view_depth_15

bcftools view A03.vcf -O b > $wkdir/12_bcftools_view_bcf/A03.bcf
bcftools view A05.vcf -O b > $wkdir/12_bcftools_view_bcf/A05.bcf
bcftools view B03.vcf -O b > $wkdir/12_bcftools_view_bcf/B03.bcf
bcftools view dec_2016_D01.vcf -O b > $wkdir/12_bcftools_view_bcf/dec_2016_D01.bcf
bcftools view dec_2016_D02.vcf -O b > $wkdir/12_bcftools_view_bcf/dec_2016_D02.bcf
bcftools view dec_2016_D03.vcf -O b > $wkdir/12_bcftools_view_bcf/dec_2016_D03.bcf
bcftools view dec_2016_D04.vcf -O b > $wkdir/12_bcftools_view_bcf/dec_2016_D04.bcf
bcftools view dec_2016_D05.vcf -O b > $wkdir/12_bcftools_view_bcf/dec_2016_D05.bcf
bcftools view dec_2016_D10.vcf -O b > $wkdir/12_bcftools_view_bcf/dec_2016_D10.bcf
bcftools view may_2017_A11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_A11.bcf
bcftools view may_2017_A12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_A12.bcf
bcftools view may_2017_B11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_B11.bcf
bcftools view may_2017_B12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_B12.bcf
bcftools view may_2017_C11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_C11.bcf
bcftools view may_2017_C12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_C12.bcf
bcftools view may_2017_D11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_D11.bcf
bcftools view may_2017_D12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_D12.bcf
bcftools view may_2017_E11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_E11.bcf
bcftools view may_2017_E12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_E12.bcf
bcftools view may_2017_F11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_F11.bcf
bcftools view may_2017_G12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_G12.bcf
bcftools view may_2017_H10.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_H10.bcf
bcftools view may_2017_H11.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_H11.bcf
bcftools view may_2017_H12.vcf -O b > $wkdir/12_bcftools_view_bcf/may_2017_H12.bcf

cd $wkdir/12_bcftools_view_bcf/

bcftools index A03.bcf
bcftools index A05.bcf
bcftools index B03.bcf
bcftools index dec_2016_D01.bcf
bcftools index dec_2016_D02.bcf
bcftools index dec_2016_D03.bcf
bcftools index dec_2016_D04.bcf
bcftools index dec_2016_D05.bcf
bcftools index dec_2016_D10.bcf
bcftools index may_2017_A11.bcf
bcftools index may_2017_A12.bcf
bcftools index may_2017_B11.bcf
bcftools index may_2017_B12.bcf
bcftools index may_2017_C11.bcf
bcftools index may_2017_C12.bcf
bcftools index may_2017_D11.bcf
bcftools index may_2017_D12.bcf
bcftools index may_2017_E11.bcf
bcftools index may_2017_E12.bcf
bcftools index may_2017_F11.bcf
bcftools index may_2017_G12.bcf
bcftools index may_2017_H10.bcf
bcftools index may_2017_H11.bcf
bcftools index may_2017_H12.bcf


#merge all bcf


mkdir $wkdir/13_bcftools_merge/

cd $wkdir/12_bcftools_view_bcf/

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
grep -v "#" inopinata_24.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 != "Sp34_ChrX"' > $wkdir/14_split_autosomes_X/autosomes/01_autosomes_no_header

#combine and sort vcf
cd $wkdir/14_split_autosomes_X/autosomes/

cat 00_header 01_autosomes_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/14_split_autosomes_X/autosomes/tmp"}' > inopinata_24_autosomes.vcf

#X chromosome

cd $wkdir/13_bcftools_merge/

#put aside header
grep "#" inopinata_24.vcf > $wkdir/14_split_autosomes_X/X/00_header


#get X
grep -v "#" inopinata_24.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 == "Sp34_ChrX"' > $wkdir/14_split_autosomes_X/X/01_X_no_header

#combine and sort vcf
cd $wkdir/14_split_autosomes_X/X/

cat 00_header 01_X_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > inopinata_24_X.vcf 


#just get sites with ~ >80% of samples having genotype calls

mkdir $wkdir/15_bcftools_view_gt_count/

cd $wkdir/14_split_autosomes_X/autosomes/

	#5/24 = 0.2083333
bcftools view -i 'COUNT(GT="mis")<5' inopinata_24_autosomes.vcf > $wkdir/15_bcftools_view_gt_count/inopinata_24_autosomes.vcf

cd $wkdir/14_split_autosomes_X/X/

	#there are four animals with unknown sex, and only using females for X stats, so the denominator is smaller for the X; 4/18 = 0.22
bcftools view -i 'COUNT(GT="mis")<4' inopinata_24_X.vcf > $wkdir/15_bcftools_view_gt_count/inopinata_24_X.vcf





#combine autosomes and X

mkdir $wkdir/16_cat/

cd $wkdir/15_bcftools_view_gt_count/
#set aside header
grep "#" inopinata_24_autosomes.vcf > $wkdir/16_cat/00_header 
#autosomes
grep -v "#" inopinata_24_autosomes.vcf >  $wkdir/16_cat/01_inopinata_24_autosomes_no_header 
#x
grep -v "#" inopinata_24_X.vcf > $wkdir/16_cat/02_inopinata_24_X_no_header 

cd $wkdir/16_cat/
#combine vcf
cat 00_header 01_inopinata_24_autosomes_no_header 02_inopinata_24_X_no_header > 03_combined_unsorted.vcf 

#sort

mkdir $wkdir/16_cat/tmp/

cat 03_combined_unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/16_cat/tmp"}' > inopinata_24_sorted.vcf 

#next, invariant sites and variant sites to get biallelic snps

#get just the invariant sites....

mkdir $wkdir/17_invariant_sites/

cd $wkdir/16_cat/

grep "#" inopinata_24_sorted.vcf > $wkdir/17_invariant_sites/header 

grep -v "#" inopinata_24_sorted.vcf > $wkdir/17_invariant_sites/inopinata_24_no_header 

cd $wkdir/17_invariant_sites/

awk 'BEGIN {FS="\t"} {OFS="\t"} $5 == "."' inopinata_24_no_header > inopinata_24_no_alt 

cat header inopinata_24_no_alt > inopinata_24_invariant_sites.vcf

#get biallelic snps with at least one copy per allele

mkdir $wkdir/18_bcftools_view_biallelic_snps_maf/

cd $wkdir/16_cat/

bcftools view -m2 -M2 -v snps --min-ac 2:minor inopinata_24_sorted.vcf > $wkdir/18_bcftools_view_biallelic_snps_maf/inopinata_24.vcf 
	#this file is also called inopinata_24_biallelic_snps.vcf and used for pca, clustering


#combine

mkdir $wkdir/19_cat_biallelic_invariant

cd $wkdir/18_bcftools_view_biallelic_snps_maf/

grep "#" inopinata_24.vcf  > $wkdir/19_cat_biallelic_invariant/00_header

grep -v "#" inopinata_24.vcf  >  $wkdir/19_cat_biallelic_invariant/01_inopinata_24_biallelic_snps_no_header

cd $wkdir/17_invariant_sites/

grep -v "#" inopinata_24_invariant_sites.vcf > $wkdir/19_cat_biallelic_invariant/02_inopinata_24_invariant_sites_no_header 

cd $wkdir/19_cat_biallelic_invariant/

cat 00_header 01_inopinata_24_biallelic_snps_no_header 02_inopinata_24_invariant_sites_no_header > 03_combined_unsorted.vcf 
#sort

mkdir $wkdir/19_cat_biallelic_invariant/tmp

cat 03_combined_unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/19_cat_biallelic_invariant/tmp"}' > inopinata_24_sorted.vcf 


#ok, split autosomes, male x and female x for pop gen stats. 

mkdir $wkdir/20_split_autosomes_X_sex

cd $wkdir/20_split_autosomes_X_sex

mkdir 00_autosomes

mkdir $wkdir/20_split_autosomes_X_sex/00_autosomes/tmp


cd $wkdir/19_cat_biallelic_invariant/

#set aside header
grep "#" inopinata_24_sorted.vcf > $wkdir/20_split_autosomes_X_sex/00_autosomes/00_header

#autosomes
grep -v "#" inopinata_24_sorted.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 != "Sp34_ChrX"' > $wkdir/20_split_autosomes_X_sex/00_autosomes/01_autosomes_no_header &

cd $wkdir/20_split_autosomes_X_sex/00_autosomes/
#sort
cat 00_header 01_autosomes_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/20_split_autosomes_X_sex/00_autosomes/tmp"}' > inopinata_24_autosomes.vcf


#cat 00_header 01_autosomes_no_header | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/20_split_autosomes_X_sex/00_autosomes/tmp"}' > inopinata_24_autosomes.vcf &


#separate males and females 

mkdir $wkdir/20_split_autosomes_X_sex/01_split_by_sex

cd $wkdir/19_cat_biallelic_invariant/

#males
bcftools view -s B03.bam,may_2017_H12.bam inopinata_24_sorted.vcf > $wkdir/20_split_autosomes_X_sex/01_split_by_sex/inopinata_males.vcf

#females
bcftools view -s A03.bam,dec_2016_D01.bam,dec_2016_D02.bam,dec_2016_D03.bam,dec_2016_D04.bam,dec_2016_D05.bam,may_2017_A11.bam,may_2017_A12.bam,may_2017_B11.bam,may_2017_C11.bam,may_2017_C12.bam,may_2017_D11.bam,may_2017_D12.bam,may_2017_E11.bam,may_2017_E12.bam,may_2017_F11.bam,may_2017_H10.bam,may_2017_H11.bam inopinata_24_sorted.vcf > $wkdir/20_split_autosomes_X_sex/01_split_by_sex/inopinata_females.vcf





#extract X for males and females
mkdir $wkdir/20_split_autosomes_X_sex/02_split_X

cd $wkdir/20_split_autosomes_X_sex/01_split_by_sex/

for i in *; do grep "#" $i > $wkdir/20_split_autosomes_X_sex/02_split_X/$i.header; done 

for i in *; do grep -v "#" $i | awk 'BEGIN {FS="\t"} {OFS="\t"} $1 == "Sp34_ChrX"' > $wkdir/20_split_autosomes_X_sex/02_split_X/${i%.vcf}.X.vcf_no_header; done

cd $wkdir/20_split_autosomes_X_sex/02_split_X/

cat inopinata_females.vcf.header inopinata_females.X.vcf_no_header  | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/20_split_autosomes_X_sex/02_split_X/tmp"}' > inopinata_females_X.vcf 

cat inopinata_males.vcf.header inopinata_males.X.vcf_no_header  | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/20_split_autosomes_X_sex/02_split_X/tmp"}' > inopinata_males_X.vcf 



mkdir $wkdir/21_popgenWindows/
mkdir $wkdir/21_popgenWindows/00_parseVCF/
mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/


#used Biopython v 1.70 for the 2021 version
#now using Biopython/1.78-foss-2020b
#module load Biopython/1.78-foss-2020b


cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/20_split_autosomes_X_sex/00_autosomes/inopinata_24_autosomes.vcf | gzip > $wkdir/21_popgenWindows/00_parseVCF/inopinata_24_autosomes.geno.gz 
	#this appeared to work

cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/21_popgenWindows/00_parseVCF/inopinata_24_autosomes.geno.gz -o  $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_autosomes_stats.csv.gz -f phased  

cd $wkdir/21_popgenWindows/01_popgenWindows_py/

gunzip inopinata_autosomes_stats.csv.gz


#now, inopinata X -- only using inopinata females for X pop gen stats 

cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf | gzip > $wkdir/21_popgenWindows/00_parseVCF/inopinata_females_X.geno.gz 

#now estimate stats

cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g /$wkdir/21_popgenWindows/00_parseVCF/inopinata_females_X.geno.gz -o  $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_females_X.csv.gz -f phased 


#combine the files
cd $wkdir/21_popgenWindows/01_popgenWindows_py/

#remove header
sed '1d' inopinata_females_X.csv > inopinata_females_X_no_header

#combine autosomes and X pi estimates
cat inopinata_autosomes_stats.csv inopinata_females_X_no_header > inopinata_autosomes_females_X_pi.csv


#inopinata FST

mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst


	#fst autosomes


cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 --analysis popPairDist -g $wkdir/21_popgenWindows/00_parseVCF/inopinata_24_autosomes.geno.gz -o $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/inopinata_autosomes_fst_pops_file.csv.gz -f phased -p ishigaki -p iriomote -p yonaguni --popsFile $wkdir/additional_files/inopinata_pops_file_islands &

	#fst female X

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 --analysis popPairDist -g $wkdir/21_popgenWindows/00_parseVCF/inopinata_females_X.geno.gz -o $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/inopinata_females_X_fst.csv.gz -f phased -p ishigaki -p iriomote -p yonaguni --popsFile $wkdir/additional_files/inopinata_pops_file_islands_females &


cd $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst

#prepare FST file for stats and figures
gunzip inopinata_autosomes_fst_pops_file.csv.gz
gunzip inopinata_females_X_fst.csv.gz

mkdir $wkdir/21_popgenWindows/01_popgenWindows_py/inopinata_fst/organize_fst_file/


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

#now using Stacks/2.59-foss-2020a

module load Stacks/2.59-foss-2020a



#FIS with stacks

# inopinata_autosomes_pop_map_stacks.txt in folder /additional_files/


mkdir $wkdir/22_stacks_populations/
mkdir $wkdir/22_stacks_populations/autosomes/


populations -V $wkdir/20_split_autosomes_X_sex/00_autosomes/inopinata_24_autosomes.vcf -O $wkdir/22_stacks_populations/autosomes/ -M $wkdir/additional_files/inopinata_autosomes_pop_map_stacks.txt  --sigma 3333 --genepop --structure --phylip


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
#now using BEDTools/2.27.1-foss-2018b
module load BEDTools/2.27.1-foss-2018b



cd $wkdir/inopinata_genome

samtools faidx caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2} ' caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa.fai | sort > inopinata_chr_sizes.txt


cd $wkdir/additional_files


bedtools makewindows -g $wkdir/inopinata_genome/inopinata_chr_sizes.txt -w 50 > inopinata.50bp.windows
bedtools makewindows -g $wkdir/inopinata_genome/inopinata_chr_sizes.txt -w 10000 > inopinata.10kb.windows


cd $wkdir/22_stacks_populations/autosomes/site_fis

bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.50bp.windows -b 01_fis_cols > 03_fis_50bp_windows



#remove missing data
awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 03_fis_50bp_windows > 04_fis_50bp_windows_no_missing_data

#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.10kb.windows -b 04_fis_50bp_windows_no_missing_data > 05_fis_10kb_windows 


#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 05_fis_10kb_windows > 06_fis_10kb_windows_no_missing_data


mkdir $wkdir/22_stacks_populations/X/


populations -V $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf -O $wkdir/22_stacks_populations/X/ -M $wkdir/additional_files/inopinata_X_females_pop_map_stacks.txt  --sigma 3333 --genepop --structure --phylip 



#prep FIS file for stats and figures
mkdir $wkdir/22_stacks_populations/X/site_fis

cd $wkdir/22_stacks_populations/X/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$3,$17}' inopinata_females_X.p.sumstats.tsv > $wkdir/22_stacks_populations/X/site_fis/01_fis_cols

cd $wkdir/22_stacks_populations/X/site_fis


#remove first two lines

sed -i '1d' 01_fis_cols

sed -i '1d' 01_fis_cols


#repair bed file (ie, replace spaces with tabs)
sed -i -e 's/ /\t/g' 01_fis_cols


#bedtools to get windowed Fis!


cd $wkdir/additional_files


grep Sp34_ChrX inopinata.50bp.windows > inopinata.50bp.windows.X
grep Sp34_ChrX inopinata.10kb.windows > inopinata.10kb.windows.X

cd $wkdir/22_stacks_populations/X/site_fis


bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.50bp.windows.X -b 01_fis_cols > 03_fis_50bp_windows 


#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 03_fis_50bp_windows > 04_fis_50bp_windows_no_missing_data


#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.10kb.windows.X -b 04_fis_50bp_windows_no_missing_data > 05_fis_10kb_windows 

#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 05_fis_10kb_windows > 06_fis_10kb_windows_no_missing_data


#combine autosomes and X chromosome

mkdir $wkdir/22_stacks_populations/combine_fis/

cat $wkdir/22_stacks_populations/autosomes/site_fis/06_fis_10kb_windows_no_missing_data $wkdir/22_stacks_populations/X/site_fis/06_fis_10kb_windows_no_missing_data > $wkdir/22_stacks_populations/combine_fis/08_fis_aut_and_x_no_missing_data


awk 'BEGIN {OFS="\t"} {print $0, "inopinata"}' $wkdir/22_stacks_populations/combine_fis/08_fis_aut_and_x_no_missing_data > $wkdir/22_stacks_populations/combine_fis/10_fis_aut_and_x_no_missing_data_inopinata


###new stuff for 2023 revision (intron/exon pi; genic/intergenic pi; codon position pi)

#using the genomics_general tools (https://github.com/simonhmartin/genomics_general/)

mkdir $wkdir/26_2022_site_annotation_genomics_general/


cd $wkdir/20_split_autosomes_X_sex/00_autosomes/

#use tabix to generate index needed for codingSiteTypes.py
#module load tabix/0.2.6-GCCcore-7.3.0

#autosomes
bgzip -c inopinata_24_autosomes.vcf > inopinata_24_autosomes.vcf.gz
tabix -p vcf inopinata_24_autosomes.vcf.gz

#X, females only
cd $wkdir/20_split_autosomes_X_sex/02_split_X/

bgzip -c inopinata_females_X.vcf > inopinata_females_X.vcf.gz
tabix -p vcf inopinata_females_X.vcf.gz

#cd to directory with scripts

cd /home/gcwoodruff/download/genomics_general/

#annotate coding site types

#autosomes
python3 codingSiteTypes.py -a $wkdir/inopinata_genome/caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations.gff3 -f gff3 -r $wkdir/inopinata_genome/caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa -v $wkdir/20_split_autosomes_X_sex/00_autosomes/inopinata_24_autosomes.vcf.gz -o $wkdir/26_2022_site_annotation_genomics_general/inopinata_24_autosomes_codingSiteTypes.out --ignoreConflicts

#X
python3 codingSiteTypes.py -a $wkdir/inopinata_genome/caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations.gff3 -f gff3 -r $wkdir/inopinata_genome/caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa -v $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf.gz -o $wkdir/26_2022_site_annotation_genomics_general/inopinata_females_X_codingSiteTypes.out --ignoreConflicts

#get sites by codon annotation

mkdir $wkdir/27_2022_site_annotation_continued/
mkdir $wkdir/27_2022_site_annotation_continued/00_bedtools_intersect

cd $wkdir/26_2022_site_annotation_genomics_general/

grep "Sp34_ChrX" inopinata_females_X_codingSiteTypes.out > inopinata_females_X_codingSiteTypes.out_grep_X
grep -v "Sp34_ChrX" inopinata_24_autosomes_codingSiteTypes.out > inopinata_24_autosomes_codingSiteTypes.out_grep_autosomes

#get relevant columns
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$2,$3,$4,$5} ' inopinata_24_autosomes_codingSiteTypes.out_grep_autosomes | sed '1d' > inopinata_24_autosomes_codingSiteTypes.out_grep_autosomes_bed

#combine autosomes and X
cat inopinata_24_autosomes_codingSiteTypes.out_grep_autosomes inopinata_females_X_codingSiteTypes.out_grep_X > inopinata_codon_site_annotations.tsv


#join coding site annotations with VCF

mkdir $wkdir/27_2022_site_annotation_continued//01_join

#sort autosomes VCF for join
sort $wkdir/20_split_autosomes_X_sex/00_autosomes/01_autosomes_no_header > $wkdir/27_2022_site_annotation_continued//01_join/00_autosomes_no_header_sort

#remove header of coding site type annotation file
sed '1d' $wkdir/26_2022_site_annotation_genomics_general/inopinata_24_autosomes_codingSiteTypes.out > $wkdir/27_2022_site_annotation_continued/00_bedtools_intersect/01_join/01_inopinata_24_autosomes_codingSiteTypes.out_no_header

mv $wkdir/27_2022_site_annotation_continued//01_join $wkdir/27_2022_site_annotation_continued/01_join

#make tmp folder for join
mkdir /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp

#join
join -j1 <(<00_autosomes_no_header_sort awk '{print $1"-"$2" "$0}' | sort -k1,1 -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/) <(<01_inopinata_24_autosomes_codingSiteTypes.out_no_header awk '{print $1"-"$2" "$0}' | sort -k1,1 -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/) > 02_join_vcf_annotations_autosomes 


#split annotated VCF by degeneracy and codon position

mkdir $wkdir/27_2022_site_annotation_continued/02_awk
mkdir $wkdir/27_2022_site_annotation_continued/02_awk/tmp
mkdir $wkdir/27_2022_site_annotation_continued/02_awk/00_degeneracy

#fold degeneracy
cd $wkdir/27_2022_site_annotation_continued/02_awk/00_degeneracy

awk '{print > $39}'  $wkdir/27_2022_site_annotation_continued/01_join/02_join_vcf_annotations_autosomes

#codon position
mkdir $wkdir/27_2022_site_annotation_continued/02_awk/01_codon_position/
cd $wkdir/27_2022_site_annotation_continued/02_awk/01_codon_position/

awk '{print > $37}'  $wkdir/27_2022_site_annotation_continued/01_join/02_join_vcf_annotations_autosomes


#X

#join coding site annotations with VCF
#sort vcf
sort $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females.X.vcf_no_header > $wkdir/27_2022_site_annotation_continued//01_join/inopinata_females.X_vcf_no_header_sort
#remove header of coding site types file
sed '1d' $wkdir/26_2022_site_annotation_genomics_general/inopinata_females_X_codingSiteTypes.out > $wkdir/27_2022_site_annotation_continued//01_join/inopinata_females_X_codingSiteTypes.out_no_header


#join
cd $wkdir/27_2022_site_annotation_continued//01_join/

join -j1 <(<inopinata_females.X_vcf_no_header_sort awk '{print $1"-"$2" "$0}' | sort -k1,1 -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/) <(<inopinata_females_X_codingSiteTypes.out_no_header awk '{print $1"-"$2" "$0}' | sort -k1,1 -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/) > join_vcf_annotations_X &


# $39 =  degeneracy
# $38 = synonymous/nonsynonymous
# $37 = codon position

#split by degeneracy and codon position

mkdir $wkdir/27_2022_site_annotation_continued/02_awk/X
mkdir $wkdir/27_2022_site_annotation_continued/02_awk/X/tmp
mkdir $wkdir/27_2022_site_annotation_continued/02_awk/X/00_degeneracy
#fold degeneracy
cd $wkdir/27_2022_site_annotation_continued/02_awk/X/00_degeneracy

awk '{print > $33}'  $wkdir/27_2022_site_annotation_continued/01_join/join_vcf_annotations_X
#codon position
mkdir $wkdir/27_2022_site_annotation_continued/02_awk/X/01_codon_position/
cd $wkdir/27_2022_site_annotation_continued/02_awk/X/01_codon_position/

awk '{print > $31}'  $wkdir/27_2022_site_annotation_continued/01_join/join_vcf_annotations_X


#preparing files for estimating pi by 10 kb windows

#autosomes
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/autosomes
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/autosomes/degneracy
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/autosomes/codon_position

#degeneracy, replace spaces with free tab, get columns right
cd $wkdir/27_2022_site_annotation_continued/02_awk/00_degeneracy/

for i in *; do sed -i -e 's/ \+/\t/g' $i; done

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34} ' $i > $wkdir/27_2022_site_annotation_continued/03_awk/autosomes/degneracy/$i; done 


cd $wkdir/27_2022_site_annotation_continued/02_awk/01_codon_position/


#codon position, replace spaces with free tab, get columns right
for i in *; do sed -i -e 's/ \+/\t/g' $i; done

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34} ' $i > $wkdir/27_2022_site_annotation_continued/03_awk/autosomes/codon_position/$i; done


#X


mkdir $wkdir/27_2022_site_annotation_continued/03_awk/
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/X
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/X/degneracy
mkdir $wkdir/27_2022_site_annotation_continued/03_awk/X/codon_position



#degeneracy, replace spaces with free tab, get columns right
cd $wkdir/27_2022_site_annotation_continued/02_awk/X/00_degeneracy/

for i in *; do sed -i -e 's/ \+/\t/g' $i; done

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28} ' $i > $wkdir/27_2022_site_annotation_continued/03_awk/X/degneracy/$i; done 


#autosomes, replace spaces with free tab, get columns right
cd $wkdir/27_2022_site_annotation_continued/02_awk/X/01_codon_position/


for i in *; do sed -i -e 's/ \+/\t/g' $i; done

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28} ' $i > $wkdir/27_2022_site_annotation_continued/03_awk/X/codon_position/$i; done

#sort


mkdir $wkdir/27_2022_site_annotation_continued/04_sort/
mkdir $wkdir/27_2022_site_annotation_continued/04_sort/autosomes
mkdir $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy
mkdir $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/codon_position

mkdir $wkdir/27_2022_site_annotation_continued/04_sort/X
mkdir $wkdir/27_2022_site_annotation_continued/04_sort/X/degneracy
mkdir $wkdir/27_2022_site_annotation_continued/04_sort/X/codon_position

#autosomes
#sort and add header
#degeneracy
cd $wkdir/27_2022_site_annotation_continued/03_awk/autosomes/degneracy/

for i in *; do cat $wkdir/20_split_autosomes_X_sex/00_autosomes/00_header $i | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/"}' > $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/$i; done 

#codon position
cd $wkdir/27_2022_site_annotation_continued/03_awk/autosomes/codon_position/

for i in *; do cat $wkdir/20_split_autosomes_X_sex/00_autosomes/00_header $i | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/"}' > $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/codon_position/$i; done 



#X

grep '#' $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf > $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X_header
#degeneracy
cd $wkdir/27_2022_site_annotation_continued/03_awk/X/degneracy/

for i in *; do cat $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X_header $i | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/"}' > $wkdir/27_2022_site_annotation_continued/04_sort/X/degneracy/$i; done 

#codon position
cd $wkdir/27_2022_site_annotation_continued/03_awk/X/codon_position/

for i in *; do cat $wkdir/20_split_autosomes_X_sex/02_split_X/inopinata_females_X_header $i | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen_4-2022/27_2022_site_annotation_continued/01_join/tmp/"}' > $wkdir/27_2022_site_annotation_continued/04_sort/X/codon_position/$i; done 


#rename files

mkdir $wkdir/27_2022_site_annotation_continued/05_popgenwindows
mkdir $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/
mkdir $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

cd $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/

mv 0 0_degen_ino_aut.vcf
mv 2 2_degen_ino_aut.vcf
mv 4 4_degen_ino_aut.vcf

cd $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/codon_position/

mv 1 First_codon_pos_ino_aut.vcf
mv 2 Second_codon_pos_ino_aut.vcf
mv 3 Third_codon_pos_ino_aut.vcf

cd $wkdir/27_2022_site_annotation_continued/04_sort/X/degneracy/

mv 0 0_degen_ino_X.vcf
mv 2 2_degen_ino_X.vcf
mv 4 4_degen_ino_X.vcf

cd $wkdir/27_2022_site_annotation_continued/04_sort/X/codon_position/

mv 1 First_codon_pos_ino_X.vcf
mv 2 Second_codon_pos_ino_X.vcf
mv 3 Third_codon_pos_ino_X.vcf

#estimate pi in 10 kb genomic windows

cd $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/


#autosomes fold-degeneracy pi in 10 kb genomic windows


cd $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/

#used Biopython v 1.70
#now using Biopython/1.78-foss-2020b
#module load Biopython/1.78-foss-2020b

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/0_degen_ino_aut.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/0_degen_ino_aut.geno.gz 

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/0_degen_ino_aut.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/0-fold_degen_ino_aut_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

gunzip 0-fold_degen_ino_aut_stats.csv.gz

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/2_degen_ino_aut.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/2_degen_ino_aut.geno.gz 

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/2_degen_ino_aut.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/2-fold_degen_ino_aut_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

gunzip 2-fold_degen_ino_aut_stats.csv.gz




#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/degneracy/4_degen_ino_aut.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/4_degen_ino_aut.geno.gz 

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/4_degen_ino_aut.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/4-fold_degen_ino_aut_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip 4-fold_degen_ino_aut_stats.csv.gz




#inopinata codon position pi 10 kb genomic windows


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/codon_position/First_codon_pos_ino_aut.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/First_codon_pos_ino_aut.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/First_codon_pos_ino_aut.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/First_codon_pos_ino_aut.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip First_codon_pos_ino_aut.csv.gz



#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/codon_position/Second_codon_pos_ino_aut.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Second_codon_pos_ino_aut.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Second_codon_pos_ino_aut.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/Second_codon_pos_ino_aut.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip Second_codon_pos_ino_aut.csv.gz




#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/autosomes/codon_position/Third_codon_pos_ino_aut.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Third_codon_pos_ino_aut.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Third_codon_pos_ino_aut.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/Third_codon_pos_ino_aut.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip Third_codon_pos_ino_aut.csv.gz




#now, inopinata X -- only using inopinata females for X pop gen stats 

#fold degeneracy X

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/X/degneracy/0_degen_ino_X.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/0_degen_ino_X.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/0_degen_ino_X.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/0-fold_degen_ino_X_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

gunzip 0-fold_degen_ino_X_stats.csv.gz


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/X/degneracy/2_degen_ino_X.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/2_degen_ino_X.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/2_degen_ino_X.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/2-fold_degen_ino_X_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

gunzip 2-fold_degen_ino_X_stats.csv.gz




#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/X/degneracy/4_degen_ino_X.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/4_degen_ino_X.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/4_degen_ino_X.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/4-fold_degen_ino_X_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip 4-fold_degen_ino_X_stats.csv.gz



#pi codon position X


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/X/codon_position/First_codon_pos_ino_X.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/First_codon_pos_ino_X.geno.gz 


#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/First_codon_pos_ino_X.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/First_codon_pos_ino_X.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip First_codon_pos_ino_X.csv.gz



#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/X/codon_position/Second_codon_pos_ino_X.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Second_codon_pos_ino_X.geno.gz 

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Second_codon_pos_ino_X.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/Second_codon_pos_ino_X.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip Second_codon_pos_ino_X.csv.gz




#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/27_2022_site_annotation_continued/04_sort/X/codon_position/Third_codon_pos_ino_X.vcf | gzip > $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Third_codon_pos_ino_X.geno.gz 
	#this appeared to work

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/27_2022_site_annotation_continued/05_popgenwindows/00_parseVCF/Third_codon_pos_ino_X.geno.gz  -o  $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/Third_codon_pos_ino_X.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip Third_codon_pos_ino_X.csv.gz




#combine pi estimates by codon position and fold degeneracy
mkdir $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/
mkdir $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk

#add column with site annotation
cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

awk 'BEGIN {FS=","} {OFS=","} {print $0,"0"} ' 0-fold_degen_ino_aut_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/0-fold_degen_ino_aut_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"0"} ' 0-fold_degen_ino_X_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/0-fold_degen_ino_X_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"2"} ' 2-fold_degen_ino_aut_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/2-fold_degen_ino_aut_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"2"} ' 2-fold_degen_ino_X_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/2-fold_degen_ino_X_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"4"} ' 4-fold_degen_ino_aut_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/4-fold_degen_ino_aut_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"4"} ' 4-fold_degen_ino_X_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/4-fold_degen_ino_X_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"1"} ' First_codon_pos_ino_aut.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/First_codon_pos_ino_aut.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"1"} ' First_codon_pos_ino_X.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/First_codon_pos_ino_X.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"2"} ' Second_codon_pos_ino_aut.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/Second_codon_pos_ino_aut.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"2"} ' Second_codon_pos_ino_X.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/Second_codon_pos_ino_X.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"3"} ' Third_codon_pos_ino_aut.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/Third_codon_pos_ino_aut.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"3"} ' Third_codon_pos_ino_X.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/Third_codon_pos_ino_X.csv

#remove headers
mkdir $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/01_sed/

cd $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/00_awk/

for i in *; do sed '1d' $i > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/01_sed/$i; done

#combine the fold degeneracy files together; the codon position files together
mkdir $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/02_cat/

cd $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/01_sed/

cat 0-fold_degen_ino_aut_stats.csv 0-fold_degen_ino_X_stats.csv 2-fold_degen_ino_aut_stats.csv 2-fold_degen_ino_X_stats.csv 4-fold_degen_ino_aut_stats.csv 4-fold_degen_ino_X_stats.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/02_cat/degeneracy_ino_pi.tmp.csv

cat First_codon_pos_ino_aut.csv First_codon_pos_ino_X.csv Second_codon_pos_ino_aut.csv Second_codon_pos_ino_X.csv Third_codon_pos_ino_aut.csv Third_codon_pos_ino_X.csv > $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/02_cat/codon_position_ino_pi.tmp.csv


#add headers
cd $wkdir/27_2022_site_annotation_continued/06_combine_popgenwindows/02_cat/

echo -e "chr,start,end,mid,sites,pi_all,fold_degeneracy" | cat - degeneracy_ino_pi.tmp.csv > degeneracy_ino_pi.csv


echo -e "chr,start,end,mid,sites,pi_all,codon_position" | cat - codon_position_ino_pi.tmp.csv > codon_position_ino_pi.csv

rm degeneracy_ino_pi.tmp.csv
rm codon_position_ino_pi.tmp.csv
#annotate with species id (for comparisons with elegans)
awk 'BEGIN {FS=","} {OFS=","} {print $0,"3","C. inopinata"} ' degeneracy_ino_pi.csv > degeneracy_ino_pi.csv.tmp
awk 'BEGIN {FS=","} {OFS=","} {print $0,"3","C. inopinata"} ' codon_position_ino_pi.csv > codon_position_ino_pi.csv.tmp


#intergenic/genic ; exon/intron pi

#get genes and exons


cd $wkdir/inopinata_genome

	#thanks https://www.biostars.org/p/112251/

awk 'OFS="\t" {print $1, "0", $2}' caenorhabditis_inopinata.PRJDB5687.WBPS16.genomic.fa.fai | sort -k1,1 -k2,2n > inopinata_chromSizes.bed

grep -v "#" caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations.gff3  > 00_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_no_pound.gff3

awk '$3 == "gene"' 00_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_no_pound.gff3 | sort -k1,1 -k4,4n > 01_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_gene.gff3

awk '$3 == "exon"' 00_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_no_pound.gff3 | sort -k1,1 -k4,4n > 02_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_exon.gff3


gff2bed < 01_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_gene.gff3 > 03__caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_gene.bed

gff2bed < 02_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_exon.gff3 > 04_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_exon.bed

#ok, get intergenic

bedtools complement -i 03__caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_gene.bed -g inopinata_chr_sizes.txt > 05_inopinata_intergenic_sorted.bed

#ok, get intronic

awk 'OFS="\t" {print $1, $2, $3, "intergenic"}' 05_inopinata_intergenic_sorted.bed > 06_inopinata_intergenic_sort.bed

cat 04_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_exon.bed 06_inopinata_intergenic_sort.bed | sort -k1,1 -k2,2n > 07_exon_intergenic_inopinata_sort.bed

sed -i -e 's/exon.*/exon/g' 07_exon_intergenic_inopinata_sort.bed

bedtools complement -i 07_exon_intergenic_inopinata_sort.bed -g inopinata_chr_sizes.txt > 08_inopinata_intron_sorted.bed

awk 'OFS="\t" {print $1, $2, $3, "intron"}' 08_inopinata_intron_sorted.bed > 09_inopinata_intron_sorted.bed

	#checked with igv, looks like it worked.


#get intersection of vcf and genomic regions

mkdir $wkdir/29_bedtools_intersect_genomic_region/ 

bedtools intersect -a $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf  -b 03__caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_gene.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_genic.vcf


bedtools intersect -a $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf  -b 06_inopinata_intergenic_sort.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intergenic.vcf


bedtools intersect -a $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf  -b 04_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_exon.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_exonic.vcf


bedtools intersect -a $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf  -b 09_inopinata_intron_sorted.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intronic.vcf


cd $wkdir/29_bedtools_intersect_genomic_region/
#add header

grep "#" $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf > $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted_header

cat $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_genic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_genic.vcf.tmp
cat $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intergenic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intergenic.vcf.tmp
cat $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_exonic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_exonic.vcf.tmp
cat $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intronic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intronic.vcf.tmp

mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_genic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_genic.vcf
mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intergenic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intergenic.vcf
mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_exonic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_exonic.vcf
mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intronic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intronic.vcf


#estimate pi of genomic regions


mkdir $wkdir/30_popgenwindows_genomic_regions/
mkdir $wkdir/30_popgenwindows_genomic_regions/00_parseVCF
mkdir $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py

module load Biopython/1.78-foss-2020b

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_genic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_genic.geno.gz 

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intergenic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_intergenic.geno.gz 

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_exonic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_exonic.geno.gz 

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_24_intronic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_intronic.geno.gz 

#cd to wherever your genomic_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_genic.geno.gz -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_24_genic_stats.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_intergenic.geno.gz  -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_24_intergenic_stats.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_exonic.geno.gz  -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_24_exonic.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_24_intronic.geno.gz  -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_24_intronic.csv.gz -f phased  

cd $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/

gunzip *

#turn into usable tsv file
mkdir $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv
cd $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py

for i in *; do awk 'BEGIN {FS=","} {OFS=","} {print $0,"inopinata"}' $i > $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv/$i; done

cd $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv

awk 'BEGIN {FS=","} {OFS=","} {print $0,"exonic"}' inopinata_24_exonic.csv > inopinata_24_exonic.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"genic"}' inopinata_24_genic_stats.csv > inopinata_24_genic_stats.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"intergenic"}' inopinata_24_intergenic_stats.csv > inopinata_24_intergenic_stats.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"intronic"}' inopinata_24_intronic.csv > inopinata_24_intronic.csv.tmp1


#combine with elegans estimates and make the tsv's for visualization
#see elegans_workflow.sh to see how these files were made
cd $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py/

for i in *; do awk 'BEGIN {FS=","} {OFS=","} {print $0,"elegans"}' $i > $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv/$i; done

cd $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv

awk 'BEGIN {FS=","} {OFS=","} {print $0,"exonic"}' elegans_24_exonic.csv > elegans_24_exonic.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"genic"}' elegans_24_genic_stats.csv > elegans_24_genic_stats.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"intergenic"}' elegans_24_intergenic_stats.csv > elegans_24_intergenic_stats.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"intronic"}' elegans_24_intronic.csv > elegans_24_intronic.csv.tmp1
#change chromosome id's for inopinata files so comparisons can be made with elegans
sed -i -e 's/Sp34_Chr1/I/g' inopinata_24_exonic.csv.tmp1
sed -i -e 's/Sp34_Chr1/I/g' inopinata_24_genic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr1/I/g' inopinata_24_intergenic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr1/I/g' inopinata_24_intronic.csv.tmp1

sed -i -e 's/Sp34_Chr2/II/g' inopinata_24_exonic.csv.tmp1
sed -i -e 's/Sp34_Chr2/II/g' inopinata_24_genic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr2/II/g' inopinata_24_intergenic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr2/II/g' inopinata_24_intronic.csv.tmp1

sed -i -e 's/Sp34_Chr3/III/g' inopinata_24_exonic.csv.tmp1
sed -i -e 's/Sp34_Chr3/III/g' inopinata_24_genic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr3/III/g' inopinata_24_intergenic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr3/III/g' inopinata_24_intronic.csv.tmp1

sed -i -e 's/Sp34_Chr4/IV/g' inopinata_24_exonic.csv.tmp1
sed -i -e 's/Sp34_Chr4/IV/g' inopinata_24_genic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr4/IV/g' inopinata_24_intergenic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr4/IV/g' inopinata_24_intronic.csv.tmp1

sed -i -e 's/Sp34_Chr5/V/g' inopinata_24_exonic.csv.tmp1
sed -i -e 's/Sp34_Chr5/V/g' inopinata_24_genic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr5/V/g' inopinata_24_intergenic_stats.csv.tmp1
sed -i -e 's/Sp34_Chr5/V/g' inopinata_24_intronic.csv.tmp1
#remove headers
sed -i '1d' elegans_24_exonic.csv.tmp1
sed -i '1d' elegans_24_genic_stats.csv.tmp1
sed -i '1d' elegans_24_intergenic_stats.csv.tmp1
sed -i '1d' elegans_24_intronic.csv.tmp1
sed -i '1d' inopinata_24_exonic.csv.tmp1
sed -i '1d' inopinata_24_genic_stats.csv.tmp1
sed -i '1d' inopinata_24_intergenic_stats.csv.tmp1
sed -i '1d' inopinata_24_intronic.csv.tmp1

#intergenic/genic ; exon/intron with the inopinata X

cd $wkdir/inopinata_genome

bedtools intersect -a $wkdir/20_split_autosomes_X_sex/02_split_X//inopinata_females_X.vcf  -b 03__caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_gene.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_genic.vcf


bedtools intersect -a $wkdir/20_split_autosomes_X_sex/02_split_X//inopinata_females_X.vcf  -b 06_inopinata_intergenic_sort.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intergenic.vcf


bedtools intersect -a $wkdir/20_split_autosomes_X_sex/02_split_X//inopinata_females_X.vcf  -b 04_caenorhabditis_inopinata.PRJDB5687.WBPS16.annotations_exon.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_exonic.vcf


bedtools intersect -a $wkdir/20_split_autosomes_X_sex/02_split_X//inopinata_females_X.vcf  -b 09_inopinata_intron_sorted.bed > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intronic.vcf


cd $wkdir/29_bedtools_intersect_genomic_region/

grep "#" $wkdir/20_split_autosomes_X_sex/02_split_X//inopinata_females_X.vcf > $wkdir/19_cat_biallelic_invariant/inopinata_X_header

cat $wkdir/19_cat_biallelic_invariant/inopinata_X_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_genic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_genic.vcf.tmp
cat $wkdir/19_cat_biallelic_invariant/inopinata_X_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intergenic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intergenic.vcf.tmp
cat $wkdir/19_cat_biallelic_invariant/inopinata_X_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_exonic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_exonic.vcf.tmp
cat $wkdir/19_cat_biallelic_invariant/inopinata_X_header $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intronic.vcf > $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intronic.vcf.tmp

mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_genic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_genic.vcf
mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intergenic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intergenic.vcf
mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_exonic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_exonic.vcf
mv $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intronic.vcf.tmp $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intronic.vcf


module load Biopython/1.78-foss-2020b

#cd to where your genomics_general files are
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_genic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_genic.geno.gz 

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intergenic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_intergenic.geno.gz 

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_exonic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_exonic.geno.gz 

python3 parseVCF.py -i  $wkdir/29_bedtools_intersect_genomic_region/inopinata_X_intronic.vcf | gzip > $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_intronic.geno.gz 


#cd to where your genomics_general files are
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_genic.geno.gz -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_X_genic_stats.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_intergenic.geno.gz  -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_X_intergenic_stats.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_exonic.geno.gz  -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_X_exonic.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/30_popgenwindows_genomic_regions/00_parseVCF/inopinata_X_intronic.geno.gz  -o  $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/inopinata_X_intronic.csv.gz -f phased  

cd $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py/

gunzip inopinata_X*

cd $wkdir/30_popgenwindows_genomic_regions/01_popgenWindows_py

awk 'BEGIN {FS=","} {OFS=","} {print $0,"inopinata"}' inopinata_X_genic_stats.csv > $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv/inopinata_X_genic_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"inopinata"}' inopinata_X_intergenic_stats.csv > $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv/inopinata_X_intergenic_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"inopinata"}' inopinata_X_exonic.csv > $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv/inopinata_X_exonic.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"inopinata"}' inopinata_X_intronic.csv > $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv/inopinata_X_intronic.csv

cd $wkdir/30_popgenwindows_genomic_regions/02_prep_tsv

awk 'BEGIN {FS=","} {OFS=","} {print $0,"exonic"}' inopinata_X_exonic.csv > inopinata_X_exonic.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"genic"}' inopinata_X_genic_stats.csv > inopinata_X_genic_stats.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"intergenic"}' inopinata_X_intergenic_stats.csv > inopinata_X_intergenic_stats.csv.tmp1
awk 'BEGIN {FS=","} {OFS=","} {print $0,"intronic"}' inopinata_X_intronic.csv > inopinata_X_intronic.csv.tmp1

sed -i -e 's/Sp34_ChrX/X/g' inopinata_X_exonic.csv.tmp1
sed -i -e 's/Sp34_ChrX/X/g' inopinata_X_genic_stats.csv.tmp1
sed -i -e 's/Sp34_ChrX/X/g' inopinata_X_intergenic_stats.csv.tmp1
sed -i -e 's/Sp34_ChrX/X/g' inopinata_X_intronic.csv.tmp1

sed -i '1d' inopinata_X_exonic.csv.tmp1
sed -i '1d'  inopinata_X_genic_stats.csv.tmp1
sed -i '1d'  inopinata_X_intergenic_stats.csv.tmp1
sed -i '1d' inopinata_X_intronic.csv.tmp1

#combine elegans and inopinata files, add headers


cat elegans_24_exonic.csv.tmp1 elegans_24_intronic.csv.tmp1 inopinata_24_exonic.csv.tmp1 inopinata_24_intronic.csv.tmp1 inopinata_X_exonic.csv.tmp1 inopinata_X_intronic.csv.tmp1 > exon_intron_pi.tmp

echo -e "scaffold,start,end,mid,sites,pi_all,species,gen_region" | cat - exon_intron_pi.tmp > exon_intron_pi.csv

cat elegans_24_genic_stats.csv.tmp1 elegans_24_intergenic_stats.csv.tmp1 inopinata_24_genic_stats.csv.tmp1 inopinata_24_intergenic_stats.csv.tmp1 inopinata_X_genic_stats.csv.tmp1 inopinata_X_intergenic_stats.csv.tmp1  > genic_intergenic_pi.tmp

echo -e "scaffold,start,end,mid,sites,pi_all,species,gen_region" | cat - genic_intergenic_pi.tmp > genic_intergenic_pi.csv








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
	#195265

	#this is the number of loci


#get locus lengths

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,$3-$2+1}' 02_inopinata_24_sorted_merge_bed > 03_inopinata_loci_lengths

echo -e "Chr\tbp_start\tbp_end\tlocus_length" | cat -  03_inopinata_loci_lengths >  inopinata_loci_lengths.tsv


#number of sites, number of snps


#number of sites

grep -v "#" $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf | wc -l 
	
	#4835565
	#this is the number of sites

#number of biallelic snps

wc -l $wkdir/19_cat_biallelic_invariant/01_inopinata_24_biallelic_snps_no_header

	#218388
	#this is the number of biallelic snps






#splitstree
#prep file for splitstree
#get 20,000 random snps for splitstree... 

mkdir $wkdir/25_splitstree


shuf -n 20000 $wkdir/19_cat_biallelic_invariant/01_inopinata_24_biallelic_snps_no_header > $wkdir/25_splitstree/inopinata_24_biallelic_snps_20k_random_snps_no_header

cat $wkdir/19_cat_biallelic_invariant/00_header $wkdir/25_splitstree/inopinata_24_biallelic_snps_20k_random_snps_no_header > $wkdir/25_splitstree/inopinata_24_biallelic_snps_20k_random_snps.vcf

#module load easybuild  GCC/6.3.0-2.27  OpenMPI/2.0.2 Python/3.7.0




###
#LD ; intrachromosomal LD

#used VCFtools/0.1.16-foss-2017b-Perl-5.26.1

#module load VCFtools/0.1.16-foss-2017b-Perl-5.26.1

mkdir $wkdir/28_VCFTools_LD/
mkdir $wkdir/28_VCFTools_LD/00_grep_chr
mkdir $wkdir/28_VCFTools_LD/01_VCFTools_LD
mkdir $wkdir/28_VCFTools_LD/02_means_R

/scratch/gcwoodruff/pop_gen/19_cat_biallelic_invariant/

#split chromosomes to estimate LD in parallel

cd $wkdir/19_cat_biallelic_invariant


grep -v "#" inopinata_24_sorted.vcf > inopinata_24_sorted_no_header

awk '$1 == "Sp34_Chr1"' inopinata_24_sorted_no_header > $wkdir/28_VCFTools_LD//00_grep_chr/Sp34_Chr1_no_header
awk '$1 == "Sp34_Chr2"' inopinata_24_sorted_no_header > $wkdir/28_VCFTools_LD//00_grep_chr/Sp34_Chr2_no_header
awk '$1 == "Sp34_Chr3"' inopinata_24_sorted_no_header > $wkdir/28_VCFTools_LD//00_grep_chr/Sp34_Chr3_no_header
awk '$1 == "Sp34_Chr4"' inopinata_24_sorted_no_header > $wkdir/28_VCFTools_LD//00_grep_chr/Sp34_Chr4_no_header
awk '$1 == "Sp34_Chr5"' inopinata_24_sorted_no_header > $wkdir/28_VCFTools_LD//00_grep_chr/Sp34_Chr5_no_header

cat inopinata_24_sorted_header $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr1_no_header > $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr1.vcf
cat inopinata_24_sorted_header $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr2_no_header > $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr2.vcf
cat inopinata_24_sorted_header $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr3_no_header > $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr3.vcf
cat inopinata_24_sorted_header $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr4_no_header > $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr4.vcf
cat inopinata_24_sorted_header $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr5_no_header > $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr5.vcf
#already split X
#link to the X vcf (your path will differ)
ln /scratch/gcwoodruff/pop_gen/20_split_autosomes_X_sex/02_split_X/inopinata_females_X.vcf /scratch/gcwoodruff/pop_gen/28_VCFTools_LD/00_grep_chr/inopinata_females_X.vcf

cd $wkdir/28_VCFTools_LD/01_VCFTools_LD/

#get pairwise LD for every site within each chromosome
vcftools --vcf $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr1.vcf --out Sp34_Chr1_out --geno-r2
vcftools --vcf $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr2.vcf --out Sp34_Chr2_out --geno-r2
vcftools --vcf $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr3.vcf --out Sp34_Chr3_out --geno-r2
vcftools --vcf $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr4.vcf --out Sp34_Chr4_out --geno-r2
vcftools --vcf $wkdir/28_VCFTools_LD/00_grep_chr/Sp34_Chr5.vcf --out Sp34_Chr5_out --geno-r2
vcftools --vcf $wkdir/28_VCFTools_LD/00_grep_chr/inopinata_females_X.vcf --out inopinata_females_X_out --geno-r2

#get mean LD for each site (script means.R in folder /additional_files/)
Rscript $wkdir/additional_files/means.R inopinata_females_X_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/inopinata_females_X_out.geno.ld
Rscript $wkdir/additional_files/means.R Sp34_Chr1_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/Sp34_Chr1_out.geno.ld
Rscript $wkdir/additional_files/means.R Sp34_Chr2_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/Sp34_Chr2_out.geno.ld
Rscript $wkdir/additional_files/means.R Sp34_Chr3_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/Sp34_Chr3_out.geno.ld
Rscript $wkdir/additional_files/means.R Sp34_Chr4_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/Sp34_Chr4_out.geno.ld
Rscript $wkdir/additional_files/means.R Sp34_Chr5_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/Sp34_Chr5_out.geno.ld

#ok, now what for LD. bedtools to get means in genomic windows.

#turn into bed
mkdir $wkdir/28_VCFTools_LD/03_bed/

cd  $wkdir/additional_files/means.R Sp34_Chr5_out.geno.ld $wkdir/28_VCFTools_LD/02_means_R/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print "X",$1-1,$1, $2}' inopinata_females_X_out.geno.ld > $wkdir/28_VCFTools_LD/03_bed/X.bed
awk 'BEGIN {FS="\t"} {OFS="\t"} {print "I",$1-1,$1, $2}' Sp34_Chr1_out.geno.ld > $wkdir/28_VCFTools_LD/03_bed/I.bed
awk 'BEGIN {FS="\t"} {OFS="\t"} {print "II",$1-1,$1, $2}' Sp34_Chr2_out.geno.ld > $wkdir/28_VCFTools_LD/03_bed/II.bed
awk 'BEGIN {FS="\t"} {OFS="\t"} {print "III",$1-1,$1, $2}' Sp34_Chr3_out.geno.ld > $wkdir/28_VCFTools_LD/03_bed/III.bed
awk 'BEGIN {FS="\t"} {OFS="\t"} {print "IV",$1-1,$1, $2}' Sp34_Chr4_out.geno.ld > $wkdir/28_VCFTools_LD/03_bed/IV.bed
awk 'BEGIN {FS="\t"} {OFS="\t"} {print "V",$1-1,$1, $2}' Sp34_Chr5_out.geno.ld > $wkdir/28_VCFTools_LD/03_bed/V.bed

cd $wkdir/28_VCFTools_LD/03_bed/

cat * > inopinata_LD.bed
#sort bed
sort -V -k1,1 -k2,2 inopinata_LD.bed > inopinata_LD_sort.bed

#bedtools windows...

mkdir $wkdir/28_VCFTools_LD/04_bedtools_windows/


cd $wkdir/28_VCFTools_LD/03_bed/

bedtools map -o mean -c 4 -a $wkdir/inopinata_genome/inopinata.50bp.windows2 -b inopinata_LD_sort.bed > $wkdir/28_VCFTools_LD/04_bedtools_windows/00_inopinata_LD_50bp_windows

cd $wkdir/28_VCFTools_LD/04_bedtools_windows/

#remove missing data
awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 00_inopinata_LD_50bp_windows > 01_inopinata_LD_50bp_windows_no_missing_data

#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/inopinata_genome/inopinata.10kb.windows2 -b 01_inopinata_LD_50bp_windows_no_missing_data > 02_inopinata_LD_10kb_windows


#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 02_inopinata_LD_10kb_windows > 03_inopinata_LD_10kb_windows_no_missing_data

#add inopinata

cd $wkdir/28_VCFTools_LD/04_bedtools_windows/


awk 'BEGIN {OFS="\t"} {print $1,$2+1,$3, "C. inopinata"}' 03_inopinata_LD_10kb_windows_no_missing_data > 04_inopinata_LD_10kb_windows_no_missing_data_species

#combine elegans and inopinata (see elegans_workflow.sh for elegans data; essentially same workflow)

cat $wkdir/elegans/30_VCFTools_LD/04_bedtools_windows/04_elegans_LD_10kb_windows_no_missing_data_species 04_inopinata_LD_10kb_windows_no_missing_data_species > 05_ele_ino_LD

#add header
echo -e "Chr\tBP\tR2\tSpecies" | cat - 05_ele_ino_LD > LD_ele_ino.tsv






cd $wkdir/25_splitstree/

#vcf2phylip.py by Edgardo M. Ortiz
#available at https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py

#vcf to nexus
python $wkdir/additional_files/vcf2phylip.py --input inopinata_24_biallelic_snps_20k_random_snps.vcf -n -f

#For splitstree figure, file "inopinata_24_biallelic_snps_20k_random_snps.min4.nexus" was loaded in SplitsTree4 (version 4.16.1, built 19 May 2020). Had to edit .nexus file-- 19507 characters (not 20000) were incorporated into the nexus. Replace "NCHAR=20000" with "NCHAR=19507" 




