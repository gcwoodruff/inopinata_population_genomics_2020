#demultiplex workflow for Woodruff et al. 2021 C. inopinata population genomics paper, "Alignment of genetic differentiation across trophic levels in a fig community"

#if you have any questions, please contact Gavin Woodruff at gcwoodruff@ou.edu

#put working directory here
wkdir="/projects/phillipslab/gavincw/pop_gen_reproduce_june_2021"

#fastq files are in folder $wkdir/00_raw_reads/

#read quality was evaluated with fastqc v. 0.11.5

mkdir $wkdir/01_FastQC

fastqc -o $wkdir/01_FastQC/ 2706_NEBNext_Index_12_S149_L004_R1_001.fastq.gz 2706_NEBNext_Index_12_S149_L004_R2_001.fastq.gz 2706_NEBNext_Index_4_S147_L004_R1_001.fastq.gz 2706_NEBNext_Index_4_S147_L004_R2_001.fastq.gz 2706_NEBNext_Index_6_S148_L004_R1_001.fastq.gz 2706_NEBNext_Index_6_S148_L004_R2_001.fastq.gz Undetermined_S0_L004_R1_001.fastq.gz Undetermined_S0_L004_R2_001.fastq.gz


#used Flip2BeRAD.py (https://github.com/tylerhether/Flip2BeRAD) to reorient reads due to BestRAD library prep
mkdir $wkdir/02_Flip2BeRAD
mkdir $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_12_S149
mkdir $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_4_S147
mkdir $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_6_S148
mkdir $wkdir/02_Flip2BeRAD/Undetermined_S0

cd $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_12_S149
python ./Flip2BeRAD.py -f $wkdir/00_raw_reads/2706_NEBNext_Index_12_S149_L004_R1_001.fastq -r $wkdir/00_raw_reads/2706_NEBNext_Index_12_S149_L004_R2_001.fastq  -b $wkdir/barcodes/barcodes_no_labels -c AATTC -m 0 -o 2
cd $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_4_S147
python ./Flip2BeRAD.py -f $wkdir/00_raw_reads/2706_NEBNext_Index_4_S147_L004_R1_001.fastq -r $wkdir/00_raw_reads/2706_NEBNext_Index_4_S147_L004_R2_001.fastq  -b $wkdir/barcodes/barcodes_no_labels -c AATTC -m 0 -o 2
cd $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_6_S148
python ./Flip2BeRAD.py -f $wkdir/00_raw_reads/2706_NEBNext_Index_6_S148_L004_R1_001.fastq -r $wkdir/00_raw_reads/2706_NEBNext_Index_6_S148_L004_R2_001.fastq  -b $wkdir/barcodes/barcodes_no_labels -c AATTC -m 0 -o 2
cd $wkdir/02_Flip2BeRAD/Undetermined_S0
python ./Flip2BeRAD.py -f $wkdir/00_raw_reads/Undetermined_S0_L004_R1_001.fastq -r $wkdir/00_raw_reads/Undetermined_S0_L004_R2_001.fastq  -b $wkdir/barcodes/barcodes_no_labels -c AATTC -m 0 -o 2

#used fastx_trimmer from FASTX-Toolkit 0.0.13 to remove first 3 bases so stacks process_radtags can find the barcodes

mkdir $wkdir/03_fastx_trimmer
mkdir $wkdir/03_fastx_trimmer/2706_NEBNext_Index_12_S149
mkdir $wkdir/03_fastx_trimmer/2706_NEBNext_Index_4_S147
mkdir $wkdir/03_fastx_trimmer/2706_NEBNext_Index_6_S148
mkdir $wkdir/03_fastx_trimmer/Undetermined_S0

fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_12_S149/filtered_reverse.fastq  -o $wkdir/03_fastx_trimmer/2706_NEBNext_Index_12_S149/filtered_reverse.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_4_S147/filtered_reverse.fastq  -o $wkdir/03_fastx_trimmer/2706_NEBNext_Index_4_S147/filtered_reverse.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_6_S148/filtered_reverse.fastq  -o $wkdir/03_fastx_trimmer/2706_NEBNext_Index_6_S148/filtered_reverse.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/Undetermined_S0/filtered_reverse.fastq  -o $wkdir/03_fastx_trimmer/Undetermined_S0/filtered_reverse.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_12_S149/filtered_forward.fastq  -o $wkdir/03_fastx_trimmer/2706_NEBNext_Index_12_S149/filtered_forward.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_4_S147/filtered_forward.fastq  -o $wkdir/03_fastx_trimmer/2706_NEBNext_Index_4_S147/filtered_forward.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/2706_NEBNext_Index_6_S148/filtered_forward.fastq  -o $wkdir/03_fastx_trimmer/2706_NEBNext_Index_6_S148/filtered_forward.fastq
fastx_trimmer -f 3 -Q 33 -i $wkdir/02_Flip2BeRAD/Undetermined_S0/filtered_forward.fastq  -o $wkdir/03_fastx_trimmer/Undetermined_S0/filtered_forward.fastq

#use stacks process_radtags (v. 2.0) to demultiplex

mkdir $wkdir/04_process_radtags
mkdir $wkdir/04_process_radtags/2706_NEBNext_Index_12_S149
mkdir $wkdir/04_process_radtags/2706_NEBNext_Index_4_S147
mkdir $wkdir/04_process_radtags/2706_NEBNext_Index_6_S148
mkdir $wkdir/04_process_radtags/Undetermined_S0

cd $wkdir/03_fastx_trimmer/2706_NEBNext_Index_12_S149/
process_radtags -1 filtered_forward.fastq -2 filtered_reverse.fastq -o $wkdir/04_process_radtags/2706_NEBNext_Index_12_S149/ -b $wkdir/barcodes/barcodes -e ecoRI -r -c -q
cd $wkdir/03_fastx_trimmer/2706_NEBNext_Index_4_S147/
process_radtags -1 filtered_forward.fastq -2 filtered_reverse.fastq -o $wkdir/04_process_radtags/2706_NEBNext_Index_4_S147/ -b $wkdir/barcodes/barcodes -e ecoRI -r -c -q
cd $wkdir/03_fastx_trimmer/2706_NEBNext_Index_6_S148/
process_radtags -1 filtered_forward.fastq -2 filtered_reverse.fastq -o $wkdir/04_process_radtags/2706_NEBNext_Index_6_S148/ -b $wkdir/barcodes/barcodes -e ecoRI -r -c -q
cd $wkdir/03_fastx_trimmer/Undetermined_S0/
process_radtags -1 filtered_forward.fastq -2 filtered_reverse.fastq -o $wkdir/04_process_radtags/Undetermined_S0/ -b $wkdir/barcodes/barcodes -e ecoRI -r -c -q


#paired-end reads corresponding to the same individual sample across different runs were concatenated and renamed with cat and mv (forward and reverse files were retained) to generate fastq files with these ids (placed in folder $wkdir/05_reads_to_align/--
#A03.1.fq 
#A03.2.fq 
#A05.1.fq 
#A05.2.fq 
#B03.1.fq 
#B03.2.fq 
#dec_2016_D01.1.fq 
#dec_2016_D01.2.fq 
#dec_2016_D02.1.fq 
#dec_2016_D02.2.fq 
#dec_2016_D03.1.fq 
#dec_2016_D03.2.fq 
#dec_2016_D04.1.fq 
#dec_2016_D04.2.fq 
#dec_2016_D05.1.fq 
#dec_2016_D05.2.fq 
#dec_2016_D10.1.fq 
#dec_2016_D10.2.fq 
#may_2017_A11.1.fq 
#may_2017_A11.2.fq 
#may_2017_A12.1.fq 
#may_2017_A12.2.fq 
#may_2017_B11.1.fq 
#may_2017_B11.2.fq 
#may_2017_B12.1.fq 
#may_2017_B12.2.fq 
#may_2017_C11.1.fq 
#may_2017_C11.2.fq 
#may_2017_C12.1.fq 
#may_2017_C12.2.fq 
#may_2017_D11.1.fq 
#may_2017_D11.2.fq 
#may_2017_D12.1.fq 
#may_2017_D12.2.fq 
#may_2017_E11.1.fq 
#may_2017_E11.2.fq 
#may_2017_E12.1.fq 
#may_2017_E12.2.fq 
#may_2017_F11.1.fq 
#may_2017_F11.2.fq 
#may_2017_G12.1.fq 
#may_2017_G12.2.fq 
#may_2017_H10.1.fq 
#may_2017_H10.2.fq 
#may_2017_H11.1.fq 
#may_2017_H11.2.fq 
#may_2017_H12.1.fq 
#may_2017_H12.2.fq 