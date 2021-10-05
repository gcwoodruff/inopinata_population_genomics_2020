#This is the workflow for generating C. elegans genotypes and population genetics stats for Woodruff et al. 2021 C. inopinata population genomics paper, "Alignment of genetic differentiation across trophic levels in a fig community"



#put working directory here
wkdir="/projects/phillipslab/gavincw/pop_gen_reproduce_june_2021"


mkdir $wkdir/elegans
cd $wkdir/elegans
mkdir 00_wget_bam
cd 00_wget_bam

#thanks CENDR for access to BAM files
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/CB4856.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/CX11314.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/DL238.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ED3017.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/EG4725.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/JT11398.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/JU258.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/JU775.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/LKC34.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/MY16.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/MY23.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/N2.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/BRC20067.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ECA246.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ECA251.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/CX11271.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/CX11276.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/CX11285.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/DL200.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/DL226.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ECA36.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ED3040.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ED3048.bam
wget https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/ED3049.bam

#define EcoRI cut sites in C. elegans genome with EMBOSS fuzznuc
#elegans genome assembly used... PRJNA13758

mkdir $wkdir/elegans/elegans_genome/

cd $wkdir/elegans/elegans_genome/

gunzip caenorhabditis_elegans.PRJNA13758.WBPS15.genomic.fa.gz

#get reference genome, thanks WormBase
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS15.genomic.fa.gz

mv caenorhabditis_elegans.PRJNA13758.WBPS15.genomic.fa elegans_genome.fa

#module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 EMBOSS/6.6.0

mkdir $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/

fuzznuc -sequence $wkdir/elegans/elegans_genome/elegans_genome.fa -pattern 'GAATTC' -complement -rformat gff -outfile $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_EcoRI_sites_out.gff



#now, make a bed file

cd $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/

grep -v '#' elegans_EcoRI_sites_out.gff | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$4,$5,"EcoRI_site"}' > elegans_EcoRI_sites_out.bed &

#remove mito

grep -v "MtDNA" elegans_EcoRI_sites_out.bed > elegans_EcoRI_sites_out_tmp.bed
mv elegans_EcoRI_sites_out_tmp.bed elegans_EcoRI_sites_out.bed

#sort the bedfile

#module load bedtools/2.25.0

bedtools sort -i elegans_EcoRI_sites_out.bed | uniq > elegans_EcoRI_sites_out_sort.bed


#get flanking bp on each side of the cut site to define regions to extract; this approximates the rad tags. 166 bp was chosen at that time because it was the average locus length in a previous implementation of the inopinata pop gen workflow.

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2-166,$3+166,$4}' elegans_EcoRI_sites_out_sort.bed > elegans_pseudo_rad_regions.bed

#repair negative value (an EcoRI site is close to the beginning of chr V)

sed -i -e 's/\t-165/\t0/g' elegans_pseudo_rad_regions.bed


#get those C. elegans alignments just in regions that flank EcoRI sites.
	#module load racs-eb SAMtools/1.7-intel-2018a

mkdir $wkdir/elegans/02_samtools_sort/

cd $wkdir/elegans/00_wget_bam

for i in *; do samtools sort $i > $wkdir/elegans/02_samtools_sort/$i; done

mkdir $wkdir/elegans/03_samtools_view_pseudo_rad

cd $wkdir/elegans/02_samtools_sort/


samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed BRC20067.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/BRC20067.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed CB4856.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/CB4856.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed CX11271.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/CX11271.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed CX11276.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/CX11276.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed CX11285.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/CX11285.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed CX11314.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/CX11314.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed DL200.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/DL200.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed DL226.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/DL226.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed DL238.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/DL238.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ECA246.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ECA246.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ECA251.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ECA251.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ECA36.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ECA36.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ED3017.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ED3017.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ED3040.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ED3040.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ED3048.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ED3048.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed ED3049.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/ED3049.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed EG4725.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/EG4725.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed JT11398.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/JT11398.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed JU258.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/JU258.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed JU775.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/JU775.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed LKC34.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/LKC34.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed MY16.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/MY16.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed MY23.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/MY23.bam
samtools view -b -L $wkdir/elegans/01_EcoRI_cut_sites_emboss_fuzznuc/elegans_pseudo_rad_regions.bed N2.bam > $wkdir/elegans/03_samtools_view_pseudo_rad/N2.bam


#do genotype calls

mkdir $wkdir/elegans/04_bcftools_mpileup_call/

cd $wkdir/elegans/03_samtools_view_pseudo_rad/

bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa BRC20067.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/BRC20067.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ECA246.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ECA246.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa CB4856.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/CB4856.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ECA251.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ECA251.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa CX11271.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/CX11271.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa CX11276.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/CX11276.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa CX11285.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/CX11285.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa CX11314.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/CX11314.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa DL200.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/DL200.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa DL226.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/DL226.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa DL238.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/DL238.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ECA36.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ECA36.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ED3017.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ED3017.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ED3040.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ED3040.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ED3048.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ED3048.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa ED3049.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/ED3049.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa EG4725.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/EG4725.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa JT11398.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/JT11398.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa JU258.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/JU258.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa JU775.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/JU775.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa LKC34.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/LKC34.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa MY16.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/MY16.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa MY23.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/MY23.bcf
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/elegans_genome.fa N2.bam | bcftools call -c -Ob -o $wkdir/elegans/04_bcftools_mpileup_call/N2.bcf


mkdir $wkdir/elegans/05_bcftools_view_depth_15

cd $wkdir/elegans/04_bcftools_mpileup_call

for i in *; do bcftools view -i 'DP>=15' $i -Ob > $wkdir/elegans/05_bcftools_view_depth_15/$i; done

mkdir $wkdir/elegans/06_bcftools_merge

cd $wkdir/elegans/05_bcftools_view_depth_15

#index bcf

for i in *; do bcftools index $i; done


#merge bcf

bcftools merge --info-rules DP:join,MQ0F:join,AF1:join,AC1:join,DP4:join,MQ:join,FQ:join -m snps BRC20067.bcf CB4853.bcf CB4856.bcf CB4858.bcf CX11271.bcf CX11276.bcf CX11285.bcf CX11314.bcf DL200.bcf DL226.bcf DL238.bcf ECA36.bcf ED3017.bcf ED3040.bcf ED3048.bcf ED3049.bcf EG4725.bcf JT11398.bcf JU258.bcf JU775.bcf LKC34.bcf MY16.bcf MY23.bcf N2.bcf -o $wkdir/elegans/06_bcftools_merge/elegans_24.vcf


#just get sites with >=80% of samples having genotype calls

mkdir $wkdir/elegans/07_bcftools_view_gt_count

cd $wkdir/elegans/06_bcftools_merge

bcftools view -i 'COUNT(GT="mis")<5' elegans_24.vcf > $wkdir/elegans/07_bcftools_view_gt_count/elegans_24.vcf &




#get just the invariant sites....

mkdir $wkdir/elegans/08_invariant_sites

cd $wkdir/elegans/07_bcftools_view_gt_count

grep "#" elegans_24.vcf > $wkdir/elegans/08_invariant_sites/header

grep -v "#" elegans_24.vcf > $wkdir/elegans/08_invariant_sites/elegans_24_no_header

cd $wkdir/elegans/08_invariant_sites/

awk 'BEGIN {FS="\t"} {OFS="\t"} $5 == "."' elegans_24_no_header > elegans_24_no_alt

cat header elegans_24_no_alt > elegans_24_invariant_sites.vcf






#get biallelic snps with at least one copy per allele

mkdir $wkdir/elegans/09_bcftools_view_biallelic_snps_maf

cd $wkdir/elegans/07_bcftools_view_gt_count

bcftools view -m2 -M2 -v snps --min-ac 2:minor elegans_24.vcf > $wkdir/elegans/09_bcftools_view_biallelic_snps_maf/elegans_24.vcf




#sort invariant sites vcf

cd $wkdir/elegans/08_invariant_sites/

cat elegans_24_invariant_sites.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > elegans_24_invariant_sites_sort.vcf


#sort biallelic snps vcf


cd $wkdir/elegans/09_bcftools_view_biallelic_snps_maf/

cat elegans_24.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > elegans_24_biallelic_snps_sort.vcf

#combine



mkdir $wkdir/elegans/10_cat

grep "#" elegans_24_biallelic_snps_sort.vcf > $wkdir/elegans/10_cat/00_header


grep -v "#" elegans_24_biallelic_snps_sort.vcf >  $wkdir/elegans/10_cat/01_elegans_24_biallelic_snps_no_header

cd $wkdir/elegans/08_invariant_sites/

grep -v "#" elegans_24_invariant_sites_sort.vcf > $wkdir/elegans/10_cat/02_elegans_24_invariant_sites_no_header

cd $wkdir/elegans/10_cat/

cat 00_header 01_elegans_24_biallelic_snps_no_header 02_elegans_24_invariant_sites_no_header > 03_combined_unsorted.vcf


#sort

cat 03_combined_unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > elegans_24_sorted.vcf




#estimate pi in 10kb windows


mkdir $wkdir/elegans/11_popgenwindows_py/

cd $wkdir/elegans/11_popgenwindows_py/

python2.7 parseVCF.py -i $wkdir/elegans/10_cat/elegans_24_sorted.vcf | gzip > elegans_24.geno.gz

python2.7 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g elegans_24.geno.gz -o elegans_stats.csv.gz -f phased  &

gunzip elegans_stats.csv.gz

awk 'BEGIN {FS=","} {OFS=","} {print $0,"C. elegans"}' elegans_stats.csv > elegans_stats.csv.tmp

mv elegans_stats.csv.tmp elegans_stats.csv

sed -i '1d' elegans_stats.csv
	#this was combined with inopinata data for the comparison in Figure 2



#Fis estimate


mkdir $wkdir/elegans/12_stacks

populations -V $wkdir/elegans/10_cat/elegans_24_sorted.vcf -O $wkdir/elegans/12_stacks/ -M $wkdir/additional_files/elegans_pop_map_stacks.txt --sigma 3333 --genepop --structure --phylip


cd $wkdir/elegans/12_stacks
mkdir site_fis

#get fis info
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$3,$17}' elegans_24_sorted.p.sumstats.tsv > $wkdir/elegans/12_stacks/site_fis/01_fis_cols

cd site_fis

#remove first two lines

sed -i '1d' 01_fis_cols

sed -i '1d' 01_fis_cols

#repair bed file (ie, replace spaces with tabs)
sed -i -e 's/ /\t/g' 01_fis_cols


#use bedtools to get fis in 10kb windows

#module load bedtools/2.25.0

#get Fis over 50bp windows
bedtools map -o mean -c 4 -a $wkdir/additional_files/elegans.50bp.windows -b 01_fis_cols > 03_fis_50bp_windows 

#exclude missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 03_fis_50bp_windows > 04_fis_50bp_windows_no_missing_data



#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/elegans.10kb.windows -b 04_fis_50bp_windows_no_missing_data > 05_fis_10kb_windows 


#exclude missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 05_fis_10kb_windows > 06_fis_10kb_windows_no_missing_data

#add header
echo -e "Chr\tBP\tFis" | cat - 06_fis_10kb_windows_no_missing_data > 07_fis.tsv

#add species id to each row
awk 'BEGIN {OFS="\t"} {print $0, "elegans"}' 06_fis_10kb_windows_no_missing_data > 08_elegans_no_header
	#this is the file used for comparison with inopinata








