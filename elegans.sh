#This is the workflow for generating C. elegans genotypes and population genetics 
#stats for Woodruff et al. 2023 C. inopinata population genomics paper, 
#"Patterns of genomic diversity in a fig-associated close relative of C. elegans"


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


wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS16.annotations.gff3.gz
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS16.canonical_geneset.gtf.gz

gunzip *


mkdir $wkdir/elegans/07_b_emboss_fuzznuc/


mkdir $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc

#emboss to find EcoRI cut sites in the C. elegans genome assembly
#module load EMBOSS/6.6.0-foss-2018b

fuzznuc -sequence $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa -pattern 'GAATTC' -complement -rformat gff -outfile $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_out.gff


#now, make a bed file

cd $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/

gff2bed < elegans_EcoRI_sites_out.gff > elegans_EcoRI_sites_out.bed

#remove mito

grep -v "MtDNA" elegans_EcoRI_sites_out.bed > elegans_EcoRI_sites_out_tmp.bed
mv elegans_EcoRI_sites_out_tmp.bed elegans_EcoRI_sites_out.bed

#sort the bedfile

#module load BEDTools/2.27.1-foss-2018b

bedtools sort -i elegans_EcoRI_sites_out.bed | uniq > elegans_EcoRI_sites_out_sort.bed

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,"EcoRI_site"}' elegans_EcoRI_sites_out_sort.bed | uniq > elegans_EcoRI_sites_out_sort.bed.tmp

#get flanking bp on each side of the cut site to define regions to extract; this approximates the rad tags. 166 bp was chosen at that time because it was the average locus length in a previous implementation of the inopinata pop gen workflow.

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2-166,$3+166,$4}' elegans_EcoRI_sites_out_sort.bed.tmp > elegans_pseudo_rad_regions.bed

#repair negative value (an EcoRI site is close to the beginning of chr V)

sed -i -e 's/\t-166/\t0/g' elegans_pseudo_rad_regions.bed

#sort alignments
mkdir $wkdir/elegans/02_samtools_sort/

cd $wkdir/elegans/00_wget_bam

for i in *; do samtools sort $i > $wkdir/elegans/02_samtools_sort/$i; done

#get those C. elegans alignments just in regions that flank EcoRI sites.
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

#call genotypes

mkdir $wkdir/elegans/12_bcftools_mpileup_call/
mkdir $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/

cd $wkdir/elegans/03_samtools_view_pseudo_rad/

bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa BRC20067.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/BRC20067.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa CB4856.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/CB4856.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa CX11271.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/CX11271.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa CX11276.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/CX11276.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa CX11285.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/CX11285.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa CX11314.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/CX11314.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa DL200.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/DL200.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa DL226.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/DL226.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa DL238.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/DL238.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ECA246.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ECA246.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ECA251.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ECA251.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ECA36.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ECA36.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ED3017.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ED3017.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ED3040.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ED3040.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ED3048.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ED3048.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa ED3049.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/ED3049.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa EG4725.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/EG4725.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa JT11398.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/JT11398.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa JU258.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/JU258.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa JU775.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/JU775.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa LKC34.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/LKC34.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa MY16.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/MY16.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa MY23.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/MY23.bam
bcftools mpileup -Ou -f $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa N2.bam | bcftools call -c -Ob -o $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/N2.bam


#filter reads with >14x coverage

mkdir $wkdir/elegans/13_bcftools_view_depth_15/
mkdir $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/
cd $wkdir/elegans/12_bcftools_mpileup_call/og_pseudo_rad/	

bcftools view -i 'DP>=15' BRC20067.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/BRC20067.vcf
bcftools view -i 'DP>=15' CB4856.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/CB4856.vcf
bcftools view -i 'DP>=15' CX11271.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/CX11271.vcf
bcftools view -i 'DP>=15' CX11276.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/CX11276.vcf
bcftools view -i 'DP>=15' CX11285.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/CX11285.vcf
bcftools view -i 'DP>=15' CX11314.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/CX11314.vcf
bcftools view -i 'DP>=15' DL200.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/DL200.vcf
bcftools view -i 'DP>=15' DL226.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/DL226.vcf
bcftools view -i 'DP>=15' DL238.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/DL238.vcf
bcftools view -i 'DP>=15' ECA246.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ECA246.vcf
bcftools view -i 'DP>=15' ECA251.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ECA251.vcf
bcftools view -i 'DP>=15' ECA36.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ECA36.vcf
bcftools view -i 'DP>=15' ED3017.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ED3017.vcf
bcftools view -i 'DP>=15' ED3040.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ED3040.vcf
bcftools view -i 'DP>=15' ED3048.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ED3048.vcf
bcftools view -i 'DP>=15' ED3049.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/ED3049.vcf
bcftools view -i 'DP>=15' EG4725.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/EG4725.vcf
bcftools view -i 'DP>=15' JT11398.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/JT11398.vcf
bcftools view -i 'DP>=15' JU258.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/JU258.vcf
bcftools view -i 'DP>=15' JU775.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/JU775.vcf
bcftools view -i 'DP>=15' LKC34.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/LKC34.vcf
bcftools view -i 'DP>=15' MY16.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/MY16.vcf
bcftools view -i 'DP>=15' MY23.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/MY23.vcf
bcftools view -i 'DP>=15' N2.bcf > $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/N2.vcf

#convert back to bcf and index


mkdir $wkdir/elegans/14_bcftools_view_bcf
mkdir $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad

cd $wkdir/elegans/13_bcftools_view_depth_15/og_pseudo_rad/

bcftools view BRC20067.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/BRC20067.bcf
bcftools view CB4856.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/CB4856.bcf
bcftools view CX11271.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/CX11271.bcf
bcftools view CX11276.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/CX11276.bcf
bcftools view CX11285.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/CX11285.bcf
bcftools view CX11314.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/CX11314.bcf
bcftools view DL200.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/DL200.bcf
bcftools view DL226.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/DL226.bcf
bcftools view DL238.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/DL238.bcf
bcftools view ECA246.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ECA246.bcf
bcftools view ECA251.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ECA251.bcf
bcftools view ECA36.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ECA36.bcf
bcftools view ED3017.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ED3017.bcf
bcftools view ED3040.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ED3040.bcf
bcftools view ED3048.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ED3048.bcf
bcftools view ED3049.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/ED3049.bcf
bcftools view EG4725.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/EG4725.bcf
bcftools view JT11398.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/JT11398.bcf
bcftools view JU258.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/JU258.bcf
bcftools view JU775.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/JU775.bcf
bcftools view LKC34.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/LKC34.bcf
bcftools view MY16.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/MY16.bcf
bcftools view MY23.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/MY23.bcf
bcftools view N2.vcf -O b > $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/N2.bcf


#merge samples


mkdir $wkdir/elegans/15_bcftools_merge/


cd $wkdir/elegans/14_bcftools_view_bcf/og_pseudo_rad/
bcftools merge --info-rules DP:join,MQ0F:join,AF1:join,AC1:join,DP4:join,MQ:join,FQ:join -m snps BRC20067.bcf CB4856.bcf CX11271.bcf CX11276.bcf CX11285.bcf CX11314.bcf DL200.bcf DL226.bcf DL238.bcf ECA246.bcf ECA251.bcf ECA36.bcf ED3017.bcf ED3040.bcf ED3048.bcf ED3049.bcf EG4725.bcf JT11398.bcf JU258.bcf JU775.bcf LKC34.bcf MY16.bcf MY23.bcf N2.bcf -o $wkdir/elegans/15_bcftools_merge/elegans_24_og_pseudo_rad.vcf

#get those sites with ~>80% samples having genotype calls

mkdir $wkdir/elegans/16_bcftools_view_gt_count/

cd $wkdir/elegans/15_bcftools_merge

bcftools view -i 'COUNT(GT="mis")<5' elegans_24_og_pseudo_rad.vcf > $wkdir/elegans/16_bcftools_view_gt_count/elegans_24_og_pseudo_rad.vcf


#sort

mkdir $wkdir/elegans/17_sort/

cd $wkdir/elegans/16_bcftools_view_gt_count

cat elegans_24_og_pseudo_rad.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen/elegans/17_sort/tmp3"}' > $wkdir/elegans/17_sort/elegans_24_og_pseudo_rad.vcf


#get invariant sites
mkdir $wkdir/elegans/18_invariant_sites/

cd $wkdir/elegans/17_sort/


grep "#" elegans_24_og_pseudo_rad.vcf > $wkdir/elegans/18_invariant_sites/elegans_24_og_pseudo_rad_header

grep -v "#" elegans_24_og_pseudo_rad.vcf > $wkdir/elegans/18_invariant_sites/elegans_24_og_pseudo_rad_no_header

cd $wkdir/elegans/18_invariant_sites/

awk 'BEGIN {FS="\t"} {OFS="\t"} $5 == "."' elegans_24_og_pseudo_rad_no_header  > elegans_24_og_pseudo_rad_no_alt

cat elegans_24_og_pseudo_rad_header elegans_24_og_pseudo_rad_no_alt > inopinata_24_og_pseudo_rad_invariant_sites.vcf

mv inopinata_24_og_pseudo_rad_invariant_sites.vcf elegans_24_og_pseudo_rad_invariant_sites.vcf

rm elegans_24_og_pseudo_rad_header
rm elegans_24_og_pseudo_rad_no_header
rm elegans_24_og_pseudo_rad_no_alt

#get biallelic snps with at least one copy per allele
mkdir $wkdir/elegans/19_bcftools_view_biallelic_snps_maf/

cd $wkdir/elegans/17_sort/

bcftools view -m2 -M2 -v snps --min-ac 2:minor elegans_24_og_pseudo_rad.vcf > $wkdir/elegans/19_bcftools_view_biallelic_snps_maf/elegans_24_og_pseudo_rad.vcf

#combine invariant and biallelic snps
mkdir $wkdir/elegans/20_cat_biallelic_invariant/
mkdir $wkdir/elegans/20_cat_biallelic_invariant/tmp
mkdir $wkdir/elegans/20_cat_biallelic_invariant/tmp2
mkdir $wkdir/elegans/20_cat_biallelic_invariant/tmp3

cd $wkdir/elegans/19_bcftools_view_biallelic_snps_maf/

grep "#" elegans_24_og_pseudo_rad.vcf > $wkdir/elegans/20_cat_biallelic_invariant/header_og_pseudo_rad
grep -v "#" elegans_24_og_pseudo_rad.vcf  >  $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_biallelic_snps_no_header_og_pseudo_rad

cd $wkdir/elegans/18_invariant_sites/

grep -v "#" elegans_24_og_pseudo_rad_invariant_sites.vcf > $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_invariant_sites_no_header

cd $wkdir/elegans/20_cat_biallelic_invariant/

cat header_og_pseudo_rad elegans_24_biallelic_snps_no_header_og_pseudo_rad elegans_24_og_pseudo_rad_invariant_sites_no_header > elegans_24_og_pseudo_rad_combined_unsorted.vcf

cat elegans_24_og_pseudo_rad_combined_unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n -T /scratch/gcwoodruff/pop_gen/elegans/20_cat_biallelic_invariant/tmp3"}' > elegans_24_og_pseudo_rad_24_sorted.vcf

rm header_og_pseudo_rad
rm elegans_24_biallelic_snps_no_header_og_pseudo_rad
rm elegans_24_og_pseudo_rad_invariant_sites_no_header
rm elegans_24_og_pseudo_rad_combined_unsorted.vcf

#get pi nucleotide diversity in 10 kb genomic windows


mkdir $wkdir/elegans/21_popgenWindows/
mkdir $wkdir/elegans/21_popgenWindows/00_parseVCF/
mkdir $wkdir/elegans/21_popgenWindows/01_popgenWindows_py/

#cd to your genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/


python3 parseVCF.py -i  $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf | gzip > $wkdir/elegans/21_popgenWindows/00_parseVCF/elegans_24_og_pseudo_rad_24.geno.gz 


#cd to your genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/21_popgenWindows/00_parseVCF/elegans_24_og_pseudo_rad_24.geno.gz -o  $wkdir/elegans/21_popgenWindows/01_popgenWindows_py/elegans_og_pseudo_rad_stats.csv.gz -f phased  

cd $wkdir/elegans/21_popgenWindows/01_popgenWindows_py/

gunzip elegans_og_pseudo_rad_stats.csv.gz


#Stacks for Fis
mkdir $wkdir/elegans/22_stacks_populations/

#module load Stacks/2.59-foss-2020a

populations -V  $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf -O $wkdir/elegans/22_stacks_populations/ -M $wkdir/additional_files/elegans_pop_map_stacks.txt  --sigma 3333 --genepop --structure --phylip 

#prep FIS file for stats and figures
mkdir $wkdir/elegans/22_stacks_populations/site_fis

cd $wkdir/elegans/22_stacks_populations/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$3,$17}' elegans_24_og_pseudo_rad_24_sorted.p.sumstats.tsv > $wkdir/elegans/22_stacks_populations/site_fis/01_fis_cols

cd $wkdir/elegans/22_stacks_populations/site_fis


#remove first two lines

sed -i '1d' 01_fis_cols

sed -i '1d' 01_fis_cols

#repair bed file (ie, replace spaces with tabs)
sed -i -e 's/ /\t/g' 01_fis_cols


#bedtools to get inbreeding coefficient Fis windows!
#now using BEDTools/2.27.1-foss-2018b
#module load BEDTools/2.27.1-foss-2018b


cd $wkdir/elegans/elegans_genome

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2} ' caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa.fai | sort > elegans_chr_sizes.txt

cd $wkdir/additional_files

bedtools makewindows -g $wkdir/elegans/elegans_genome/elegans_chr_sizes.txt -w 50 > elegans.50bp.windows
bedtools makewindows -g $wkdir/elegans/elegans_genome/elegans_chr_sizes.txt -w 10000 > elegans.10kb.windows

cd $wkdir/elegans/22_stacks_populations/site_fis

bedtools map -o mean -c 4 -a $wkdir/additional_files/elegans.50bp.windows -b 01_fis_cols > 03_fis_50bp_windows


#remove missing data
awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 03_fis_50bp_windows > 04_fis_50bp_windows_no_missing_data

#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/elegans.10kb.windows -b 04_fis_50bp_windows_no_missing_data > 05_fis_10kb_windows


#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 05_fis_10kb_windows > 06_fis_10kb_windows_no_missing_data

awk 'BEGIN {OFS="\t"} {print $0, "elegans"}' 06_fis_10kb_windows_no_missing_data > 07_fis_10kb_windows_no_missing_data_elegans



#coding site types!!!
cd $wkdir/elegans/20_cat_biallelic_invariant/

#module load tabix/0.2.6-GCCcore-7.3.0

#index vcf
bgzip -c elegans_24_og_pseudo_rad_24_sorted.vcf > elegans_24_og_pseudo_rad_24_sorted.vcf.gz
tabix -p vcf elegans_24_og_pseudo_rad_24_sorted.vcf.gz

#cd to your genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

mkdir  $wkdir/elegans/23_site_annotation_genomics_general

#module load Biopython/1.78-foss-2020b

#get coding site types

python3 codingSiteTypes.py -a $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.canonical_geneset.gtf -f gtf -r $wkdir/elegans/elegans_genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa -v $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf.gz -o $wkdir/elegans/23_site_annotation_genomics_general/elegans_24_og_pseudo_rad_codingSiteTypes.out --ignoreConflicts

mkdir  $wkdir/elegans/23_site_annotation_genomics_general/00_bedtools_intersect
cd $wkdir/elegans/23_site_annotation_genomics_general/


#get intersection of site types with site genotypes

mv elegans_24_og_pseudo_rad_codingSiteTypes.out 00_elegans_24_og_pseudo_rad_codingSiteTypes.out

grep -v "#" $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf > $wkdir/elegans/23_site_annotation_genomics_general/01_elegans_24_og_pseudo_rad_24_sorted_no_header_vcf

cd $wkdir/elegans/23_site_annotation_genomics_general/

sed '1d' 00_elegans_24_og_pseudo_rad_codingSiteTypes.out > 02_elegans_24_og_pseudo_rad_codingSiteTypes.out_no_header

#remove mitochondrial genome
grep -v "MtDNA" 02_elegans_24_og_pseudo_rad_codingSiteTypes.out_no_header > 03_elegans_24_og_pseudo_rad_codingSiteTypes.out_no_header_no_mito

mkdir $wkdir/elegans/23_site_annotation_genomics_general/tmp

mkdir $wkdir/elegans/23_site_annotation_genomics_general/tmp2

join -j1 <(<01_elegans_24_og_pseudo_rad_24_sorted_no_header_vcf awk '{print $1"-"$2" "$0}' | sort -k1,1 -T /scratch/gcwoodruff/pop_gen/elegans/23_site_annotation_genomics_general/tmp/) <(<03_elegans_24_og_pseudo_rad_codingSiteTypes.out_no_header_no_mito awk '{print $1"-"$2" "$0}' | sort -k1,1 -T /scratch/gcwoodruff/pop_gen/elegans/23_site_annotation_genomics_general/tmp2/) > 04_join_vcf_annotations_autosomes


#separate by degeneracy and codon position

mkdir $wkdir/elegans/24_awk_codon_sites

mkdir $wkdir/elegans/24_awk_codon_sites/00_degeneracy
mkdir $wkdir/elegans/24_awk_codon_sites/01_codon_position

cd $wkdir/elegans/24_awk_codon_sites/00_degeneracy

awk '{print > $39}'  $wkdir/elegans/23_site_annotation_genomics_general/04_join_vcf_annotations_autosomes


cd $wkdir/elegans/24_awk_codon_sites/01_codon_position

awk '{print > $37}'  $wkdir/elegans/23_site_annotation_genomics_general/04_join_vcf_annotations_autosomes


#clean up vcf's

mkdir $wkdir/elegans/25_awk_codon_sites

mkdir $wkdir/elegans/25_awk_codon_sites/00_degeneracy
mkdir $wkdir/elegans/25_awk_codon_sites/01_codon_position


cd $wkdir/elegans/24_awk_codon_sites/00_degeneracy

for i in *; do sed -i -e 's/ \+/\t/g' $i; done

#head -1 0 | sed 's/\t/\n/g'

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34} ' $i > $wkdir/elegans/25_awk_codon_sites/00_degeneracy/$i; done &


cd $wkdir/elegans/24_awk_codon_sites/01_codon_position

for i in *; do sed -i -e 's/ \+/\t/g' $i; done

#head -1 0 | sed 's/\t/\n/g'

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34} ' $i > $wkdir/elegans/25_awk_codon_sites/01_codon_position/$i; done &

#add vcf headers

mkdir $wkdir/elegans/26_cat_codon_sites

mkdir $wkdir/elegans/26_cat_codon_sites/00_degeneracy
mkdir $wkdir/elegans/26_cat_codon_sites/01_codon_position

cd $wkdir/elegans/20_cat_biallelic_invariant/

#grep "#" elegans_24_og_pseudo_rad_24_sorted.vcf > $wkdir/elegans/20_cat_biallelic_invariant/header_wgs
#I know the header is 47 lines

head -47 elegans_24_og_pseudo_rad_24_sorted.vcf > elegans_24_og_pseudo_rad_24_sorted_header

cd $wkdir/elegans/25_awk_codon_sites/00_degeneracy/

for i in *; do cat $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted_header $i >  $wkdir/elegans/26_cat_codon_sites/00_degeneracy/$i; done


cd $wkdir/elegans/25_awk_codon_sites/01_codon_position/

for i in *; do cat $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted_header $i >  $wkdir/elegans/26_cat_codon_sites/01_codon_position/$i; done

#rename

cd $wkdir/elegans/26_cat_codon_sites/00_degeneracy/

mv 0 0_degen_elegans.vcf
mv 2 2_degen_elegans.vcf
mv 4 4_degen_elegans.vcf

cd $wkdir/elegans/26_cat_codon_sites/01_codon_position/

mv 1 First_codon_pos_elegans.vcf
mv 2 Second_codon_pos_elegans.vcf
mv 3 Third_codon_pos_elegans.vcf

mkdir $wkdir/elegans/27_popgenwindows_codon_sites
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/01_popgenWindows_py
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py

#okay, we forgot to sort vcf files
cd $wkdir/elegans/26_cat_codon_sites/00_degeneracy/


cat 0_degen_elegans.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > 0_degen_elegans.vcf.tmp
mv 0_degen_elegans.vcf.tmp 0_degen_elegans.vcf

cat 2_degen_elegans.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > 2_degen_elegans.vcf.tmp
mv 2_degen_elegans.vcf.tmp 2_degen_elegans.vcf

cat 4_degen_elegans.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > 4_degen_elegans.vcf.tmp
mv 4_degen_elegans.vcf.tmp 4_degen_elegans.vcf

cd $wkdir/elegans/26_cat_codon_sites/01_codon_position/

cat First_codon_pos_elegans.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > First_codon_pos_elegans.vcf.tmp
mv First_codon_pos_elegans.vcf.tmp First_codon_pos_elegans.vcf

cat Second_codon_pos_elegans.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Second_codon_pos_elegans.vcf.tmp
mv Second_codon_pos_elegans.vcf.tmp Second_codon_pos_elegans.vcf

cat Third_codon_pos_elegans.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > Third_codon_pos_elegans.vcf.tmp
mv Third_codon_pos_elegans.vcf.tmp Third_codon_pos_elegans.vcf



#get pi in 10kb genomic windows across codon site types
#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/26_cat_codon_sites/00_degeneracy/0_degen_elegans.vcf | gzip > $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF/0_degen_elegans.geno.gz 

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF/0_degen_elegans.geno.gz  -o  $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/01_popgenWindows_py/0-fold_degen_elegans_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

gunzip 0-fold_degen_elegans_stats.csv.gz

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/26_cat_codon_sites/00_degeneracy/2_degen_elegans.vcf | gzip > $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF/2_degen_elegans.geno.gz 

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF/2_degen_elegans.geno.gz  -o  $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/01_popgenWindows_py/2-fold_degen_elegans_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/

gunzip 2-fold_degen_elegans_stats.csv.gz


#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/26_cat_codon_sites/00_degeneracy/4_degen_elegans.vcf | gzip > $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF/4_degen_elegans.geno.gz 

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/00_parseVCF/4_degen_elegans.geno.gz  -o  $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/01_popgenWindows_py/4-fold_degen_elegans_stats.csv.gz -f phased  

cd $wkdir/27_2022_site_annotation_continued/05_popgenwindows/01_popgenWindows_py/
gunzip 4-fold_degen_elegans_stats.csv.gz



#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/26_cat_codon_sites/01_codon_position/First_codon_pos_elegans.vcf | gzip > $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF/First_codon_pos_elegans.geno.gz 

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF/First_codon_pos_elegans.geno.gz  -o  $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py/First_codon_pos_elegans.csv.gz -f phased  

cd $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py/
gunzip First_codon_pos_elegans.csv.gz


#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/26_cat_codon_sites/01_codon_position/Second_codon_pos_elegans.vcf | gzip > $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF/Second_codon_pos_elegans.geno.gz 

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF/Second_codon_pos_elegans.geno.gz  -o  $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py/Second_codon_pos_elegans.csv.gz -f phased  

cd $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py/
gunzip Second_codon_pos_elegans.csv.gz



#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/26_cat_codon_sites/01_codon_position/Third_codon_pos_elegans.vcf | gzip > $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF/Third_codon_pos_elegans.geno.gz 

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/00_parseVCF/Third_codon_pos_elegans.geno.gz  -o  $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py/Third_codon_pos_elegans.csv.gz -f phased  

cd $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py/
gunzip Third_codon_pos_elegans.csv.gz
 
#add site type names to csv's


cd $wkdir/elegans/27_popgenwindows_codon_sites/00_degeneracy/01_popgenWindows_py
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk
mkdir $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/01_sed

awk 'BEGIN {FS=","} {OFS=","} {print $0,"0"} ' 0-fold_degen_elegans_stats.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/0-fold_degen_elegans_stats.csv

awk 'BEGIN {FS=","} {OFS=","} {print $0,"2"} ' 2-fold_degen_elegans_stats.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/2-fold_degen_elegans_stats.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"4"} ' 4-fold_degen_elegans_stats.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/4-fold_degen_elegans_stats.csv


cd $wkdir/elegans/27_popgenwindows_codon_sites/01_codon_position/01_popgenWindows_py
awk 'BEGIN {FS=","} {OFS=","} {print $0,"1"} ' First_codon_pos_elegans.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/First_codon_pos_elegans.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"2"} ' Second_codon_pos_elegans.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/Second_codon_pos_elegans.csv
awk 'BEGIN {FS=","} {OFS=","} {print $0,"3"} ' Third_codon_pos_elegans.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/Third_codon_pos_elegans.csv


#clean up and combine csv's with codon site types in comparable categories (fold-degeneracy, codon position)


mkdir $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/01_sed/

cd $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/00_awk/

for i in *; do sed '1d' $i > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/01_sed/$i; done

mkdir $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/02_cat/

cd $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/01_sed/

cat 0-fold_degen_elegans_stats.csv  2-fold_degen_elegans_stats.csv  4-fold_degen_elegans_stats.csv  > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/02_cat/degeneracy_elegans_pi.tmp.csv

cat First_codon_pos_elegans.csv  Second_codon_pos_elegans.csv  Third_codon_pos_elegans.csv > $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/02_cat/codon_position_elegans_pi.tmp.csv

cd $wkdir/elegans/27_popgenwindows_codon_sites/02_combine_popgenwindows/02_cat/

echo -e "chr,start,end,mid,sites,pi_all,fold_degeneracy" | cat - degeneracy_elegans_pi.tmp.csv > degeneracy_elegans_pi.csv


echo -e "chr,start,end,mid,sites,pi_all,codon_position" | cat - codon_position_elegans_pi.tmp.csv > codon_position_elegans_pi.csv


###ok, genic, intergenic, exonic, intronic [also see workflow for this in inopinata]....

cd $wkdir/elegans/elegans_genome

grep -v "#" caenorhabditis_elegans.PRJNA13758.WBPS16.annotations.gff3  > 00_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_no_pound.gff3

awk '$3 == "gene"' 00_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_no_pound.gff3 | sort -k1,1 -k4,4n > 01_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene.gff3

grep "biotype=protein_coding" 01_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene.gff3 > 01b_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene_protein-coding.gff3

awk '$3 == "exon"' 00_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_no_pound.gff3 | sort -k1,1 -k4,4n > 02_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_exon.gff3

#gff to bed

export PATH=$PATH:/home/gcwoodruff/download/bin

gff2bed < 01b_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene_protein-coding.gff3 > 03_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene.bed

gff2bed < 02_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_exon.gff3 > 04_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_exon.bed

#ok, get intergenic

bedtools complement -i 03_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene.bed -g elegans_chr_sizes.txt > 05_elegans_intergenic_sorted.bed

#ok, get intronic

awk 'OFS="\t" {print $1, $2, $3, "intergenic"}' 05_elegans_intergenic_sorted.bed > 06_elegans_intergenic_sort.bed

cat 04_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_exon.bed 06_elegans_intergenic_sort.bed | sort -k1,1 -k2,2n > 07_exon_intergenic_elegans_sort.bed

sed -i -e 's/exon.*/exon/g' 07_exon_intergenic_elegans_sort.bed

bedtools complement -i 07_exon_intergenic_elegans_sort.bed -g elegans_chr_sizes.txt > 08_elegans_intron_sorted.bed

awk 'OFS="\t" {print $1, $2, $3, "intron"}' 08_elegans_intron_sorted.bed > 09_elegans_intron_sorted.bed


#get intersection of vcf and genomic regions

mkdir $wkdir/elegans/28_bedtools_intersect_genomic_region/ 



bedtools intersect -a $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf  -b 03_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_gene.bed > $wkdir/elegans/28_bedtools_intersect_genomic_region//elegans_24_genic.vcf


bedtools intersect -a $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf  -b 06_elegans_intergenic_sort.bed > $wkdir/elegans/28_bedtools_intersect_genomic_region//elegans_24_intergenic.vcf


bedtools intersect -a $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf  -b 04_caenorhabditis_elegans.PRJNA13758.WBPS16.annotations_exon.bed > $wkdir/elegans/28_bedtools_intersect_genomic_region//elegans_24_exonic.vcf


bedtools intersect -a $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf  -b 09_elegans_intron_sorted.bed > $wkdir/elegans/28_bedtools_intersect_genomic_region//elegans_24_intronic.vcf


cd $wkdir/elegans/28_bedtools_intersect_genomic_region/

cat $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_genic.vcf > $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_genic.vcf.tmp
cat $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intergenic.vcf > $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intergenic.vcf.tmp
cat $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_exonic.vcf > $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_exonic.vcf.tmp
cat $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intronic.vcf > $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intronic.vcf.tmp

mv $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_genic.vcf.tmp $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_genic.vcf
mv $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intergenic.vcf.tmp $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intergenic.vcf
mv $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_exonic.vcf.tmp $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_exonic.vcf
mv $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intronic.vcf.tmp $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intronic.vcf


#estimate pi of genomic regions


mkdir $wkdir/elegans/29_popgenwindows_genomic_regions/
mkdir $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF
mkdir $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py

#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/VCF_processing/


python3 parseVCF.py -i  $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_genic.vcf | gzip > $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_genic.geno.gz 

python3 parseVCF.py -i  $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intergenic.vcf | gzip > $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_intergenic.geno.gz 

python3 parseVCF.py -i  $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_exonic.vcf | gzip > $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_exonic.geno.gz 

python3 parseVCF.py -i  $wkdir/elegans/28_bedtools_intersect_genomic_region/elegans_24_intronic.vcf | gzip > $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_intronic.geno.gz 


#cd to genomics_general folder
cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_genic.geno.gz -o  $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py/elegans_24_genic_stats.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_intergenic.geno.gz  -o  $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py/elegans_24_intergenic_stats.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_exonic.geno.gz  -o  $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py/elegans_24_exonic.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/29_popgenwindows_genomic_regions/00_parseVCF/elegans_24_intronic.geno.gz  -o  $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py/elegans_24_intronic.csv.gz -f phased  

cd $wkdir/elegans/29_popgenwindows_genomic_regions/01_popgenWindows_py/

gunzip *

	#files added to inopinata estimates in align_genotype_pop_gen.sh



#LD


#estimate LD

mkdir $wkdir/elegans/30_VCFTools_LD/
mkdir $wkdir/elegans/30_VCFTools_LD/00_grep_chr
mkdir $wkdir/elegans/30_VCFTools_LD/01_VCFTools_LD
mkdir $wkdir/elegans/30_VCFTools_LD/02_means_R

#split chromosomes to parallelize



cd $wkdir/elegans/20_cat_biallelic_invariant

grep -v "#" elegans_24_og_pseudo_rad_24_sorted.vcf > elegans_24_og_pseudo_rad_24_sorted_no_header

awk '$1 == "I"' elegans_24_og_pseudo_rad_24_sorted_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/I_no_header
awk '$1 == "II"' elegans_24_og_pseudo_rad_24_sorted_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/II_no_header
awk '$1 == "III"' elegans_24_og_pseudo_rad_24_sorted_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/III_no_header
awk '$1 == "IV"' elegans_24_og_pseudo_rad_24_sorted_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/IV_no_header
awk '$1 == "V"' elegans_24_og_pseudo_rad_24_sorted_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/V_no_header
awk '$1 == "X"' elegans_24_og_pseudo_rad_24_sorted_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/X_no_header


cat elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/30_VCFTools_LD/00_grep_chr/I_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/I.vcf
cat elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/30_VCFTools_LD/00_grep_chr/II_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/II.vcf
cat elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/30_VCFTools_LD/00_grep_chr/III_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/III.vcf
cat elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/30_VCFTools_LD/00_grep_chr/IV_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/IV.vcf
cat elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/30_VCFTools_LD/00_grep_chr/V_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/V.vcf
cat elegans_24_og_pseudo_rad_24_sorted_header $wkdir/elegans/30_VCFTools_LD/00_grep_chr/X_no_header > $wkdir/elegans/30_VCFTools_LD/00_grep_chr/X.vcf




module load VCFtools/0.1.16-foss-2017b-Perl-5.26.1


cd $wkdir/elegans/30_VCFTools_LD/01_VCFTools_LD/


vcftools --vcf $wkdir/elegans/30_VCFTools_LD/00_grep_chr/I.vcf --out I_out --geno-r2
vcftools --vcf $wkdir/elegans/30_VCFTools_LD/00_grep_chr/II.vcf --out II_out --geno-r2
vcftools --vcf $wkdir/elegans/30_VCFTools_LD/00_grep_chr/III.vcf --out III_out --geno-r2
vcftools --vcf $wkdir/elegans/30_VCFTools_LD/00_grep_chr/IV.vcf --out IV_out --geno-r2
vcftools --vcf $wkdir/elegans/30_VCFTools_LD/00_grep_chr/V.vcf --out V_out --geno-r2
vcftools --vcf $wkdir/elegans/30_VCFTools_LD/00_grep_chr/X.vcf --out X_out --geno-r2


#turn into bed
mkdir $wkdir/elegans/30_VCFTools_LD/03_bed/

cd $wkdir/elegans/30_VCFTools_LD/02_means_R/


for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $3,$1-1,$1, $2}' $i > $wkdir/elegans/30_VCFTools_LD/03_bed/$i; done


cd $wkdir/elegans/30_VCFTools_LD/03_bed/

cat * > elegans_LD.bed
#sort bed
sort -V -k1,1 -k2,2 elegans_LD.bed > elegans_LD_sort.bed

#bedtools windows...

mkdir $wkdir/elegans/30_VCFTools_LD/04_bedtools_windows/


cd $wkdir/elegans/30_VCFTools_LD/03_bed/

bedtools map -o mean -c 4 -a $wkdir/additional_files/elegans.50bp.windows -b elegans_LD_sort.bed > $wkdir/elegans/30_VCFTools_LD/04_bedtools_windows/00_elegans_LD_50bp_windows

cd $wkdir/elegans/30_VCFTools_LD/04_bedtools_windows/

#remove missing data
awk 'BEGIN {OFS="\t"} $4 != "." {print $0}' 00_elegans_LD_50bp_windows > 01_elegans_LD_50bp_windows_no_missing_data

#now 10kb windows

bedtools map -o mean -c 4 -a $wkdir/additional_files/elegans.10kb.windows -b 01_elegans_LD_50bp_windows_no_missing_data > 02_elegans_LD_10kb_windows


#remove missing data

awk 'BEGIN {OFS="\t"} $4 != "." {print $1,$2,$4}' 02_elegans_LD_10kb_windows > 03_elegans_LD_10kb_windows_no_missing_data


#add elegans

awk 'BEGIN {OFS="\t"} {print $1,$2+1,$3, "C. elegans"}' 03_elegans_LD_10kb_windows_no_missing_data > 04_elegans_LD_10kb_windows_no_missing_data_species
