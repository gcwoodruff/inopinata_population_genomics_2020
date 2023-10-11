cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/

mkdir pop_gen_revisions_2023

cd pop_gen_revisions_2023

mkdir elegans

wkdir='/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023'

cp -r /scratch/gcwoodruff/pop_gen_transfer_may_2023/pop_gen_oscer_transfer_4-28-22/pop_gen_4-2022/18_bcftools_view_biallelic_snps_maf/ /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/


#ok, vcf to bed to eventually generate 1kb SNP windows

export PATH=$PATH:/home/gcwoodruff/download/bin/

mkdir $wkdir/revisions_hyperdivergent

mkdir $wkdir/revisions_hyperdivergent/00_vcf2bed

cd $wkdir/18_bcftools_view_biallelic_snps_maf/

vcf2bed < inopinata_24.vcf > $wkdir/revisions_hyperdivergent/00_vcf2bed/inopinata_24_billalelic_snps.bed

cd $wkdir/elegans/19_bcftools_view_biallelic_snps_maf/


vcf2bed < elegans_24_og_pseudo_rad.vcf > $wkdir/revisions_hyperdivergent/00_vcf2bed/elegans_24_og_pseudo_rad_billalelic_snps.bed


vcf2bed < elegans_24_wgs.vcf > $wkdir/revisions_hyperdivergent/00_vcf2bed/elegans_24_wgs_billalelic_snps.bed

#okay, reduce and add a 1, get ready for some bedtools

mkdir $wkdir/revisions_hyperdivergent/01_awk

cd $wkdir/revisions_hyperdivergent/00_vcf2bed/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,"1", "C. inopinata RAD"}' inopinata_24_billalelic_snps.bed > $wkdir/revisions_hyperdivergent/01_awk/inopinata_24_billalelic_snps.bed


awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,"1", "C. elegans og_pseudo_RAD"}' elegans_24_og_pseudo_rad_billalelic_snps.bed > $wkdir/revisions_hyperdivergent/01_awk/elegans_24_og_pseudo_rad_billalelic_snps.bed

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3,"1", "C. elegans WGS"}' elegans_24_wgs_billalelic_snps.bed > $wkdir/revisions_hyperdivergent/01_awk/elegans_24_wgs_billalelic_snps.bed

#ok, bedtools....

mkdir $wkdir/additional_files

#make 1 kb windows
module load BEDTools/2.27.1-foss-2018b

cd $wkdir/additional_files/inopinata_genome/

bedtools makewindows -g inopinata_chr_sizes.txt -w 1000 > inopinata_1kb_windows.bed

cd $wkdir/additional_files/elegans_genome/

bedtools makewindows -g elegans_chr_sizes.txt -w 1000 > elegans_1kb_windows.bed

#okay, now, let's do the adding up of snps across 1 kb windows

mkdir $wkdir/revisions_hyperdivergent/02_bedtools_map/

cd $wkdir/revisions_hyperdivergent/01_awk/

bedtools map -a $wkdir/additional_files/inopinata_genome/inopinata_1kb_windows.bed -b inopinata_24_billalelic_snps.bed -c 4 -o sum > $wkdir/revisions_hyperdivergent/02_bedtools_map/inopinata_24_num_SNPs_1kb_windows.bed


bedtools map -a $wkdir/additional_files/elegans_genome/elegans_1kb_windows.bed -b elegans_24_og_pseudo_rad_billalelic_snps.bed -c 4 -o sum > $wkdir/revisions_hyperdivergent/02_bedtools_map/elegans_24_og_pseudo_rad_num_SNPs_1kb_windows.bed


bedtools map -a $wkdir/additional_files/elegans_genome/elegans_1kb_windows.bed -b elegans_24_wgs_billalelic_snps.bed -c 4 -o sum > $wkdir/revisions_hyperdivergent/02_bedtools_map/elegans_24_wgs_num_SNPs_1kb_windows.bed


#replace "." with NA

mkdir $wkdir/revisions_hyperdivergent/03_sed/

cd  $wkdir/revisions_hyperdivergent/02_bedtools_map/

for i in *; do sed 's/\./NA/g' $i > $wkdir/revisions_hyperdivergent/03_sed/$i; done

scp -r gcwoodruff@dtn2.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/revisions_hyperdivergent/03_sed/ /Users/gavin/genome/pop_gen_revisions_9-2023/














