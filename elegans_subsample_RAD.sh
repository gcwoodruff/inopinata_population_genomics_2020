
#elegans-- how does not getting all RAD loci impact things?

#start with 30 bp loci

wkdir="/scratch/gcwoodruff/pop_gen/"

#get the sites, and get just the unique ones

cd $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3}' elegans_EcoRI_sites_out.bed > elegans_EcoRI_sites_out_col_1-3.bed

sort elegans_EcoRI_sites_out_col_1-3.bed | uniq | sort -V -k1,1 -k2,2 > elegans_EcoRI_sites_out_col_1-3_uniq.bed

#ok, get 30 bp loci

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2-12,$3+12}' elegans_EcoRI_sites_out_col_1-3_uniq.bed > elegans_EcoRI_sites_30_bp_loci.bed


#are they 30 bp

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,$3-$2}' elegans_EcoRI_sites_30_bp_loci.bed | head -50
#cool, looks like it worked.
#okay, let's get some loci from the vcfs......

#extract loci

mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/

mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites

cp $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/100_sites.bed

wc -l $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed

#45935

#0.1	4594
#0.2	9187
#0.3	13781
#0.4	18374
#0.5	22968
#0.6	27561
#0.7	32154
#0.8	36748
#0.9	41342

shuf -n 4594 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/10_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/10_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/10_sites.bed

wc -l $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/10_sites.bed
#4594
wc -l $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/10_sites.bed.tmp
#4594

#cool it seems to work


shuf -n 9187 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/20_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/20_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/20_sites.bed




shuf -n 13781 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/30_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/30_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/30_sites.bed



shuf -n 18374 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/40_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/40_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/40_sites.bed


shuf -n 22968 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/50_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/50_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/50_sites.bed



shuf -n 27561 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/60_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/60_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/60_sites.bed



shuf -n 32154 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/70_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/70_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/70_sites.bed




shuf -n 36748 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/80_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/80_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/80_sites.bed



shuf -n 41342 $wkdir/elegans/07_b_emboss_fuzznuc/00_emboss_fuzznuc/elegans_EcoRI_sites_30_bp_loci.bed > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/90_sites.bed.tmp

sort $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/90_sites.bed.tmp  | uniq | sort -V -k1,1 -k2,2 > $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/90_sites.bed


cd $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/

rm *.bed.tmp


#yay, we got some sites, okay.
#get overlapping sites in VCF!

mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/


module load VCFtools/0.1.16-foss-2017b-Perl-5.26.1

cd $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/

vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/100_sites.bed --out 100_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/90_sites.bed --out 90_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/80_sites.bed --out 80_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/70_sites.bed --out 70_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/60_sites.bed --out 60_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/50_sites.bed --out 50_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/40_sites.bed --out 40_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/30_sites.bed --out 30_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/20_sites.bed --out 20_sites --recode --keep-INFO-all
vcftools --vcf $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf --bed $wkdir/elegans/31_get_fraction_EcoRI_regions/00_get_sites/10_sites.bed --out 10_sites --recode --keep-INFO-all

wkdir="/scratch/gcwoodruff/pop_gen/"

mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/

mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf
mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py


module load Biopython/1.78-foss-2020b

cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/100_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/100_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/90_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/90_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/80_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/80_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/70_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/70_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/60_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/60_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/50_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/50_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/40_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/40_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/30_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/30_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/20_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/20_sites.geno.gz
python3 parseVCF.py -i  $wkdir/elegans/31_get_fraction_EcoRI_regions/01_vcftools_get_sites/10_sites.recode.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/10_sites.geno.gz


cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/100_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/100_sites.csv.gz -f phased  

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/90_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/90_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/80_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/80_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/70_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/700_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/60_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/60_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/50_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/50_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/40_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/40_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/30_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/30_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/20_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/20_sites.csv.gz -f phased  
python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/10_sites.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/10_sites.csv.gz -f phased  

#also wgs


cd /home/gcwoodruff/download/genomics_general/VCF_processing/

python3 parseVCF.py -i  $wkdir/elegans/20_cat_biallelic_invariant/elegans_24_wgs_24_sorted.vcf | gzip > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/elegans_24_wgs.geno.gz

cd /home/gcwoodruff/download/genomics_general/

python3 popgenWindows.py --windType coordinate -w 10000 -s 10000 -m 160 -g $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/00_parse_vcf/elegans_24_wgs.geno.gz -o  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/elegans_24_wgs.csv.gz -f phased  


#okay, prep for R

#add 

mkdir $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/

cd $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/01_popgenWindows_py/

gunzip *

awk 'BEGIN {OFS=","} {print $0, "100%"}' 100_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/100_sites.csv
awk 'BEGIN {OFS=","} {print $0, "10%"}' 10_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/10_sites.csv
awk 'BEGIN {OFS=","} {print $0, "20%"}' 20_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/20_sites.csv
awk 'BEGIN {OFS=","} {print $0, "30%"}' 30_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/30_sites.csv
awk 'BEGIN {OFS=","} {print $0, "40%"}' 40_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/40_sites.csv
awk 'BEGIN {OFS=","} {print $0, "50%"}' 50_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/50_sites.csv
awk 'BEGIN {OFS=","} {print $0, "60%"}' 60_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/60_sites.csv
awk 'BEGIN {OFS=","} {print $0, "70%"}' 70_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/70_sites.csv
awk 'BEGIN {OFS=","} {print $0, "80%"}' 80_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/80_sites.csv
awk 'BEGIN {OFS=","} {print $0, "90%"}' 90_sites.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/90_sites.csv
awk 'BEGIN {OFS=","} {print $0, "WGS"}' elegans_24_wgs.csv > $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/elegans_24_wgs.csv

cd $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/02_awk/

for i in *; do sed -i '1d' $i; done

cat * > elegans_random_rad_pi.tmp

echo -e "Chr,start,end,mid,sites,Pi,perc_RAD" | cat -  elegans_random_rad_pi.tmp >  $wkdir/elegans/31_get_fraction_EcoRI_regions/02_popgenwindows/elegans_random_rad_pi.csv
