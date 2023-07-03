#getting coverage info for all sites, all samples, before and after genotyping/filtering

#put working directory here
wkdir="/projects/phillipslab/gavincw/pop_gen_reproduce_june_2021"

mkdir $wkdir/24_coverage

#coverage after alignment, and after filtering for unique alignments, but before genotype filtering

mkdir $wkdir/24_coverage/00_before_genotype_filtering

mkdir $wkdir/24_coverage/00_before_genotype_filtering/00_samtools_depth

cd $wkdir/08_samtools_sort

for i in *; do samtools depth $i > $wkdir/24_coverage/00_before_genotype_filtering/00_samtools_depth/$i; done &


#make bed file

mkdir $wkdir/24_coverage/00_before_genotype_filtering/01_awk_bed

cd $wkdir/24_coverage/00_before_genotype_filtering/00_samtools_depth/

for i in *; do awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$2,$3}' $i > $wkdir/24_coverage/00_before_genotype_filtering/01_awk_bed/$i; done &



#ok, sort bed file

mkdir $wkdir/24_coverage/00_before_genotype_filtering/02_depth_bedtools_sort

#module load bedtools/2.25.0

cd $wkdir/24_coverage/00_before_genotype_filtering/01_awk_bed

for i in *; do bedtools sort -i $i  > $wkdir/24_coverage/00_before_genotype_filtering/02_depth_bedtools_sort/$i; done &

#get average coverages for windows

mkdir $wkdir/24_coverage/00_before_genotype_filtering/03_depth_bedtools_map

cd $wkdir/24_coverage/00_before_genotype_filtering/02_depth_bedtools_sort/

for i in *; do bedtools map -o mean -c 4 -a $wkdir/additional_files/inopinata.10kb.windows -b $i > $wkdir/24_coverage/00_before_genotype_filtering/03_depth_bedtools_map/$i; done &

#replace "." with "0" ; ie, missing data is zero coverage

cd $wkdir/24_coverage/00_before_genotype_filtering/03_depth_bedtools_map/

for i in *; do sed -i -e 's/\.$/0/g' $i; done

#label rows with bam file ID

mkdir $wkdir/24_coverage/00_before_genotype_filtering/04_depth_awk
cd $wkdir/24_coverage/00_before_genotype_filtering/03_depth_bedtools_map/

for i in *; do awk 'BEGIN{FS="\t";OFS="\t"} {print $0,FILENAME}' $i > $wkdir/24_coverage/00_before_genotype_filtering/04_depth_awk/$i; done

#combine


mkdir $wkdir/24_coverage/00_before_genotype_filtering/05_depth_cat

cd $wkdir/24_coverage/00_before_genotype_filtering/04_depth_awk_div

cat * > $wkdir/24_coverage/00_before_genotype_filtering/05_depth_cat/all_depth

cd $wkdir/24_coverage/00_before_genotype_filtering/05_depth_cat/

echo -e "Chr\tBP_start\tBP_end\tCov\tsample" | cat - all_depth > all_depth.tmp

mv all_depth.tmp all_depth

mv all_depth all_depth.tsv



#get EcoRI cut sites in inopinata genome assembly

mkdir $wkdir/24_coverage/01_inopinata_EcoRI_cut_sites/

fuzznuc -sequence $wkdir/inopinata_genome/inopinata_genome.fa -pattern 'GAATTC' -rformat gff -outfile $wkdir/24_coverage/01_inopinata_EcoRI_cut_sites/inopinata_EcoRI_sites_out.gff

#gff to bed

cd $wkdir/24_coverage/01_inopinata_EcoRI_cut_sites/inopinata_EcoRI_sites_out.gff

grep -v "#" inopinata_EcoRI_sites_out.gff | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$4,$5,1}' > inopinata_EcoRI_sites.bed

#number of cut sites per 10kb window

bedtools map -o sum -c 4 -a $wkdir/additional_files/inopinata.10kb.windows -b inopinata_EcoRI_sites.bed > inopinata_EcoRI_sites.10kb_windows.bed

echo -e "Chr\tBP_start\tBP_end\tEcoRI_cut_sites" | cat - inopinata_EcoRI_sites.10kb_windows.bed > inopinata_EcoRI_sites.tsv



#coverage after genotype filtering


mkdir $wkdir/24_coverage/02_after_genotype_filtering

mkdir $wkdir/24_coverage/02_after_genotype_filtering/00_bcftools_query_DP

#getting the proper file, the MAF filtered, properly called X, with invariant sites, including males and females, $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf

cd $wkdir/19_cat_biallelic_invariant

#use bedtools to get the needed features from the vcf

#get site and depth info. #i am using this file and not the vcf files split by sex and chromosome because I want to get site coverage of all males and females and chromosomes following coverage and MAF filters.

bcftools query -f '%CHROM\t%POS\t%DP\n' $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf > $wkdir/24_coverage/02_after_genotype_filtering/00_bcftools_query_DP/inopinata_24_sorted.tsv &


#extract site and genotype info
mkdir $wkdir/24_coverage/02_after_genotype_filtering/01_awk/

cd $wkdir/19_cat_biallelic_invariant/

grep -v "#" $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33}' > $wkdir/24_coverage/02_after_genotype_filtering/01_awk/inopinata_24_sorted.tsv &



#extract genotype info
mkdir $wkdir/24_coverage/02_after_genotype_filtering/02_awk/


grep -v "#" $wkdir/19_cat_biallelic_invariant/inopinata_24_sorted.vcf | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33}' > $wkdir/24_coverage/02_after_genotype_filtering/02_awk/inopinata_24_sorted.tsv &


#get just depth
mkdir $wkdir/24_coverage/02_after_genotype_filtering/03_awk/

cd $wkdir/24_coverage/02_after_genotype_filtering/00_bcftools_query_DP/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $3}' inopinata_24_sorted.tsv > $wkdir/24_coverage/02_after_genotype_filtering/03_awk/inopinata_24_sorted.tsv




#replace commas with tabs
cd $wkdir/24_coverage/02_after_genotype_filtering/03_awk/

sed -i -e 's/,/\t/g' inopinata_24_sorted.tsv
#add header
cd  $wkdir/24_coverage/02_after_genotype_filtering/02_awk/

echo -e "C1\tC2\tC3\tC4\tC5\tC6\tC7\tC8\tC9\tC10\tC11\tC12\tC13\tC14\tC15\tC16\tC17\tC18\tC19\tC20\tC21\tC22\tC23\tC24" | cat - inopinata_24_sorted.tsv > inopinata_24_sorted_gt.tsv

#add header
cd  $wkdir/24_coverage/02_after_genotype_filtering/03_awk/

echo -e "C1\tC2\tC3\tC4\tC5\tC6\tC7\tC8\tC9\tC10\tC11\tC12\tC13\tC14\tC15\tC16\tC17\tC18\tC19\tC20\tC21\tC22\tC23\tC24" | cat - inopinata_24_sorted.tsv > inopinata_24_sorted_dp.tsv




#split to parallelize


mkdir $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt

mkdir $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp


#genotypes
cd $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/

split -l 100000 -d $wkdir/24_coverage/02_after_genotype_filtering/02_awk/inopinata_24_sorted.tsv gt_

for i in *; do echo -e "C1\tC2\tC3\tC4\tC5\tC6\tC7\tC8\tC9\tC10\tC11\tC12\tC13\tC14\tC15\tC16\tC17\tC18\tC19\tC20\tC21\tC22\tC23\tC24" | cat - $i > $i.tmp; done

for i in *.tmp; do mv -- "$i" "${i%.tmp}"; done

#depth
cd $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/

split -l 100000 -d $wkdir/24_coverage/02_after_genotype_filtering/03_awk/inopinata_24_sorted.tsv dp_

for i in *; do echo -e "C1\tC2\tC3\tC4\tC5\tC6\tC7\tC8\tC9\tC10\tC11\tC12\tC13\tC14\tC15\tC16\tC17\tC18\tC19\tC20\tC21\tC22\tC23\tC24" | cat - $i > $i.tmp; done

for i in *.tmp; do mv -- "$i" "${i%.tmp}"; done


mkdir $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/

#use script fix_vcf_coverage_args.R in parallel to extract usable site coverage. 
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_00 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_00 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_00
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_01 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_01 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_01
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_02 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_02 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_02
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_03 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_03 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_03
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_04 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_04 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_04
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_05 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_05 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_05
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_06 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_06 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_06
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_07 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_07 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_07
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_08 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_08 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_08
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_09 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_09 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_09
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_10 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_10 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_10
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_11 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_11 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_11
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_12 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_12 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_12
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_13 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_13 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_13
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_14 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_14 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_14
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_15 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_15 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_15
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_16 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_16 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_16
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_17 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_17 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_17
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_18 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_18 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_18
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_19 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_19 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_19
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_20 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_20 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_20
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_21 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_21 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_21
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_22 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_22 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_22
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_23 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_23 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_23
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_24 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_24 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_24
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_25 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_25 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_25
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_26 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_26 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_26
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_27 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_27 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_27
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_28 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_28 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_28
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_29 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_29 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_29
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_30 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_30 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_30
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_31 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_31 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_31
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_32 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_32 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_32
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_33 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_33 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_33
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_34 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_34 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_34
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_35 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_35 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_35
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_36 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_36 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_36
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_37 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_37 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_37
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_38 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_38 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_38
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_39 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_39 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_39
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_40 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_40 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_40
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_41 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_41 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_41
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_42 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_42 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_42
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_43 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_43 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_43
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_44 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_44 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_44
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_45 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_45 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_45
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_46 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_46 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_46
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_47 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_47 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_47
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_48 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_48 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_48
Rscript $wkdir/additional_files/fix_vcf_coverage_args.R $wkdir/24_coverage/02_after_genotype_filtering/04_split_gt/gt_49 $wkdir/24_coverage/02_after_genotype_filtering/05_split_dp/dp_49 $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/out_49



#remove headers
cd $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/

for i in *; do sed -i '1d' $i; done &

#combine 
mkdir $wkdir/24_coverage/02_after_genotype_filtering/07_cat/

cd $wkdir/24_coverage/02_after_genotype_filtering/06_fix_vcf_coverage/

cat * > $wkdir/24_coverage/02_after_genotype_filtering/07_cat/site_coverage.tsv


#get site positions
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2}' $wkdir/24_coverage/02_after_genotype_filtering/00_bcftools_query_DP/inopinata_24_sorted.tsv > $wkdir/24_coverage/02_after_genotype_filtering/07_cat/sites.tsv


#paste site positions with depth
cd $wkdir/24_coverage/02_after_genotype_filtering/07_cat/

paste sites.tsv site_coverage.tsv > site_coverage.tsv.tmp

mv site_coverage.tsv.tmp site_coverage.tsv


#add header
echo -e "Chr\tBP\tA03\tA05\tB03\tdec_2016_D01\tdec_2016_D02\tdec_2016_D03\tdec_2016_D04\tdec_2016_D05\tdec_2016_D10\tmay_2017_A11\tmay_2017_A12\tmay_2017_B11\tmay_2017_B12\tmay_2017_C11\tmay_2017_C12\tmay_2017_D11\tmay_2017_D12\tmay_2017_E11\tmay_2017_E12\tmay_2017_F11\tmay_2017_G12\tmay_2017_H10\tmay_2017_H11\tmay_2017_H12" | cat - site_coverage.tsv > site_coverage.tsv.tmp

mv site_coverage.tsv.tmp site_coverage.tsv













