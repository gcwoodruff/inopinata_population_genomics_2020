
#here, getting the number of sites across various genomic regions, chromosomes, and chromosome arms


wkdir='/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023'

mkdir $wkdir/revisions_number_of_sites_regions


mkdir $wkdir/revisions_number_of_sites_regions/00_vcf/

module load BCFtools/1.11-GCC-10.2.0


#inopinata all sites, /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_June_2023/19_cat_biallelic_invariant/inopinata_24_sorted.vcf

#just get the biallelic inopinata sites
bcftools view --types snps -m 2 -M 2 /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_June_2023/19_cat_biallelic_invariant/inopinata_24_sorted.vcf > $wkdir/revisions_number_of_sites_regions/00_vcf/inopinata_24_sorted_biallelic.vcf

#elegans all sites, /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf

#just get the biallelic inopinata sites
bcftools view --types snps -m 2 -M 2 /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf > $wkdir/revisions_number_of_sites_regions/00_vcf/elegans_24_og_pseudo_rad_24_sorted_biallelic.vcf

#ok, let's remove some headers (biallelic sites)

mkdir $wkdir/revisions_number_of_sites_regions/01_grep/

cd $wkdir/revisions_number_of_sites_regions/00_vcf/

grep -v "#" inopinata_24_sorted_biallelic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/inopinata_24_sorted_biallelic.txt

grep -v "#" elegans_24_og_pseudo_rad_24_sorted_biallelic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/elegans_24_og_pseudo_rad_24_sorted_biallelic.txt

#all sites, remove headers


grep -v "#" /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_June_2023/19_cat_biallelic_invariant/inopinata_24_sorted.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/inopinata_24_sorted_all_sites.txt

grep -v "#" /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/elegans_24_og_pseudo_rad_24_sorted_all_sites.txt

#inopinata codon sites, /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_27_site_annotation_continued/27_2022_site_annotation_continued/04_sort/


cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_27_site_annotation_continued/27_2022_site_annotation_continued/04_sort/autosomes/codon_position

grep -v "#" First_codon_pos_ino_aut.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/First_codon_pos_ino_aut.txt

grep -v "#" Second_codon_pos_ino_aut.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/Second_codon_pos_ino_aut.txt

grep -v "#" Third_codon_pos_ino_aut.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/Third_codon_pos_ino_aut.txt

cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_27_site_annotation_continued/27_2022_site_annotation_continued/04_sort/autosomes/degneracy

grep -v "#" 0_degen_ino_aut.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/0_degen_ino_aut.txt

grep -v "#" 2_degen_ino_aut.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/2_degen_ino_aut.txt

grep -v "#" 4_degen_ino_aut.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/4_degen_ino_aut.txt


cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_27_site_annotation_continued/27_2022_site_annotation_continued/04_sort/X/codon_position

grep -v "#" First_codon_pos_ino_X.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/First_codon_pos_ino_X.txt

grep -v "#" Second_codon_pos_ino_X.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/Second_codon_pos_ino_X.txt

grep -v "#" Third_codon_pos_ino_X.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/Third_codon_pos_ino_X.txt


cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen_27_site_annotation_continued/27_2022_site_annotation_continued/04_sort/X/degneracy


grep -v "#" 0_degen_ino_X.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/0_degen_ino_X.txt

grep -v "#" 2_degen_ino_X.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/2_degen_ino_X.txt

grep -v "#" 4_degen_ino_X.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/4_degen_ino_X.txt


#inopinata genomic regions sites, /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/29_bedtools_intersect_genomic_region

cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/29_bedtools_intersect_genomic_region

grep -v "#" inopinata_24_exonic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/inopinata_24_exonic.txt

grep -v "#" inopinata_24_genic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/inopinata_24_genic.txt

grep -v "#" inopinata_24_intergenic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/inopinata_24_intergenic.txt

grep -v "#" inopinata_24_intronic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/inopinata_24_intronic.txt

#elegans pseudo-rad codon sites, /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/26_cat_codon_sites

cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/26_cat_codon_sites/00_degeneracy

grep -v "#" 0_degen_elegans.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/0_degen_elegans.txt

grep -v "#" 2_degen_elegans.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/2_degen_elegans.txt

grep -v "#" 4_degen_elegans.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/4_degen_elegans.txt


cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/26_cat_codon_sites/01_codon_position

grep -v "#" First_codon_pos_elegans.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/First_codon_pos_elegans.txt

grep -v "#" Second_codon_pos_elegans.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/Second_codon_pos_elegans.txt

grep -v "#" Third_codon_pos_elegans.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/Third_codon_pos_elegans.txt



#elegans genomic regions sites, /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/28_bedtools_intersect_genomic_region

cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/SCRATCH_10-2023/gcwoodruff/pop_gen/elegans/28_bedtools_intersect_genomic_region

grep -v "#" elegans_24_exonic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/elegans_24_exonic.txt

grep -v "#" elegans_24_genic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/elegans_24_genic.txt

grep -v "#" elegans_24_intergenic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/elegans_24_intergenic.txt

grep -v "#" elegans_24_intronic.vcf > $wkdir/revisions_number_of_sites_regions/01_grep/elegans_24_intronic.txt

#ok, put some autosomes and X's together

cd $wkdir/revisions_number_of_sites_regions/01_grep/

cat First_codon_pos_ino_aut.txt First_codon_pos_ino_X.txt > First_codon_pos_ino.txt
cat Second_codon_pos_ino_aut.txt Second_codon_pos_ino_X.txt > Second_codon_pos_ino.txt
cat Third_codon_pos_ino_aut.txt Third_codon_pos_ino_X.txt > Third_codon_pos_ino.txt
cat 0_degen_ino_aut.txt 0_degen_ino_X.txt > 0_degen_ino.txt
cat 2_degen_ino_aut.txt 2_degen_ino_X.txt > 2_degen_ino.txt
cat 4_degen_ino_aut.txt 4_degen_ino_X.txt > 4_degen_ino.txt

rm First_codon_pos_ino_aut.txt
rm First_codon_pos_ino_X.txt
rm Second_codon_pos_ino_aut.txt
rm Second_codon_pos_ino_X.txt
rm Third_codon_pos_ino_aut.txt
rm Third_codon_pos_ino_X.txt
rm 0_degen_ino_aut.txt
rm 0_degen_ino_X.txt
rm 2_degen_ino_aut.txt
rm 2_degen_ino_X.txt
rm 4_degen_ino_aut.txt
rm 4_degen_ino_X.txt

#ok, get just the invariant sites

awk '$5 == "."' elegans_24_og_pseudo_rad_24_sorted_all_sites.txt > elegans_24_og_pseudo_rad_24_sorted_invariant_sites.txt

awk '$5 == "."' inopinata_24_sorted_all_sites.txt > inopinata_24_sorted_invariant_sites.txt

#ok, add species columns and data type columns

mkdir $wkdir/revisions_number_of_sites_regions/02_awk_species/

cd $wkdir/revisions_number_of_sites_regions/01_grep/


awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "0-fold degenerate"}' 0_degen_elegans.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/0_degen_elegans.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "0-fold degenerate"}' 0_degen_ino.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/0_degen_ino.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "2-fold degenerate"}' 2_degen_elegans.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/2_degen_elegans.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "2-fold degenerate"}' 2_degen_ino.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/2_degen_ino.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "4-fold degenerate"}' 4_degen_elegans.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/4_degen_elegans.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "4-fold degenerate"}' 4_degen_ino.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/4_degen_ino.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "exon"}' elegans_24_exonic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_exonic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "genic"}' elegans_24_genic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_genic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "intergenic"}' elegans_24_intergenic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_intergenic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "inron"}' elegans_24_intronic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_intronic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "all sites"}' elegans_24_og_pseudo_rad_24_sorted_all_sites.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_og_pseudo_rad_24_sorted_all_sites.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "variant sites"}' elegans_24_og_pseudo_rad_24_sorted_biallelic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_og_pseudo_rad_24_sorted_biallelic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "invariant sites"}' elegans_24_og_pseudo_rad_24_sorted_invariant_sites.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/elegans_24_og_pseudo_rad_24_sorted_invariant_sites.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "first codon position"}' First_codon_pos_elegans.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/First_codon_pos_elegans.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "first codon position"}' First_codon_pos_ino.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/First_codon_pos_ino.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "exon"}' inopinata_24_exonic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_exonic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "genic"}' inopinata_24_genic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_genic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "intergenic"}' inopinata_24_intergenic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_intergenic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "intron"}' inopinata_24_intronic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_intronic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "all sites"}' inopinata_24_sorted_all_sites.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_sorted_all_sites.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "variant sites"}' inopinata_24_sorted_biallelic.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_sorted_biallelic.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "invariant sites"}' inopinata_24_sorted_invariant_sites.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/inopinata_24_sorted_invariant_sites.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "second codon position"}' Second_codon_pos_elegans.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/Second_codon_pos_elegans.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "second codon position"}' Second_codon_pos_ino.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/Second_codon_pos_ino.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. elegans", "third codon position"}' Third_codon_pos_elegans.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/Third_codon_pos_elegans.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $2, "C. inopinata", "third codon position"}' Third_codon_pos_ino.txt > $wkdir/revisions_number_of_sites_regions/02_awk_species/Third_codon_pos_ino.txt

#sort unique

mkdir $wkdir/revisions_number_of_sites_regions/03_sort_uniq/


cd $wkdir/revisions_number_of_sites_regions/02_awk_species/



#cat

mkdir $wkdir/revisions_number_of_sites_regions/03_cat/

cd $wkdir/revisions_number_of_sites_regions/02_awk_species/

cat * > $wkdir/revisions_number_of_sites_regions/03_cat/site_positions.tsv


#okay, R next, see file "number_of_sites.R"


