#convert vcf to bed

vcf2bed < /Volumes/BIGG2/oscer_transfer_11-16-2022/gcwoodruff/pop_gen/elegans/20_cat_biallelic_invariant/elegans_24_og_pseudo_rad_24_sorted.vcf > /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted.bed

#just get the relevant columns
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3}' 20220216_c_elegans_divergent_regions_strain.bed > 20220216_c_elegans_divergent_regions_strain_first_three_columns.bed

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2,$3}' elegans_24_og_pseudo_rad_24_sorted.bed > elegans_24_og_pseudo_rad_24_sorted_first_three_columns.bed

#sort uniq

sort 20220216_c_elegans_divergent_regions_strain_first_three_columns.bed | uniq > 20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq.bed

bedtools sort -i 20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq.bed > 20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq_bedtools_sort.bed



#merge

bedtools merge -i elegans_24_og_pseudo_rad_24_sorted_first_three_columns.bed > elegans_24_og_pseudo_rad_24_sorted_first_three_columns_merged.bed

bedtools merge -i 20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq_bedtools_sort.bed > 20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq_bedtools_sort_merge.bed

#get overlap of my vcf with the reported hyperdivergent regions with bed
# -wa hyperdivergent -b my variants

#contiguous regions
bedtools intersect -wa -a /Users/gavin/genome/pop_gen_revisions_9-2023/20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq_bedtools_sort_merge.bed -b /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted_first_three_columns_merged.bed > /Users/gavin/genome/pop_gen_revisions_9-2023/hyperdivergent_regions_with_elegans_rad_sites.bed

#sites
bedtools intersect -a /Users/gavin/genome/pop_gen_revisions_9-2023/20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq_bedtools_sort_merge.bed -b /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted_first_three_columns.bed > /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted_first_three_columns_number_of_sites_in_regions.bed

#remove duplicate rows

sort /Users/gavin/genome/pop_gen_revisions_9-2023/hyperdivergent_regions_with_elegans_rad_sites.bed | uniq > /Users/gavin/genome/pop_gen_revisions_9-2023/hyperdivergent_regions_with_elegans_rad_sites_sort_uniq.bed

sort /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted_first_three_columns_number_of_sites_in_regions.bed | uniq > /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted_first_three_columns_number_of_sites_in_regions_sort_uniq.bed


#how many hyperdivergent regions?

wc -l 20220216_c_elegans_divergent_regions_strain_first_three_columns_sort_uniq_bedtools_sort_merge.bed
#312 

#how many overlap with pseudo-rad sites?
wc -l hyperdivergent_regions_with_elegans_rad_sites_sort_uniq.bed
#299

#299/312 = 0.9583333

#how many elegans pseudo-rad sites are in this region?

wc -l /Users/gavin/genome/pop_gen_revisions_9-2023/elegans_24_og_pseudo_rad_24_sorted_first_three_columns_number_of_sites_in_regions_sort_uniq.bed

#2105793
#that's a lot!

#how many sites do we have again in the pseudo-rad data??
wc -l elegans_24_og_pseudo_rad_24_sorted.bed
#14169530

#2105793/14169530 (0.1486142 , 15%)


