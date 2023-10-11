
#get data in there

dat <- read.table("site_positions.tsv",header=FALSE,sep="\t")

ele <- dat[dat$V3 == "C. elegans",]
ino <- dat[dat$V3 == "C. inopinata",]


ele_I <- ele[ele$V1 == "I",]
ele_II <- ele[ele$V1 == "II",]
ele_III <- ele[ele$V1 == "III",]
ele_IV <- ele[ele$V1 == "IV",]
ele_V <- ele[ele$V1 == "V",]
ele_X <- ele[ele$V1 == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$V2)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$V2)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$V2)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$V2)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$V2)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$V2)/8859471)/2

ino_I <- ino[ino$V1 == "Sp34_Chr1",]
ino_II <- ino[ino$V1 == "Sp34_Chr2",]
ino_III <- ino[ino$V1 == "Sp34_Chr3",]
ino_IV <- ino[ino$V1 == "Sp34_Chr4",]
ino_V <- ino[ino$V1 == "Sp34_Chr5",]
ino_X <- ino[ino$V1 == "Sp34_ChrX",]

ino_I$norm_dist_center <- abs((10297276-ino_I$V2)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$V2)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$V2)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$V2)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$V2)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$V2)/9095254)/2



all_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)


#define chromosome arms and centers
all_dat$chr_str_type <- as.factor(ifelse(all_dat$norm_dist_center >= 0.25,"Arms", "Centers"))

all_dat$category <- as.factor(paste(all_dat$V3,all_dat$V4, all_dat$chr_str_type))

number_of_sites_by_category <- as.data.frame(table(all_dat$category))

write.table(number_of_sites_by_category,"number_of_sites_by_category.tsv",row.names=FALSE, col.names=FALSE,quote=FALSE,sep='\t')

number_of_sites_by_category
#
#
#
# 
#                                         Var1    Freq
#1           C. elegans 0-fold degenerate Arms 1009858
#2        C. elegans 0-fold degenerate Centers 1786338
#3           C. elegans 2-fold degenerate Arms  454803
#4        C. elegans 2-fold degenerate Centers  803971
#5           C. elegans 4-fold degenerate Arms  247472
#6        C. elegans 4-fold degenerate Centers  434948
#7                   C. elegans all sites Arms 5730396
#8                C. elegans all sites Centers 8439007
#9                        C. elegans exon Arms 2242327
#10                    C. elegans exon Centers 3754036
#11       C. elegans first codon position Arms  584840
#12    C. elegans first codon position Centers 1015495
#13                      C. elegans genic Arms 3716973
#14                   C. elegans genic Centers 5835402
#15                      C. elegans inron Arms 1692637
#16                   C. elegans inron Centers 2260424
#17                 C. elegans intergenic Arms 2013423
#18              C. elegans intergenic Centers 2603605
#19            C. elegans invariant sites Arms 5654537
#20         C. elegans invariant sites Centers 8395704
#21      C. elegans second codon position Arms  584986
#22   C. elegans second codon position Centers 1015621
#23       C. elegans third codon position Arms  584285
#24    C. elegans third codon position Centers 1015290
#25              C. elegans variant sites Arms   75860
#26           C. elegans variant sites Centers   43304
#27        C. inopinata 0-fold degenerate Arms  361876
#28     C. inopinata 0-fold degenerate Centers  386747
#29        C. inopinata 2-fold degenerate Arms  167512
#30     C. inopinata 2-fold degenerate Centers  177744
#31        C. inopinata 4-fold degenerate Arms   87903
#32     C. inopinata 4-fold degenerate Centers   93097
#33                C. inopinata all sites Arms 2205424
#34             C. inopinata all sites Centers 2630024
#35                     C. inopinata exon Arms  656932
#36                  C. inopinata exon Centers  691067
#37     C. inopinata first codon position Arms  217746
#38  C. inopinata first codon position Centers  228871
#39                    C. inopinata genic Arms 1328288
#40                 C. inopinata genic Centers 1511527
#41               C. inopinata intergenic Arms  877136
#42            C. inopinata intergenic Centers 1118497
#43                   C. inopinata intron Arms  671356
#44                C. inopinata intron Centers  820460
#45          C. inopinata invariant sites Arms 2094662
#46       C. inopinata invariant sites Centers 2522407
#47    C. inopinata second codon position Arms  218252
#48 C. inopinata second codon position Centers  229358
#49     C. inopinata third codon position Arms  214954
#50  C. inopinata third codon position Centers  226684
#51            C. inopinata variant sites Arms  110765
#52         C. inopinata variant sites Centers  107623







#old, wrong data
                                         #Var1    Freq
#1           C. elegans 0-fold degenerate Arms 1009858
#2        C. elegans 0-fold degenerate Centers 1786339
#3           C. elegans 2-fold degenerate Arms  454803
#4        C. elegans 2-fold degenerate Centers  803972
#5           C. elegans 4-fold degenerate Arms  247472
#6        C. elegans 4-fold degenerate Centers  434948
#7                   C. elegans all sites Arms 5730446
#8                C. elegans all sites Centers 8439084
#9                        C. elegans exon Arms 3707270
#10                    C. elegans exon Centers 7034991
#11       C. elegans first codon position Arms  584840
#12    C. elegans first codon position Centers 1015496
#13                      C. elegans genic Arms 3829012
#14                   C. elegans genic Centers 6003937
#15                      C. elegans inron Arms 1692655
#16                   C. elegans inron Centers 2260467
#17                 C. elegans intergenic Arms 2013450
#18              C. elegans intergenic Centers 2603628
#19            C. elegans invariant sites Arms 5654586
#20         C. elegans invariant sites Centers 8395780
#21      C. elegans second codon position Arms  584986
#22   C. elegans second codon position Centers 1015622
#23       C. elegans third codon position Arms  584285
#24    C. elegans third codon position Centers 1015290
#25              C. elegans variant sites Arms   75860
#26           C. elegans variant sites Centers   43304
#27        C. inopinata 0-fold degenerate Arms  361879
#28     C. inopinata 0-fold degenerate Centers  386752
#29        C. inopinata 2-fold degenerate Arms  167513
#30     C. inopinata 2-fold degenerate Centers  177749
#31        C. inopinata 4-fold degenerate Arms   87905
#32     C. inopinata 4-fold degenerate Centers   93098
#33                C. inopinata all sites Arms 2205470
#34             C. inopinata all sites Centers 2630095
#35                     C. inopinata exon Arms  657852
#36                  C. inopinata exon Centers  691171
#37     C. inopinata first codon position Arms  217749
#38  C. inopinata first codon position Centers  228875
#39                    C. inopinata genic Arms 1339196
#40                 C. inopinata genic Centers 1525931
#41               C. inopinata intergenic Arms  877150
#42            C. inopinata intergenic Centers 1118534
#43                   C. inopinata intron Arms  671381
#44                C. inopinata intron Centers  820483
#45          C. inopinata invariant sites Arms 2094705
#46       C. inopinata invariant sites Centers 2522472
#47    C. inopinata second codon position Arms  218252
#48 C. inopinata second codon position Centers  229361
#49     C. inopinata third codon position Arms  214957
#50  C. inopinata third codon position Centers  226688
#51            C. inopinata variant sites Arms  110765
#52         C. inopinata variant sites Centers  107623
#
#
