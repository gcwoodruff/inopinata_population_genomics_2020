library(effsize)
library(adegenet)
library(ggplot2)
library(vcfR)
library(reshape2)
library(RColorBrewer)
library(lemon)
library(rstatix)


###nucleotide diversity and its genomic landscape

dat <- read.csv("ino_elg_pi.csv", header=TRUE)

dat$MB <- dat$start/1000000



#alright, pi summary stats...

summary(dat[dat$species == "C. elegans",]$pi_all)

#> summary(dat[dat$species == "C. elegans",]$pi_all)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.000500 0.000900 0.002169 0.001900 0.059600


summary(dat[dat$species == "C. inopinata",]$pi_all)
#> summary(dat[dat$species == "C. inopinata",]$pi_all)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.00590 0.00960 0.01129 0.01470 0.08500


var(dat[dat$species == "C. elegans",]$pi_all)

#[1] 1.956366e-05


var(dat[dat$species == "C. inopinata",]$pi_all)

#> var(dat[dat$species == "C. inopinata",]$pi_all)
#[1] 6.008132e-05

var(dat[dat$species == "C. inopinata",]$pi_all)/var(dat[dat$species == "C. elegans",]$pi_all)

#[1] 3.071067

#split by species, chromosome; get normalized chromosome positions (0 is end of chr, 0.5 is center of chr)
ele <- dat[dat$species == "C. elegans",]
ino <- dat[dat$species == "C. inopinata",]

ele_I <- ele[ele$scaffold == "I",]
ele_II <- ele[ele$scaffold == "II",]
ele_III <- ele[ele$scaffold == "III",]
ele_IV <- ele[ele$scaffold == "IV",]
ele_V <- ele[ele$scaffold == "V",]
ele_X <- ele[ele$scaffold == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$start)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$start)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$start)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$start)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$start)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$start)/8859471)/2

ino_I <- ino[ino$scaffold == "I",]
ino_II <- ino[ino$scaffold == "II",]
ino_III <- ino[ino$scaffold == "III",]
ino_IV <- ino[ino$scaffold == "IV",]
ino_V <- ino[ino$scaffold == "V",]
ino_X <- ino[ino$scaffold == "X",]

ino_I$norm_dist_center <- abs((10297276-ino_I$start)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$start)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$start)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$start)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$start)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$start)/9095254)/2


#put back together
pi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)
#classify position by being on chromosome arm or chromosome center
pi_dat$chr_str_type <- ifelse(pi_dat$norm_dist_center >= 0.25,"arms", "centers")


ele <- pi_dat[pi_dat$species == "C. elegans",]
ino <- pi_dat[pi_dat$species == "C. inopinata",]

#elegans: mean pi in chromosome arms, mean pi in chromosome centers, mean arms-mean centers
c(mean(ele[ele$chr_str_type == "arms",]$pi_all),mean(ele[ele$chr_str_type == "centers",]$pi_all), mean(ele[ele$chr_str_type == "arms",]$pi_all)-mean(ele[ele$chr_str_type == "centers",]$pi_all))

#[1] 0.0035030511 0.0009121709 0.0025908802

#inopinata: mean pi in chromosome arms, mean pi in chromosome centers, mean arms-mean centers


c(mean(ino[ino$chr_str_type == "arms",]$pi_all),mean(ino[ino$chr_str_type == "centers",]$pi_all), mean(ino[ino$chr_str_type == "arms",]$pi_all)-mean(ino[ino$chr_str_type == "centers",]$pi_all))

#[1] 0.012385479 0.010361194 0.002024285


#elegans, summary stats for chromosome arms and centers

summary(ele[ele$chr_str_type =="arms",]$pi_all)

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.000900 0.001700 0.003503 0.003300 0.059600

summary(ele[ele$chr_str_type =="centers",]$pi_all)

#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000000 0.0003000 0.0006000 0.0009122 0.0010000 0.0330000


#elegans, cohen's d effect size of chromosome position (arms or centers) on pi
cohen.d(ele[ele$chr_str_type =="arms",]$pi_all,ele[ele$chr_str_type =="centers",]$pi_all)


#Cohen's d
#
#d estimate: 0.6125738 (medium)
#95 percent confidence interval:
#    lower     upper
#0.5715993 0.6535483




#inopinata, summary stats for chromosome arms and centers



summary(ino[ino$chr_str_type =="arms",]$pi_all)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.00640 0.01060 0.01239 0.01630 0.07500

summary(ino[ino$chr_str_type =="centers",]$pi_all)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.00550 0.00890 0.01036 0.01340 0.08500

cohen.d(ino[ino$chr_str_type =="arms",]$pi_all,ino[ino$chr_str_type =="centers",]$pi_all)


#Cohen's d
#
#d estimate: 0.2633785 (small)
#95 percent confidence interval:
#    lower     upper
#0.2167175 0.3100394



#what about arm-center difference in pi on the X??

#elegans
ele_x <- ele[ele$scaffold == "X",]
ele_aut <- ele[ele$scaffold != "X",]

#elegans X
cohen.d(ele_x[ele_x$chr_str_type =="arms",]$pi_all,ele_x[ele_x$chr_str_type =="centers",]$pi_all)

#Cohen's d
#
#d estimate: 0.0839734 (negligible)
#95 percent confidence interval:
#      lower       upper
#-0.01063521  0.17858201

#elegans autosomes
cohen.d(ele_aut[ele_aut$chr_str_type =="arms",]$pi_all,ele_aut[ele_aut$chr_str_type =="centers",]$pi_all)


#Cohen's d
#
#d estimate: 0.6901529 (medium)
#95 percent confidence interval:
#    lower     upper
#0.6446425 0.7356633


#inopinata x v autosomes
ino_x <- ino[ino$scaffold == "X",]
ino_aut <- ino[ino$scaffold != "X",]


cohen.d(ino_x[ino_x$chr_str_type =="arms",]$pi_all,ino_x[ino_x$chr_str_type =="centers",]$pi_all)

#Cohen's d
#
#d estimate: 0.2457338 (small)
#95 percent confidence interval:
#    lower     upper
#0.1016421 0.3898256



cohen.d(ino_aut[ino_aut$chr_str_type =="arms",]$pi_all,ino_aut[ino_aut$chr_str_type =="centers",]$pi_all)

#Cohen's d
#
#d estimate: 0.2736956 (small)
#95 percent confidence interval:
#    lower     upper
#0.2243528 0.3230385




####
#genic intergenic pi
####

genexdat <- read.csv("genic_intergenic_pi.csv", header=TRUE,stringsAsFactors = TRUE)


genexdat$startMB <- genexdat$start/1000000



levels(genexdat$species)[levels(genexdat$species)=='elegans'] <- 'C. elegans'
levels(genexdat$species)[levels(genexdat$species)=='inopinata'] <- 'C. inopinata'
levels(genexdat$gen_region)[levels(genexdat$gen_region)=='genic'] <- 'Genic'
levels(genexdat$gen_region)[levels(genexdat$gen_region)=='intergenic'] <- 'Intergenic'


ele <- genexdat[genexdat$species == "C. elegans",]
ino <- genexdat[genexdat$species == "C. inopinata",]

ele_I <- ele[ele$scaffold == "I",]
ele_II <- ele[ele$scaffold == "II",]
ele_III <- ele[ele$scaffold == "III",]
ele_IV <- ele[ele$scaffold == "IV",]
ele_V <- ele[ele$scaffold == "V",]
ele_X <- ele[ele$scaffold == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$start)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$start)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$start)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$start)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$start)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$start)/8859471)/2

ino_I <- ino[ino$scaffold == "I",]
ino_II <- ino[ino$scaffold == "II",]
ino_III <- ino[ino$scaffold == "III",]
ino_IV <- ino[ino$scaffold == "IV",]
ino_V <- ino[ino$scaffold == "V",]
ino_X <- ino[ino$scaffold == "X",]

ino_I$norm_dist_center <- abs((10297276-ino_I$start)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$start)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$start)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$start)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$start)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$start)/9095254)/2



genexdat_b <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)




#define chromosome arms and centers
genexdat_b$chr_str_type <- ifelse(genexdat_b$norm_dist_center >= 0.25,"Arms", "Centers")

genexdat_b$species.gen_region <- paste(genexdat_b$species,genexdat_b$gen_region)

#cohen d effect size
as.data.frame(genexdat_b %>% cohens_d(pi_all ~ species.gen_region,ci = TRUE))

#     .y.                group1                  group2    effsize   n1   n2
#1 pi_all      C. elegans Genic   C. elegans Intergenic -0.0470822 8227 5808
#2 pi_all      C. elegans Genic      C. inopinata Genic -1.3706648 8227 4828
#3 pi_all      C. elegans Genic C. inopinata Intergenic -1.6378532 8227 3892
#4 pi_all C. elegans Intergenic      C. inopinata Genic -1.3488767 5808 4828
#5 pi_all C. elegans Intergenic C. inopinata Intergenic -1.6191998 5808 3892
#6 pi_all    C. inopinata Genic C. inopinata Intergenic -0.3703356 4828 3892
#  conf.low conf.high  magnitude
#1    -0.08     -0.01 negligible
#2    -1.42     -1.33      large
#3    -1.69     -1.59      large
#4    -1.40     -1.30      large
#5    -1.68     -1.57      large
#6    -0.41     -0.33      small


#write.table(as.data.frame(genexdat_b %>% cohens_d(pi_all ~ species.gen_region,ci = TRUE)),"cohens_d_pi_means_genic-intergenic.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#percent differences
mean_df <- aggregate(genexdat_b$pi_all,list(genexdat_b$species.gen_region), mean)
#get pairwise percent differences among those means
mean_perc_diffm <- ((outer(mean_df$x, mean_df$x, `-`))/mean_df$x)*-100
#round
mean_perc_diffm <- round(mean_perc_diffm)
#get column and row ids right
colnames(mean_perc_diffm) <- mean_df$Group.1
rownames(mean_perc_diffm) <- mean_df$Group.1
mean_perc_diffm

#                        C. elegans Genic C. elegans Intergenic
#C. elegans Genic                       0                     9
#C. elegans Intergenic                 -9                     0
#C. inopinata Genic                   -81                   -79
#C. inopinata Intergenic              -85                   -84
#                        C. inopinata Genic C. inopinata Intergenic
#C. elegans Genic                       429                     590
#C. elegans Intergenic                  384                     531
#C. inopinata Genic                       0                      30
#C. inopinata Intergenic                -23                       0

#write.table(mean_perc_diffm,"percent_differences_pi_means_genic-intergenic.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#wilcoxon rank sum test
as.data.frame(genexdat_b %>% wilcox_test(pi_all ~ species.gen_region,p.adjust.method="BH"))
#     .y.                group1                  group2   n1   n2 statistic
#1 pi_all      C. elegans Genic   C. elegans Intergenic 8227 5808  21472738
#2 pi_all      C. elegans Genic      C. inopinata Genic 8227 4828   3258034
#3 pi_all      C. elegans Genic C. inopinata Intergenic 8227 3892   1723822
#4 pi_all C. elegans Intergenic      C. inopinata Genic 5808 4828   2530960
#5 pi_all C. elegans Intergenic C. inopinata Intergenic 5808 3892   1315282
#6 pi_all    C. inopinata Genic C. inopinata Intergenic 4828 3892   7295671
#         p    p.adj p.adj.signif
#1 1.33e-24 1.33e-24         ****
#2 0.00e+00 0.00e+00         ****
#3 0.00e+00 0.00e+00         ****
#4 0.00e+00 0.00e+00         ****
#5 0.00e+00 0.00e+00         ****
#6 3.51e-72 4.21e-72         ****


#write.table(as.data.frame(genexdat_b %>% wilcox_test(pi_all ~ species.gen_region,p.adjust.method="BH")),"wilcox_test_genic-intergenic.tsv",sep="\t", quote=FALSE,row.names=TRUE)


summarystats <- aggregate(pi_all ~ species.gen_region,FUN=summary, data=genexdat_b)

variancestats <- aggregate(pi_all ~ species.gen_region,FUN=var, data=genexdat_b)

samplesizestats <- aggregate(pi_all ~ species.gen_region,FUN=length, data=genexdat_b)

summarystats <- cbind(summarystats,variance = variancestats$pi_all,sample_size = samplesizestats$pi_all)
summarystats

#       species.gen_region pi_all.Min. pi_all.1st Qu. pi_all.Median pi_all.Mean
#1        C. elegans Genic 0.000000000    0.000400000   0.000800000 0.001907172
#2   C. elegans Intergenic 0.000000000    0.000500000   0.001100000 0.002085537
#3      C. inopinata Genic 0.000000000    0.004900000   0.008500000 0.010096417
#4 C. inopinata Intergenic 0.000000000    0.006800000   0.011300000 0.013151413
#  pi_all.3rd Qu. pi_all.Max.     variance sample_size
#1    0.001800000 0.054600000 1.477751e-05        8227
#2    0.002200000 0.054700000 1.392627e-05        5808
#3    0.013400000 0.061000000 5.661542e-05        4828
#4    0.017400000 0.086100000 7.948531e-05        3892
#write.table(summarystats,"summary_stats_genic-intergenic.tsv",sep="\t", quote=FALSE,row.names=FALSE)


#genic/intergenic arms and centers

genexdat_b$species.gen_region.chr_str_type <- paste(genexdat_b$species,genexdat_b$gen_region,genexdat_b$chr_str_type)


#cohen d effect size
as.data.frame(genexdat_b %>% cohens_d(pi_all ~ species.gen_region.chr_str_type,ci = TRUE))
#      .y.                        group1                          group2
#1  pi_all         C. elegans Genic Arms        C. elegans Genic Centers
#2  pi_all         C. elegans Genic Arms      C. elegans Intergenic Arms
#3  pi_all         C. elegans Genic Arms   C. elegans Intergenic Centers
#4  pi_all         C. elegans Genic Arms         C. inopinata Genic Arms
#5  pi_all         C. elegans Genic Arms      C. inopinata Genic Centers
#6  pi_all         C. elegans Genic Arms    C. inopinata Intergenic Arms
#7  pi_all         C. elegans Genic Arms C. inopinata Intergenic Centers
#8  pi_all      C. elegans Genic Centers      C. elegans Intergenic Arms
#9  pi_all      C. elegans Genic Centers   C. elegans Intergenic Centers
#10 pi_all      C. elegans Genic Centers         C. inopinata Genic Arms
#11 pi_all      C. elegans Genic Centers      C. inopinata Genic Centers
#12 pi_all      C. elegans Genic Centers    C. inopinata Intergenic Arms
#13 pi_all      C. elegans Genic Centers C. inopinata Intergenic Centers
#14 pi_all    C. elegans Intergenic Arms   C. elegans Intergenic Centers
#15 pi_all    C. elegans Intergenic Arms         C. inopinata Genic Arms
#16 pi_all    C. elegans Intergenic Arms      C. inopinata Genic Centers
#17 pi_all    C. elegans Intergenic Arms    C. inopinata Intergenic Arms
#18 pi_all    C. elegans Intergenic Arms C. inopinata Intergenic Centers
#19 pi_all C. elegans Intergenic Centers         C. inopinata Genic Arms
#20 pi_all C. elegans Intergenic Centers      C. inopinata Genic Centers
#21 pi_all C. elegans Intergenic Centers    C. inopinata Intergenic Arms
#22 pi_all C. elegans Intergenic Centers C. inopinata Intergenic Centers
#23 pi_all       C. inopinata Genic Arms      C. inopinata Genic Centers
#24 pi_all       C. inopinata Genic Arms    C. inopinata Intergenic Arms
#25 pi_all       C. inopinata Genic Arms C. inopinata Intergenic Centers
#26 pi_all    C. inopinata Genic Centers    C. inopinata Intergenic Arms
#27 pi_all    C. inopinata Genic Centers C. inopinata Intergenic Centers
#28 pi_all  C. inopinata Intergenic Arms C. inopinata Intergenic Centers
#       effsize   n1   n2 conf.low conf.high  magnitude
#1   0.61865169 3812 4415     0.59      0.65   moderate
#2  -0.01399409 3812 2678    -0.06      0.03 negligible
#3   0.53250081 3812 3130     0.50      0.57   moderate
#4  -1.18095662 3812 2246    -1.24     -1.12      large
#5  -0.99328177 3812 2582    -1.06     -0.93      large
#6  -1.48511421 3812 1725    -1.56     -1.41      large
#7  -1.29776952 3812 2167    -1.37     -1.23      large
#8  -0.66991423 4415 2678    -0.71     -0.63   moderate
#9  -0.16742584 4415 3130    -0.21     -0.12 negligible
#10 -1.75763722 4415 2246    -1.83     -1.69      large
#11 -1.72450668 4415 2582    -1.80     -1.66      large
#12 -2.00284499 4415 1725    -2.09     -1.92      large
#13 -1.91488897 4415 2167    -2.00     -1.82      large
#14  0.57759670 2678 3130     0.54      0.62   moderate
#15 -1.18785268 2678 2246    -1.26     -1.12      large
#16 -1.00125865 2678 2582    -1.07     -0.93      large
#17 -1.49318832 2678 1725    -1.57     -1.42      large
#18 -1.30696120 2678 2167    -1.38     -1.23      large
#19 -1.69545894 3130 2246    -1.77     -1.63      large
#20 -1.64201175 3130 2582    -1.72     -1.57      large
#21 -1.94895082 3130 1725    -2.03     -1.87      large
#22 -1.84903735 3130 2167    -1.94     -1.77      large
#23  0.31027111 2246 2582     0.25      0.37      small
#24 -0.36638106 2246 1725    -0.43     -0.31      small
#25 -0.07592101 2246 2167    -0.13     -0.02 negligible
#26 -0.68252717 2582 1725    -0.74     -0.62   moderate
#27 -0.40049678 2582 2167    -0.46     -0.35      small
#28  0.30082916 1725 2167     0.24      0.37      small

#write.table(as.data.frame(genexdat_b %>% cohens_d(pi_all ~ species.gen_region.chr_str_type,ci = TRUE)),"cohens_d_pi_means_genic-intergenic_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#percent differences
mean_df <- aggregate(genexdat_b$pi_all,list(genexdat_b$species.gen_region.chr_str_type), mean)
#get pairwise percent differences among those means
mean_perc_diffm <- ((outer(mean_df$x, mean_df$x, `-`))/mean_df$x)*-100
#round
mean_perc_diffm <- round(mean_perc_diffm)
#get column and row ids right
colnames(mean_perc_diffm) <- mean_df$Group.1
rownames(mean_perc_diffm) <- mean_df$Group.1
mean_perc_diffm

#                                C. elegans Genic Arms C. elegans Genic Centers
#C. elegans Genic Arms                               0                      -74
#C. elegans Genic Centers                          284                        0
#C. elegans Intergenic Arms                         -2                      -75
#C. elegans Intergenic Centers                     187                      -25
#C. inopinata Genic Arms                           -72                      -93
#C. inopinata Genic Centers                        -65                      -91
#C. inopinata Intergenic Arms                      -78                      -94
#C. inopinata Intergenic Centers                   -74                      -93
#                                C. elegans Intergenic Arms
#C. elegans Genic Arms                                    2
#C. elegans Genic Centers                               293
#C. elegans Intergenic Arms                               0
#C. elegans Intergenic Centers                          193
#C. inopinata Genic Arms                                -71
#C. inopinata Genic Centers                             -64
#C. inopinata Intergenic Arms                           -78
#C. inopinata Intergenic Centers                        -73
#                                C. elegans Intergenic Centers
#C. elegans Genic Arms                                     -65
#C. elegans Genic Centers                                   34
#C. elegans Intergenic Arms                                -66
#C. elegans Intergenic Centers                               0
#C. inopinata Genic Arms                                   -90
#C. inopinata Genic Centers                                -88
#C. inopinata Intergenic Arms                              -92
#C. inopinata Intergenic Centers                           -91
#                                C. inopinata Genic Arms
#C. elegans Genic Arms                               259
#C. elegans Genic Centers                           1278
#C. elegans Intergenic Arms                          251
#C. elegans Intergenic Centers                       927
#C. inopinata Genic Arms                               0
#C. inopinata Genic Centers                           26
#C. inopinata Intergenic Arms                        -23
#C. inopinata Intergenic Centers                      -5
#                                C. inopinata Genic Centers
#C. elegans Genic Arms                                  185
#C. elegans Genic Centers                               995
#C. elegans Intergenic Arms                             179
#C. elegans Intergenic Centers                          717
#C. inopinata Genic Arms                                -21
#C. inopinata Genic Centers                               0
#C. inopinata Intergenic Arms                           -38
#C. inopinata Intergenic Centers                        -25
#                                C. inopinata Intergenic Arms
#C. elegans Genic Arms                                    363
#C. elegans Genic Centers                                1679
#C. elegans Intergenic Arms                               353
#C. elegans Intergenic Centers                           1227
#C. inopinata Genic Arms                                   29
#C. inopinata Genic Centers                                62
#C. inopinata Intergenic Arms                               0
#C. inopinata Intergenic Centers                           22
#                                C. inopinata Intergenic Centers
#C. elegans Genic Arms                                       278
#C. elegans Genic Centers                                   1354
#C. elegans Intergenic Arms                                  270
#C. elegans Intergenic Centers                               984
#C. inopinata Genic Arms                                       5
#C. inopinata Genic Centers                                   33
#C. inopinata Intergenic Arms                                -18
#C. inopinata Intergenic Centers                               0

#write.table(mean_perc_diffm,"percent_differences_pi_means_genic-intergenic_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#wilcoxon rank sum test
as.data.frame(genexdat_b %>% wilcox_test(pi_all ~ species.gen_region.chr_str_type,p.adjust.method="BH"))

#      .y.                        group1                          group2   n1
#1  pi_all         C. elegans Genic Arms        C. elegans Genic Centers 3812
#2  pi_all         C. elegans Genic Arms      C. elegans Intergenic Arms 3812
#3  pi_all         C. elegans Genic Arms   C. elegans Intergenic Centers 3812
#4  pi_all         C. elegans Genic Arms         C. inopinata Genic Arms 3812
#5  pi_all         C. elegans Genic Arms      C. inopinata Genic Centers 3812
#6  pi_all         C. elegans Genic Arms    C. inopinata Intergenic Arms 3812
#7  pi_all         C. elegans Genic Arms C. inopinata Intergenic Centers 3812
#8  pi_all      C. elegans Genic Centers      C. elegans Intergenic Arms 4415
#9  pi_all      C. elegans Genic Centers   C. elegans Intergenic Centers 4415
#10 pi_all      C. elegans Genic Centers         C. inopinata Genic Arms 4415
#11 pi_all      C. elegans Genic Centers      C. inopinata Genic Centers 4415
#12 pi_all      C. elegans Genic Centers    C. inopinata Intergenic Arms 4415
#13 pi_all      C. elegans Genic Centers C. inopinata Intergenic Centers 4415
#14 pi_all    C. elegans Intergenic Arms   C. elegans Intergenic Centers 2678
#15 pi_all    C. elegans Intergenic Arms         C. inopinata Genic Arms 2678
#16 pi_all    C. elegans Intergenic Arms      C. inopinata Genic Centers 2678
#17 pi_all    C. elegans Intergenic Arms    C. inopinata Intergenic Arms 2678
#18 pi_all    C. elegans Intergenic Arms C. inopinata Intergenic Centers 2678
#19 pi_all C. elegans Intergenic Centers         C. inopinata Genic Arms 3130
#20 pi_all C. elegans Intergenic Centers      C. inopinata Genic Centers 3130
#21 pi_all C. elegans Intergenic Centers    C. inopinata Intergenic Arms 3130
#22 pi_all C. elegans Intergenic Centers C. inopinata Intergenic Centers 3130
#23 pi_all       C. inopinata Genic Arms      C. inopinata Genic Centers 2246
#24 pi_all       C. inopinata Genic Arms    C. inopinata Intergenic Arms 2246
#25 pi_all       C. inopinata Genic Arms C. inopinata Intergenic Centers 2246
#26 pi_all    C. inopinata Genic Centers    C. inopinata Intergenic Arms 2582
#27 pi_all    C. inopinata Genic Centers C. inopinata Intergenic Centers 2582
#28 pi_all  C. inopinata Intergenic Arms C. inopinata Intergenic Centers 1725
#     n2  statistic         p     p.adj p.adj.signif
#1  4415 13006648.5  0.00e+00  0.00e+00         ****
#2  2678  4796222.0  3.37e-05  3.49e-05         ****
#3  3130  8471356.0 4.80e-200 6.72e-200         ****
#4  2246  1070095.5  0.00e+00  0.00e+00         ****
#5  2582  1501715.5  0.00e+00  0.00e+00         ****
#6  1725   550066.0  0.00e+00  0.00e+00         ****
#7  2167   878166.5  0.00e+00  0.00e+00         ****
#8  2678  2374924.5  0.00e+00  0.00e+00         ****
#9  3130  5830235.5  4.00e-31  4.67e-31         ****
#10 2246   289848.0  0.00e+00  0.00e+00         ****
#11 2582   396375.0  0.00e+00  0.00e+00         ****
#12 1725   107664.5  0.00e+00  0.00e+00         ****
#13 2167   187924.5  0.00e+00  0.00e+00         ****
#14 3130  6206221.0 5.78e-220 8.52e-220         ****
#15 2246   773950.0  0.00e+00  0.00e+00         ****
#16 2582  1099089.5  0.00e+00  0.00e+00         ****
#17 1725   389785.5  0.00e+00  0.00e+00         ****
#18 2167   630388.0  0.00e+00  0.00e+00         ****
#19 2246   272735.0  0.00e+00  0.00e+00         ****
#20 2582   385186.0  0.00e+00  0.00e+00         ****
#21 1725   107964.5  0.00e+00  0.00e+00         ****
#22 2167   187143.5  0.00e+00  0.00e+00         ****
#23 2582  3377786.0  4.20e-23  4.70e-23         ****
#24 1725  1507650.0  3.80e-33  4.63e-33         ****
#25 2167  2277280.5  2.22e-04  2.22e-04          ***
#26 1725  1371758.5 1.74e-101 2.32e-101         ****
#27 2167  2138982.0  1.67e-44  2.13e-44         ****
#28 2167  2179836.5  4.46e-19  4.80e-19         ****

#write.table(as.data.frame(genexdat_b %>% wilcox_test(pi_all ~ species.gen_region.chr_str_type,p.adjust.method="BH")),"wilcox_test_genic-intergenic_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#summary stats
summarystats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=summary, data=genexdat_b)

variancestats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=var, data=genexdat_b)

samplesizestats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=length, data=genexdat_b)

summarystats <- cbind(summarystats,variance = variancestats$pi_all,sample_size = samplesizestats$pi_all)
summarystats

#  species.gen_region.chr_str_type  pi_all.Min. pi_all.1st Qu. pi_all.Median
#1           C. elegans Genic Arms 0.0000000000   0.0007000000  0.0015000000
#2        C. elegans Genic Centers 0.0000000000   0.0002000000  0.0005000000
#3      C. elegans Intergenic Arms 0.0000000000   0.0009000000  0.0018000000
#4   C. elegans Intergenic Centers 0.0000000000   0.0003000000  0.0007000000
#5         C. inopinata Genic Arms 0.0000000000   0.0054250000  0.0096000000
#6      C. inopinata Genic Centers 0.0000000000   0.0044000000  0.0076000000
#7    C. inopinata Intergenic Arms 0.0000000000   0.0076000000  0.0127000000
#8 C. inopinata Intergenic Centers 0.0000000000   0.0063000000  0.0104000000
#   pi_all.Mean pi_all.3rd Qu.  pi_all.Max.     variance sample_size
#1 0.0031629591   0.0031250000 0.0546000000 2.644806e-05        3812
#2 0.0008228992   0.0009000000 0.0245000000 2.166795e-06        4415
#3 0.0032330471   0.0034000000 0.0547000000 2.372004e-05        2678
#4 0.0011037380   0.0013000000 0.0322000000 3.460494e-06        3130
#5 0.0113403829   0.0152000000 0.0602000000 6.944676e-05        2246
#6 0.0090143300   0.0120000000 0.0610000000 4.295847e-05        2582
#7 0.0146429565   0.0197000000 0.0790000000 9.305909e-05        1725
#8 0.0119640978   0.0159000000 0.0861000000 6.553600e-05        2167

#write.table(summarystats,"summary_stats_genic-intergenic_arms-centers.tsv",sep="\t", quote=FALSE,row.names=FALSE)



####
#exon intron pi
####


exindat <- read.csv("exon_intron_pi.csv", header=TRUE,stringsAsFactors = TRUE)


exindat$startMB <- exindat$start/1000000

levels(exindat$species)[levels(exindat$species)=='elegans'] <- 'C. elegans'
levels(exindat$species)[levels(exindat$species)=='inopinata'] <- 'C. inopinata'
levels(exindat$gen_region)[levels(exindat$gen_region)=='exonic'] <- 'Exon'
levels(exindat$gen_region)[levels(exindat$gen_region)=='intronic'] <- 'Intron'


ele <- exindat[exindat$species == "C. elegans",]
ino <- exindat[exindat$species == "C. inopinata",]

ele_I <- ele[ele$scaffold == "I",]
ele_II <- ele[ele$scaffold == "II",]
ele_III <- ele[ele$scaffold == "III",]
ele_IV <- ele[ele$scaffold == "IV",]
ele_V <- ele[ele$scaffold == "V",]
ele_X <- ele[ele$scaffold == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$start)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$start)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$start)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$start)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$start)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$start)/8859471)/2

ino_I <- ino[ino$scaffold == "I",]
ino_II <- ino[ino$scaffold == "II",]
ino_III <- ino[ino$scaffold == "III",]
ino_IV <- ino[ino$scaffold == "IV",]
ino_V <- ino[ino$scaffold == "V",]
ino_X <- ino[ino$scaffold == "X",]

ino_I$norm_dist_center <- abs((10297276-ino_I$start)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$start)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$start)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$start)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$start)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$start)/9095254)/2


exinpi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)


#define chromosome arms and centers
exinpi_dat$chr_str_type <- ifelse(exinpi_dat$norm_dist_center >= 0.25,"Arms", "Centers")


exinpi_dat$gen_region <- factor(exinpi_dat$gen_region , levels=c('Intron','Exon'))


exinpi_dat$species.gen_region <- paste(exinpi_dat$species,exinpi_dat$gen_region)


#cohen d eff size

as.data.frame(exinpi_dat %>% cohens_d(pi_all ~ species.gen_region,ci = TRUE))

#     .y.            group1              group2     effsize   n1   n2 conf.low conf.high  magnitude
#1 pi_all   C. elegans Exon   C. elegans Intron -0.03956498 7720 6342    -0.07   -0.0087 negligible
#2 pi_all   C. elegans Exon   C. inopinata Exon -0.87714255 7720 2873    -0.93   -0.8300      large
#3 pi_all   C. elegans Exon C. inopinata Intron -1.72241826 7720 3016    -1.79   -1.6600      large
#4 pi_all C. elegans Intron   C. inopinata Exon -0.93046330 6342 2873    -0.98   -0.8800      large
#5 pi_all C. elegans Intron C. inopinata Intron -1.78830066 6342 3016    -1.85   -1.7300      large
#6 pi_all C. inopinata Exon C. inopinata Intron -0.97211137 2873 3016    -1.03   -0.9100      large


#write.table(as.data.frame(exinpi_dat %>% cohens_d(pi_all ~ species.gen_region,ci = TRUE)),"cohens_d_pi_means_exon-intron.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#percent differences
mean_df <- aggregate(exinpi_dat$pi_all,list(exinpi_dat$species.gen_region), mean)
#get pairwise percent differences among those means
mean_perc_diffm <- ((outer(mean_df$x, mean_df$x, `-`))/mean_df$x)*-100
#round
mean_perc_diffm <- round(mean_perc_diffm)
#get column and row ids right
colnames(mean_perc_diffm) <- mean_df$Group.1
rownames(mean_perc_diffm) <- mean_df$Group.1
mean_perc_diffm

#                    C. elegans Exon C. elegans Intron C. inopinata Exon C. inopinata Intron
#C. elegans Exon                   0                 8               263                 685
#C. elegans Intron                -8                 0               235                 625
#C. inopinata Exon               -72               -70                 0                 116
#C. inopinata Intron             -87               -86               -54                   0


#write.table(mean_perc_diffm,"percent_differences_pi_means_exon-intron.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#wilcoxon rank-sum test
as.data.frame(exinpi_dat %>% wilcox_test(pi_all ~ species.gen_region,p.adjust.method="BH"))

     #.y.            group1              group2   n1   n2  statistic         p       p.adj p.adj.signif
#1 pi_all   C. elegans Exon   C. elegans Intron 7720 6342 19218040.0 1.91e-107 1 1.91e-107         ****
#2 pi_all   C. elegans Exon   C. inopinata Exon 7720 2873  3109539.5  0.00e+00 2  0.00e+00         ****
#3 pi_all   C. elegans Exon C. inopinata Intron 7720 3016  1140656.5  0.00e+00 3  0.00e+00         ****
#4 pi_all C. elegans Intron   C. inopinata Exon 6342 2873  3183794.5  0.00e+00 4  0.00e+00         ****
#5 pi_all C. elegans Intron C. inopinata Intron 6342 3016   925266.5  0.00e+00 5  0.00e+00         ****
#6 pi_all C. inopinata Exon C. inopinata Intron 2873 3016  1829834.0  0.00e+00 6  0.00e+00         ****


#write.table(as.data.frame(exinpi_dat %>% wilcox_test(pi_all ~ species.gen_region,p.adjust.method="BH")),"wilcox_test_exon-intron.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#summary stats
summarystats <- aggregate(pi_all ~ species.gen_region,FUN=summary, data=exinpi_dat)

variancestats <- aggregate(pi_all ~ species.gen_region,FUN=var, data=exinpi_dat)

samplesizestats <- aggregate(pi_all ~ species.gen_region,FUN=length, data=exinpi_dat)

summarystats <- cbind(summarystats,variance = variancestats$pi_all,sample_size = samplesizestats$pi_all)
summarystats
#   species.gen_region pi_all.Min. pi_all.1st Qu. pi_all.Median pi_all.Mean
#1     C. elegans Exon 0.000000000    0.000200000   0.000500000 0.001697811
#2   C. elegans Intron 0.000000000    0.000400000   0.001000000 0.001839577
#3   C. inopinata Exon 0.000000000    0.002400000   0.004700000 0.006163105
#4 C. inopinata Intron 0.000000000    0.007200000   0.011700000 0.013330968
#  pi_all.3rd Qu. pi_all.Max.     variance sample_size
#1    0.001400000 0.054600000 1.716308e-05        7720
#2    0.002200000 0.036900000 8.514579e-06        6342
#3    0.008100000 0.069400000 3.466793e-05        2873
#4    0.017900000 0.081900000 7.406908e-05        3016


#write.table(summarystats,"summary_stats_exon-intron.tsv",sep="\t", quote=FALSE,row.names=FALSE)

#exon intron; chromsome arms and centers

exinpi_dat$species.gen_region.chr_str_type <- paste(exinpi_dat$species,exinpi_dat$gen_region,exinpi_dat$chr_str_type)


#cohen d effect size
as.data.frame(exinpi_dat %>% cohens_d(pi_all ~ species.gen_region.chr_str_type,ci = TRUE))
#      .y.                    group1                      group2      effsize
#1  pi_all      C. elegans Exon Arms     C. elegans Exon Centers  0.537607629
#2  pi_all      C. elegans Exon Arms      C. elegans Intron Arms -0.004511046
#3  pi_all      C. elegans Exon Arms   C. elegans Intron Centers  0.460118040
#4  pi_all      C. elegans Exon Arms      C. inopinata Exon Arms -0.673997227
#5  pi_all      C. elegans Exon Arms   C. inopinata Exon Centers -0.457616828
#6  pi_all      C. elegans Exon Arms    C. inopinata Intron Arms -1.572508673
#7  pi_all      C. elegans Exon Arms C. inopinata Intron Centers -1.341200617
#8  pi_all   C. elegans Exon Centers      C. elegans Intron Arms -0.805969195
#9  pi_all   C. elegans Exon Centers   C. elegans Intron Centers -0.188444531
#10 pi_all   C. elegans Exon Centers      C. inopinata Exon Arms -1.327822453
#11 pi_all   C. elegans Exon Centers   C. inopinata Exon Centers -1.275701129
#12 pi_all   C. elegans Exon Centers    C. inopinata Intron Arms -2.136211146
#13 pi_all   C. elegans Exon Centers C. inopinata Intron Centers -2.076337244
#14 pi_all    C. elegans Intron Arms   C. elegans Intron Centers  0.681856435
#15 pi_all    C. elegans Intron Arms      C. inopinata Exon Arms -0.771932396
#16 pi_all    C. elegans Intron Arms   C. inopinata Exon Centers -0.554240561
#17 pi_all    C. elegans Intron Arms    C. inopinata Intron Arms -1.705360738
#18 pi_all    C. elegans Intron Arms C. inopinata Intron Centers -1.511049467
#19 pi_all C. elegans Intron Centers      C. inopinata Exon Arms -1.255167656
#20 pi_all C. elegans Intron Centers   C. inopinata Exon Centers -1.176901371
#21 pi_all C. elegans Intron Centers    C. inopinata Intron Arms -2.083562998
#22 pi_all C. elegans Intron Centers C. inopinata Intron Centers -2.006508230
#23 pi_all    C. inopinata Exon Arms   C. inopinata Exon Centers  0.294420941
#24 pi_all    C. inopinata Exon Arms    C. inopinata Intron Arms -0.993418879
#25 pi_all    C. inopinata Exon Arms C. inopinata Intron Centers -0.664440192
#26 pi_all C. inopinata Exon Centers    C. inopinata Intron Arms -1.302119615
#27 pi_all C. inopinata Exon Centers C. inopinata Intron Centers -1.016448098
#28 pi_all  C. inopinata Intron Arms C. inopinata Intron Centers  0.408995216
#     n1   n2 conf.low conf.high  magnitude
#1  3540 4180     0.51      0.57   moderate
#2  3540 2809    -0.06      0.04 negligible
#3  3540 3533     0.43      0.49      small
#4  3540 1362    -0.74     -0.61   moderate
#5  3540 1511    -0.52     -0.39      small
#6  3540 1371    -1.66     -1.49      large
#7  3540 1645    -1.42     -1.26      large
#8  4180 2809    -0.86     -0.76      large
#9  4180 3533    -0.23     -0.14 negligible
#10 4180 1362    -1.44     -1.22      large
#11 4180 1511    -1.35     -1.21      large
#12 4180 1371    -2.24     -2.03      large
#13 4180 1645    -2.17     -1.98      large
#14 2809 3533     0.63      0.73   moderate
#15 2809 1362    -0.84     -0.71   moderate
#16 2809 1511    -0.62     -0.50   moderate
#17 2809 1371    -1.80     -1.63      large
#18 2809 1645    -1.59     -1.44      large
#19 3533 1362    -1.37     -1.16      large
#20 3533 1511    -1.25     -1.11      large
#21 3533 1371    -2.18     -1.99      large
#22 3533 1645    -2.10     -1.92      large
#23 1362 1511     0.22      0.37      small
#24 1362 1371    -1.08     -0.91      large
#25 1362 1645    -0.75     -0.58   moderate
#26 1511 1371    -1.39     -1.23      large
#27 1511 1645    -1.10     -0.93      large
#28 1371 1645     0.34      0.47      small

#write.table(as.data.frame(exinpi_dat %>% cohens_d(pi_all ~ species.gen_region.chr_str_type,ci = TRUE)),"cohens_d_pi_means_exon-intron_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#percent differences
mean_df <- aggregate(exinpi_dat$pi_all,list(exinpi_dat$species.gen_region.chr_str_type), mean)
#get pairwise percent differences among those means
mean_perc_diffm <- ((outer(mean_df$x, mean_df$x, `-`))/mean_df$x)*-100
#round
mean_perc_diffm <- round(mean_perc_diffm)
#get column and row ids right
colnames(mean_perc_diffm) <- mean_df$Group.1
rownames(mean_perc_diffm) <- mean_df$Group.1
mean_perc_diffm
#                            C. elegans Exon Arms C. elegans Exon Centers
#C. elegans Exon Arms                           0                     -77
#C. elegans Exon Centers                      330                       0
#C. elegans Intron Arms                        -1                     -77
#C. elegans Intron Centers                    198                     -31
#C. inopinata Exon Arms                       -59                     -90
#C. inopinata Exon Centers                    -46                     -87
#C. inopinata Intron Arms                     -81                     -96
#C. inopinata Intron Centers                  -75                     -94
#                            C. elegans Intron Arms C. elegans Intron Centers
#C. elegans Exon Arms                             1                       -66
#C. elegans Exon Centers                        333                        44
#C. elegans Intron Arms                           0                       -67
#C. elegans Intron Centers                      200                         0
#C. inopinata Exon Arms                         -59                       -86
#C. inopinata Exon Centers                      -45                       -82
#C. inopinata Intron Arms                       -81                       -94
#C. inopinata Intron Centers                    -75                       -92
#                            C. inopinata Exon Arms C. inopinata Exon Centers
#C. elegans Exon Arms                           143                        84
#C. elegans Exon Centers                        946                       690
#C. elegans Intron Arms                         142                        83
#C. elegans Intron Centers                      625                       448
#C. inopinata Exon Arms                           0                       -24
#C. inopinata Exon Centers                       32                         0
#C. inopinata Intron Arms                       -54                       -65
#C. inopinata Intron Centers                    -40                       -55
#                            C. inopinata Intron Arms
#C. elegans Exon Arms                             425
#C. elegans Exon Centers                         2153
#C. elegans Intron Arms                           421
#C. elegans Intron Centers                       1461
#C. inopinata Exon Arms                           115
#C. inopinata Exon Centers                        185
#C. inopinata Intron Arms                           0
#C. inopinata Intron Centers                       30
#                            C. inopinata Intron Centers
#C. elegans Exon Arms                                304
#C. elegans Exon Centers                            1637
#C. elegans Intron Arms                              301
#C. elegans Intron Centers                          1104
#C. inopinata Exon Arms                               66
#C. inopinata Exon Centers                           120
#C. inopinata Intron Arms                            -23
#C. inopinata Intron Centers                           0

#write.table(mean_perc_diffm,"percent_differences_pi_means_exon-intron_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#wilcoxon rank sum test
as.data.frame(exinpi_dat %>% wilcox_test(pi_all ~ species.gen_region.chr_str_type,p.adjust.method="BH"))
#      .y.                    group1                      group2   n1   n2
#1  pi_all      C. elegans Exon Arms     C. elegans Exon Centers 3540 4180
#2  pi_all      C. elegans Exon Arms      C. elegans Intron Arms 3540 2809
#3  pi_all      C. elegans Exon Arms   C. elegans Intron Centers 3540 3533
#4  pi_all      C. elegans Exon Arms      C. inopinata Exon Arms 3540 1362
#5  pi_all      C. elegans Exon Arms   C. inopinata Exon Centers 3540 1511
#6  pi_all      C. elegans Exon Arms    C. inopinata Intron Arms 3540 1371
#7  pi_all      C. elegans Exon Arms C. inopinata Intron Centers 3540 1645
#8  pi_all   C. elegans Exon Centers      C. elegans Intron Arms 4180 2809
#9  pi_all   C. elegans Exon Centers   C. elegans Intron Centers 4180 3533
#10 pi_all   C. elegans Exon Centers      C. inopinata Exon Arms 4180 1362
#11 pi_all   C. elegans Exon Centers   C. inopinata Exon Centers 4180 1511
#12 pi_all   C. elegans Exon Centers    C. inopinata Intron Arms 4180 1371
#13 pi_all   C. elegans Exon Centers C. inopinata Intron Centers 4180 1645
#14 pi_all    C. elegans Intron Arms   C. elegans Intron Centers 2809 3533
#15 pi_all    C. elegans Intron Arms      C. inopinata Exon Arms 2809 1362
#16 pi_all    C. elegans Intron Arms   C. inopinata Exon Centers 2809 1511
#17 pi_all    C. elegans Intron Arms    C. inopinata Intron Arms 2809 1371
#18 pi_all    C. elegans Intron Arms C. inopinata Intron Centers 2809 1645
#19 pi_all C. elegans Intron Centers      C. inopinata Exon Arms 3533 1362
#20 pi_all C. elegans Intron Centers   C. inopinata Exon Centers 3533 1511
#21 pi_all C. elegans Intron Centers    C. inopinata Intron Arms 3533 1371
#22 pi_all C. elegans Intron Centers C. inopinata Intron Centers 3533 1645
#23 pi_all    C. inopinata Exon Arms   C. inopinata Exon Centers 1362 1511
#24 pi_all    C. inopinata Exon Arms    C. inopinata Intron Arms 1362 1371
#25 pi_all    C. inopinata Exon Arms C. inopinata Intron Centers 1362 1645
#26 pi_all C. inopinata Exon Centers    C. inopinata Intron Arms 1511 1371
#27 pi_all C. inopinata Exon Centers C. inopinata Intron Centers 1511 1645
#28 pi_all  C. inopinata Intron Arms C. inopinata Intron Centers 1371 1645
#    statistic         p     p.adj p.adj.signif
#1  10353994.0 5.16e-203 8.50e-203         ****
#2   3828047.0  4.12e-56  4.44e-56         ****
#3   7703373.5  2.92e-64  3.27e-64         ****
#4    922649.0 7.70e-247 1.44e-246         ****
#5   1311271.0 8.80e-182 1.37e-181         ****
#6    375584.5  0.00e+00  0.00e+00         ****
#7    567189.5  0.00e+00  0.00e+00         ****
#8   2019521.5  0.00e+00  0.00e+00         ****
#9   5667098.0  4.08e-70  4.76e-70         ****
#10   329328.0  0.00e+00  0.00e+00         ****
#11   546291.5  0.00e+00  0.00e+00         ****
#12    84441.0  0.00e+00  0.00e+00         ****
#13   113441.5  0.00e+00  0.00e+00         ****
#14  7546361.0 1.67e-279 3.34e-279         ****
#15   878333.0 4.48e-177 6.27e-177         ****
#16  1298199.5  1.18e-98  1.50e-98         ****
#17   267234.0  0.00e+00  0.00e+00         ****
#18   438884.0  0.00e+00  0.00e+00         ****
#19   375838.0  0.00e+00  0.00e+00         ****
#20   631424.0  0.00e+00  0.00e+00         ****
#21    90029.5  0.00e+00  0.00e+00         ****
#22   129119.0  0.00e+00  0.00e+00         ****
#23  1215119.5  5.11e-17  5.11e-17         ****
#24   388088.5 3.60e-154 4.80e-154         ****
#25   623076.5  1.01e-97  1.23e-97         ****
#26   309715.5 2.33e-232 4.08e-232         ****
#27   508953.5 4.04e-181 5.95e-181         ****
#28  1378201.0  6.83e-26  7.08e-26         ****

#write.table(as.data.frame(exinpi_dat %>% wilcox_test(pi_all ~ species.gen_region.chr_str_type,p.adjust.method="BH")),"wilcox_test_genic-intergenic_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#summary stats
summarystats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=summary, data=exinpi_dat)

variancestats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=var, data=exinpi_dat)

samplesizestats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=length, data=exinpi_dat)

summarystats <- cbind(summarystats,variance = variancestats$pi_all,sample_size = samplesizestats$pi_all)
summarystats
#  species.gen_region.chr_str_type  pi_all.Min. pi_all.1st Qu. pi_all.Median
#1            C. elegans Exon Arms 0.0000000000   0.0003000000  0.0010000000
#2         C. elegans Exon Centers 0.0000000000   0.0001000000  0.0003000000
#3          C. elegans Intron Arms 0.0000000000   0.0009000000  0.0019000000
#4       C. elegans Intron Centers 0.0000000000   0.0002000000  0.0006000000
#5          C. inopinata Exon Arms 0.0000000000   0.0030000000  0.0053000000
#6       C. inopinata Exon Centers 0.0000000000   0.0020000000  0.0040000000
#7        C. inopinata Intron Arms 0.0000000000   0.0084000000  0.0135000000
#8     C. inopinata Intron Centers 0.0000000000   0.0065000000  0.0103000000
#   pi_all.Mean pi_all.3rd Qu.  pi_all.Max.     variance sample_size
#1 0.0029042938   0.0026000000 0.0546000000 3.221175e-05        3540
#2 0.0006760526   0.0007000000 0.0247000000 2.145909e-06        4180
#3 0.0029258455   0.0034000000 0.0369000000 1.343806e-05        2809
#4 0.0009759128   0.0012000000 0.0343000000 2.918171e-06        3533
#5 0.0070718796   0.0091000000 0.0694000000 4.425681e-05        1362
#6 0.0053439444   0.0072000000 0.0376000000 2.463182e-05        1511
#7 0.0152331875   0.0204000000 0.0819000000 9.072794e-05        1371
#8 0.0117455927   0.0157000000 0.0610000000 5.469923e-05        1645

#write.table(summarystats,"summary_stats_exon-intron_arms-centers.tsv",sep="\t", quote=FALSE,row.names=FALSE)




##
##codon positions
##

cod_pos_ele <- read.csv("codon_position_elegans_pi.csv", header=TRUE)

cod_pos_ino <- read.csv("codon_position_ino_pi.csv", header=TRUE)

cod_pos_ino$chr <- as.factor(cod_pos_ino$chr)

levels(cod_pos_ino$chr)[levels(cod_pos_ino$chr)=='Sp34_Chr1'] <- 'I'
levels(cod_pos_ino$chr)[levels(cod_pos_ino$chr)=='Sp34_Chr2'] <- 'II'
levels(cod_pos_ino$chr)[levels(cod_pos_ino$chr)=='Sp34_Chr3'] <- 'III'
levels(cod_pos_ino$chr)[levels(cod_pos_ino$chr)=='Sp34_Chr4'] <- 'IV'
levels(cod_pos_ino$chr)[levels(cod_pos_ino$chr)=='Sp34_Chr5'] <- 'V'
levels(cod_pos_ino$chr)[levels(cod_pos_ino$chr)=='Sp34_ChrX'] <- 'X'

cod_pos_ino$species <- "C. inopinata"
cod_pos_ele$species <- "C. elegans"

cod_pos_all <- rbind(cod_pos_ele,cod_pos_ino)

cod_pos_all$codon_position <- as.factor(cod_pos_all$codon_position)

levels(cod_pos_all$codon_position)[levels(cod_pos_all$codon_position)=='1'] <- 'One'
levels(cod_pos_all$codon_position)[levels(cod_pos_all$codon_position)=='2'] <- 'Two'
levels(cod_pos_all$codon_position)[levels(cod_pos_all$codon_position)=='3'] <- 'Three'


cod_pos_all$species <- as.factor(cod_pos_all$species)


cod_pos_all$startMB <- cod_pos_all$start/1000000

names(cod_pos_all)[names(cod_pos_all) == "codon_position"] <- "gen_region"

ele <- cod_pos_all[cod_pos_all$species == "C. elegans",]
ino <- cod_pos_all[cod_pos_all$species == "C. inopinata",]


ele_I <- ele[ele$chr == "I",]
ele_II <- ele[ele$chr == "II",]
ele_III <- ele[ele$chr == "III",]
ele_IV <- ele[ele$chr == "IV",]
ele_V <- ele[ele$chr == "V",]
ele_X <- ele[ele$chr == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$start)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$start)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$start)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$start)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$start)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$start)/8859471)/2

ino_I <- ino[ino$chr == "I",]
ino_II <- ino[ino$chr == "II",]
ino_III <- ino[ino$chr == "III",]
ino_IV <- ino[ino$chr == "IV",]
ino_V <- ino[ino$chr == "V",]
ino_X <- ino[ino$chr == "X",]

ino_I$norm_dist_center <- abs((10297276-ino_I$start)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$start)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$start)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$start)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$start)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$start)/9095254)/2

cod_pos_pi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)

#define chromosome arms and centers
cod_pos_pi_dat$chr_str_type <- ifelse(cod_pos_pi_dat$norm_dist_center >= 0.25,"Arms", "Centers")


cod_pos_pi_dat$species.gen_region <- paste(cod_pos_pi_dat$species,cod_pos_pi_dat$gen_region)

#cohen d effect size
cdes_cod_pos <- as.data.frame(cod_pos_pi_dat %>% cohens_d(pi_all ~ species.gen_region,ci = TRUE))
cdes_cod_pos

#      .y.             group1             group2       effsize   n1   n2
#1  pi_all     C. elegans One   C. elegans Three -0.2787857816 3920 3918
#2  pi_all     C. elegans One     C. elegans Two  0.0888572531 3920 3922
#3  pi_all     C. elegans One   C. inopinata One -0.6407700760 3920  709
#4  pi_all     C. elegans One C. inopinata Three -1.7247108159 3920  692
#5  pi_all     C. elegans One   C. inopinata Two -0.4534404354 3920  722
#6  pi_all   C. elegans Three     C. elegans Two  0.3293903207 3918 3922
#7  pi_all   C. elegans Three   C. inopinata One -0.1708301428 3918  709
#8  pi_all   C. elegans Three C. inopinata Three -1.2477953438 3918  692
#9  pi_all   C. elegans Three   C. inopinata Two  0.0009683785 3918  722
#10 pi_all     C. elegans Two   C. inopinata One -0.7307156767 3922  709
#11 pi_all     C. elegans Two C. inopinata Three -1.7843505267 3922  692
#12 pi_all     C. elegans Two   C. inopinata Two -0.5584333219 3922  722
#13 pi_all   C. inopinata One C. inopinata Three -1.2536877586  709  692
#14 pi_all   C. inopinata One   C. inopinata Two  0.2419723849  709  722
#15 pi_all C. inopinata Three   C. inopinata Two  1.4615224477  692  722
#   conf.low conf.high  magnitude
#1     -0.31     -0.25      small
#2      0.05      0.13 negligible
#3     -0.72     -0.58   moderate
#4     -1.83     -1.63      large
#5     -0.53     -0.38      small
#6      0.30      0.36      small
#7     -0.24     -0.11 negligible
#8     -1.35     -1.14      large
#9     -0.06      0.06 negligible
#10    -0.81     -0.66   moderate
#11    -1.89     -1.69      large
#12    -0.63     -0.49   moderate
#13    -1.37     -1.15      large
#14     0.14      0.34      small
#15     1.37      1.57      large

#write.table(cdes_cod_pos,"cohens_d_pi_means_codon_positions.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#percent differences
mean_df <- aggregate(cod_pos_pi_dat$pi_all,list(cod_pos_pi_dat$species.gen_region), mean)
#get pairwise percent differences among those means
mean_perc_diffm <- ((outer(mean_df$x, mean_df$x, `-`))/mean_df$x)*-100
#round
mean_perc_diffm <- round(mean_perc_diffm)
#get column and row ids right
colnames(mean_perc_diffm) <- mean_df$Group.1
rownames(mean_perc_diffm) <- mean_df$Group.1
mean_perc_diffm

#                   C. elegans One C. elegans Three C. elegans Two
#C. elegans One                  0              138            -21
#C. elegans Three              -58                0            -67
#C. elegans Two                 27              203              0
#C. inopinata One              -70              -29            -76
#C. inopinata Three            -92              -80            -93
#C. inopinata Two              -58                0            -67
#                   C. inopinata One C. inopinata Three C. inopinata Two
#C. elegans One                  233               1079              138
#C. elegans Three                 40                395                0
#C. elegans Two                  324               1399              202
#C. inopinata One                  0                254              -29
#C. inopinata Three              -72                  0              -80
#C. inopinata Two                 40                396                0

#write.table(mean_perc_diffm,"percent_differences_pi_means_codon_positions.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#wilcoxon rank sum test
as.data.frame(cod_pos_pi_dat %>% wilcox_test(pi_all ~ species.gen_region,p.adjust.method="BH"))

#      .y.             group1             group2   n1   n2  statistic         p
#1  pi_all     C. elegans One   C. elegans Three 3920 3918  5801268.5  2.25e-82
#2  pi_all     C. elegans One     C. elegans Two 3920 3922  8176150.5  2.09e-07
#3  pi_all     C. elegans One   C. inopinata One 3920  709   687682.5 7.84e-110
#4  pi_all     C. elegans One C. inopinata Three 3920  692   103846.0  0.00e+00
#5  pi_all     C. elegans One   C. inopinata Two 3920  722   938858.0  5.09e-51
#6  pi_all   C. elegans Three     C. elegans Two 3918 3922 10041584.5 3.08e-130
#7  pi_all   C. elegans Three   C. inopinata One 3918  709   948129.0  5.11e-42
#8  pi_all   C. elegans Three C. inopinata Three 3918  692   219313.5 6.99e-275
#9  pi_all   C. elegans Three   C. inopinata Two 3918  722  1221955.5  4.14e-09
#10 pi_all     C. elegans Two   C. inopinata One 3922  709   621068.0 3.40e-135
#11 pi_all     C. elegans Two C. inopinata Three 3922  692    80817.5  0.00e+00
#12 pi_all     C. elegans Two   C. inopinata Two 3922  722   865655.5  1.62e-69
#13 pi_all   C. inopinata One C. inopinata Three  709  692    68781.0 2.40e-120
#14 pi_all   C. inopinata One   C. inopinata Two  709  722   299233.0  2.50e-08
#15 pi_all C. inopinata Three   C. inopinata Two  692  722   449848.5 3.05e-150
#       p.adj p.adj.signif
#1   3.75e-82         ****
#2   2.09e-07         ****
#3  1.47e-109         ****
#4   0.00e+00         ****
#5   6.94e-51         ****
#6  7.70e-130         ****
#7   6.39e-42         ****
#8  3.49e-274         ****
#9   4.78e-09         ****
#10 1.02e-134         ****
#11  0.00e+00         ****
#12  2.43e-69         ****
#13 5.14e-120         ****
#14  2.68e-08         ****
#15 1.14e-149         ****

#write.table(as.data.frame(cod_pos_pi_dat %>% wilcox_test(pi_all ~ species.gen_region,p.adjust.method="BH")),"wilcox_test_codon_positions.tsv",sep="\t", quote=FALSE,row.names=TRUE)


summarystats <- aggregate(pi_all ~ species.gen_region,FUN=summary, data=cod_pos_pi_dat)

variancestats <- aggregate(pi_all ~ species.gen_region,FUN=var, data=cod_pos_pi_dat)

samplesizestats <- aggregate(pi_all ~ species.gen_region,FUN=length, data=cod_pos_pi_dat)

summarystats <- cbind(summarystats,variance = variancestats$pi_all,sample_size = samplesizestats$pi_all)
summarystats
#  species.gen_region  pi_all.Min. pi_all.1st Qu. pi_all.Median  pi_all.Mean
#1     C. elegans One 0.0000000000   0.0000000000  0.0002000000 0.0009696939
#2   C. elegans Three 0.0000000000   0.0000000000  0.0005000000 0.0023104135
#3     C. elegans Two 0.0000000000   0.0000000000  0.0000000000 0.0007624936
#4   C. inopinata One 0.0000000000   0.0005000000  0.0019000000 0.0032321580
#5 C. inopinata Three 0.0000000000   0.0052750000  0.0099000000 0.0114332370
#6   C. inopinata Two 0.0000000000   0.0000000000  0.0011000000 0.0023055402
#  pi_all.3rd Qu.  pi_all.Max.     variance sample_size
#1   0.0009000000 0.0466000000 6.481388e-06        3920
#2   0.0017000000 0.0723000000 3.977425e-05        3918
#3   0.0007000000 0.0449000000 4.393500e-06        3922
#4   0.0044000000 0.0361000000 1.845244e-05         709
#5   0.0155250000 0.0481000000 6.713169e-05         692
#6   0.0029000000 0.0248000000 1.087672e-05         722

#write.table(summarystats,"summary_stats_codon_positions.tsv",sep="\t", quote=FALSE,row.names=FALSE)

#genic/intergenic arms and centers

cod_pos_pi_dat$species.gen_region.chr_str_type <- paste(cod_pos_pi_dat$species,cod_pos_pi_dat$gen_region,cod_pos_pi_dat$chr_str_type)


#cohen d effect size
cohd <- as.data.frame(cod_pos_pi_dat %>% cohens_d(pi_all ~ species.gen_region.chr_str_type,ci = TRUE))
cohd





#write.table(cohd,"cohens_d_pi_means_codon_positions_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#codon positions, chromosome arms and centers


#genic/intergenic arms and centers

cod_pos_pi_dat$species.gen_region.chr_str_type <- paste(cod_pos_pi_dat$species,cod_pos_pi_dat$gen_region,cod_pos_pi_dat$chr_str_type)


#cohen d effect size
cohd <- as.data.frame(cod_pos_pi_dat %>% cohens_d(pi_all ~ species.gen_region.chr_str_type,ci = TRUE))
cohd

write.table(cohd,"cohens_d_pi_means_codon_positions_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#percent differences
mean_df <- aggregate(cod_pos_pi_dat$pi_all,list(cod_pos_pi_dat$species.gen_region.chr_str_type), mean)
#get pairwise percent differences among those means
mean_perc_diffm <- ((outer(mean_df$x, mean_df$x, `-`))/mean_df$x)*-100
#round
mean_perc_diffm <- round(mean_perc_diffm)
#get column and row ids right
colnames(mean_perc_diffm) <- mean_df$Group.1
rownames(mean_perc_diffm) <- mean_df$Group.1
mean_perc_diffm

#write.table(mean_perc_diffm,"percent_differences_pi_means_codon_positions_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)


#wilcoxon rank sum test
wilcod <- as.data.frame(cod_pos_pi_dat %>% wilcox_test(pi_all ~ species.gen_region.chr_str_type,p.adjust.method="BH"))
wilcod
#write.table(wilcod,"wilcox_test_codon_positions_arms-centers.tsv",sep="\t", quote=FALSE,row.names=TRUE)

#summary stats
summarystats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=summary, data=cod_pos_pi_dat)

variancestats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=var, data=cod_pos_pi_dat)

samplesizestats <- aggregate(pi_all ~ species.gen_region.chr_str_type,FUN=length, data=cod_pos_pi_dat)

summarystats <- cbind(summarystats,variance = variancestats$pi_all,sample_size = samplesizestats$pi_all)
summarystats

#write.table(summarystats,"summary_stats_codon_positions_arms-centers.tsv",sep="\t", quote=FALSE,row.names=FALSE)


###FST
dat <- read.csv("all_fst_combined.csv", header=TRUE)


#inopinata FST summary stats

summary(dat$FST)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#-0.438600 -0.009825  0.010200  0.016015  0.035100  0.749600


sd(dat$FST)

#[1] 0.04859959


dat$MB <- dat$BP_start/1000000

dat$pop_pair <- factor(dat$pop_pair, levels = c("ishi-irio","irio-yona","ishi-yona"))

#FST summary statistics  per pop pair
tapply(dat$FST,dat$pop_pair,summary)

#$`ishi-irio`
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#-0.369400 -0.006225  0.010000  0.015784  0.032400  0.629000
#
#$`irio-yona`
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.43860 -0.01920  0.00680  0.01308  0.03608  0.74960
#
#$`ishi-yona`
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#-0.381800 -0.006075  0.013300  0.019176  0.037200  0.556300


tapply(dat$FST,dat$pop_pair,sd)

# ishi-irio  irio-yona  ishi-yona
#0.03938055 0.05867758 0.04556242

tapply(dat$FST,dat$pop_pair,IQR)

#ishi-irio irio-yona ishi-yona
# 0.038625  0.055275  0.043275


###FIS


new_dat <- read.table("fis_elegans_and_inopinata.tsv", header = TRUE, sep = "\t")

tapply(new_dat$Fis,new_dat$species,summary)

#$elegans
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-0.3692  0.9618  1.0000  0.9360  1.0000  1.0000
#
#$inopinata
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.73214  0.01519  0.11752  0.14531  0.23510  1.00000



tapply(new_dat$Fis,new_dat$species,sd)

#  elegans inopinata
#0.1579432 0.2034266




###PCA



vcf <- read.vcfR("inopinata_24_biallelic_snps.vcf", verbose = FALSE)

x <- vcfR2genlight(vcf)

pca <- glPca(x , nf=20)

pca_df <- as.data.frame(pca$scores)

pca_df$Island <- c("Yonaguni", "Yonaguni", "Yonaguni", "Iriomote", "Iriomote", "Ishigaki", "Ishigaki", "Ishigaki", "Yonaguni", "Ishigaki", "Iriomote", "Ishigaki", "Ishigaki", "Iriomote", "Ishigaki", "Iriomote", "Ishigaki", "Iriomote", "Ishigaki", "Yonaguni", "Ishigaki", "Ishigaki", "Yonaguni", "Ishigaki")


pca$eig

# [1] 1.730104e+03 1.168930e+03 9.756727e+02 9.064308e+02 8.079734e+02
# [6] 7.410195e+02 6.105813e+02 5.811965e+02 5.210590e+02 4.881713e+02
#[11] 4.713327e+02 4.660548e+02 4.303158e+02 4.144107e+02 3.973199e+02
#[16] 3.438410e+02 3.384153e+02 2.848523e+02 2.415126e+02 1.750465e+02
#[21] 1.470126e+02 1.087125e+02 7.739248e+01 6.089259e-10

#variance explained
pca_tot_eig <- sum(pca$eig)

(pca$eig)/sum(pca$eig)

# [1] 1.392174e-01 9.406102e-02 7.851007e-02 7.293834e-02 6.501570e-02
# [6] 5.962808e-02 4.913203e-02 4.676751e-02 4.192838e-02 3.928198e-02
#[11] 3.792702e-02 3.750232e-02 3.462649e-02 3.334665e-02 3.197139e-02
#[16] 2.766807e-02 2.723148e-02 2.292139e-02 1.943395e-02 1.408558e-02
#[21] 1.182976e-02 8.747838e-03 6.227589e-03 4.899882e-14

#PC1 explains 14% of the variance; PC2 explains 9% of the variance

write.table(pca_df, file = "24_inopinata_worms_pca.tsv", quote = FALSE, sep = "\t")


#plots


ggplot(pca_df, aes(x = PC1, y = PC2)) + geom_point(aes(colour=Island)) + scale_colour_manual(values=c(c("#111A1A","#FF0000","#005AB3"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=14, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=16), legend.title=element_text(size = 15), legend.text = element_text(size = 14), legend.key = element_blank()) + xlab("PC1 (14%)") + ylab("PC2 (9%)") + scale_x_continuous(limits=c(-110,50),breaks=c(-100,-75,-50,-25,0,25,50))

ggsave("pca_plot_with_legend_6-30-21.png", width=6, height=6.5,units = "in")  


ggplot(pca_df, aes(x = PC1, y = PC2)) + geom_point(aes(colour=Island)) + scale_colour_manual(values=c(c("#111A1A","#FF0000","#005AB3"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=14, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=16), legend.position = "none") + xlab("PC1 (14%)") + ylab("PC2 (9%)") + scale_x_continuous(limits=c(-110,50),breaks=c(-100,-75,-50,-25,0,25,50))


ggsave("pca_plot_6-30-21.pdf", width=3.5, height=3,units = "in",useDingbats=FALSE)  



#DAPC clustering, #https://grunwaldlab.github.io/Population_Genetics_in_R/clustering_plot.html



#get data in there
vcf <- read.vcfR("inopinata_24_biallelic_snps.vcf", verbose = FALSE)

genl <- vcfR2genlight(vcf)


#get clusters
maxK <- 20
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(genl, n.pca = 23, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}



#plot BIC for k clusters

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

#get best k that maximizes BIC
my_df[which.max(my_df$BIC),]$K

#[1] 12

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)") 
p1

summary(my_df$BIC)

ggplot(my_df, aes(x = K, y = BIC)) + geom_dotplot(binaxis="y",stackdir="center",alpha=1,colour="#0570b0",fill="#0570b0",binwidth=.33) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5,size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14)) + xlab("Number of groups (K)") + ylab("BIC") + scale_y_continuous(limits=c(210,240),breaks=c(210,220,230,240))

ggsave("k_means_clustering_BIC.png", width=8, height=6,units = "in")  

ggsave("k_means_clustering_BIC.pdf", width=8, height=6,units = "in",useDingbats=FALSE)  


#DAPC
my_k <- 3:16

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(genl, n.pca = 23, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(genl, pop = grp_l[[i]]$grp, n.pca = 23, n.da = my_k[i])
#  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}


#plot classification of clusters for each individual

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]

samp_loc <- read.table("sample_locations.txt", header = TRUE, sep = "\t")

samp_loc$tree_id <- as.factor(samp_loc$tree_id)

tmp$sample_id <- rownames(tmp)
tmp$island <- samp_loc$island
tmp$tree_id <- samp_loc$tree_id
tmp$fig_id <- samp_loc$fig_id

tmp <- melt(tmp, id = c("sample_id", "K","island","tree_id","fig_id"))
names(tmp)[6:7] <- c("Group", "Posterior")
my_df <- tmp


for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$sample_id <- rownames(tmp)
  tmp$island <- samp_loc$island
  tmp$tree_id <- samp_loc$tree_id
  tmp$fig_id <- samp_loc$fig_id
  tmp <- melt(tmp, id = c("sample_id", "K","island","tree_id","fig_id"))
  names(tmp)[6:7] <- c("Group", "Posterior")
  my_df <- rbind(my_df, tmp)
}


#number of islands = 3
k3 <- my_df[my_df$K == "3",]
#number of plants = 5
k5 <- my_df[my_df$K == "5",]
#number of figs = 9
k9 <- my_df[my_df$K == "9",]
#best K by BIC = 12
k12 <- my_df[my_df$K == "12",]

plotdf <- rbind(k3,k5,k9,k12)

#write.table(my_df,file="dapc_clusters_k_3-16.tsv",sep='\t',quote=FALSE,row.names = FALSE)

plotdf$fig_id <- factor(plotdf$fig_id ,levels = c("89-4","89-3","89-2","87-3","72-17","72-16","44-2","44-1","22-2"))

plotdf$sample_id <-factor(plotdf$sample_id, levels = c("may_2017_H11.bam","may_2017_F11.bam","A05.bam","dec_2016_D10.bam","B03.bam","A03.bam","may_2017_D11.bam","may_2017_C11.bam","may_2017_A12.bam","may_2017_E11.bam","dec_2016_D02.bam","dec_2016_D01.bam","may_2017_H10.bam","may_2017_E12.bam","may_2017_D12.bam","may_2017_C12.bam","dec_2016_D05.bam","dec_2016_D04.bam","dec_2016_D03.bam","may_2017_H12.bam","may_2017_G12.bam","may_2017_B12.bam","may_2017_B11.bam","may_2017_A11.bam"))

colourCount = length(unique(k12$K))
getPalette = colorRampPalette(brewer.pal(12, "Set1"))

write.table(plotdf,file="dapc_clusters_for_figure.tsv",sep='\t',quote=FALSE,row.names = FALSE)


ggplot(plotdf, aes(x = sample_id, y = Posterior, fill = Group)) + geom_bar(stat = "identity") + scale_fill_manual(values=getPalette(12)) + facet_rep_wrap( ~ K, ncol=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=16), axis.text.x=element_blank(),axis.text.y=element_text(colour="black", size=15),strip.text.x = element_text(size=15), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Sample") + ylab("Posterior") +scale_y_continuous(limits=c(0,1),breaks=c(0,1)) 

#ggsave("color_plot_no_labels.pdf",height=6,width=8,units="in",useDingbats=FALSE)


