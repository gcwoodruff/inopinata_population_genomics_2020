library(effsize)
library(adegenet)
library(ggplot2)
library(vcfR)
library(reshape2)
library(RColorBrewer)
library(lemon)


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


