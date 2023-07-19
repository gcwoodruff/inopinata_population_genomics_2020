
#figures for Woodruff et al. 2023 C. inopinata population genomics paper.


library(ggplot2)
library(lemon)
library(ggforce)
library(patchwork)
library(RColorBrewer)
library(rstatix)
library(cowplot)
library(reshape2)

#pi, nucleotide diversity figures

dat <- read.csv("ino_elg_pi.csv", header=TRUE)


levels(dat$species)[levels(dat$species)=="elegans"] <- "C. elegans"
levels(dat$species)[levels(dat$species)=="inopinata"] <- "C. inopinata"


dat$species <-factor(dat$species, levels = c("C. elegans", "C. inopinata"))

dat$MB <- dat$start/1000000

#Figure 2A
ggplot(dat, aes(x = MB, y = pi_all)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ scaffold) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("π (Nucleotide diversity)") + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15))




#ggsave("inopinata_elegans_pi_genome.pdf", width=6.5, height=4,units = "in", useDingbats=FALSE)  
#ggsave("inopinata_elegans_pi_genome.png", width=6.5, height=4,units = "in")  


#normalize by chromosome position

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


pi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)

#define chromosome arms and centers
pi_dat$chr_str_type <- ifelse(pi_dat$norm_dist_center >= 0.25,"Arms", "Centers")


#Figure 2B
ggplot(pi_dat, aes(x = species, y = pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.2, alpha=0.2,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12, face="italic"), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=14), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=12, colour="black"),legend.key=element_blank(),legend.position = "none") + xlab("Species") + ylab("π (Nucleotide diversity)") + labs(colour = "Chromosome\nregion",size=12) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) + scale_y_continuous(limits=c(0,0.08),breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08))



#ggsave("inopinata_elegans_pi_sina_no_legend_8-7-20.pdf", width=4, height=3,units = "in", useDingbats=FALSE)  
#ggsave("inopinata_elegans_pi_sina_no_legend_8-7-20.png", width=4, height=3,units = "in")  



#genomic landscapes of nucleotide diversity, genic/intergenic & exon/intron, figure 3


exindat <- read.csv("exon_intron_pi.csv", header=TRUE,stringsAsFactors = TRUE)


exindat$startMB <- exindat$start/1000000

levels(exindat$species)[levels(exindat$species)=='elegans'] <- 'C. elegans'
levels(exindat$species)[levels(exindat$species)=='inopinata'] <- 'C. inopinata'
levels(exindat$gen_region)[levels(exindat$gen_region)=='exonic'] <- 'Exon'
levels(exindat$gen_region)[levels(exindat$gen_region)=='intronic'] <- 'Intron'


genexdat <- read.csv("genic_intergenic_pi.csv", header=TRUE,stringsAsFactors = TRUE)


genexdat$startMB <- genexdat$start/1000000



levels(genexdat$species)[levels(genexdat$species)=='elegans'] <- 'C. elegans'
levels(genexdat$species)[levels(genexdat$species)=='inopinata'] <- 'C. inopinata'
levels(genexdat$gen_region)[levels(genexdat$gen_region)=='genic'] <- 'Genic'
levels(genexdat$gen_region)[levels(genexdat$gen_region)=='intergenic'] <- 'Intergenic'





exinchrplot <- ggplot(exindat, aes(x = startMB, y = pi_all)) + geom_point(alpha=0.25, size=0.25,aes(colour=gen_region))  + stat_smooth(size=0.75, se=FALSE, aes(colour=gen_region)) + facet_rep_grid(species ~ scaffold) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title.y=element_text(size=14), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.key=element_blank()) + scale_colour_manual(values=c("#fc8d62","#8da0cb")) + xlab("Position (MB)") + ylab("π (Nucleotide diversity)") + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15)) + scale_y_continuous(limits=c(0,0.08),breaks=c(0,0.02,0.04,0.06,0.08)) + labs(colour="Gene\nregion")


genexchrplot <- ggplot(genexdat, aes(x = startMB, y = pi_all)) + geom_point(alpha=0.25, size=0.25,aes(colour=gen_region))  + stat_smooth(size=0.5, se=FALSE, aes(colour=gen_region)) + facet_rep_grid(species ~ scaffold) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_blank(),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.key=element_blank()) + scale_colour_manual(values=c("#4daf4a","#377eb8")) + xlab("Position (MB)") + ylab("π (Nucleotide diversity)") + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15)) + scale_y_continuous(limits=c(0,0.08),breaks=c(0,0.02,0.04,0.06,0.08)) + labs(colour="Genomic\nregion")


#this is Figure 3
exinchrplot/genexchrplot


#Figure 4


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


a <- ggplot(cod_pos_all[cod_pos_all$species=="C. elegans",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=gen_region),size=0.5, alpha=0.5,scale="width") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank()) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)) +ylab("π")

b <- ggplot(cod_pos_all[cod_pos_all$species=="C. inopinata",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=gen_region),size=0.5, alpha=0.5,scale="width") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3")) + ylim(0,0.05)


c <- a+b


d <- ggplot(exindat[exindat$species=="C. elegans",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=gen_region),size=0.5, alpha=0.5,scale="width") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank()) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#fc8d62","#8da0cb"))+ scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)) +ylab("π (Nucleotide diversity)")

e <- ggplot(exindat[exindat$species=="C. inopinata",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=gen_region),size=0.5, alpha=0.5,scale="width") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#fc8d62","#8da0cb")) + ylim(0,0.05)

f <- d+e



g <- ggplot(genexdat[genexdat$species=="C. elegans",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=gen_region),size=0.5, alpha=0.5,scale="width") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank()) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#4daf4a","#377eb8"))+ scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)) +ylab("π")

h <- ggplot(genexdat[genexdat$species=="C. inopinata",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=gen_region),size=0.5, alpha=0.5,scale="width") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#4daf4a","#377eb8")) + ylim(0,0.05)

i <- g+h


#this is figure 4
c/f/i

#figure 5, now with arm-center differences

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



genintgen_pi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)


#define chromosome arms and centers
genintgen_pi_dat$chr_str_type <- ifelse(genintgen_pi_dat$norm_dist_center >= 0.25,"Arms", "Centers")


#codon positions
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


j <- ggplot(cod_pos_pi_dat[cod_pos_pi_dat$species=="C. elegans",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank()) + scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)) + ylab("π")



k <- ggplot(cod_pos_pi_dat[cod_pos_pi_dat$species=="C. inopinata",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()) + ylim(0,0.05)

l <- j+k


m <- ggplot(exinpi_dat[exinpi_dat$species=="C. elegans",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank()) + scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)) + ylab("π (Nucleotide diversity)")



n <- ggplot(exinpi_dat[exinpi_dat$species=="C. inopinata",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()) + ylim(0,0.05)

o <- m+n

p <- ggplot(genintgen_pi_dat[genintgen_pi_dat$species=="C. elegans",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank()) + scale_y_continuous(limits=c(0,0.015),breaks=c(0,0.0025,0.005,0.0075,0.01,0.0125,0.015)) + ylab("π")



q <- ggplot(genintgen_pi_dat[genintgen_pi_dat$species=="C. inopinata",], aes(x=gen_region,y=pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()) + ylim(0,0.05)

r <- p+q

#this is figure 5
l/o/r

#figure 6



func_gen_eff <- read.table("functional_genomic_regions_effect_sizes.tsv", header=TRUE,sep="\t")


AC_func_gen_eff <- read.table("functional_genomic_regions_arm-center_effect_sizes.tsv", header=TRUE,sep="\t")


func_gen_eff$x_axis <- as.factor(paste(func_gen_eff$species, func_gen_eff$group1, func_gen_eff$group2))

func_gen_eff$x_axis <- factor(func_gen_eff$x_axis, levels=c("C. elegans Genic Intergenic", "C. elegans Exon Intron", "C. elegans One Two","C. elegans One Three","C. elegans Two Three","C. inopinata Genic Intergenic", "C. inopinata Exon Intron", "C. inopinata One Two","C. inopinata One Three","C. inopinata Two Three"))


func_gen_eff$x_axis2 <- as.factor(paste(func_gen_eff$group1, func_gen_eff$group2))


func_gen_eff$x_axis2 <- factor(func_gen_eff$x_axis2, levels=c("Genic Intergenic", "Exon Intron", "One Two","One Three","Two Three"))


#changing valence of one value to have everything positive, then will go and edit x axis accordingly
func_gen_eff[5, "effsize2"] <- func_gen_eff[5, "effsize2"]*-1
func_gen_eff[5, "conf_high"] <- func_gen_eff[5, "conf_high"]*-1
func_gen_eff[5, "conf_low"] <- func_gen_eff[5, "conf_low"]*-1

func_gen_eff[8, "effsize2"] <- func_gen_eff[8, "effsize2"]*-1
func_gen_eff[8, "conf_high"] <- func_gen_eff[8, "conf_high"]*-1
func_gen_eff[8, "conf_low"] <- func_gen_eff[8, "conf_low"]*-1

levels(func_gen_eff$x_axis2)[levels(func_gen_eff$x_axis2)=="Genic Intergenic"] <- "Intergenic-Genic"
levels(func_gen_eff$x_axis2)[levels(func_gen_eff$x_axis2)=="Exon Intron"] <- "Intron-Exon"
levels(func_gen_eff$x_axis2)[levels(func_gen_eff$x_axis2)=="One Two"] <- "One-Two"
levels(func_gen_eff$x_axis2)[levels(func_gen_eff$x_axis2)=="One Three"] <- "Three-One"
levels(func_gen_eff$x_axis2)[levels(func_gen_eff$x_axis2)=="Two Three"] <- "Three-Two"

a <- ggplot(func_gen_eff,aes(x=x_axis2, y=effsize2,fill=species)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin=conf_low, ymax=conf_high), colour="black", width=.1, position=position_dodge(0.9)) + theme_cowplot() + scale_fill_brewer(palette = "Set1")  + ylab("Pi among functional\ncategory effect size") + xlab ("Comparison") 

#ok, other plot


AC_func_gen_eff$gen_region <- as.factor(AC_func_gen_eff$gen_region)

AC_func_gen_eff$gen_region <- factor(AC_func_gen_eff$gen_region, levels=c("Genic","Intergenic","Exon","Intron","One","Two","Three"))
AC_func_gen_eff[is.na(AC_func_gen_eff)] <- "Intergenic"

b <- ggplot(AC_func_gen_eff,aes(x=gen_region, y=effsize,fill=species)) + geom_bar(stat="identity",position="dodge") + geom_errorbar(aes(ymin=conf.low, ymax=conf.high), colour="black", width=.1, position=position_dodge(0.9)) + theme_cowplot() + scale_fill_brewer(palette = "Set1")  + ylab("Pi arm-center\neffect size") + xlab ("Genomic functional category") 

#this is figure 6
a/b

#fstat figures (figure 7)



dat <- read.csv("all_fst_combined.csv", header=TRUE)


dat$MB <- dat$BP_start/1000000


dat$pop_pair <- factor(dat$pop_pair, levels = c("ishi-irio","irio-yona","ishi-yona"))





#fis
new_dat <- read.table("fis_elegans_and_inopinata.tsv", header = TRUE, sep = "\t")


new_dat$MB <- new_dat$BP/1000000

new_dat$species <- as.factor(new_dat$species)

levels(new_dat$species)[levels(new_dat$species)=="elegans"] <- "C. elegans"
levels(new_dat$species)[levels(new_dat$species)=="inopinata"] <- "C. inopinata"


ino_fis_dat <- new_dat[new_dat$species == "C. inopinata",]

ele_fis_dat <- new_dat[new_dat$species == "C. elegans",]

#this is supplemental figure 8
ggplot(new_dat, aes(x = MB, y = Fis)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab(expression(italic(F[IS]))) + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15)) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))



#ggsave("inopinata_elegans_fis_genome.pdf", width=6.5, height=4,units = "in", useDingbats=FALSE)  

#ggsave("inopinata_elegans_fis_genome.png", width=6.5, height=4,units = "in")  

#prep data for composite plot


ino_fis_prep <- data.frame(Chr=ino_fis_dat$Chr,BP=ino_fis_dat$BP,F=ino_fis_dat$Fis,MB=ino_fis_dat$MB)
ino_fis_prep$pop_pair <- "fis"

ino_fst_prep <- data.frame(Chr=dat$Chr,BP=dat$BP_start,F=dat$FST,MB=dat$MB,pop_pair=dat$pop_pair)

f_stats <- rbind(ino_fst_prep,ino_fis_prep)



x_labels <- c(expression(paste("Ir-Is ",italic(F[ST]))),expression(paste("Ir-Yo ",italic(F[ST]))),expression(paste("Is-Yo ",italic(F[ST]))),expression(italic(F[IS])))

#this is figure 7a
ggplot(f_stats, aes(x = MB, y = F)) + geom_point(aes(colour=pop_pair),alpha=0.25, size=0.15) + geom_smooth(aes(colour=pop_pair), size=0.75,se=FALSE) + facet_rep_wrap( ~ Chr, nrow=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.title=element_text(size = 12), legend.text = element_text(size = 11), legend.key = element_blank(),legend.position = "none") + xlab("Position (MB)") + ylab(expression(italic(F)))  + scale_colour_manual(values=c("#a6bddb","#2b8cbe","#045a8d","red")) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + scale_x_continuous(breaks=c(0,10,20))



#ggsave("f_genome_plot.png", width=7.25, height=4.75, units = "in")


#composite plot




a <- ggplot(f_stats, aes(x = MB, y = F)) + geom_point(aes(colour=pop_pair),alpha=0.25, size=0.15) + geom_smooth(aes(colour=pop_pair), size=0.75,se=FALSE) + facet_rep_wrap( ~ Chr, nrow=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.title=element_text(size = 12), legend.text = element_text(size = 11), legend.key = element_blank(),legend.position = "none") + xlab("Position (MB)") + ylab(expression(italic(F)))  + scale_colour_manual(values=c("#a6bddb","#2b8cbe","#045a8d","red")) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + scale_x_continuous(breaks=c(0,10,20))



b <- ggplot(f_stats, aes(x = pop_pair, y = F)) + geom_sina(size=0.15, alpha=0.25,scale="width",aes(colour=pop_pair)) + scale_colour_manual(values=c("#a6bddb","#2b8cbe","#045a8d","red")) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13), axis.text.y = element_text(colour="black", size=13),legend.position = "none",axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black")) + xlab(expression(paste(italic(F)," statistic"))) + ylab(expression(italic(F))) + scale_x_discrete(labels= x_labels) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))

#this is figure 7
a/b

#ggsave("f_composite_plot.png", width=7, height=5.5, units = "in")

#ggsave("f_composite_plot.pdf", width=7, height=5.5, units = "in", useDingbats=FALSE)


#natural history worm histograms


#get data in there; data from 2018 woodruff and phillips 2018 bmc ecology paper
ecol_dat <- read.table("okinawa_septica_stage_worm_presence_data_R_edit_5-3-18.csv", sep="\t", header=T)

#getting data types right
ecol_dat$Latitude <- as.numeric(as.matrix(ecol_dat$Latitude))
ecol_dat$Longitude <- as.numeric(as.matrix(ecol_dat$Longitude))
ecol_dat$Fig_stage <- as.factor(ecol_dat$Fig_stage)
ecol_dat$Plant <- as.factor(ecol_dat$Plant)
ecol_dat$worm_presence <- as.factor(ecol_dat$worm_presence)

#set aside pollinated and unpollinated

ecol_dat$pollination[ecol_dat$Fig_stage != "1"] <- "yes"
ecol_dat$pollination <- as.factor(ecol_dat$pollination)
levels(ecol_dat$pollination) <- c("yes", "no")
ecol_dat$pollination[ecol_dat$Fig_stage != "1"] <- "yes"
ecol_dat$pollination[ecol_dat$Fig_stage == "1"] <- "no"

ecol_dat[,c(8,18)]

pollinated <- ecol_dat[ecol_dat$pollination == "yes",]
unpollinated <- ecol_dat[ecol_dat$pollination == "no",]


#hm, some missing data here, let's assume that if the thing is unpollinated, then there are no foundresses...

ecol_dat_b <- ecol_dat

ecol_dat_b$foundresses[ecol_dat_b$pollination == "no"] <- 0


levels(ecol_dat_b$pollination) <- c("Pollinated figs","Unpollinated figs")

pollinated_figs_foundress_data <- ecol_dat_b[ecol_dat_b$pollination == "Pollinated figs",]$foundresses

worms_per_wasp_data <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,3,3,5,6)
	#this is actually the data, see 2018 woodruff and phillips 2018 bmc ecology paper


a <- ggplot() + aes(pollinated_figs_foundress_data) + geom_histogram(binwidth=1, fill="#111A1A",colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.position="none", strip.text = element_text(colour="black", size=12), strip.background = element_rect(colour="white", fill="white"),axis.title=element_text(size=13), axis.ticks = element_line(colour = "black")) + xlab("Number of wasps per fig") + ylab("Frequency") + scale_x_continuous(breaks=c(0:11)) + scale_y_continuous(limits=c(0,60),breaks=c(0,10,20,30,40,50,60))


b <- ggplot() + aes(worms_per_wasp_data) + geom_histogram(binwidth=1, fill="#FF1A00",colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.position="none", strip.text = element_text(colour="black", size=12), strip.background = element_rect(colour="white", fill="white"),axis.title=element_text(size=13), axis.ticks = element_line(colour = "black")) + xlab("Number of worms per wasp") + ylab("Frequency") + scale_x_continuous(breaks=c(0:6)) + scale_y_continuous(limits=c(0,20),breaks=c(0,5,10,15,20))


#figure 8a & figure 8b
a/b

#ggsave("natural_history_histograms.png",height=5,width=4,units="in")
#ggsave("natural_history_histograms.pdf",height=5,width=4,units="in")


#comparison of worm F stats with fig wasp F stats


dat <- read.table("bisulcatus_data.tsv", header = TRUE, sep = "\t")

#figure 8c
ggplot(dat, aes(x = f_stat, y = value)) + geom_hline(yintercept=0.14531,size=0.8,linetype="dashed") + geom_hline(yintercept=0.016015,size=0.8,linetype="dotted") + geom_sina(aes(colour=paper),size=1, alpha=1,scale="width") + stat_summary(aes(group=paper),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13, face="italic"), axis.text.y = element_text(colour="black", size=13),legend.position = "none",axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black"),legend.key=element_blank()) + scale_y_continuous(limits=c(-0.1,0.5),breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5)) + ylab(expression(italic(F))) + xlab("Pollinating fig wasp F statistic")

#ggsave("bisculatus_sina_no_legend.pdf",height=5,width=4,useDingbats=FALSE)

#supplemental figure , coverage


before_dat <- read.table("all_depth_before_genotype_filter.tsv",sep="\t",header=TRUE)

EcoRI_sites <- read.table("inopinata_EcoRI_sites.tsv",sep="\t",header=TRUE)

EcoRI_sites$Chr <- as.factor(EcoRI_sites$Chr)

levels(EcoRI_sites$Chr)[levels(EcoRI_sites$Chr)=="I"] <- "Sp34_Chr1"
levels(EcoRI_sites$Chr)[levels(EcoRI_sites$Chr)=="II"] <- "Sp34_Chr2"
levels(EcoRI_sites$Chr)[levels(EcoRI_sites$Chr)=="III"] <- "Sp34_Chr3"
levels(EcoRI_sites$Chr)[levels(EcoRI_sites$Chr)=="IV"] <- "Sp34_Chr4"
levels(EcoRI_sites$Chr)[levels(EcoRI_sites$Chr)=="V"] <- "Sp34_Chr5"
levels(EcoRI_sites$Chr)[levels(EcoRI_sites$Chr)=="X"] <- "Sp34_ChrX"


EcoRI_sites$BP_start.plusone <- EcoRI_sites$BP_start+1

EcoRI_sites$chr.bp_start <- as.factor(paste(EcoRI_sites$Chr,EcoRI_sites$BP_start.plusone))

EcoRI_sites_ii <- EcoRI_sites[EcoRI_sites$Chr != "CSP34.Sp34_mitochondrion",]

EcoRI_sites_ii$EcoRI_cut_sites <- as.numeric(gsub('\\.', 0, EcoRI_sites_ii$EcoRI_cut_sites))

summary(EcoRI_sites_ii$EcoRI_cut_sites)

#figure

before_dat$BP_start.plusone <- before_dat$BP_start+1

before_dat$chr.bp_start <- as.factor(paste(before_dat$Chr,before_dat$BP_start.plusone))

summary_before_dat <- before_dat %>% group_by(chr.bp_start) %>% get_summary_stats(Cov, type="common")

mean_before_dat <- data.frame(chr.bp_start = summary_before_dat$chr.bp_start, cov= summary_before_dat$mean)

out <- do.call(rbind,strsplit(as.character(mean_before_dat$chr.bp_start),' '))

mean_before_dat$Chr <- as.factor(out[,1])
mean_before_dat$BP <- as.numeric(out[,2])


mean_before_dat$MB <- mean_before_dat$BP/1000000

mean_before_dat_ii <- mean_before_dat[mean_before_dat$Chr != "CSP34.Sp34_mitochondrion",]


fig_df <- merge(x = mean_before_dat_ii, y = EcoRI_sites_ii, by = "chr.bp_start", all = TRUE)

fig_df_ii <- data.frame(Chr = fig_df$Chr.x, MB = fig_df$MB, Cov= fig_df$cov, EcoRI_cut_sites= fig_df$EcoRI_cut_sites)

fig_df_ii$coverage_per_cut_site = fig_df$cov/fig_df$EcoRI_cut_sites


fig_df_melt <- melt(fig_df_ii, id=c("Chr","MB"))

levels(fig_df_melt$variable)[levels(fig_df_melt$variable)=="Cov"] <- "Coverage"
levels(fig_df_melt$variable)[levels(fig_df_melt$variable)=="EcoRI_cut_sites"] <- "EcoRI cut sites"
levels(fig_df_melt$variable)[levels(fig_df_melt$variable)=="coverage_per_cut_site"] <- "Coverage per cut site"


fig_df_melt$variable <-factor(fig_df_melt$variable, levels = c("Coverage", "EcoRI cut sites", "Coverage per cut site"))

fig_df_melt[is.na(fig_df_melt$value),]

fig_df_ii[is.infinite(fig_df_ii$coverage_per_cut_site),]



#supplemental figure 1
ggplot(fig_df_melt, aes(x = MB, y = value)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(variable ~ Chr,scales="free_y") + theme_cowplot() + theme(strip.background = element_blank()) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15))


#ggsave("coverage.png", width=7, height=6,units = "in")  



#supplemental figure, nucleotide diversity across Caenorhabditis genus


dat<- read.table("pi_estimates_caenorhabditis_with_inopinata_for_presentation_10-16-20.tsv", sep="\t", header=T)


species_order <- c("C. brenneri","C. sinica","C. remanei","C. latens","C. japonica","C. inopinata","C. briggsae","C. elegans","C. tropicalis")

dat$species <- factor(dat$species, levels = c("C. brenneri","C. sinica","C. remanei","C. latens","C. japonica","C. inopinata","C. briggsae","C. elegans","C. tropicalis"))

#supplemental figure 2
ggplot(dat,aes(x=reorder(species,-PI), y=PI)) + geom_bar(stat="identity",aes(fill=reproductive_mode)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=14,face="italic",angle = 45, hjust =1), axis.text.y = element_text(colour="black", size=15),legend.position="none", strip.text.x = element_text(size=15),axis.title=element_text(size=17), strip.background = element_rect(colour="white", fill="white"),axis.ticks = element_line(colour = "black")) + xlab("Species") + ylab("π (nucleotide diversity)") + scale_fill_brewer(breaks=c("gon", "self"), labels=c("Male/Female", "Selfer")) + annotate("segment", x=0.5, xend=6.5, y=.17, yend=.17, arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) + annotate("segment", x=6.6, xend=9.6, y=.17, yend=.17, arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) + annotate("text", x=3.5, y=0.177, label="Obligate outcrossers", size=6) + annotate("text", x="C. briggsae", y=0.176, label="Selfers", size=6)


#ggsave("caeno_pi_WITH_inopinata.png", width=7, height=6, units = "in")




#supplemental figure 3 , codon fold degeneracy


#Fold degeneracy
fol_deg_ele <- read.csv("degeneracy_elegans_pi.csv", header=TRUE)

fol_deg_ino <- read.csv("degeneracy_ino_pi.csv", header=TRUE)

fol_deg_ino$chr <- as.factor(fol_deg_ino$chr)

levels(fol_deg_ino$chr)[levels(fol_deg_ino$chr)=='Sp34_Chr1'] <- 'I'
levels(fol_deg_ino$chr)[levels(fol_deg_ino$chr)=='Sp34_Chr2'] <- 'II'
levels(fol_deg_ino$chr)[levels(fol_deg_ino$chr)=='Sp34_Chr3'] <- 'III'
levels(fol_deg_ino$chr)[levels(fol_deg_ino$chr)=='Sp34_Chr4'] <- 'IV'
levels(fol_deg_ino$chr)[levels(fol_deg_ino$chr)=='Sp34_Chr5'] <- 'V'
levels(fol_deg_ino$chr)[levels(fol_deg_ino$chr)=='Sp34_ChrX'] <- 'X'

fol_deg_ino$species <- "C. inopinata"
fol_deg_ele$species <- "C. elegans"

fol_deg_all <- rbind(fol_deg_ele,fol_deg_ino)

fol_deg_all$fold_degeneracy <- as.factor(fol_deg_all$fold_degeneracy)

levels(fol_deg_all$fold_degeneracy)[levels(fol_deg_all$fold_degeneracy)=='0'] <- 'Zero'
levels(fol_deg_all$fold_degeneracy)[levels(fol_deg_all$fold_degeneracy)=='2'] <- 'Two'
levels(fol_deg_all$fold_degeneracy)[levels(fol_deg_all$fold_degeneracy)=='4'] <- 'Four'


#not enough for genomic landscape, but what about all together?


ele <- fol_deg_all[fol_deg_all$species == "C. elegans",]
ino <- fol_deg_all[fol_deg_all$species == "C. inopinata",]


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


fol_deg_pi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)

#define chromosome arms and centers
fol_deg_pi_dat$chr_str_type <- ifelse(fol_deg_pi_dat$norm_dist_center >= 0.25,"Arms", "Centers")

#this is supplemental figure 3
ggplot(fol_deg_pi_dat, aes(x = fold_degeneracy, y = pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + facet_rep_wrap( ~ species,nrow=1) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=14), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=12, colour="black",face="italic"),legend.key=element_blank()) + xlab("Codon Site Fold-degeneracy") + ylab("π (Nucleotide diversity)") + labs(colour = "Chromosome\nregion",size=12) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) + scale_y_continuous(limits=c(0,0.09),breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09))


#supplemental figure 4, LD


dat <- read.table("LD_ele_ino.tsv", header=TRUE,sep="\t")



dat$MB <- dat$BP/1000000

#this is supplemental figure 4
ggplot(dat, aes(x = MB, y = R2)) + geom_point(alpha=0.5, size=0.5,colour="black")  + stat_smooth(size=0.5, se=FALSE,colour="blue") + facet_rep_grid(Species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Intrachromosomal R2") + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15))


#supplemental figure 5, LD sina


ele <- dat[dat$Species == "C. elegans",]
ino <- dat[dat$Species == "C. inopinata",]


ele_I <- ele[ele$Chr == "I",]
ele_II <- ele[ele$Chr == "II",]
ele_III <- ele[ele$Chr == "III",]
ele_IV <- ele[ele$Chr == "IV",]
ele_V <- ele[ele$Chr == "V",]
ele_X <- ele[ele$Chr == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$BP)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$BP)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$BP)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$BP)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$BP)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$BP)/8859471)/2

ino_I <- ino[ino$Chr == "I",]
ino_II <- ino[ino$Chr == "II",]
ino_III <- ino[ino$Chr == "III",]
ino_IV <- ino[ino$Chr == "IV",]
ino_V <- ino[ino$Chr == "V",]
ino_X <- ino[ino$Chr == "X",]

ino_I$norm_dist_center <- abs((10297276-ino_I$BP)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$BP)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$BP)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$BP)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$BP)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$BP)/9095254)/2


R2_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)

#define Chromosome arms and centers
R2_dat$chr_str_type <- ifelse(R2_dat$norm_dist_center >= 0.25,"Arms", "Centers")

#this is supplemental figure 5
ggplot(R2_dat, aes(x = Species, y = R2)) + geom_sina(aes(colour=chr_str_type),size=0.5, alpha=0.5,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=14), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=12, colour="black",face="italic"),legend.key=element_blank()) + xlab("Species") + labs(colour = "Chromosome\nregion",size=12) + ylab("Intrachromosomal R2")  + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) 



#Supplemental figure, PCA


pca_df <- read.table("24_inopinata_worms_pca.tsv", header = TRUE, sep = "\t")


#plots


ggplot(pca_df, aes(x = PC1, y = PC2)) + geom_point(aes(colour=Island)) + scale_colour_manual(values=c(c("#111A1A","#FF0000","#005AB3"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=14, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=16), legend.position = "none") + xlab("PC1 (14%)") + ylab("PC2 (9%)") + scale_x_continuous(limits=c(-110,50),breaks=c(-100,-75,-50,-25,0,25,50))


#ggsave("pca_plot_6-30-21.pdf", width=3.5, height=3,units = "in",useDingbats=FALSE)  




#Supplemental figure, DAPC clustering


plotdf <- read.table("dapc_clusters_for_figure.tsv", header = TRUE, sep = "\t")

plotdf$K <- as.factor(plotdf$K)
plotdf$Group <- as.factor(plotdf$Group)

#best K by BIC = 12
k12 <- plotdf[plotdf$K == "12",]


colourCount = length(unique(k12$K))
getPalette = colorRampPalette(brewer.pal(12, "Set1"))


plotdf$fig_id <- factor(plotdf$fig_id ,levels = c("89-4","89-3","89-2","87-3","72-17","72-16","44-2","44-1","22-2"))

plotdf$sample_id <-factor(plotdf$sample_id, levels = c("may_2017_H11.bam","may_2017_F11.bam","A05.bam","dec_2016_D10.bam","B03.bam","A03.bam","may_2017_D11.bam","may_2017_C11.bam","may_2017_A12.bam","may_2017_E11.bam","dec_2016_D02.bam","dec_2016_D01.bam","may_2017_H10.bam","may_2017_E12.bam","may_2017_D12.bam","may_2017_C12.bam","dec_2016_D05.bam","dec_2016_D04.bam","dec_2016_D03.bam","may_2017_H12.bam","may_2017_G12.bam","may_2017_B12.bam","may_2017_B11.bam","may_2017_A11.bam"))


ggplot(plotdf, aes(x = sample_id, y = Posterior, fill = Group)) + geom_bar(stat = "identity") + scale_fill_manual(values=getPalette(12)) + facet_rep_wrap( ~ K, ncol=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=16), axis.text.x=element_blank(),axis.text.y=element_text(colour="black", size=15),strip.text.x = element_text(size=15), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Sample") + ylab("Posterior") +scale_y_continuous(limits=c(0,1),breaks=c(0,1)) 

#ggsave("DAPC_cluster_plot_no_labels.pdf",height=6,width=8,units="in",useDingbats=FALSE)


#supplemental figure for BIC plot following DAPC in file "statistics.R"

#supplemental figure 8


new_dat <- read.table("fis_elegans_and_inopinata.tsv", header = TRUE, sep = "\t")


new_dat$MB <- new_dat$BP/1000000

new_dat$species <- as.factor(new_dat$species)

levels(new_dat$species)[levels(new_dat$species)=="elegans"] <- "C. elegans"
levels(new_dat$species)[levels(new_dat$species)=="inopinata"] <- "C. inopinata"


ino_fis_dat <- new_dat[new_dat$species == "C. inopinata",]

ele_fis_dat <- new_dat[new_dat$species == "C. elegans",]

#this is supplemental figure 8
ggplot(new_dat, aes(x = MB, y = Fis)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab(expression(italic(F[IS]))) + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15)) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))


#supplemental figure 9


dat <- read.csv("elegans_random_rad_pi.csv",header=T)

dat$perc_RAD <- as.factor(dat$perc_RAD)

dat$perc_RAD <- factor(dat$perc_RAD, levels=c("20%","30%","40%","50%","60%","70%","80%","90%","100%","WGS"))


aggregate(Pi ~ perc_RAD, data=dat,FUN="mean")

#   perc_RAD          Pi
#1       20% 0.000250000
#2       30% 0.002697222
#3       40% 0.001910843
#4       50% 0.001370192
#5       60% 0.001344363
#6       70% 0.001231000
#7       80% 0.001244982
#8       90% 0.001308101
#9      100% 0.001244584
#10      WGS 0.002534557



	#mean C. elegans pi from https://elifesciences.org/articles/62587 , https://cdn.elifesciences.org/articles/62587/elife-62587-fig2-data1-v2.tsv.zip

#this is supplemental figure 9
ggplot(dat, aes(x = perc_RAD, y = Pi)) + geom_sina(scale="width",alpha=0.25,size=0.75) +geom_hline(yintercept=0.005318837452971,colour="blue",linetype="dotted") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9))  + scale_colour_brewer(palette="Set1") + theme_cowplot()  + ylab("Pi")+ ylim(0,0.05) + ylab("π (Nucleotide diversity)") + xlab("% N2 EcoRI sites retained")
