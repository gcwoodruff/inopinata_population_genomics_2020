#set that working directory

setwd("/Users/gavin/genome/pop_gen_revisions_9-2023/03_sed")

#get data in there

inop <- read.table("inopinata_24_num_SNPs_1kb_windows.bed",sep="\t")

#is NA?

inop$V5 <- is.na(inop$V4)

#counts of conecutive NA and non-NA

inop$V6 <- sequence(rle(as.character(inop$V5))$lengths)

#alright, get just the non-NA

inop_snp_windows <- inop[inop$V5 == FALSE,]

#get places where the counter is > 8

inop_snp_windows_nine <- inop_snp_windows[inop_snp_windows$V6 > 8,]

#not many, but something!

#             V1       V2       V3 V4    V5 V6
#3678  Sp34_Chr1  3677000  3678000 21 FALSE  9
#3679  Sp34_Chr1  3678000  3679000  8 FALSE 10
#3680  Sp34_Chr1  3679000  3680000  3 FALSE 11
#3681  Sp34_Chr1  3680000  3681000 30 FALSE 12
#5346  Sp34_Chr1  5345000  5346000  8 FALSE  9 #5
#10621 Sp34_Chr1 10620000 10621000  4 FALSE  9 #6
#37329 Sp34_Chr2 16733000 16734000  3 FALSE  9 #7
#42074 Sp34_Chr3  1361000  1362000 10 FALSE  9 #8
#45818 Sp34_Chr3  5105000  5106000 44 FALSE  9 #9
#58825 Sp34_Chr3 18112000 18113000 56 FALSE  9 #10
#58826 Sp34_Chr3 18113000 18114000  4 FALSE 10 #11
#64837 Sp34_Chr4  4687000  4688000 16 FALSE  9 #12
#72751 Sp34_Chr4 12601000 12602000  9 FALSE  9 #13
#72752 Sp34_Chr4 12602000 12603000  1 FALSE 10 #14
#91740 Sp34_Chr5 10572000 10573000  6 FALSE  9 #15
#91741 Sp34_Chr5 10573000 10574000  3 FALSE 10 #16
#91742 Sp34_Chr5 10574000 10575000  5 FALSE 11 #17
#94176 Sp34_Chr5 13008000 13009000  1 FALSE  9 #18
#94177 Sp34_Chr5 13009000 13010000 13 FALSE 10 #19

nrow(inop_snp_windows_nine[inop_snp_windows_nine$V6 == 9,])
#[1] 11
#11 windows

inop_snp_windows_nine$window <- seq(1,nrow(inop_snp_windows_nine),1)


#Define some ranges to extract
upper<- as.integer(rownames(inop_snp_windows_nine)) - 8
lower<- as.integer(rownames(inop_snp_windows_nine))

inop_consec_hyper <- NULL
some_rows <- NULL
number_rows_greater_15_snps <- NULL
newdf <- NULL

for (row in 1:nrow(inop_snp_windows_nine)) {
 some_rows <- as.character(upper[row]:lower[row])
 the_inop_some_rows <- inop[some_rows,]
 number_rows_greater_15_snps <- nrow(the_inop_some_rows[the_inop_some_rows$V4 > 15,])
 newdf <- data.frame(upper_row_name = upper[row],lower_row_name=lower[row],number_rows_greater_15_snps=number_rows_greater_15_snps)
 inop_consec_hyper <- rbind(inop_consec_hyper,newdf)
}

inop_consec_hyper

#   upper_row_name lower_row_name number_rows_greater_15_snps
#1            3670           3678                           3
#2            3671           3679                           3
#3            3672           3680                           3
#4            3673           3681                           4
#5            5338           5346                           2
#6           10613          10621                           1
#7           37321          37329                           2
#8           42066          42074                           4
#9           45810          45818                           2
#10          58817          58825                           6
#11          58818          58826                           5
#12          64829          64837                           7
#13          72743          72751                           2
#14          72744          72752                           2
#15          91732          91740                           0
#16          91733          91741                           0
#17          91734          91742                           0
#18          94168          94176                           0
#19          94169          94177                           0

#ZERO hyperdivergent regions by the Lee et al. 2021 definition. Unsurprising, considering the definition assumes WGS.

#and to be sure,

inop_consec_hyper[inop_consec_hyper$number_rows_greater_15_snps>8,]

#let's look at the elegans pseudo-rad data.




elerad <- read.table("elegans_24_og_pseudo_rad_num_SNPs_1kb_windows.bed",sep="\t")

#is NA?

elerad$V5 <- is.na(elerad$V4)

#counts of conecutive NA and non-NA

elerad$V6 <- sequence(rle(as.character(elerad$V5))$lengths)

#alright, get just the non-NA

elerad_snp_windows <- elerad[elerad$V5 == FALSE,]

#get places where the counter is > 8

elerad_snp_windows_nine <- elerad_snp_windows[elerad_snp_windows$V6 > 8,]


#       V1       V2       V3 V4    V5 V6
#10275   I 10274000 10275000 24 FALSE  9
#10276   I 10275000 10276000 22 FALSE 10
#10277   I 10276000 10277000 28 FALSE 11
#17538  II  2464000  2465000  5 FALSE  9
#17539  II  2465000  2466000 13 FALSE 10
#17540  II  2466000  2467000 58 FALSE 11
#17586  II  2512000  2513000  5 FALSE  9
#17587  II  2513000  2514000  6 FALSE 10
#17588  II  2514000  2515000 10 FALSE 11
#17589  II  2515000  2516000 20 FALSE 12
#17590  II  2516000  2517000  2 FALSE 13
#18704  II  3630000  3631000 13 FALSE  9
#18705  II  3631000  3632000  2 FALSE 10
#18706  II  3632000  3633000  8 FALSE 11
#18707  II  3633000  3634000  6 FALSE 12
#26401  II 11327000 11328000  4 FALSE  9
#26402  II 11328000 11329000  3 FALSE 10
#27384  II 12310000 12311000  8 FALSE  9
#29293  II 14219000 14220000  7 FALSE  9
#38309 III  7955000  7956000  1 FALSE  9
#38310 III  7956000  7957000  1 FALSE 10
#40308 III  9954000  9955000  4 FALSE  9
#46758  IV  2620000  2621000  5 FALSE  9
#50754  IV  6616000  6617000  7 FALSE  9
#50755  IV  6617000  6618000  9 FALSE 10
#50756  IV  6618000  6619000 10 FALSE 11
#50997  IV  6859000  6860000  1 FALSE  9
#50998  IV  6860000  6861000  1 FALSE 10
#65129   V  3483000  3484000  5 FALSE  9
#65782   V  4136000  4137000  2 FALSE  9
#65783   V  4137000  4138000 19 FALSE 10
#65784   V  4138000  4139000  3 FALSE 11
#67309   V  5663000  5664000  2 FALSE  9
#69416   V  7770000  7771000  3 FALSE  9
#69417   V  7771000  7772000  3 FALSE 10
#89493   X  6922000  6923000  4 FALSE  9
#89494   X  6923000  6924000 11 FALSE 10
#89495   X  6924000  6925000  6 FALSE 11
#96856   X 14285000 14286000  3 FALSE  9


nrow(elerad_snp_windows[elerad_snp_windows$V6 == 9,])
#18 windows

elerad_snp_windows_nine$window <- seq(1,nrow(elerad_snp_windows_nine),1)

upper<- as.integer(rownames(elerad_snp_windows_nine)) - 8
lower<- as.integer(rownames(elerad_snp_windows_nine))


elerad_consec_hyper <- NULL


some_rows <- as.character(upper[1]:lower[1])
the_elerad_some_rows <- elerad[some_rows,]
number_rows_greater_15_snps <- nrow(the_elerad_some_rows[the_elerad_some_rows$V4 > 15,])
newdf <- data.frame(upper_row_name = upper[1],lower_row_name=lower[1],number_rows_greater_15_snps=number_rows_greater_15_snps)
elerad_consec_hyper <- rbind(elerad_consec_hyper,newdf)


#[1] 8
#not 9!

#okay, let's keep going, but this time, for loop


elerad_consec_hyper <- NULL
some_rows <- NULL
number_rows_greater_15_snps <- NULL
newdf <- NULL

for (row in 1:nrow(elerad_snp_windows_nine)) {
 some_rows <- as.character(upper[row]:lower[row])
 the_elerad_some_rows <- elerad[some_rows,]
 number_rows_greater_15_snps <- nrow(the_elerad_some_rows[the_elerad_some_rows$V4 > 15,])
 newdf <- data.frame(upper_row_name = upper[row],lower_row_name=lower[row],number_rows_greater_15_snps=number_rows_greater_15_snps)
 elerad_consec_hyper <- rbind(elerad_consec_hyper,newdf)
}

elerad_consec_hyper

#   upper_row_name lower_row_name number_rows_greater_15_snps
#1           10267          10275                           8
#2           10268          10276                           8
#3           10269          10277                           8
#4           17530          17538                           3
#5           17531          17539                           2
#6           17532          17540                           3
#7           17578          17586                           0
#8           17579          17587                           0
#9           17580          17588                           0
#10          17581          17589                           1
#11          17582          17590                           1
#12          18696          18704                           0
#13          18697          18705                           0
#14          18698          18706                           0
#15          18699          18707                           0
#16          26393          26401                           0
#17          26394          26402                           0
#18          27376          27384                           0
#19          29285          29293                           0
#20          38301          38309                           0
#21          38302          38310                           0
#22          40300          40308                           4
#23          46750          46758                           1
#24          50746          50754                           2
#25          50747          50755                           2
#26          50748          50756                           2
#27          50989          50997                           0
#28          50990          50998                           0
#29          65121          65129                           0
#30          65774          65782                           2
#31          65775          65783                           3
#32          65776          65784                           3
#33          67301          67309                           0
#34          69408          69416                           0
#35          69409          69417                           0
#36          89485          89493                           0
#37          89486          89494                           0
#38          89487          89495                           0
#39          96848          96856                           0

#none!
#and to be sure,

elerad_consec_hyper[elerad_consec_hyper$number_rows_greater_15_snps>8,]
#[1] upper_row_name              lower_row_name
#[3] number_rows_greater_15_snps
#<0 rows> (or 0-length row.names)

#yep, none.




#let's look at the elegans WGS data.



elewgs <- read.table("elegans_24_wgs_num_SNPs_1kb_windows.bed",sep="\t")

#is NA?

elewgs$V5 <- is.na(elewgs$V4)

#counts of conecutive NA and non-NA

elewgs$V6 <- sequence(rle(as.character(elewgs$V5))$lengths)

#alright, get just the non-NA

elewgs_snp_windows <- elewgs[elewgs$V5 == FALSE,]

#get places where the counter is > 8

elewgs_snp_windows_nine <- elewgs_snp_windows[elewgs_snp_windows$V6 > 8,]

nrow(elewgs_snp_windows[elewgs_snp_windows$V6 == 9,])
#[1] 2394
#way more windows

elewgs_snp_windows_nine$window <- seq(1,nrow(elewgs_snp_windows_nine),1)

upper<- as.integer(rownames(elewgs_snp_windows_nine)) - 8
lower<- as.integer(rownames(elewgs_snp_windows_nine))


elewgs_consec_hyper <- NULL

#okay, let's keep going, but this time, for loop


elewgs_consec_hyper <- NULL
some_rows <- NULL
number_rows_greater_15_snps <- NULL
newdf <- NULL

for (row in 1:nrow(elewgs_snp_windows_nine)) {
 some_rows <- as.character(upper[row]:lower[row])
 the_elewgs_some_rows <- elewgs[some_rows,]
 number_rows_greater_15_snps <- nrow(the_elewgs_some_rows[the_elewgs_some_rows$V4 > 15,])
 newdf <- data.frame(upper_row_name = upper[row],lower_row_name=lower[row],number_rows_greater_15_snps=number_rows_greater_15_snps)
 elewgs_consec_hyper <- rbind(elewgs_consec_hyper,newdf)
}

#yep, way more windows.... this takes a while to run


elewgs_consec_hyper
#yep, more rows

#and to be sure,

elewgs_consec_hyper[elewgs_consec_hyper$number_rows_greater_15_snps>8,]
#okay, lots of windows, some of which are surely overlapping

nrow(elewgs_consec_hyper[elewgs_consec_hyper$number_rows_greater_15_snps>8,])
#[1] 880
#880 rows representing potentially overlapping sets of 9 or more 1kb windows that harbor at least 16 SNP's consecutively

#okay, do the consecutive window thing here to get the non-overlapping windows

#subset the regions that have more than 8 consecutive 1kb windows with greater than 15 SNP's

elewgs_consec_hyper_grt_eight <- elewgs_consec_hyper[elewgs_consec_hyper$number_rows_greater_15_snps>8,]
#make a new, duplicate df with a fake row of data to extend the df by one
elewgs_consec_hyper_grt_eight_b <- rbind(elewgs_consec_hyper_grt_eight, (c(1:3)))
#basically, get the position difference between rows by generating a nearly duplicate column that is offset by one row so you are finding the position difference between rows. Here, if the difference is 1, then you know those windows are consecutive and should be condensed into the same hyperdivergent region.
elewgs_consec_hyper_grt_eight_b$diff <- elewgs_consec_hyper_grt_eight_b$upper_row_name - c(0,elewgs_consec_hyper_grt_eight$upper_row_name)

#if the rows are consecutive, keep counting, if not, re-start the count
elewgs_consec_hyper_grt_eight_b$seq <- sequence(rle(as.character(elewgs_consec_hyper_grt_eight_b$diff))$lengths)
#great, how many ones are there (minus the fake row I added to do the subtraction)? this is the number

nrow(elewgs_consec_hyper_grt_eight_b[elewgs_consec_hyper_grt_eight_b$seq == 1,])-1
#[1] 248
#248 hyperdivergent regions!

#ok, making figures to illustrate points.

figure_data <- data.frame(data=c("C. inopinata RAD", "C. elegans pseduo-RAD","C. elegans WGS"),number_of_consecutive_1kb_windows_grt_8=c(11,18,2394),number_of_hyperdivergent_regions=c(0,0,248))
library(ggplot2)
library(cowplot)
library(lemon)
library(ggrepel)

ggplot(figure_data,aes(x=number_of_consecutive_1kb_windows_grt_8,y=number_of_hyperdivergent_regions,label=data)) + geom_point() + geom_text_repel() + theme_cowplot() + xlab("Number of consecutive 1 KB windows greater than 8\n(with segregating sites)") + ylab("Number of hyperdivergent regions") + theme(plot.margin = margin(, 0.5, , , "cm"))


ggsave("hyperdivergent_regions_scatterplot.png",height=8,width=10,units="in",bg="white")

#one more figure

inop$data <- "C. inopinata RAD"

elerad$data <- "C. elegans pseduo-RAD"

elewgs$data <- "C. elegans WGS"

all_seg_sites <- rbind(inop,elerad,elewgs)

all_seg_sites_na_omit <- na.omit(all_seg_sites)

library(ggforce)

ggplot(all_seg_sites_na_omit, aes(x = data, y = V4)) + geom_sina(size=1, alpha=1,scale="width") + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.5, colour="red",) + theme_cowplot() + ylab("Segregating sites per 1kb window")


ggsave("segregating_sites_sina_plot.png",height=8,width=10,units="in",bg="white")


#from my first attempt at this (10/3), just going through each window one at a time. for elegans, I wrote a for loop that I re-used for inopinata.




#there are so few of these, I can just go through them by hand for now.
#the first window appears to be a stretch of 12 1kb windows

some_rows <- as.character(upper[1]:lower[4])

inop[some_rows,]

#Window 1

#            V1      V2      V3 V4    V5 V6
#3669 Sp34_Chr1 3668000 3669000 NA  TRUE  4
#3670 Sp34_Chr1 3669000 3670000  1 FALSE  1
#3671 Sp34_Chr1 3670000 3671000  2 FALSE  2
#3672 Sp34_Chr1 3671000 3672000  4 FALSE  3
#3673 Sp34_Chr1 3672000 3673000 25 FALSE  4
#3674 Sp34_Chr1 3673000 3674000 23 FALSE  5
#3675 Sp34_Chr1 3674000 3675000  1 FALSE  6
#3676 Sp34_Chr1 3675000 3676000  4 FALSE  7
#3677 Sp34_Chr1 3676000 3677000  4 FALSE  8
#3678 Sp34_Chr1 3677000 3678000 21 FALSE  9
#3679 Sp34_Chr1 3678000 3679000  8 FALSE 10
#3680 Sp34_Chr1 3679000 3680000  3 FALSE 11
#3681 Sp34_Chr1 3680000 3681000 30 FALSE 12
#3682 Sp34_Chr1 3681000 3682000 NA  TRUE  1

#a bunch of these windows have <16 SNP's!

#ok, next window

some_rows <- as.character(upper[5]:lower[5])

inop[some_rows,]

#Window 2

#            V1      V2      V3 V4    V5 V6
#5337 Sp34_Chr1 5336000 5337000 NA  TRUE  1
#5338 Sp34_Chr1 5337000 5338000  1 FALSE  1
#5339 Sp34_Chr1 5338000 5339000  3 FALSE  2
#5340 Sp34_Chr1 5339000 5340000  2 FALSE  3
#5341 Sp34_Chr1 5340000 5341000 22 FALSE  4
#5342 Sp34_Chr1 5341000 5342000 12 FALSE  5
#5343 Sp34_Chr1 5342000 5343000 18 FALSE  6
#5344 Sp34_Chr1 5343000 5344000  4 FALSE  7
#5345 Sp34_Chr1 5344000 5345000 13 FALSE  8
#5346 Sp34_Chr1 5345000 5346000  8 FALSE  9
#5347 Sp34_Chr1 5346000 5347000 NA  TRUE  1

#again, a bunch of sites with number of SNPs < 16

#keep going


some_rows <- as.character(upper[6]:lower[6])

inop[some_rows,]

#> inop[some_rows,]
#             V1       V2       V3 V4    V5 V6
#10612 Sp34_Chr1 10611000 10612000 NA  TRUE  2
#10613 Sp34_Chr1 10612000 10613000 16 FALSE  1
#10614 Sp34_Chr1 10613000 10614000  9 FALSE  2
#10615 Sp34_Chr1 10614000 10615000 10 FALSE  3
#10616 Sp34_Chr1 10615000 10616000  7 FALSE  4
#10617 Sp34_Chr1 10616000 10617000 12 FALSE  5
#10618 Sp34_Chr1 10617000 10618000  2 FALSE  6
#10619 Sp34_Chr1 10618000 10619000  1 FALSE  7
#10620 Sp34_Chr1 10619000 10620000  6 FALSE  8
#10621 Sp34_Chr1 10620000 10621000  4 FALSE  9
#10622 Sp34_Chr1 10621000 10622000 NA  TRUE  1

#again, a bunch of sites with number of SNPs < 16


some_rows <- as.character(upper[7]:lower[7])

inop[some_rows,]

#             V1       V2       V3 V4    V5 V6
#37320 Sp34_Chr2 16724000 16725000 NA  TRUE  2
#37321 Sp34_Chr2 16725000 16726000  5 FALSE  1
#37322 Sp34_Chr2 16726000 16727000  5 FALSE  2
#37323 Sp34_Chr2 16727000 16728000 14 FALSE  3
#37324 Sp34_Chr2 16728000 16729000 13 FALSE  4
#37325 Sp34_Chr2 16729000 16730000  7 FALSE  5
#37326 Sp34_Chr2 16730000 16731000  5 FALSE  6
#37327 Sp34_Chr2 16731000 16732000 21 FALSE  7
#37328 Sp34_Chr2 16732000 16733000 22 FALSE  8
#37329 Sp34_Chr2 16733000 16734000  3 FALSE  9
#37330 Sp34_Chr2 16734000 16735000 NA  TRUE  1


#again, a bunch of sites with number of SNPs < 16


some_rows <- as.character(upper[8]:lower[8])

inop[some_rows,]

#             V1      V2      V3 V4    V5 V6
#42065 Sp34_Chr3 1352000 1353000 NA  TRUE  2
#42066 Sp34_Chr3 1353000 1354000  7 FALSE  1
#42067 Sp34_Chr3 1354000 1355000  4 FALSE  2
#42068 Sp34_Chr3 1355000 1356000 17 FALSE  3
#42069 Sp34_Chr3 1356000 1357000 31 FALSE  4
#42070 Sp34_Chr3 1357000 1358000 12 FALSE  5
#42071 Sp34_Chr3 1358000 1359000 26 FALSE  6
#42072 Sp34_Chr3 1359000 1360000 26 FALSE  7
#42073 Sp34_Chr3 1360000 1361000 14 FALSE  8
#42074 Sp34_Chr3 1361000 1362000 10 FALSE  9
#42075 Sp34_Chr3 1362000 1363000 NA  TRUE  1


some_rows <- as.character(upper[9]:lower[9])

inop[some_rows,]

#             V1      V2      V3 V4    V5 V6
#45809 Sp34_Chr3 5096000 5097000 NA  TRUE 17
#45810 Sp34_Chr3 5097000 5098000 12 FALSE  1
#45811 Sp34_Chr3 5098000 5099000  3 FALSE  2
#45812 Sp34_Chr3 5099000 5100000  7 FALSE  3
#45813 Sp34_Chr3 5100000 5101000 16 FALSE  4
#45814 Sp34_Chr3 5101000 5102000  9 FALSE  5
#45815 Sp34_Chr3 5102000 5103000 12 FALSE  6
#45816 Sp34_Chr3 5103000 5104000 12 FALSE  7
#45817 Sp34_Chr3 5104000 5105000 14 FALSE  8
#45818 Sp34_Chr3 5105000 5106000 44 FALSE  9
#45819 Sp34_Chr3 5106000 5107000 NA  TRUE  1

#nope!


some_rows <- as.character(upper[10]:lower[11])

inop[some_rows,]
#             V1       V2       V3 V4    V5 V6
#58816 Sp34_Chr3 18103000 18104000 NA  TRUE 11
#58817 Sp34_Chr3 18104000 18105000 33 FALSE  1
#58818 Sp34_Chr3 18105000 18106000 46 FALSE  2
#58819 Sp34_Chr3 18106000 18107000 14 FALSE  3
#58820 Sp34_Chr3 18107000 18108000 10 FALSE  4
#58821 Sp34_Chr3 18108000 18109000  2 FALSE  5
#58822 Sp34_Chr3 18109000 18110000 33 FALSE  6
#58823 Sp34_Chr3 18110000 18111000 20 FALSE  7
#58824 Sp34_Chr3 18111000 18112000 17 FALSE  8
#58825 Sp34_Chr3 18112000 18113000 56 FALSE  9
#58826 Sp34_Chr3 18113000 18114000  4 FALSE 10
#58827 Sp34_Chr3 18114000 18115000 NA  TRUE  1

#nope, but closer!!


some_rows <- as.character(upper[12]:lower[12])

inop[some_rows,]
#             V1      V2      V3 V4    V5 V6
#64828 Sp34_Chr4 4678000 4679000 NA  TRUE  3
#64829 Sp34_Chr4 4679000 4680000 17 FALSE  1
#64830 Sp34_Chr4 4680000 4681000 20 FALSE  2
#64831 Sp34_Chr4 4681000 4682000  5 FALSE  3
#64832 Sp34_Chr4 4682000 4683000 13 FALSE  4
#64833 Sp34_Chr4 4683000 4684000 34 FALSE  5
#64834 Sp34_Chr4 4684000 4685000 19 FALSE  6
#64835 Sp34_Chr4 4685000 4686000 42 FALSE  7
#64836 Sp34_Chr4 4686000 4687000 26 FALSE  8
#64837 Sp34_Chr4 4687000 4688000 16 FALSE  9
#64838 Sp34_Chr4 4688000 4689000 NA  TRUE  1

#nope



some_rows <- as.character(upper[13]:lower[14])

inop[some_rows,]
#             V1       V2       V3 V4    V5 V6
#72742 Sp34_Chr4 12592000 12593000 NA  TRUE  1
#72743 Sp34_Chr4 12593000 12594000  8 FALSE  1
#72744 Sp34_Chr4 12594000 12595000  5 FALSE  2
#72745 Sp34_Chr4 12595000 12596000 22 FALSE  3
#72746 Sp34_Chr4 12596000 12597000  2 FALSE  4
#72747 Sp34_Chr4 12597000 12598000  4 FALSE  5
#72748 Sp34_Chr4 12598000 12599000 23 FALSE  6
#72749 Sp34_Chr4 12599000 12600000  5 FALSE  7
#72750 Sp34_Chr4 12600000 12601000  6 FALSE  8
#72751 Sp34_Chr4 12601000 12602000  9 FALSE  9
#72752 Sp34_Chr4 12602000 12603000  1 FALSE 10
#72753 Sp34_Chr4 12603000 12604000 NA  TRUE  1

#nope!

some_rows <- as.character(upper[15]:lower[17])

inop[some_rows,]

#             V1       V2       V3 V4    V5 V6
#91731 Sp34_Chr5 10563000 10564000 NA  TRUE  2
#91732 Sp34_Chr5 10564000 10565000  1 FALSE  1
#91733 Sp34_Chr5 10565000 10566000 10 FALSE  2
#91734 Sp34_Chr5 10566000 10567000  3 FALSE  3
#91735 Sp34_Chr5 10567000 10568000  2 FALSE  4
#91736 Sp34_Chr5 10568000 10569000 11 FALSE  5
#91737 Sp34_Chr5 10569000 10570000  3 FALSE  6
#91738 Sp34_Chr5 10570000 10571000  1 FALSE  7
#91739 Sp34_Chr5 10571000 10572000  2 FALSE  8
#91740 Sp34_Chr5 10572000 10573000  6 FALSE  9
#91741 Sp34_Chr5 10573000 10574000  3 FALSE 10
#91742 Sp34_Chr5 10574000 10575000  5 FALSE 11
#91743 Sp34_Chr5 10575000 10576000 NA  TRUE  1

#nope!

some_rows <- as.character(upper[18]:lower[19])

inop[some_rows,]

#             V1       V2       V3 V4    V5 V6
#94167 Sp34_Chr5 12999000 13000000 NA  TRUE  8
#94168 Sp34_Chr5 13000000 13001000  4 FALSE  1
#94169 Sp34_Chr5 13001000 13002000  6 FALSE  2
#94170 Sp34_Chr5 13002000 13003000  5 FALSE  3
#94171 Sp34_Chr5 13003000 13004000  7 FALSE  4
#94172 Sp34_Chr5 13004000 13005000  7 FALSE  5
#94173 Sp34_Chr5 13005000 13006000 13 FALSE  6
#94174 Sp34_Chr5 13006000 13007000  2 FALSE  7
#94175 Sp34_Chr5 13007000 13008000  5 FALSE  8
#94176 Sp34_Chr5 13008000 13009000  1 FALSE  9
#94177 Sp34_Chr5 13009000 13010000 13 FALSE 10
#94178 Sp34_Chr5 13010000 13011000 NA  TRUE  1

#nope


