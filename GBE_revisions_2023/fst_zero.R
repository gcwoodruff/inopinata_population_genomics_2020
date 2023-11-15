#Is mean Fst significantly different from zero?

setwd("/Users/gavin/genome/pop_gen_revisions_9-2023/")

dat <- read.csv("all_fst_combined.csv", header=T,stringsAsFactors=T)

ir_yo <- dat[dat$pop_pair == "irio-yona",]
ish_ir <- dat[dat$pop_pair == "ishi-irio",]
ish_yo <- dat[dat$pop_pair == "ishi-yona",]

#linear model approach
summary(lm(ir_yo$FST ~ 1))

#Call:
#lm(formula = ir_yo$FST ~ 1)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.45168 -0.03228 -0.00628  0.02299  0.73652
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0130801  0.0006951   18.82   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.05868 on 7125 degrees of freedom

summary(lm(ish_ir$FST ~ 1))

#Call:
#lm(formula = ish_ir$FST ~ 1)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.38518 -0.02201 -0.00578  0.01662  0.61322
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0157836  0.0004661   33.87   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.03938 on 7139 degrees of freedom

summary(lm(ish_yo$FST ~ 1))

#Call:
#lm(formula = ish_yo$FST ~ 1)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.40098 -0.02525 -0.00588  0.01802  0.53712
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept) 0.0191759  0.0005393   35.56   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.04556 on 7137 degrees of freedom


#wilcox test

wilcox.test(ir_yo$FST)
#	Wilcoxon signed rank test with continuity correction
#
#data:  ir_yo$FST
#V = 15469938, p-value < 2.2e-16
#alternative hypothesis: true location is not equal to 0

wilcox.test(ish_ir$FST)

#	Wilcoxon signed rank test with continuity correction
#
#data:  ish_ir$FST
#V = 18554136, p-value < 2.2e-16
#alternative hypothesis: true location is not equal to 0

wilcox.test(ish_yo$FST)

#	Wilcoxon signed rank test with continuity correction
#
#data:  ish_yo$FST
#V = 18958978, p-value < 2.2e-16
#alternative hypothesis: true location is not equal to 0

