###############################################################################################################
########################### Example NHP univaraite risk factor analysis ############################
###############################################################################################################

rm(list=ls())


library(tidyverse)
library(ggplot2)
library(stringr)
library(summarytools)
library(writexl)
library(ggpubr)
library(lsmeans)
library(lme4)

 
x <- read.csv("inputfile.csv")

## Null model
m0 <- glmer(cbind(monkey, 365 - monkey) ~ + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m0)

## Age
x <- x %>%
  mutate(ageCat = cut(age3, breaks=c(-Inf, 5, 15, 30, 55, Inf), labels = c("0", "1", "2","3", "4")))
t <- data.frame(table(x$ageCat, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
x$ageCat <- factor(x$ageCat)
m2 <- glmer(monkey ~ ageCat + (1 | houseID), data = x, family = gaussian)
summary(m2)
se <- sqrt(diag(vcov(m2)))
cof <- cbind(Est = fixef(m2), LL = fixef(m2) - 1.96 * se, UL = fixef(m2) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m2)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$Var <- revalue(t$Var, c("0"="Under 5", "1"="5-15", "2"="15-30", "3"="30-55", "4"="Over 55"))
t$q <- "Age category"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- t

## Gender, significantly impacts model ***
t <- data.frame(table(x$male, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ gender + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "gender"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Ethnicity ***
x$eth <- revalue(x$Bangsa, c("Cina"="other", "Filipino"="other", "Indonesian"="other", "Kadazan"="Dusun", "Melayu"="other", 
                               "Kegayan"="other", "Tambunuo"="other", "Ubian"="other", "Lain-lain"="other"))
x$eth[x$eth==""] <- "other"
x$eth <- factor(x$eth)
t <- data.frame(table(x$eth, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ eth + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "ethnicity"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Treatment seeking behavior / illness

## Previous malaria ***
x$hadMalaria[x$hadMalaria==""] <- "N"
x$hadMalaria <- factor(x$hadMalaria)
t <- data.frame(table(x$hadMalaria, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ hadMalaria + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "previously diagnosed malaria"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Take medication for malaria ***
x$malariaMedic[x$malariaMedic==""] <- "N"
x$malariaMedic <- factor(x$malariaMedic)
t <- data.frame(table(x$malariaMedic, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ malariaMedic + (1 | houseID), data = x, family = gaussian)
summary(m1)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "anti-malaria medication"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Medicine from clinic *
x$med <- factor(str_detect(x$feverMonth, "umonkey dari klinik"))
t <- data.frame(table(x$med, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ med + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "treatment seeking - clinic meds"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)


## Traditional medicine *
x$trad <- factor(str_detect(x$feverMonth, "umonkey tradisional"))
t <- data.frame(table(x$trad, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ trad + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "treatment seeking - traditional med"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Meet doctor in hospital **
x$doc <- factor(str_detect(x$feverMonth, "doktor hospital"))
t <- data.frame(table(x$doc, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ doc + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "treatment seeking - hospital"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Do nothing when sick ***
x$notx <- factor(str_detect(x$feverMonth, "Tiada"))
t <- data.frame(table(x$notx, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ notx + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "no treatment seeking"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Meet traditional healer
x$heal <- factor(str_detect(x$feverMonth, "bantuan pawang"))
t <- data.frame(table(x$heal, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ heal + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "treatment seeking - healer"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Occupation fishing, office/shop, student ***, rubber * --> hunting not included 
summary(x$occupation)
x$occ <- revalue(x$occupation, c("Suri rumah"="none", "Tiada"="none", "Pesara "="none", "Pekerja pejamonkey swasta / kerajaan"="office/shop", 
                                   "Pekedai"="office/shop", "Lain-lain"="other", "Kontraktor"="other", "Pemandu"="other", "Nelayan"="fishing",
                                   "Pekerja ladang getah"="rubber", "Pekerja ladang kelapa sawit"="palmoil", "Petani"="farmer", "Pelajar"="student"))
x$occ[x$occ==""] <- "none"
x$occ[x$occ=="Pesara"] <- "none"
x$occ <- factor(x$occ)
t <- data.frame(table(x$occ, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ occ + (1 | houseID), data = x, family = gaussian)
summary(m1)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "occupation"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Individual farmwork ***
x$farmWork[x$farmWork ==""] <- "N"
x$farmWork <- factor(x$farmWork)
t <- data.frame(table(x$farmWork, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ farmWork + (1 | houseID), data = x, family = gaussian)
summary(m1)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "farm work"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Occupation place ***
summary(x$occPlace)
x$occPlace <- revalue(x$occPlace, c("Dalam atau dekat kampung"= "kampung", "Dalam atau dekat kampung lain dalam daerah yang sama"=
                                        "district", "Dalam atau sekeliling rumah"="house", "Di luar daerah, dalam negeri Sabah"=
                                        "outside district", "Di luar Sabah"="outside district", "Tidak kerja atau sekolah"="house"))
x$occPlace[x$occPlace==""] <- "house"
x$occPlace <- factor(x$occPlace)
t <- data.frame(table(x$occPlace, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ occPlace + (1 | houseID), data = x, family = gaussian)
summary(m1)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "occupation place"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Time travelling to/from work in early morning (11pm - 6am) ***
x$goAM <- factor(str_detect(x$goTime, "awal pagi"))
x$backAM <- factor(str_detect(x$backTime, "awal pagi"))
x$AM <- factor(ifelse(x$goAM == "TRUE" | x$backAM == "TRUE", "Y", "N"))
t <- data.frame(table(x$AM, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ AM + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "travel to/from work 11pm-6am"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)
summary(m1)

## Time travelling to/from work in evening (5pm - 10pm) ***
x$goPM <- factor(str_detect(x$goTime, "10pm"))
x$backPM <- factor(str_detect(x$backTime, "10pm"))
x$PM <- factor(ifelse(x$goPM == "TRUE" | x$backPM == "TRUE", "Y", "N"))
t <- data.frame(table(x$PM, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ PM + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "travel to/from work 5-10pm"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)
summary(m1)
###########################################################################################################################################################################
## Walk to work
x$walk <- factor(str_detect(x$goWork, "kaki"))
summary(x$walk)
t <- data.frame(table(x$walk, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ walk + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "walk to work"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Walk to work through forest 
summary(x$goWork)
x$forestPath <- factor(str_detect(x$goWork, "hutan"))
t <- data.frame(table(x$forestPath, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ forestPath + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "walk to work through forest"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Frequency of farm work (none, less than 15 days, over 15 days per month)
x$goToFarm <- findInterval(x$goToFarm, c(0, 15))
x$goToFarm[is.na(x$goToFarm)] <- 0
x$goToFarm <- factor(x$goToFarm)
t <- data.frame(table(x$goToFarm, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ goToFarm + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$Var <- revalue(t$Var, c("1"="under 15 days", "2"="over 15 days"))
t$q <- "frequency of farm work per month"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Go to forest
x$goForest[x$goForest == ""] <- "N"
x$goForest <- factor(x$goForest)
t <- data.frame(table(x$goForest, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ goForest + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "go to forest"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Time travelling to/from forest in early morning (11pm - 6am)
x$forGoAM <- factor(str_detect(x$goToForest, "awal pagi"))
x$forBackAM <- factor(str_detect(x$backForest, "awal pagi"))
x$forAM <- factor(ifelse(x$forGoAM == "TRUE" | x$forBackAM == "TRUE", "Y", "N"))
x$forAM <- factor(x$forAM)
t <- data.frame(table(x$forAM, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ forAM + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "go to forest 11pm - 6am"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

# Time travelling to/from forest in evening (5pm - 10pm)
x$forGoPM <- factor(str_detect(x$goToForest, "10pm"))
x$forBackPM <- factor(str_detect(x$backForest, "10pm"))
x$forPM <- factor(ifelse(x$forGoPM == "TRUE" | x$forBackPM == "TRUE", "Y", "N"))
x$forPM <- factor(x$forPM)
t <- data.frame(table(x$forPM, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ forPM + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "go to forest 5-10pm"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Any forest travel at night
x$forNight <- factor(ifelse(x$forAM == "Y" | x$forPM == "Y", "Y", "N"))
summary(x$forNight)
x$forNight <- factor(x$forNight)
t <- data.frame(table(x$forNight, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ forNight + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "go to forest 5pm - 6am"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Activities in forest - hunting
x$hunt <- factor(str_detect(x$actForest, "Memburu"))
x$hunt <- factor(x$hunt)
t <- data.frame(table(x$hunt, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ hunt + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "hunting"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Activities in forest - collecting wood
x$wood <- factor(str_detect(x$actForest, "kayu"))
x$wood <- factor(x$wood)
t <- data.frame(table(x$wood, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ wood + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "collect wood in forest"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Clearing land
x$clear <- x$DalamBeberapaTahunIniAdakahKamuMelakukanAktivitiihkanKawa
x$clear[x$clear==""] <- "N"
x$clear <- factor(x$clear)
t <- data.frame(table(x$clear, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ clear + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "clear land"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Construction activities
x$const <- x$PernahkahKamuTerlimonkeyDalamAktivitiPembinaanSepertiRumahJalanRaya
x$const[x$const==""] <- "N"
x$const <- factor(x$const)
t <- data.frame(table(x$const, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ const + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "construction"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Other activities
x$activity <- revalue(x$activity, c("Bersukan"="sport", "Lain-lain"="other", "Lepak-lepak di dalam rumah"
                                      ="none", "Lepak-lepak di luar rumah"="outside house", "Memancing ikan"=
                                        "fishing", "Memburu"="other", "Tiada"="none"))
x$activity[x$activity==""]<- "none"
x$activity <- factor(x$activity)
t <- data.frame(table(x$activity, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ activity + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "activity"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Time of activities - early morning
x$actAM <- factor(str_detect(x$whenAct, "awal pagi"))
x$actAM <- factor(x$actAM)
t <- data.frame(table(x$actAM, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ actAM + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "early morning activities"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Time of activities - early morning
x$actPM <- factor(str_detect(x$whenAct, "10pm"))
x$actPM <- factor(x$actPM)
t <- data.frame(table(x$actPM, x$monkey))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ actPM + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "evening activities"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Bath outside
x$bathOut <- factor(str_detect(x$bathPlace, "di luar"))
x$bathOut <- factor(x$bathOut)
t <- data.frame(table(x$bathOut, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey ~ bathOut + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "bathe outside"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Bath at the river
x$bathRiv <- factor(str_detect(x$bathPlace, "Sungai"))
x$bathRiv <- factor(x$bathRiv)
t <- data.frame(table(x$bathRiv, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ bathRiv + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "bathe at river"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Bathe inside
x$bathIn <- factor(ifelse(x$bathOut == "FALSE" & x$bathRiv == "FALSE", "Y", "N"))
x$bathIn <- factor(x$bathIn)
t <- data.frame(table(x$bathIn, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ bathIn + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "bathe inside"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Time of bathing
x$bathT <- x$BilaMasaBiasanyaKamuMandi
x$bathAM <- factor(str_detect(x$bathT, "awal pagi"))
x$bathPM <- factor(str_detect(x$bathT, "10pm"))
x$bathAMout <- factor(ifelse(x$bathIn == "N" & x$bathAM == "TRUE", "Y", "N"))
x$bathPMout <- factor(ifelse(x$bathIn == "N" & x$bathPM == "TRUE", "Y", "N"))
x$bathNite <- factor(ifelse(x$bathAMout=="Y" | x$bathPMout =="Y", "Y", "N"))
t <- data.frame(table(x$bathNite, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ bathNite + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "bathe outside at night"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Time outside house in evening
x$TimeOut <- x$BerapaLamaSelalunyaAndaBeradaDiLuarRumahPadaWaktuMalam6PetangHin
x$TimeOut[x$TimeOut==""]<- "Tidak tahu"
x$TimeOut <- factor(x$TimeOut)
t <- data.frame(table(x$TimeOut, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ TimeOut + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "time outside house at night"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Have stayed outside kampung in past month
summary(x$outsideKmpg)
x$outsideKmpg[x$outsideKmpg==""]<- "N"
x$outsideKmpg <- factor(x$outsideKmpg)
t <- data.frame(table(x$outsideKmpg, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ outsideKmpg + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "stay outside kampung in past month"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Location stayed outside kampung
x$stayOut <- revalue(x$stayOut, c("Dalam atau dekat kampung lain dalam daerah sama"="district", 
                                    "Di luar daerah dalam Sabah"= "outside", "Di luar Sabah"="outside",
                                    "Ladang/ estet dekat kampung ini"="plantation", "Hutan dekat kampung ini"
                                    ="forest"))
x$stayOut <- factor(x$stayOut)
t <- data.frame(table(x$stayOut, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ stayOut + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "where stay outside kampung"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Sleep outside walls
x$outsideWall[x$outsideWall==""]<- "N"
x$outsideWall <- factor(x$outsideWall)
t <- data.frame(table(x$outsideWall, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ outsideWall + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "slept outside walls"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Mosquito prevention - bednet 
x$bednet <- factor(str_detect(x$preventMosq, "Kelambu"))
t <- data.frame(table(x$bednet, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ bednet + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "use bednet"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Mosquito prevention - insecticide
x$insect <- factor(str_detect(x$preventMosq, "Ubat"))
t <- data.frame(table(x$insect, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ insect + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "use insecticide"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Mosquito prevention - fan
x$fan <- factor(str_detect(x$preventMosq, "Kipas"))
t <- data.frame(table(x$fan, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ fan + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "use fan to prevent mosquitoes"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Mosquito prevention - smoke
x$smoke <- factor(str_detect(x$preventMosq, "Asap"))
t <- data.frame(table(x$smoke, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ smoke + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "use smoke to prevent mosquitoes"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Mosquito prevention - screens
x$screen <- factor(str_detect(x$preventMosq, "Jaring"))
t <- data.frame(table(x$screen, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ screen + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "use screen to prevent mosquitoes"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Don't use any mosquito prevention practices
x$noMos <- factor(str_detect(x$preventMosq, "Tiada"))
t <- data.frame(table(x$noMos, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ noMos + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "don't use any mosquito prevention"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Contact with monkeys
x$sawMonkey[x$sawMonkey==""]<- "N"
x$sawMonkey <- factor(x$sawMonkey)
t <- data.frame(table(x$sawMonkey, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ sawMonkey + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "contact with monkeys"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Type of monkey - long tailed macaque
x$ltmac <- factor(str_detect(x$monkeyType, "panjang"))
t <- data.frame(table(x$ltmac, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ ltmac + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "long tailed macaques"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Type of monkey- pig tailed macaque
x$stmac <- factor(str_detect(x$monkeyType, "pendek"))
x$stmac2 <- factor(str_detect(x$monkeyType, "Beruk"))
x$stmac <- factor(ifelse(x$stmac =="TRUE"|x$stmac2 =="TRUE", "Y", "N"))
t <- data.frame(table(x$stmac, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ stmac + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "pig tailed macaques"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Places monkeys seen - around the house
x$mhous <- factor(str_detect(x$placeSawMonkey, "rumah"))
t <- data.frame(table(x$mhous, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ mhous + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "monkeys seen around the house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Places monkeys seen - around the kampung
x$mkam <- factor(str_detect(x$placeSawMonkey, "kampung"))
t <- data.frame(table(x$mkam, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ mkam + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "monkeys seen around kampung"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Places monkeys seen - around farms/ plantations
x$mgar <- factor(str_detect(x$placeSawMonkey, "kebun"))
x$mplant <- factor(str_detect(x$placeSawMonkey, "ladang"))
x$mfarm <- factor(ifelse(x$mgar =="TRUE"|x$mplant =="TRUE", "Y", "N"))
t <- data.frame(table(x$mfarm, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ mfarm + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "monkeys seen around farm or plantation"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Frequency of monkey sightings
x$totalSawMonkey <- revalue(x$totalSawMonkey, c("Beberapa kali dalam satu bulan"="monthly",
                                                  "Beberapa kali dalam satu tahun"="yearly", "Lain-lain (Contoh ketika berkebun)"
                                                  = "yearly", "Minggu-minggu"="weekly", "Tiap-tiap hari"="daily"))
x$totalSawMonkey <- factor(x$totalSawMonkey)
t <- data.frame(table(x$totalSawMonkey, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ totalSawMonkey + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "frequency of monkey sightings"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

### Write individual level results
#write.csv(df.sum, "ind_univariate_p50.csv")
#write.csv(df.sum, "ind_univariate_p75.csv")

################################ Household level factors ##################################

## SES
x$wealth <- factor(x$wealth)
t <- data.frame(table(x$wealth, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ wealth + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "socioeconomic status"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Length of time at house
x$stayPeriod <- factor(x$stayPeriod)
t <- data.frame(table(x$stayPeriod, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ stayPeriod + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "length of stay at house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## House age
x$houseAge <- factor(x$houseAge)
t <- data.frame(table(x$houseAge, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ houseAge + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house age"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household head education
x$hhed <- factor(x$hhed)
t <- data.frame(table(x$hhed, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ hhed + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household head education"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Roof type - zink
x$roof <- factor(x$roof)
t <- data.frame(table(x$roof, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ roof + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "corrugated iron roof"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Floor
x$floor <- factor(x$floor)
t <- data.frame(table(x$floor, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ floor + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "concrete or tile floor"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Wall type - wood or bamboo walls
x$wall <- factor(x$wall)
t <- data.frame(table(x$wall, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ wall + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "wood or bamboo walls"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## House height
x$hgt <- factor(x$hgt)
t <- data.frame(table(x$hgt, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ hgt + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house height"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Gaps in house
x$bdgap <- factor(x$bdgap)
t <- data.frame(table(x$bdgap, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ bdgap + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house gaps"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Windows that can close
x$wind <- factor(x$wind)
t <- data.frame(table(x$wind, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ wind + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "windows in house can close"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Insect screens in house windows
x$insc <- factor(x$insc)
t <- data.frame(table(x$insc, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ insc + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "insect screens in house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Kitchen outside house
x$kitch <- factor(x$kitch)
t <- data.frame(table(x$kitch, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ kitch + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "kitchen outside house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Toilet
x$toilet <- factor(x$toilet)
t <- data.frame(table(x$toilet, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ toilet + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house has toilet"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Toilet inside house
x$toiletIn <- factor(x$toiletIn)
t <- data.frame(table(x$toiletIn, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ toiletIn + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "toilet is inside house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Piped water inside house
x$pipe <- factor(x$pipe)
t <- data.frame(table(x$pipe, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ pipe + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "piped water inside house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Cattle
x$cow <- factor(x$cow)
t <- data.frame(table(x$cow, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ cow + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "cattle"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Buffalo
x$buff <- factor(x$buff)
t <- data.frame(table(x$buff, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ buff + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "buffalo"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Goats
x$goat <- factor(x$goat)
t <- data.frame(table(x$goat, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ goat + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "goat"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Pigs
x$pig <- factor(x$pig)
t <- data.frame(table(x$pig, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ pig + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "pigs"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Pet monkey
x$petm <- factor(x$petm)
t <- data.frame(table(x$petm, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ petm + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house has pet monkey"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Pet monkey in kampung
x$petkg <- factor(x$petkg)
t <- data.frame(table(x$petkg, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ petkg + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "pet monkey in kampung"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Monkey raid household crops
x$cropm <- factor(x$cropm)
t <- data.frame(table(x$cropm, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ cropm + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "monkeys raid household crops"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## River near house
x$river <- factor(x$river)
t <- data.frame(table(x$river, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ river + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "river near house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Pond near house
x$pond <- factor(x$pond)
t <- data.frame(table(x$pond, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ pond + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "pond near house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Well near house
x$well <- factor(x$well)
t <- data.frame(table(x$well, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ well + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "well near house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Plastic containers of water near house
x$plast <- factor(x$plast)
t <- data.frame(table(x$plast, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ plast + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "plastic containers of water near house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Lake near house
x$lake <- factor(x$lake)
t <- data.frame(table(x$lake, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ lake + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "lake near house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Sea
x$sea <- factor(x$sea)
t <- data.frame(table(x$sea, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ sea + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "Sea near house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Uncovered water containers seen
x$cont <- factor(x$cont)
t <- data.frame(table(x$cont, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ cont + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "uncovered water containers seen"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household farms fruit
x$fruit <- factor(x$fruit)
t <- data.frame(table(x$fruit, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ fruit + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "farming - fruit"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household farming rubber
x$rub <- factor(x$rub)
t <- data.frame(table(x$rub, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ rub + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household farming rubber"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household farming corn
x$corn <- factor(x$corn)
t <- data.frame(table(x$corn, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ corn + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household farming corn"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household has livestock
x$livestk <- factor(x$livestk)
t <- data.frame(table(x$livestk, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ livestk + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household has livestock"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household farms vegetables
x$veg <- factor(x$veg)
t <- data.frame(table(x$veg, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ veg + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household farms vegetables"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household has rice paddies
x$paddy <- factor(x$paddy)
t <- data.frame(table(x$paddy, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ paddy + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household has rice paddies"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Household farms palm oil
x$po <- factor(x$po)
t <- data.frame(table(x$po, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ po + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household farms palm oil"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Amount of land farmed
x$farmNow <- factor(x$farmNow)
t <- data.frame(table(x$farmNow, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ farmNow + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "household land farmed"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Distance of farm from house
x$farmDist <- factor(x$farmDist)
t <- data.frame(table(x$farmDist, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ farmDist + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "distance of farm from house"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Use insecticide for farming
x$farmIns <- factor(x$farmIns)
t <- data.frame(table(x$farmIns, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ farmIns + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "use insecticide for farming"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Swidden farming
x$swid <- factor(x$swid)
t <- data.frame(table(x$swid, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ swid + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "swidden farming"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Collect wood from forest
x$collectWood <- factor(x$collectWood)
t <- data.frame(table(x$collectWood, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ collectWood + (1 | houseID), data = x, family = gaussian, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house collects wood from forest"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Collect food from forest
x$forestFood <- factor(x$forestFood)
t <- data.frame(table(x$forestFood, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ forestFood + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house collects food from forest"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Collect medicines from forest
x$forestMed <- factor(x$forestMed)
t <- data.frame(table(x$forestMed, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ forestMed + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "house collects medicine from forest"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Distance from clinic (by quartile)
x$clinc <- findInterval(x$clin, c(0,3.7, 6.8, 10))
x$clinc <- factor(x$clinc)
t <- data.frame(table(x$clinc, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ clinc + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "distance to nearest clinic"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Distance from hospital (by quartile)
x$hospc <- findInterval(x$hosp, c(0,9.7, 18, 27.7))
x$hospc <- factor(x$hospc)
t <- data.frame(table(x$hospc, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ hospc + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "distance to nearest hospital"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
df.sum <- rbind(df.sum, t)

## Elevation 
x$demc <- findInterval(x$dem, c(0,50, 250, 500))
x$demc <- factor(x$demc)
t <- data.frame(table(x$demc, x$bat))
tp <- subset(t, Var2==1)
tn <- subset(t, Var2==0)
t <- cbind(tp, tn)
names(t) <- c("Var", "Var2", "pos", "Var3", "Var4", "Freq")
t$total <- t$pos + t$Freq
t <- t[c("Var", "pos", "total")]
m1 <- glmer(monkey~ demc + (1 | houseID), data = x, family = gaussian)
se <- sqrt(diag(vcov(m1)))
cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
cof <- exp(cof)
a <- anova(m0, m1)
pv <- a$`Pr(>Chisq)`[2]
t <- cbind(t, cof)
t$pvalue <- pv
t$q <- "elevation"
t[1,4] <- "ref"; t[1,5] <- "ref"; t[1,6] <- "ref"
t$Var <- revalue(t$Var, c("1"="Under 50 MSL", "2"="50-250 MSL", "3"="250-500 MSL", "4"="Over 500 MSL"))
df.sum <- rbind(df.sum, t)

write.csv(df.sum, "initial_univariate.csv")
write.csv(x, "prep_for_multivariate.csv")
