###############################################################################################################
########################### Example NHP multivaraite risk factor analysis ############################
###############################################################################################################
########################### Test all factors with p values < or = 0.2     ########################### 

rm(list=ls())

install.packages('see')
library(performance)
library(ggplot2)
library(ggpubr)
library(lme4)
library(plyr)
library(stringr)
library(dplyr)
library(lmtest)

x <- read.csv("inputfile.csv")

## Null model 
m0 <- glmer(cbind(monkey, 365 - monkey) ~ + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m0)

#test random effect
m0.glm <- glm(cbind(monkey, 365 - monkey)~1, data=x, family=binomial)
AIC(logLik(m0.glm))
AIC(logLik(m0))

### sea near house
x$sea <- factor(x$sea)
m1 <- glmer(cbind(monkey, 365 - monkey) ~ sea + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m1)
lrtest(m0, m1)

## Distance from hospital (by quartile)
x$hospc <- findInterval(x$hosp, c(0,9.7, 18, 27.7))
x$hospc <- factor(x$hospc)
m2 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m2)
lrtest(m1, m2)

## Elevation
x$demc <- findInterval(x$dem, c(0,50, 250, 500))
x$demc <- factor(x$demc) 
m3 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m3)
lrtest(m2, m3)

## Occupation
x$occ <- factor(x$occ, levels = c("none", "office/shop", "other","fishing",
                                     "rubber", "palmoil", "farmer", "student"))
m4 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m4)
lrtest(m4, m3)

## House height
x$hgt <- factor(x$hgt)
x$hgt <- revalue(x$hgt, c("less1"="ground"))
m5 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m5)
lrtest(m4, m5)

## Age
m6 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ age3 +(1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m6)
lrtest(m5, m6)

## Age category
x$ageCat <- as.factor(x$ageCat)
m7 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat +(1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m7)
lrtest(m5, m6, m7)

## Ethnicity
x$eth <- factor(x$eth)
m8 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ age3 + eth + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m8)
lrtest(m7, m8)
sjPlot::tab_model(m8)

## Socioeconomic status
x$wealth <- factor(x$wealth)
m9 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m9)
lrtest(m7, m9)

## Piped water inside house
x$pipe <- factor(x$pipe)
m10 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + pipe + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m10)
lrtest(m9, m10)


## Gaps in house
x$bdgap <- factor(x$bdgap)
m11 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + bdgap + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m11)
lrtest(m9, m11)

## Toilet
x$toilet <- factor(x$toilet)
m12 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + toilet + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m12)
lrtest(m12, m9)

## Cattle
x$cow <- factor(x$cow)
m13 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + cow + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m13)
lrtest(m9, m13)

## Household farming corn
x$corn <- factor(x$corn)
m14 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + corn + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m14)
lrtest(m14, m9)

## Individual farm work
m15 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + farmWork + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m15)
lrtest(m15, m9)

## Activities in forest - collecting wood
x$wood <- factor(str_detect(x$actForest, "kayu"))
x$wood <- factor(x$wood)
m16 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m16)
lrtest(m9, m16)

## Clearing land
x$clear <- x$DalamBeberapaTahunIniAdakahKamuMelakukanAktivitiihkanKawa
x$clear[x$clear==""] <- "N"
x$clear <- factor(x$clear)
m17 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + clear + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m17)
lrtest(m17, m16)

## Construction activities
x$const <- x$PernahkahKamuTerlibatDalamAktivitiPembinaanSepertiRumahJalanRaya
x$const[x$const==""] <- "N"
x$const <- factor(x$const)
m18 <- glmer(monkey ~ sea + hospc + demc + occ + hgt+ age3 + wealth + pipe + bdgap + cow + corn + farmWork + wood + const + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m18)
lrtest(m16, m18)  

## Household head education
x$hhed <- factor(x$hhed)
m19 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + hhed + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m19)
lrtest(m16, m19)

## Other activities, outside the house
x$activity <- revalue(x$activity, c("Bersukan"="sport", "Lain-lain"="other", "Lepak-lepak di dalam rumah"
                                    ="none", "Lepak-lepak di luar rumah"="outside house", "Memancing ikan"=
                                      "fishing", "Memburu"="other", "Tiada"="none"))
x$activity[x$activity==""]<- "none"
x$activity <- factor(x$activity, levels = c("none", "outside house", "sport", "fishing", "other"))
m20 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m20)
lrtest(m16, m20)

## Bathe outside
x$bathOut <- factor(str_detect(x$bathPlace, "di luar"))
x$bathOut <- factor(x$bathOut)
m21 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + bathOut + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m21)
lrtest(m21, m20)

## Go to forest
x$goForest[x$goForest == ""] <- "N"
x$goForest <- factor(x$goForest)
m22 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + goForest + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m22)
lrtest(m20, m22)

## Walk to work through forest
m23 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + forestPath + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m23)
lrtest(m20, m23)

## Frequency of farm work (none, less than 15 days, over 15 days per month)
x$goToFarm <- findInterval(x$goToFarm, c(0, 15))
x$goToFarm[is.na(x$goToFarm)] <- 0
x$goToFarm <- factor(x$goToFarm)
m24 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + goToFarm + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m24)
lrtest(m20, m24)

## River near house
x$river <- factor(x$river)
m25 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + river + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m25)
lrtest(m20, m25)

## Floor
x$floor <- factor(x$floor)
m26 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + floor + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m26)
lrtest(m20, m26)

## Well near house
x$well <- factor(x$well)
m27 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m27)
lrtest(m20, m27)

## Pond near house
x$pond <- factor(x$pond)
m28 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + pond + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m28)
lrtest(m27, m28)

## Distance from clinic (by quartile)
x$clinc <- findInterval(x$clin, c(0,3.7, 6.8, 10))
x$clinc <- factor(x$clinc)
m29 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m29)
lrtest(m27, m29) 

## Walk to work
m30 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + walk + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m30)
lrtest(m29, m30) 

## Household farming rubber
x$rub <- factor(x$rub)
m31 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m31)
lrtest(m29, m31) 

## Medicine from clinic
x$med <- factor(str_detect(x$feverMonth, "umonkey dari klinik"))
m32 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m32)
lrtest(m31, m32) 
sjPlot::tab_model(m32)

## Mosquito prevention - insecticide
x$insect <- factor(str_detect(x$preventMosq, "Ubat"))
m33 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + insect + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m33)
lrtest(m33, m32) 

## Time travelling to/from work in early morning (11pm - 6am)
x$goAM <- factor(str_detect(x$goTime, "awal pagi"))
x$backAM <- factor(str_detect(x$backTime, "awal pagi"))
x$AM <- factor(ifelse(x$goAM == "TRUE" | x$backAM == "TRUE", "Y", "N"))
m34 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + AM + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m34)
lrtest(m34, m32) 

## Previous malaria
x$hadMalaria[x$hadMalaria==""] <- "N"
x$hadMalaria <- factor(x$hadMalaria)
m35 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + hadMalaria + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m35)
lrtest(m35, m32)

## Time travelling to/from work in evening (5pm - 10pm)
x$goPM <- factor(str_detect(x$goTime, "10pm"))
x$backPM <- factor(str_detect(x$backTime, "10pm"))
x$PM <- factor(ifelse(x$goPM == "TRUE" | x$backPM == "TRUE", "Y", "N"))
m36 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + PM + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m36)
lrtest(m36, m32)

## Any forest travel at night
x$forGoAM <- factor(str_detect(x$goToForest, "awal pagi"))
x$forBackAM <- factor(str_detect(x$backForest, "awal pagi"))
x$forAM <- factor(ifelse(x$forGoAM == "TRUE" | x$forBackAM == "TRUE", "Y", "N"))
x$forAM <- factor(x$forAM)
x$forGoPM <- factor(str_detect(x$goToForest, "10pm"))
x$forBackPM <- factor(str_detect(x$backForest, "10pm"))
x$forPM <- factor(ifelse(x$forGoPM == "TRUE" | x$forBackPM == "TRUE", "Y", "N"))
x$forPM <- factor(x$forPM)
x$forNight <- factor(ifelse(x$forAM == "Y" | x$forPM == "Y", "Y", "N"))
x$forNight <- factor(x$forNight)
m37 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + forNight + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m37)
lrtest(m37, m32)

## Distance of farm from house
x$farmDist <- factor(x$farmDist, levels = c("none", "house", "kampung", "outside"))
m38 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + farmDist + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m38)
lrtest(m38, m32)

# Time travelling to/from forest in evening (5pm - 10pm)
x$forGoPM <- factor(str_detect(x$goToForest, "10pm"))
x$forBackPM <- factor(str_detect(x$backForest, "10pm"))
x$forPM <- factor(ifelse(x$forGoPM == "TRUE" | x$forBackPM == "TRUE", "Y", "N"))
x$forPM <- factor(x$forPM)
m39 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + forPM + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m39)
lrtest(m39, m32)

## Occupation place
m40 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + occPlace + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m40)
lrtest(m32, m40)

## Household has livestock
x$livestk <- factor(x$livestk)
m41 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + livestk + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m41)
lrtest(m32, m41)

## Wall type - wood or bamboo walls
x$wall <- factor(x$wall)
m42 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + wall + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m42)
lrtest(m32, m42)

## Buffalo
x$buff <- factor(x$buff)
m43 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + buff + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m43)
lrtest(m32, m43)

## Location stayed outside kampung
x$stayOut <- revalue(x$stayOut, c("Dalam atau dekat kampung lain dalam daerah sama"="district", 
                                  "Di luar daerah dalam Sabah"= "outside", "Di luar Sabah"="outside",
                                  "Ladang/ estet dekat kampung ini"="plantation", "Hutan dekat kampung ini"
                                  ="forest"))
x$stayOut <- factor(x$stayOut)
m44 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + stayOut + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m44)
lrtest(m32, m44)

## House age
x$houseAge <- factor(x$houseAge)
m45 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + houseAge + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m45)
lrtest(m32, m45)

#household farms palm oil
x$po <- factor(x$po)
m46 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m46)
lrtest(m32, m46)


## Amount of land farmed
x$farmNow <- factor(x$farmNow)
m47 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m47)
lrtest(m46, m47)

## Kitchen outside house
x$kitch <- factor(x$kitch)
m48 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m48)
lrtest(m46, m48)

## Have stayed outside kampung in past month
x$outsideKmpg[x$outsideKmpg==""]<- "N"
x$outsideKmpg <- factor(x$outsideKmpg)
m49 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + outsideKmpg + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m49)
lrtest(m48, m49)

## Gender
m50 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + gender + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m50)
lrtest(m48, m50)

## Lake near house
x$lake <- factor(x$lake)
m51 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + lake + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m51)
lrtest(m48, m51)

## Activities in forest - hunting
x$hunt <- factor(str_detect(x$actForest, "Memburu"))
x$hunt <- factor(x$hunt)
m52 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + hunt + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m52)
lrtest(m48, m52)

## Meet doctor in hospital
x$doc <- factor(str_detect(x$feverMonth, "doktor hospital"))
m53 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m53)
lrtest(m48, m53)

## Goats
x$goat <- factor(x$goat)
m54 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m54)
lrtest(m48, m54)

## Meet traditional healer
x$heal <- factor(str_detect(x$feverMonth, "bantuan pawang"))
m55 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + heal + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m55)
lrtest(m54, m55)

## Length of time at house
x$stayPeriod <- factor(x$stayPeriod)
m56 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + stayPeriod + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m56)
lrtest(m54, m56)

## Household farms fruit
x$fruit <- factor(x$fruit)
m57 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + stayPeriod  + fruit + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m57)
lrtest(m54, m57)

## Time traveling to/from forest in early morning (11pm - 6am)
x$forGoAM <- factor(str_detect(x$goToForest, "awal pagi"))
x$forBackAM <- factor(str_detect(x$backForest, "awal pagi"))
x$forAM <- factor(ifelse(x$forGoAM == "TRUE" | x$forBackAM == "TRUE", "Y", "N"))
x$forAM <- factor(x$forAM)
m58 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + stayPeriod  + fruit + forAM + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m58)
lrtest(m58, m57)

## Pigs
x$pig <- factor(x$pig)
m59 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + stayPeriod + fruit + pig + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m59)
lrtest(m57, m59)


## Swiddening
x$swid <- factor(x$swid)
m60 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + stayPeriod  + fruit + swid + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m60)
lrtest(m57, m60)

## Collect wood from forest
x$collectWood <- factor(x$collectWood)
m61 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch  + doc + goat + stayPeriod  + fruit + collectWood + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m61)
lrtest(m57, m61)

## Windows that can close
x$wind <- factor(x$wind)
m62 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + doc + goat + stayPeriod + fruit + collectWood + wind + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m62)
lrtest(m57, m62)

## Bathe at the river
x$bathRiv <- factor(str_detect(x$bathPlace, "Sungai"))
x$bathRiv <- factor(x$bathRiv)
m63 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + doc + goat + stayPeriod + fruit + collectWood + wind + bathRiv + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m63)
lrtest(m62, m63)

## Bathe inside
x$bathIn <- factor(ifelse(x$bathOut == "FALSE" & x$bathRiv == "FALSE", "Y", "N"))
x$bathIn <- factor(x$bathIn)
m64 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + doc + goat + stayPeriod + fruit + collectWood + wind + bathRiv + bathIn + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m64)
lrtest(m64, m62)

## Time of bathing
x$bathT <- x$BilaMasaBiasanyaKamuMandi
x$bathAM <- factor(str_detect(x$bathT, "awal pagi"))
x$bathPM <- factor(str_detect(x$bathT, "10pm"))
x$bathAMout <- factor(ifelse(x$bathIn == "N" & x$bathAM == "TRUE", "Y", "N"))
x$bathPMout <- factor(ifelse(x$bathIn == "N" & x$bathPM == "TRUE", "Y", "N"))
x$bathNite <- factor(ifelse(x$bathAMout=="Y" | x$bathPMout =="Y", "Y", "N"))
m65 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + doc + goat + stayPeriod + fruit + collectWood + wind + bathNite + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m65)
lrtest(m65, m62)

## Time outside house in evening
x$TimeOut <- x$BerapaLamaSelalunyaAndaBeradaDiLuarRumahPadaWaktuMalam6PetangHin
x$TimeOut[x$TimeOut==""]<- "Tidak tahu"
x$TimeOut <- factor(x$TimeOut)
m66 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + doc + goat + stayPeriod + fruit + collectWood + wind + TimeOut + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m66)
lrtest(m66, m62)

## Plastic containers of water near house
x$plast <- factor(x$plast)
m67 <- glmer(cbind(monkey, 365 - monkey) ~ sea + hospc + demc + occ + hgt+ ageCat + wealth + wood + activity + well + clinc + rub + med + po + farmNow + kitch + doc + goat + stayPeriod + fruit + collectWood + wind + plast + (1 | houseID), data = x, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
summary(m67)
lrtest(m67, m62)

summary(m62) #final model
sjPlot::tab_model(m62)
plot(m62, pch = 20, col = "black", lty = "dotted")


plot_model(m62, sort.est = TRUE) + 
  theme(panel.background = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black", size = 14),
        text = element_text(size = 14),
        plot.title = element_text(size=14)) + ggtitle("Odds ratios and 95% confidence 
intervals for NHP sightings") +
  labs(x = "Household and individual level risk factors")+
  scale_x_discrete(labels=c("ageCat" = "Age group (every 15 years)",
                            "hospc" = "Travel time to nearest hospital (in quartiles)",
                            "sea" = "House near sea",
                            "farmWorkY" = "Engages in farm work",
                            "corn" = "Household farms corn", 
                            "bdgap" = "Gaps in eaves of house",
                            "fruit" = "Household farms fruit",
                            "demc" = "Household elevation (in  
meters above sea level)", 
                            "doc" = "Treatment-seeking behavior during fever: 
go to hospital",
                            "wind" = "Number of windows in house that can close",
                            "woodTRUE" = "Household collects wood from forest",
                            "activityother" = "Other evening activities outside the house",
                            "activityoutside house" ="Spends time outside the house in the evening",
                            "activitysport" = "Plays sports in the evening",
                            "activityfishing" = "Fishing in the evening")) 
