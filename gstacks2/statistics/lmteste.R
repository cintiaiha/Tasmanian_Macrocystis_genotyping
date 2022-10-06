library(tidyverse)
library(vcfR)

rm(list = ls())
setwd("~/Cintia/Genotyping/gstacks2/statistics/")

#import vcf file
vcf <- read.vcfR("populations.nomiss.recode.vcf", verbose = FALSE)
show(vcf)

#Extract only genotyping and transform as numeric and transpose
gt <- as.data.frame(t(extract.gt(vcf, element = "GT", as.numeric = TRUE)))

#Select SNPs that are NOT iqual in all samples
gt <- select_if(gt, !(apply(gt, 2, function(a) length(unique(a))==TRUE)))

gt <- rownames_to_column(gt, var = "id")
rownames(gt)
head(gt)


#Add phenotype
pop <- read.table("phenotype.txt", col.names = c("id", "population", "region", "phenotype"))
#pop <- pop %>% mutate(phenotype = case_when(
#  phenotype == "NO_TOLERANT" ~ 0,
#  phenotype == "THERM_TOLERANT" ~ 1
#))
head(pop)

#compare if samples are in the same order
gt$id == pop$id

#Add phenotype to gt
data <- merge(pop, gt, by = "id")
head(data)
write.table(data, file = "snpdata.txt", sep = "\t")

therm <- data %>% filter(phenotype == "THERM_TOLERANT") %>% arrange(region)




####################################################
# test 1 SNP
test <- data[,1:5]
head(test)

mod1 <- lm(formula = phenotype~`37332:85:+`, data)
layout(matrix(1:4, 2))
plot(mod1)
summary(mod1)

####################################################

depVarlist <- setdiff(colnames(data), c("id", "population", "region", "phenotype"))

allmodels <- lapply(depVarlist, function(x){
  lm(formula = paste0("`", x, "` ~ phenotype"), data = data, na.action = na.omit)
})

names(allmodels) <- depVarlist

allcoef <- data %>% select(-phenotype, -population, -region)

allcoef[,-1] <- sapply(2:ncol(allcoef), function(x){
  coef <- allcoef[,x]
  coef[!is.na(coef)] = allmodels[[x-1]]$coef
  coef
})




allResiduals <- data %>% select(-phenotype, -population, -region)

allResiduals[,-1] = sapply(2:ncol(allResiduals), function(x){
  residuals = allResiduals[,x]
  residuals[!is.na(residuals)] = allmodels[[x-1]]$residuals
  residuals
})

sapply(allmodels, coef)

summaries <- lapply(allmodels, summary)
coef <- lapply(summaries, function(x) x$coefficients[, c(1,4)])


?sapply()
?all




allResiduals$

plot(allResiduals)




#Random sampling

no_tolerant_d <- data %>% filter(phenotype == "NO_TOLERANT")

no_tolerant_d2 <- no_tolerant_d %>% filter(!(id %in% group1))

no_tol_id2 <- no_tolerant_d2$id


group1 <- sample(no_tol_id, size = 7)
group2 <- sample(no_tol_id2, size = 7)



################ DRAFT ################

#convert vcfR to tidy format
tidy <- vcfR2tidy(vcf)
head(tidy)
names(tidy)
tidy$meta
tidy$gt$ChromKey
tidy$fix$ID