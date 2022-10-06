#install.packages("RcppArmadillo")
#install.packages("CompQuadForm")
#install.packages("foreach")
#install.packages("Matrix")
#install.packages("testthat")
#BiocManager::install(c("SeqArray", "SeqVarTools"))
#devtools::install_github("hanchenphd/GMMAT")


library(vcfR)
library(tidyverse)
library("GWASTools")
library("SNPRelate")
library(GMMAT)



#clean environments
rm(list = ls())

setwd("~/Cintia/Genotyping/gstacks2/statistics/")


#Import phenotype and info sample
pop_map <- read_delim("phenotype.txt", delim = "\t", col_names = c("SampleID", "Population", "Region", "Phenotype"))
head(pop_map)

pop_map %>% mutate(pheno = case_when(
  Phenotype == "NO_TOLERANT" ~ 0,
  Phenotype == "THERM_TOLERANT" ~ 1
))

colnames(pop_map$SampleID)

#Convert to GDS
vcf.fn <- "C:/Users/iha002/Documents/Cintia/Genotyping/gstacks2/filtering/populations.filtered.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "snps.gds", verbose = TRUE)
snpgdsSummary("snps.gds")
gds <- GdsGenotypeReader("snps.gds")
#snpID <- getSnpID(gds)
#chromosome <- as.integer(getChromosome(gds))
#position <- getPosition(gds)
#alleleA <- getAlleleA(gds)
#alleleB <- getAlleleB(gds)
#rsID <- getVariable(gds, "snp.rs.id")
#filter <- getVariable(gds, "snp.annot/filter")
#snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position, rsID, alleleA, alleleB, filter, stringsAsFactors = FALSE))
#genoData <- GenotypeData(gds, snpAnnot = snpAnnot)

#Matrices of covariance structure
#Using Standardized GRM from GEMMA
grm <- as.matrix(read.table("Vmatrix2.sXX.txt", check.names = FALSE))
head(grm)

model0 <- glmmkin(Phenotype ~ )

