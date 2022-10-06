library(tidyverse)

setwd("~/CSIRO-ANACC/Genotyping/Genotyping/gstacks2/")

fst <- read_delim("populations.fst_AI-FMC.tsv", delim = "\t", col_names = TRUE)
names(fst) <- gsub(" ", "_", names(fst))

x <- boxplot(fst$AMOVA_Fst)

sort(x$out)

ggplot(fst, aes(`# Locus ID`, ``))