library(vcfR)
library(adegenet)
library(poppr)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ape)


##########################################################
### Trying my own pipeline -> maybe stupid thing to do ###
########################################################## 
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
#tutorial dapc
#https://www.nature.com/articles/s41437-020-0348-2

rm(list = ls())
setwd("~/CSIRO-ANACC/Genotyping/Genotyping/gstacks2/")


####################
### Filtered VCF ###
####################


#Import vcf
vcf <- read.vcfR("populations.filtered.recode.renamed.vcf", verbose = FALSE)
show(vcf)
## ***** Object of Class vcfR *****
## 43 samples
## 1505 CHROMs
## 5,911 variants
## Object size: 19.7 Mb
## 0 percent missing data
## *****        *****         ****


######################################
#####   Genotyping analysis     ######
######################################

#Import data attributes
pop.data <- read_delim("pop_map2.renamed", delim = "\t", col_names = c("SampleID", "Population", "Region", "Phenotype"))
#head(pop.data)

pop.data$Population <- factor(pop.data$Population, levels=c("A", "B", "C", "D", "E", "F"))

#check that all the samples in the VCF and the population data frame are included
#colnames(gstacks_vcf@gt)[-1]
#pop.data$SampleID
all(colnames(vcf@gt)[-1] == pop.data$SampleID)
## TRUE

#Converting to genlight object
gl <- vcfR2genlight(vcf)
pop(gl) <- pop.data$Population
#gl@pop
#gl@ind.names

rm(vcf, pop.data)

cols <- brewer.pal(n = nPop(gl), name = "Dark2")

###PCA

pca <- glPca(gl, nf=40)
pca$scores

barplot(100*pca$eig/sum(pca$eig), col = heat.colors(50), main = "PCA Eigenvalues")
head(100*pca$eig/sum(pca$eig))
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

pca.scores <- as.data.frame(pca$scores)
pca.scores$pop <- pop(gl)

set.seed(9)
ggplot(pca.scores, aes(x=PC1, y=PC2, colour=pop)) + 
  geom_point(size=2) + 
  stat_ellipse(level = 0.9, size = 0.5) + 
  scale_color_manual(name="Populations", values = cols) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  labs(x="PC1 (32.7%)", y="PC2 (6.2%)") +
  theme_bw()


### DAPC

#Infer number of clusters K
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(gl, n.pca=40, choose.n.clust = FALSE, max.n.clust = maxK, criterion = "diffNgroup")
  myMat[i,] <- grp$Kstat
}

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)

grp$grp

Kplot <- ggplot(my_df, aes(K, BIC)) + geom_boxplot() + theme_bw() + xlab("Number of groups (K)")
table(pop(gl), grp$grp)

table.value(table(pop(gl), grp$grp), col.labels = paste("Inferred cluster", 1:6), row.labels = levels(pop(gl)))


dapc1 <- dapc(gl, n.pca = 2, n.da = 2)
scatter(dapc1, scree.da=TRUE, posi.da = "bottomright", col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "topleft", cleg = 0.85, scree.pca = TRUE)

dapc2 <- dapc(gl,n.pca = 2, pop = gl$pop, n.da = 2)
scatter(dapc2, scree.da=TRUE, posi.da = "bottomright", col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "topleft", cleg = 0.85, scree.pca = TRUE)
optim <- optim.a.score(x = dapc2,n.pca = 1:ncol(dapc2$tab),plot = F,n.sim = 100)
dapc.optim<-dapc(gl, n.pca = optim$best, pop = gl$pop, n.da = 2)
scatter(dapc.optim, scree.da=TRUE, posi.da = "bottomright", col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "topleft", cleg = 0.85, scree.pca = TRUE)
scatter(dapc.optim, scree.da=TRUE, posi.da = "bottomright", col = cols, cex = 3, posi.leg = "topleft", cleg = 0.85, scree.pca = TRUE)

compoplot(dapc1, col = cols, posi = 'top')

dapc1$posterior



dapc.results <- as.data.frame(dapc1$posterior)
dapc.results$pop <- pop(gl)
dapc.results$indNames <- rownames(dapc.results)

library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)


dapc.results$name <- factor(dapc.results$name, levels = c("A", "B", "C", "D", "E", "F"))


colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")


t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) + geom_bar(stat = 'identity') +
  scale_fill_manual(name="Populations", values = cols) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(y="Posterior membership probability", x=NULL)






## Adegenet viginnet

x <- find.clusters(gl, n.pca=40, choose.n.clust = FALSE, max.n.clust = 10, criterion = "diffNgroup")
plot(x$Kstat)
x$stat
x$grp

dapc.x <- dapc(gl, var.contrib = TRUE, scale = FALSE, n.pca = 40, pop = x$grp, n.da = 10)
dapc.x$n.da
#scatter(dapc.x, scree.da=FALSE, col = cols, cex = 3, legend = TRUE, clabel = F, posi.leg = "bottomleft", posi.pca = "topleft", cleg = 0.85)

optim.x <- optim.a.score(x = dapc.x, n.pca = 1:ncol(dapc.x$tab), plot = F, n.sim = 100)
optim.x$best


optim.dapc.x <- dapc(gl, var.contrib = TRUE, scale = FALSE, n.pca = 2, pop = x$grp, n.da = 10)
scatter(optim.dapc.x, scree.da=FALSE, col = cols, cex = 3, legend = TRUE, clabel = F, posi.leg = "bottomleft", posi.pca = "topleft", cleg = 0.85)



#find.grp <- find.clusters(x = genID, max.n.clust=10,n.pca = 500, choose.n.clust = F, criterion = "diffNgroup")
#dapc.with.groups<-dapc(genID, var.contrib = TRUE, scale = FALSE, n.pca = 30, pop = find.grp$grp, n.da = 5)
#optim<-optim.a.score(x = dapc.with.groups,n.pca = 1:ncol(dapc.with.groups$tab),plot = F,n.sim = 100)
#optim.dapc.with.groups<-dapc(genID, var.contrib = TRUE, scale = FALSE, n.pca = optim$best, pop = find.grp$grp, n.da = 5)


plot(x$Kstat)

table(pop(gl), grp$grp)
table.value(table(pop(gl), grp$grp), col.labels = paste("Infered_cluster", 1:6), row.labels = levels(pop(gl)))



#test many values of K

my_k <- 2:6

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  grp_l[[i]] <- find.clusters(gl, n.pca=40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")

ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) + 
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))


tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Population <- pop.data$Population
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Population <- pop.data$Population
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Population, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_color_brewer(palette="Dark2")
p3 <- p3 + scale_fill_manual(values=c(my_pal))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3



adegenetTutorial("dapc")
###



dapc <- dapc(gl, n.pca = 3, n.da = 2)
scatter(optim.dapc.x, scree.da=FALSE, col = cols, cex = 3, legend = TRUE, clabel = F, posi.leg = "topleft", cleg = 0.85)

compoplot(optim.dapc.x, col = cols, posi = 'top')


dapc.results <- as.data.frame(optim.dapc.x$posterior)
dapc.results$pop <- pop(gl)
dapc.results$indNames <- rownames(dapc.results)



library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
head(dapc.results, n = 6)


dapc.results$name <- factor(dapc.results$name, levels = c("AI", "PI", "HI", "FMC", "SC", "SH"))


colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")



t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) + geom_bar(stat = 'identity') +
  scale_fill_manual(name="Populations", values = cols) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(y="Posterior membership probability", x=NULL)

library(ggpubr)
library(reshape2)
ggarrange(dapcp, t, ncol=1, nrow =2)


## Genetic differentiation

pop <- factor(pop.data$SampleID)

myDiff <- genetic_diff(gstacks_vcf, pops = pop, method = 'nei')
dpf <- melt(myDiff[,c(46, 90, 92, 93)], varnames=c('Index', 'Sample'), value.name = "Depth", na.rm=TRUE)

ggplot(dpf, aes(x=variable, y=Depth)) + geom_violin(fill='blue', adjust = 1.2)


################################################################################
##Try linkage desequilibrium

#This test is useful to determine if populations are clonal (where significant disequilibrium is expected due to linkage among loci)
#or sexual (where linkage among loci is not expected). The null hypothesis tested is that alleles observed at different loci are not
#linked if populations are sexual while alleles recombine freely into new genotypes during the process of sexual reproduction.
#In molecular ecology we typically use the index of association or related indices to test this phenomenon.

library("SNPRelate")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/radiator")
library(radiator)
require(SeqArray)

setwd("C:/Users/iha002/Genotyping_test")

summary_strata("strata.txt")

data <- read_vcf(data = "populations.filtered.recode.vcf", strata = "strata.txt", 
                 parallel.core = 1L)


d1 <- filter_ld(data, filter.short.ld = "mac", filter.long.ld = NULL, parallel.core = 1L)




n############################################################################################
colnames(gstacks_vcf@gt)[-1]
pop <- as.factor(c(colnames(gstacks_vcf@gt)[-1]))

myDiff <- genetic_diff(gstacks_vcf, pops=pop, method = 'nei')
head(myDiff)

library(reshape2)
head(myDiff)

dpf <- melt(myDiff[,c(3:8,19)], varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
