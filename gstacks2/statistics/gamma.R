library(tidyverse)

rm(list = ls())

setwd("~/Cintia/Genotyping/gstacks2/statistics/")

#http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf

hyp.params <- read.table("bslmm2.hyp.txt", header = T)
head(hyp.params)

# Get mean, median, and 95% ETPI of hyperparameters
# h-> approximation to proportion of phenotypic variance explained by variants (PVE)

h <- c("h", mean(hyp.params$h), quantile(hyp.params$h, probs = c(0.5,0.025,0.975)))

pve <- c("PVE", mean(hyp.params$pve), quantile(hyp.params$pve, probs = c(0.5,0.025,0.975)))

# rho-> approximation to proportion of genetic variance explained by variants with major effect (PGE)
# rho=0 -> pure LMM, highly polygenic basis
# rho=1 -> pure BVSR, few major effect loci
rho <- c("rho", mean(hyp.params$rho), quantile(hyp.params$rho, probs = c(0.5,0.025,0.975)))

pge<- c("PGE", mean(hyp.params$pge), quantile(hyp.params$pge, probs = c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi <- c("pi", mean(hyp.params$pi), quantile(hyp.params$pi, probs = c(0.5,0.025,0.975)))

#n.gamma -> number of variants with major effect
n.gamma <- c("n.gamma", mean(hyp.params$n_gamma), quantile(hyp.params$n_gamma, probs = c(0.5,0.025,0.975)))

# get table of hyperparameters and save it to a file

hyp.params.table<-as.data.frame(rbind(h,pve,rho,pge,pi,n.gamma),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
hyp.params.table

##  hyperparam              mean      median           2.5%        97.5%
##1          h  0.90261842134608   0.9970537   0.0699391025    0.9997154
##2        PVE  0.85538094504468    0.954801    0.031513958    0.9973944
##3        rho  0.93467606268984   0.9990107   0.2157856275    0.9997411
##4        PGE 0.877008001674224    0.990784   0.0238821455     0.997158
##5         pi  0.00966776157842 0.004249402 0.003024403075 0.0932391225
##6    n.gamma           14.5059           7              4          119



# plot traces and distributions of hyperparameters
# set up layout
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))


pdf(file="hyperparameters2.pdf", width=8.3, height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))
#h
plot(hyp.params$h, type="l", ylab="h", main="h")
hist(hyp.params$h, main="", xlab="h")
plot(density(hyp.params$h), main="", xlab="h")


# PVE
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE")
hist(hyp.params$pve, main="", xlab="PVE")
plot(density(hyp.params$pve), main="", xlab="PVE")

# rho
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(hyp.params$rho, type="l", ylab="rho", main="rho")
hist(hyp.params$rho, main="", xlab="rho")
plot(density(hyp.params$rho), main="", xlab="rho")

# PGE
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE")
hist(hyp.params$pge, main="", xlab="PGE")
plot(density(hyp.params$pge), main="", xlab="PGE")


# pi
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="", xlab="pi")
plot(density(hyp.params$pi), main="", xlab="pi")

#No gamma
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(hyp.params$n_gamma, type="l", ylab="No gamma", main="No gamma")
hist(hyp.params$n_gamma, main="No gamma", xlab="No gamma")
plot(density(hyp.params$n_gamma), main="No gamma", xlab="No gamma")

dev.off()

# Load parameters
params <- read.table("bslmm.param.txt", header = T)

# Get variants with sparse effect size on phenotypes add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)
head(params)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]
head(params.effects)

# sort by decreasing effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10)

# variants with the highest sparse effects
# top 1% variants (above 99% quantile)
top1<- params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]

# top 0.1% variants (above 99.9% quantile)
top01<- params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]

# top 0.01% variants (above 99.99% quantile)
top001<- params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]

##### PIP #####

# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC
# (the number of times it appears in the MCMC with a non-zero sparse effect)
# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]


# plot variants PIPs across linkage groups/chromosomes
# Prepare data
# sort by linkage group and position
chr<-c(gsub("C","100",params$chr))
chr<-c(gsub("scaffold","200",chr))
head(chr)
head(params$chr)
tail(chr)
tail(params$chr)
#change names chr
params["chr"]<-chr

params.sort<-params[order(as.numeric(params$chr), params$rs),]
head(params.sort)

# get list of linkage groups/chromosomes
chrs<-sort(as.numeric(unique(params$chr)))

# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,1),ylab="PIP",xlab="linkage group", xaxt="n")

# plot grey bands for chromosome/linkage groups
start<-1
lab.pos<-vector()
for (ch in chrs){
  size<-nrow(params.sort[params.sort$chr==ch,])
  cat ("CH: ", ch, "\n")
  colour<-"light grey"
  if (ch%%2 > 0){
    polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour,
            border=colour)
  }
  cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
  txtpos<-start+size/2
  lab.pos<-c(lab.pos, txtpos)
  start<-start+size
}

# Add variants outside linkage groups
chrs<-c(chrs,"NA")
size<-nrow(params.sort[params.sort$chr=="NA",])
lab.pos<-c(lab.pos, start+size/2)

# Add x axis labels
axis(side=1,at=lab.pos,labels=chrs,tick=F)

# plot PIP for all variants
# rank of variants across linkage groups
x<-seq(1,length(params.sort$gamma),1)
# PIP
y<-params.sort$gamma
# sparse effect size, used for dot size
z<-params.sort$eff
# log-transform to enhance visibility
z[z==0]<-0.00000000001
z<-1/abs(log(z))
# plot
symbols(x,y,circles=z, bg="black",inches=1/5, fg=NULL,add=T)

# highligh high PIP variants (PIP>=0.25)
# plot threshold line
abline(h=0.25,lty=3,col="dark grey")
# rank of high PIP variants across linkage groups
x<-match(params.sort$gamma[params.sort$gamma>=0.25],params.sort$gamma)
# PIP
y<-params.sort$gamma[params.sort$gamma>=0.25]
# sparse effect size, used for dot size
z<-params.sort$eff[params.sort$gamma>=0.25]
z<-1/abs(log(z))
symbols(x,y,circles=z, bg="red",inches=1/5,fg=NULL,add=T)
# ------------------------------------------------------------------------------
# add label high PIP variants
text(x,y,labels=params.sort$rs[params.sort$gamma>=0.25], adj=c(0,0), cex=0.5)
# ------------------------------------------------------------------------------
params.sort$rs[params.sort$gamma>=0.25]



############################ Linear model #####################################

rm(list = ls())

lmm <- read.table("lmm_all.assoc.txt", header = TRUE)
head(lmm)
ggplot(lmm, aes(rs, p_score)) + geom_point()
ggplot(lmm, aes(p_wald)) + geom_density()

lmm0.05pscore <- lmm %>% filter(p_score < 0.05)

plot(lmm$l_mle)

plot(lmm0.05pscore$p_score)
summary(lmm$p_score)
