# C	 assim. Esses data tem dados de SNPs de cada individuo e informaC'C5es a mais. Os SNPs sC#o valores de 0 e 1 mesmo.
# Eu quero testar se cada SNP tem relaC'C#o com o phenotype
# phenotype: Termo tolerante = 1; nC#o tolerante = 0

library(data.table) # para facilitar algumas operaC'C5es e para preservar os nomes dos SNPs
library(ggplot2)
setwd("~/Cintia/Genotyping/gstacks2/statistics/")

## Carregar e checar dados =====================================================

dat <- fread("snpdata.txt", drop = c(1,3,4))
# Vai dar um aviso, porque a primeira coluna estC! sem nome.
# Quando usar write.table e similares, se o data.frame nC#o tiver nomes para as linhas,
# use a opC'C#o `row.names = FALSE`, ou vai gerar uma coluna sem nome com os C-ndices.

# Cortar os SNPs que sC#o 0 ou 1 para todos, jC! que nC#o sC#o informativos:
colsel <- apply(dat, 2, function(x) all(x == 0) | all(x == 1))
dat <- dat[, !colsel, with = FALSE]
rm(colsel)

# Nomes das colunas com SNPs
snps <- names(dat)[-(1:2)]


## Ver se tem algum SNP que C) 1 para todos de um fenC3tipo e 0 para todos do outro

# Se todas as linhagens sC#o 0 para cada SNP por fenC3tipo
ntol <- dat[, lapply(.SD, function(x) all(x == 0)), by = phenotype, .SDcols = snps]

# Se todas as linhagens sC#o 1 para cada SNP por fenC3tipo
tol <- dat[, lapply(.SD, function(x) all(x == 1)), by = phenotype, .SDcols = snps]

# Agora um truque. Queremos saber qual C) TRUE para tol em um fenC3tipo e TRUE para ntol no outro fenC3tipo
# R codifica FALSE como 0 e TRUE como 1, entC#o ao soma fenC3tipo 0 de uma condiC'C#o com fenC3tipo 1 da outra,
# se houver algum TRUE para ambos o resultado serC! 2. Vamos checar se algum bate com essa condiC'C#o:
any(ntol + tol[2:1] > 1)
# Deu FALSO, entC#o nenhum SNP bate com todos.

rm(ntol, tol)


## Coeficiente de similaridade =================================================
# Essa C) a abordagem que acho mais adequada. Existem dezenas de mC)tricas de
# similaridade/dissimilaridade, algumas prC3prias para comparar variC!veis binC!rias.

source("binary_similarity_coefficients.R")

# Porque estC! desbalanceado, penso que o coeficientes de Consonni-Todeschini ou de Michael serC#o melhores
sims <- sapply(snps, function(x) similarity(dat$phenotype, dat[[x]], "michael"))

## GrC!fico bonitinho

p <- ggplot(as.data.table(sims, TRUE)[abs(sims) > .1],
            aes(reorder(rn, abs(sims)), sims)) +
       geom_col(aes(fill = as.factor(ifelse(sims > 0, "p", "n")))) +
       scale_fill_discrete(type = c('#533100', '#003B2F'), guide = "none") +
       theme_minimal() +
       coord_flip() +
       labs(x = NULL, y = "correlation to tolerance",
            caption = "only SNPs with correlation greater than ??0.1")

ggplot(as.data.table(sims, TRUE)[abs(sims) > .1],
            aes(reorder(rn, abs(sims)), sims)) +
  geom_col(aes(fill = as.factor(ifelse(sims > 0, "p", "n")))) +
  scale_fill_discrete(type = c('#533100', '#003B2F'), guide = "none") +
  theme_minimal() +
  coord_flip() +
  labs(x = NULL, y = "correlation to tolerance",
  caption = "only SNPs with correlation greater than ??0.1")

ggsave("snps.png", p, width = 6, height = 9, dpi = 100, bg = "white")


## CorrelaC'C#o ==================================================================
# NC#o C) ideal para variC!veis binC!rias. Mas o coeficiente phi para dados binC!rios
# se aproxima da correlaC'C#o de Pearson para valores baixos (caso dos seus dados).
# Se alguC)m achar ruim de usar dissimilaridade, pode usar correlaC'C#o. O resultado
# serC! o mesmo de usar similarity(..., method = "phi")

cors <- sapply(snps, function(x) cor(dat$phenotype, dat[[x]]))

# Verificando a diferenC'a:
plot(sims, cors)
abline(0,1)


## Modelo linear (lm) ==========================================================
# NC#o C) recomendado para variC!veis binC!rias, mas se a galera da C!rea costuma
# usar, vamos fazer, sC3 para comparar.

lms <- lapply(snps, function(snp) lm(dat[[snp]] ~ dat$phenotype))
names(lms) <- snps

# Uma vez que os resultados estC#o em uma lista, pode usar lapply para obter
# os coeficientes, resC-duos, etc, de todos. Pode obter de um sC3 usando os
# nomes, por exemplo coef(lms$`24659:130:+`)
# do.call("rbind", lista) vai unir as linhas de todos os resultados
coefs <- do.call("rbind", lapply(lms, coef))

# Comparando com similaridade/correlaC'C#o:
plot(sims, coefs[,2])
abline(0, 1)

# Parecido. Mas vai na similaridade, que C) mais adequada.
