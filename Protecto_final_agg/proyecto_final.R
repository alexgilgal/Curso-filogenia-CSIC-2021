#Proyecto Alejandro Gil

require(phytools)
require(geiger)
require(lmtest)
require(diversitree)
require(nlme)

# Cargamos los datos

# Dataframe:

transposons <- read.delim("Tarea2_Transposable_elements.tsv", dec = ",")
transposons$Species <- gsub(" ", "_", transposons$Species)

tree <- read.newick("Ensembl_tree.newick") # this is a tree with 1992 species of teleost fishes

plot(tree, type = "fan")
# Necesitamos ajustar el árbol a los datos

tree$tip.label # necesitamos usar grep buscar los nombres del Dataframe en nuestro árbol

selected_species <- NULL

for (i in transposons$Species) {
  
  selected_species <- c(selected_species, grep(i, tree$tip.label, value = T, ignore.case = T))
  
}

selected_species # we end up with 26 species with filogenetic info

tree_pruned <- drop.tip(tree,tree$tip.label[! (tree$tip.label %in% selected_species)] )

# tenemos que poner en el mismo orden los datos de transposons, para que haya correspondencia entre nombres
tree_pruned$tip.label

transposons$Species # aquí tenemos las 30 especies
transposons_filtered <- transposons
keep <- NULL

for (i in transposons_filtered$Species){
  if (any(grepl(i, tree_pruned$tip.label, ignore.case = T))) {
    keep <- c(keep, TRUE)
    transposons_filtered$Species[transposons_filtered$Species == i] <- grep(i, tree_pruned$tip.label, value = T, ignore.case = T)
  } else {
    keep <- c(keep, FALSE)
  }
}

transposons_filtered <- transposons_filtered[keep,]

rownames(transposons_filtered) <- tolower(transposons_filtered$Species)

# sólo nos falta ponerlos en el mismo orden

transposons_filtered <- transposons_filtered[tree_pruned$tip.label,]

transposons_filtered$Species
tree_pruned$tip.label

plot.phylo(tree_pruned, type = "fan")

transposon.K<-transposon.l <- transposon.K.pval <- transposon.l.pval <- vector()

for (i in colnames(transposons_filtered[, -c(1, 27)])) {
  
  message(i)
  
  k <- phylosig(tree = tree_pruned,
                x =transposons_filtered[,i] , method="K",nsim=100, test=T)
  transposon.K[i] <- k$K
  
  transposon.K.pval[i] <- k$P
  
  l <- phylosig(tree = tree_pruned,
                x =transposons_filtered[,i], method="lambda",nsim=100, test=T)
  
  transposon.l[i] <- l$lambda
  
  transposon.l.pval[i] <- l$P
}

Transposon_signal <- data.frame(Variables = colnames(transposons_filtered[, -c(1, 27)]),
                                K = transposon.K,
                                K_pvalue = transposon.K.pval,
                                Lambda = transposon.l,
                                Lambda_pvalue = transposon.l.pval)

Transposon_signal # Este dataset puede ser util para el proyecto final

write.table( Transposon_signal[c(5:8,29),], "~/Dropbox/CABD/2021/curso filogenia/signal.txt", row.names = F, quote = F, sep = "\t")

plot(tree, show.tip.label = T, type = "fan", cex = 0.4,no.margin = T)

plot(tree_pruned,  show.tip.label = T, type = "fan", cex = 0.8,no.margin = F)

trans_matrix <- as.matrix(transposons_filtered[,c("genome_length","repeat_length","CLASSI", "CLASSII")])

colnames(trans_matrix) <- c("Tamaño \n genoma", "Tamaño \n repeticiones", "Clase I", "Clase II")
require(viridisLite)
phylo.heatmap(tree_pruned, trans_matrix, standardize = T, colors = viridis(100), mar=c(4,4,4,4))

## RCA

# CR1 superfamily

cr1 <- transposons_filtered$CR1

names(cr1) <- transposons_filtered$Species

# Vamos a ver el fenograma

phenogram(tree_pruned, cr1) # parece que sigue un modelo browniano con tendencia

cr1_ancestral <- anc.ML(tree_pruned, cr1, maxit=2000,model = "BM")

plot(tree_pruned, show.tip.label=T, no.margin = T)#   
tiplabels(round(cr1,2), frame = "none", adj = 2)# valores actuales 
nodelabels(round(cr1_ancestral$ace,2), cex=.7) # valores reconstruidos

fancyTree(tree_pruned,type="phenogram95",x = cr1,spread.cost=c(1,0))
## PPCA

sapply(transposons_filtered, class)

transposons_filtered2 <- transposons_filtered[,c("genome_length","repeat_length","CLASSI", "CLASSII")]

log1=function(x){return(log(x+1))}
require(plyr)
scaler=function(x){s=scale(x, center=T);attributes(s)=NULL;return(s) }
df_scaled=colwise(log1)(transposons_filtered2)
df_scaled=colwise(scaler)(df_scaled)
rownames(df_scaled)=rownames(transposons_filtered2)

colMeans(df_scaled)

PCA_phylo <- phyl.pca(tree_pruned,df_scaled,method="BM", mode="corr")

summary(PCA_phylo)

PCA_phylo$Evec


## Regresión

# we have to use log 10 to normalize the data, but it is still interpretable

x <- log10(transposons_filtered$repeat_length)

y <- log10(transposons_filtered$genome_length)

df <- data.frame(species = transposons_filtered$Species,
                 x = x,
                 y = y)

# Vamos a generar diferentes modelos, usando ML para comparar

m1 <- gls(y~x,method="ML", corBrownian(1, tree_pruned,form = ~species),data =  df)
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m1)) #  Por los pelos, pero siguen una normal
plot(residuals(m1)~x, main = "Residuals Brownian model")

m2 <- (gls(y~x, method="ML",  corMartins(1, tree_pruned,form = ~species,fixed = T),data =  df))
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m2)) # no lo son
plot(residuals(m2)~x, main = "Residuals OU Model")

m3 <- (gls(y~x, method="ML",  corBlomberg(1, tree_pruned ,form = ~species, fixed = T),data =  df))
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m3)) # Por los pelos, pero siguen una normal
plot(residuals(m3)~x, main = "Residuals Blomberg model (ACDC)")

m4 <- (gls(y~x, method="ML",  corPagel(0.5, tree_pruned ,form = ~species, fixed = F),data =  df))
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m4)) # No lo son
plot(residuals(m4)~x, main = "Residuals Pagel model (Lambda)")


AIC(m1,m2,m3,m4)
BIC(m1,m2,m3,m4)

# vamos a verlos en una tabla

indices <- data.frame(df = AIC(m1,m2,m3,m4)$df, AIC = AIC(m1,m2,m3,m4)$AIC, BIC =BIC(m1,m2,m3,m4)$BIC )

indices

plot(df$y ~ df$x, main = "Regresión bajo modelo BM", ylab = "log10(Genome length)", xlab = "log10(Repeats)")

abline(m1, lty = 4, lwd = 2, col = "red")

summary(m1)$tTable

## regresión múltiple


y <- log10(transposons_filtered$genome_length +1 )

x1 <- log10(transposons_filtered$CLASSI +1)

x2 <- log10(transposons_filtered$CLASSII +1)

df <- data.frame(species = transposons_filtered$Species,
                 x1 = x1,
                 y = y,
                 x2 = x2)

plot(y ~ x1)

plot(y ~ x2)
require(nortest)

m1 <- gls(y~x1+x2,method="ML", corBrownian(1, tree_pruned,form = ~species),data =  df)
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m1)) #  No lo son
plot(residuals(m1)~x, main = "Residuals Brownian model")

m2 <- (gls(y~x1+x2, method="ML",  corMartins(1, tree_pruned,form = ~species,fixed = T),data =  df))
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m2)) # no lo son
plot(residuals(m2)~x, main = "Residuals OU Model")

m3 <- (gls(y~x1+x2, method="ML",  corBlomberg(1, tree_pruned ,form = ~species, fixed = T),data =  df))
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m3)) # Por los pelos, pero siguen una normal
plot(residuals(m3)~x, main = "Residuals Blomberg model (ACDC)")

m4 <- (gls(y~x1+x2, method="ML",  corPagel(0.5, tree_pruned ,form = ~species, fixed = F),data =  df))
shapiro.test(chol(solve(vcv(tree_pruned)))%*%residuals(m4)) # No lo son
plot(residuals(m4)~x, main = "Residuals Pagel model (Lambda)")


indices <- data.frame(df = AIC(m1,m2,m3,m4)$df, AIC = AIC(m1,m2,m3,m4)$AIC, BIC =BIC(m1,m2,m3,m4)$BIC )

indices

summary(m1)$tTable
