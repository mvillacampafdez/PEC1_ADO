## LIBRARIES
library(GEOquery)
library(affy)
library(oligo)
library(arrayQualityMetrics)
library(ggplot2)
library(ggrepel)
library("pd.clariom.s.human")
library(limma)
library(hgu133plus2cdf)
require(genefilter)
library("hgu133plus2.db")
library(gplots)
library(ReactomePA)
require(org.Hs.eg.db)

# Lo primero es leer los ficheros para guardar las lecturas
if(!file.exists("geo_downloads")) dir.create("geo_downloads")
if(!file.exists("results"))  dir.create("results", recursive=TRUE)
GEOgse <- getGEO(GEO="GSE106151", filename=NULL, destdir="./geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
class(GEOgse) # 'list'
# como se obtiene un objeto de tipo lista, se descarta el resto de objetos para obtener solo el objeto tipo ExpressionSet
GEOgse <- GEOgse[[1]]
class(GEOgse) # 'EspressionSet'

# Se comprueba que se han cargado de forma correcta los datos
head(exprs(GEOgse))

if(!file.exists('./geo_downloads/GSE106151'))
  getGEOSuppFiles('GSE106151', makeDirectory=T, baseDir="geo_downloads")


untar("geo_downloads/GSE106151/GSE106151_RAW.tar", exdir="geo_downloads/GSE106151/CEL")
list.files("geo_downloads/GSE106151/CEL")

# se obtiene el listado de los ficheros de las 18 muestras de partida
my.cels <- sort(list.files('./geo_downloads/GSE106151/CEL'))
my.cels

# Se preparan los datos para guardar los datos fenotípicos:
phenoD <- as.data.frame(pData(GEOgse), stringsAsFactors=F)[, c("title", "geo_accession", "description")]
phenoD <- phenoD[order(rownames(phenoD)), ]

# Los nombres de las filas deben de corresponder con los nombres de los ficheros (en los datos fenotípicos):
rownames(phenoD) <- my.cels

write.table(phenoD, file=paste0('./geo_downloads/GSE106151/GSE106151_SelectPhenoData.txt'), sep="\t", quote=F)

## Se obtienen los datos de partida para estudiar: 'rawData'
rawData <- ReadAffy(celfile.path="geo_downloads/GSE106151/CEL", phenoData="geo_downloads/GSE106151/GSE106151_SelectPhenoData.txt")
head(rawData)

### 1. Identificar que grupos hay y a qué grupo pertenece cada muestra.

# A partir de los nombres de los ficheros de las muestras, se pueden crear dos nuevas columnas
# con la numeración de las muestras y el grupo al que pertenecen:
my.cels
# Erlot 1 -> tratamiento con Erlotinib con concentración 1 microMolar
# Erlot100 -> tratamiento con Erlotinib con concentración 100 nanoMolar
# NoTreat -> Sin tratamiento
# DMSO -> tratamiento con DMSO
pData(rawData)$sample.levels <- c(rep("Erlot1", 5), rep("Erlot100", 5), 
                                  rep("NoTreat", 4), rep("DMSO", 4))
pData(rawData)$sample.labels <- c(paste("Erlot1", 1:5, sep="."), paste("Erlot100", 1:5, sep="."), 
                                  paste("NoTreat", 1:4, sep="."), paste("DMSO", 1:4, sep="."))
head(pData(rawData))

# Obtener el resumen estadístico de las intensidades en cada muestra:
rawDaata <- rawData
colnames(exprs(rawDaata)) <- pData(rawData)$sample.labels
apply(exprs(rawDaata),2,summary)


### 2. Control de calidad de los datos crudos
arrayQualityMetrics(rawData, outdir = file.path("./results", "QCDir.Raw"), force=TRUE)

plotPCA3 <- function(datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
   data <- prcomp(t(datos),scale=scale)
   # plot adjustments
   dataDf <- data.frame(data$x)
   Group <- factor
   loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
   # main plot
   p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
     theme_classic() +
     geom_hline(yintercept = 0, color = "gray70") +
     geom_vline(xintercept = 0, color = "gray70") +
     geom_point(aes(color = Group), alpha = 0.55, size = 3) +
     coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
     scale_fill_discrete(name = "Group")
   # avoiding labels superposition
   p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
     labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
     ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
     theme(plot.title = element_text(hjust = 0.5)) +
     scale_color_manual(values=colores)
}

# PCA a partir de los datos crudos
plotPCA3(exprs(rawData), labels = pData(rawData)$sample.labels, pData(rawData)$sample.levels, 
          title="Datos sin transformar", scale = FALSE, size = 3, 
          colores = c("red", "blue", "green", "yellow"))
# boxplot de los datos crudos
boxplot(exprs(rawData),names=pData(rawData)$sample.labels, 
          col=c("red", "blue", "green", "yellow"), main="Datos sin transformar")
legend("topright", legend=c('Erlot1','Erlot100','NoTreat','DMSO'), fill=c("red", "blue", "green", "yellow"))
# density plot de los datos crudos
plotDensities(exprs(rawData), legend=F, col=c("red", "blue", "green", "yellow"), main="Datos sin transformar")
legend("topright", legend=c('Erlot1','Erlot100','NoTreat','DMSO'), fill=c("red", "blue", "green", "yellow"))


### 3. Normalización
eset_rma <- affy::rma(rawData)

### 4. [Control de calidad de los datos normalizados] (opcional)
arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"), force=TRUE)

# PCA a partir de las expresiones normalizadas
plotPCA3(exprs(eset_rma), labels = pData(eset_rma)$sample.labels, pData(eset_rma)$sample.levels, 
         title="Datos normalizados", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))
# boxplot de las expresiones normalizadas
boxplot(exprs(eset_rma),names=pData(eset_rma)$sample.labels, 
        col=c(rep("red",5), rep("blue",5), rep("green",4), rep("yellow",4)), main="Datos normalizados")
legend("topright", legend=c('Erlot1','Erlot100','NoTreat','DMSO'), fill=c("red", "blue", "green", "yellow"))
# density plot de las expresiones normalizadas
plotDensities(exprs(eset_rma), legend=F, col=c("red", "blue", "green", "yellow"), main="Datos normalizados")
legend("topright", legend=c('Erlot1','Erlot100','NoTreat','DMSO'), fill=c("red", "blue", "green", "yellow"))

### 5. Filtraje no específico [opcional]
filtered <- nsFilter(eset_rma, require.entrez=TRUE,
                     remove.dupEntrez=TRUE, var.func=IQR,
                     var.cutoff=0.5, var.filter=TRUE,
                     filterByQuantile=TRUE, feature.exclude="^AFFX")
names(filtered)
class(filtered$eset)
dim(exprs(filtered$eset))
print(filtered$filter.log)
eset_filtered <-filtered$eset
resultsDir <- "C:/Users/MARINA/Desktop/MásterUOC/Análisis de datos ómicos/PEC1/GSE106151/results"
save(eset_rma, eset_filtered, file=file.path(resultsDir, "eset_NormFilt.Rda"))

# PCA a partir de los datos filtrados
plotPCA3(exprs(eset_filtered), labels = pData(eset_filtered)$sample.labels, pData(eset_filtered)$sample.levels, 
         title="Datos filtrados", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))
# boxplot de los datos filtrados
boxplot(exprs(eset_filtered),names=pData(eset_filtered)$sample.labels, 
        col=c(rep("red",5), rep("blue",5), rep("green",4), rep("yellow",4)), main="Datos filtrados")
legend("topright", legend=c('Erlot1','Erlot100','NoTreat','DMSO'), fill=c("red", "blue", "green", "yellow"))
# density plot de los datos filtrados
plotDensities(exprs(eset_filtered), legend=F, col=c("red", "blue", "green", "yellow"), main="Datos filtrados")
legend("topright", legend=c('Erlot1','Erlot100','NoTreat','DMSO'), fill=c("red", "blue", "green", "yellow"))

### 6. Identificación de genes diferencialmente expresados
if (!exists("eset_filtered")) load (file="C:/Users/MARINA/Desktop/MásterUOC/Análisis de datos ómicos/PEC1/GSE106151/results/eset_NormFilt.Rda")

# A partir del paquete limma
# Hay que obtener la matriz de diseño de las muestras en primer lugar:
design <- model.matrix(~0 + sample.levels, pData(eset_filtered))
rownames(design) <- pData(eset_filtered)$sample.labels
colnames(design) <- levels(as.factor(pData(eset_filtered)$sample.levels))
design

#Obtener la matriz de contrastes 
cont.matrix <- makeContrasts(Er1 =(Erlot1-NoTreat), 
                             Er100=(Erlot100-NoTreat), 
                             DMSO=(DMSO-NoTreat), 
                             ErDif =(Erlot100-Erlot1),
                             levels=design)
cont.matrix
fit.main <- lmFit(eset_filtered, design)
fit.main <- contrasts.fit(fit.main,cont.matrix)
fit.main <- eBayes(fit.main)
class(fit.main)

# Listas DEG:
topTab_Er1 <- topTable(fit.main, number=nrow(fit.main), coef='Er1',adjust='fdr')
head(topTab_Er1)
dim(topTab_Er1[topTab_Er1$adj.P.Val<0.05,])[1]
dim(topTab_Er1[topTab_Er1$adj.P.Val<0.01,])[1]

topTab_Er100 <- topTable(fit.main, number=nrow(fit.main), coef='Er100',adjust='fdr')
head(topTab_Er100)
dim(topTab_Er100[topTab_Er100$adj.P.Val<0.05,])[1]
dim(topTab_Er100[topTab_Er100$adj.P.Val<0.01,])[1]

topTab_DMSO <- topTable(fit.main, number=nrow(fit.main), coef='DMSO',adjust='fdr')
head(topTab_DMSO)
dim(topTab_DMSO[topTab_DMSO$adj.P.Val<0.05,])[1]
dim(topTab_DMSO[topTab_DMSO$adj.P.Val<0.01,])[1]

topTab_ErDif <- topTable(fit.main, number=nrow(fit.main), coef='ErDif',adjust='fdr')
head(topTab_ErDif)
dim(topTab_ErDif[topTab_ErDif$adj.P.Val<0.05,])[1]
dim(topTab_ErDif[topTab_ErDif$adj.P.Val<0.01,])[1]

### 7. Anotación de los resultados

annotatedTopTable <- function(topTab, anotPackage)
 {
   topTab <- cbind(PROBEID=rownames(topTab), topTab)
   myProbes <- rownames(topTab)
   thePackage <- eval(parse(text = anotPackage))
   geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
   annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
 return(annotatedTopTab)
 }

topAnnotated_Er1 <- annotatedTopTable(topTab_Er1,
 anotPackage="hgu133plus2.db")
write.csv(topAnnotated_Er1, file="./results/topAnnotated_Er1.csv")

topAnnotated_Er100 <- annotatedTopTable(topTab_Er100,
 anotPackage="hgu133plus2.db")
write.csv(topAnnotated_Er100, file="./results/topAnnotated_Er100.csv")

topAnnotated_DMSO <- annotatedTopTable(topTab_DMSO,
 anotPackage="hgu133plus2.db")
write.csv(topAnnotated_DMSO, file="./results/topAnnotated_DMSO.csv")

topAnnotated_ErDif <- annotatedTopTable(topTab_ErDif,
 anotPackage="hgu133plus2.db")
write.csv(topAnnotated_ErDif, file="./results/topAnnotated_ErDif.csv")

# Visualizar expresión diferencial obtenida:
geneSymbols <- select(hgu133plus2.db, rownames(fit.main), c("SYMBOL")) 
SYMBOLS<- geneSymbols$SYMBOL

# 'Er1'
volcanoplot(fit.main, coef=1, highlight=4, names=SYMBOLS, 
  main=paste("Genes diferencialmente expresados", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))
abline(h=c(2))

# 'Er100'
volcanoplot(fit.main, coef=2, highlight=4, names=SYMBOLS, 
            main=paste("Genes diferencialmente expresados", colnames(cont.matrix)[2], sep="\n"))
abline(v=c(-1,1))
abline(h=c(2))

# 'DMSO'
volcanoplot(fit.main, coef=3, highlight=4, names=SYMBOLS, 
            main=paste("Genes diferencialmente expresados", colnames(cont.matrix)[3], sep="\n"))
abline(v=c(-1,1))
abline(h=c(2))

# 'ErDif'
volcanoplot(fit.main, coef=4, highlight=4, names=SYMBOLS, 
            main=paste("Genes diferencialmente expresados", colnames(cont.matrix)[4], sep="\n"))
abline(v=c(-1,1))
abline(h=c(2))

### 8. Comparación entre distintas comparaciones 
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.01, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,]
summary(res)
vennDiagram (res.selected[,1:3], cex=0.9)
vennDiagram (res.selected[,2:4], cex=0.9)
vennDiagram (res.selected[,c(1,3,4)], cex=0.9)
vennDiagram (res.selected[,c(1,2,4)], cex=0.9)

# se puede visualizar también gracias a un HeatMap:
probesInHeatmap <- rownames(res.selected)
HMdata <- exprs(eset_filtered)[rownames(exprs(eset_filtered)) %in% probesInHeatmap,]
geneSymbols <- select(hgu133plus2.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
write.csv(HMdata, file = file.path("./results/data4Heatmap.csv"))
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)

# HeatMap sin agrupar, ordenados según columnas
heatmap.2(HMdata,
        Rowv = FALSE,
        Colv = FALSE,
        main = "Genes diferencialmente expresados \n FDR < 0,1, logFC >=1",
        scale = "row",
        col = my_palette,
        sepcolor = "white",
        sepwidth = c(0.05,0.05),
        cexRow = 0.5,
        cexCol = 0.9,
        key = TRUE,
        keysize = 1.5,
        density.info = "histogram",
        ColSideColors = c(rep("red",5),rep("blue",5), rep("green",4), rep("yellow",4)),
        tracecol = NULL,
        dendrogram = "none",
        srtCol = 30)

# Heatmap agrupados por similaridad
heatmap.2(HMdata,
        Rowv = TRUE,
        Colv = TRUE,
        dendrogram = "both",
        main = "Genes diferencialmente expresados \n FDR < 0,1, logFC >=1",
        scale = "row",
        col = my_palette,
        sepcolor = "white",
        sepwidth = c(0.05,0.05),
        cexRow = 0.5,
        cexCol = 0.9,
        key = TRUE,
        keysize = 1.5,
        density.info = "histogram",
        ColSideColors = c(rep("red",5),rep("blue",5), rep("green",4), rep("yellow",4)),
        tracecol = NULL,
        srtCol = 30)

### 9. Análisis de significación biológica ("Gene Enrichment Analysis")
listOfTables <- list(Er1 = topTab_Er1, 
    Er100 = topTab_Er100, DMSO = topTab_DMSO,
    ErDif = topTab_ErDif)
listOfSelected <- list()

for (i in 1:length(listOfTables)){
   # select the toptable
   topTab <- listOfTables[[i]]
   # select the genes to be included in the analysis
   whichGenes<-topTab["adj.P.Val"]<0.01
   selectedIDs <- rownames(topTab)[whichGenes]
   # convert the ID to Entrez
   EntrezIDs<- select(hgu133plus2.db, selectedIDs, c("ENTREZID"))
   EntrezIDs <- EntrezIDs$ENTREZID
   listOfSelected[[i]] <- EntrezIDs
   names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)

mapped_genes2GO <- mappedkeys(org.Hs.egGO)
mapped_genes2KEGG <- mappedkeys(org.Hs.egPATH)
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)

# library(ReactomeAP)
listOfData <- listOfSelected[1:4]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)){
   genesIn <- listOfData[[i]]
   comparison <- comparisonsNames[i]
   enrich.result <- enrichPathway(gene = genesIn,
                                  pvalueCutoff = 0.05,
                                  readable = T,
                                  pAdjustMethod = "BH",
                                  organism = "human",
                                  universe = universe)
   
   cat("##################################")
   cat("\nComparison: ", comparison,"\n")
   print(head(enrich.result))
 
   if (length(rownames(enrich.result@result)) != 0) {
   write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
              row.names = FALSE)
   
   pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
     print(barplot(enrich.result, showCategory = 15, font.size = 4, 
             title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
   dev.off()
   
   pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
     print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
          vertex.label.cex = 0.75))
   dev.off()
   }
}
funcBiolEr1 <- read.csv('results/ReactomePA.Results.Er1.csv')
head(funcBiolEr1,n=1)
