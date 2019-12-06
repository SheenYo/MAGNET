MS<-as.numeric(gsub("\\D", "", Module)) #list of significant modules

GL_all=list()

for(i in 1:29){
  GL_all[[i]]=Kang_genes$symbol[Kang_genes$Module==paste0("M",i)]
}

Genes<-list()
for(i in 1:29)
  Genes[[i]]<-Kang_genes$symbol[Kang_genes$Module==paste0("M",i)]

#GL_all is same as Genes here

Genes_expression<-list()
pcatest<-list()
for (i in 1:29){
  Genes_expression[[i]]<-matriz[,which(colnames(matriz) %in% Genes[[i]])]
  pcatest=prcomp(t(as.matrix(Genes_expression[[i]])),retx=TRUE)
}


#PCA test
PCA<-data.frame(pcatest$rotation)
PCA$donor_name<-rownames(PCA)
PC1<-data.frame(PCA[,c(1,6)])

#Combining the age with expression data
list <- strsplit(sampleInfo$age, " ")
library("plyr")
df <- ldply(list)
colnames(df) <- c("Age", "time")

sampleInfo<-cbind(sampleInfo[,1:9],df)
sampleInfo$Age<-as.numeric(sampleInfo$Age)

sampleInfo$period<-ifelse(sampleInfo$time=="pcw",sampleInfo$Age*7,ifelse(sampleInfo$time=="yrs",sampleInfo$Age*365+270,ifelse(sampleInfo$time=="mos",sampleInfo$Age*30+270,NA)))

#We need it just for the donor names
source("/home/afsheenyousaf/Desktop/ShinyApp_EnrichmentTesting/MAGNET/RefData/mergewithorder.R",local=TRUE)
PCA_matrix<-merge.with.order(PC1,sampleInfo,by.y="SampleID",by.x="donor_name",keep_order=1)

#Select which have phenotype info present 
matriz2<-matriz[which(rownames(matriz) %in% PCA_matrix$donor_name),]
FactorGenes_expression<-list()
#Factors here mean modules
for (i in 1:29){
  FactorGenes_expression[[i]]<-matriz2[,which(colnames(matriz2) %in% Genes[[i]])]
}


FactorseGE<-list()
for (i in 1:29){
  FactorseGE[[i]]<-FactorGenes_expression[[i]]
}

allModgenes=NULL
colors=vector()
for ( i in 1:29){
  allModgenes=cbind(allModgenes,FactorseGE[[i]])
  colors=c(colors, rep(i, ncol(FactorseGE[[i]])))
}

lengths=unlist(lapply(FactorGenes_expression, ncol), use.names = F)

MEorig=moduleEigengenes(allModgenes, colors)


PCA_matrixfreeze=PCA_matrix

index=!PCA_matrix$structure_acronym %in% c("URL", "DTH", "CGE","LGE", "MGE",  "Ocx", "PCx", "M1C-S1C","DIE", "TCx", "CB")
PCA_matrix=PCA_matrix[index,]
ME = MEorig$eigengenes[index,]
timepoints=seq(56,15000, length.out=1000)
matrix(c("CB", "THA", "CBC", "MD"), ncol=2 ) -> cnm
for (j in MS) {
  MEmod=ME[,j]
  toplot=data.frame(matrix(NA, nrow=length(table(PCA_matrix$structure_acronym)), ncol=998))
  rownames(toplot)=unique(PCA_matrix$structure_acronym)
  target <- c("OFC", "DFC", "VFC", "MFC","M1C","S1C","IPC","A1C","STC","ITC","V1C","HIP","AMY","STR","MD","CBC")
  toplot<-toplot[c(6,2,8,5,11,12,10,9,7,4,14,3,1,13,16,15),]
  
  
  for ( i in unique(PCA_matrix$structure_acronym)){
    index=PCA_matrix$structure_acronym==i
    LOESS=loess(MEmod[index]~PCA_matrix$period[index])
    toplot[i,]=predict(LOESS,newdata = round(exp(seq(log(56),log(15000), length.out=998)),2))
    colnames(toplot)[c(1,77,282,392,640,803,996)]<- c("1pcw","21pcw","Birth","1.3years","5.4years","13.6years","40.7years")
  }
  
  
  cols=colorRampPalette(c(  "dodgerblue3","white","darkorange"))(100)
  labvec <- c(rep(NA, 1000))
  
  labvec[c(1,77,282,392,640,803,996)] <- c("1pcw","21pcw","Birth","1.3years","5.4years","13.6years","40.7years")
  
  
  #pdf(paste0("Module",1,"Heatmap.pdf"), paper="a4r")
  toplot<-toplot[,1:998]
  date<-c(1:998)
  dateY<-paste0(round(date/365,2),"_Years")
  
  names(toplot)<-dateY
  
  
  cerebroScale2(as.matrix(toplot),clamp=NULL,divData = FALSE,center_zero=FALSE) -> ex1_scaled
  par(xpd=FALSE) 
  heatmap.2(as.matrix(toplot), col = cols, 
            trace = "none", 
            na.color = "grey",
            Colv = F, Rowv = F,labCol = labvec,
            breaks = seq(-0.1,0.1, length.out=101),
            symkey =F,
            key.title = "",dendrogram = "none",
            key.xlab = "eigengene",
            density.info = "none",
            #main=paste("Module",1),
            srtCol=90,tracecol = "none", cexRow = 1,
            add.expr=eval.parent(abline(v=282),axis(1,at=c(1,77,282,392,640,803,996),labels =FALSE)),cexCol = 1)
  #print(Heaty)
  #dev.off()
}
