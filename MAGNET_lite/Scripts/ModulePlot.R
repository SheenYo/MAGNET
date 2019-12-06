pheno="Phenotype"
GOI1 <- GOI_1()
GOI_Pheno <- GOI_P()

GL_all=list()
for(i in 1:29){
  GL_all[[i]]=Kang_genes$symbol[Kang_genes$Module==paste0("M",i)]
}
Sys.sleep(1)
names(GL_all)=paste0("Module ",1:29)

GL_all[["Module all"]]=Kang_genes$symbol[!Kang_genes$Module==""]

Sys.sleep(1)
UniversalGeneset=KangUni_Final$Symbol 
pval<-list()
Testframe<-data.frame()
Resultsall=list()
Resultsall[["beta"]]=data.frame(matrix(NA, length(GL_all), length(GOI1)))
Resultsall[["SE"]]=data.frame(matrix(NA, length(GL_all), length(GOI1)))
Resultsall[["Pval"]]=data.frame(matrix(NA, length(GL_all), length(GOI1)))
Resultsall[["OR"]]=data.frame(matrix(NA, length(GL_all), length(GOI1)))
Resultsall[["ORL"]]=data.frame(matrix(NA, length(GL_all), length(GOI1)))
Resultsall[["ORU"]]=data.frame(matrix(NA, length(GL_all), length(GOI1)))

for(i in 1:length(GL_all)){
  for(j in 1:length(GOI1)){
    Modulegene=GL_all[[i]]
    Factorgene=GOI1[[j]]
    Testframe<-fisher.test(table(factor(UniversalGeneset %in% Factorgene,levels=c("TRUE","FALSE")),
                                 factor(UniversalGeneset %in% Modulegene,levels=c("TRUE","FALSE"))))
    beta=log(Testframe$estimate)
    Resultsall[["beta"]][i,j]=beta
    Resultsall[["SE"]][i,j]=abs(beta-log(Testframe$conf.int[1]))/1.96
    Resultsall[["Pval"]][i,j]=Testframe$p.value
    Resultsall[["OR"]][i,j]=(Testframe$estimate)
    Resultsall[["ORL"]][i,j]=(Testframe$conf.int[1])
    Resultsall[["ORU"]][i,j]=(Testframe$conf.int[2])
  }
}
Sys.sleep(1)

df <- data.frame(matrix(unlist(Resultsall), nrow=30, byrow=F))
colnames(df)<-c("beta","SE","P","OR","ORL","ORU")
df$Padj=p.adjust(df$"P", method="bonferroni")
df$Module=c(paste0("Module ", 1:29), "Module all")
Sys.sleep(2)
write.table(df,file=paste0(pheno,"Enrichment_output.xls"),sep="\t",quote=FALSE,row.names=FALSE)


convertpval=function(x){
  sapply(x, function(x){ifelse(x<=0.01, "**", ifelse(x<=0.05, "*", ""))})
}

setwd(Home)


load("RefData/DataPreprocessing.RData") #Load the Kang expression data of all genes 
datExprPlot=matriz #Expression data of Kang loaded as Rdata object DataPreprocessing.RData

setwd(Home)
setwd("RefData")

Sys.sleep(1)

KangModules_v1<-Kang_genes[is.na(Kang_genes$Module)==FALSE,] #Selecting dataset without NA modules

AllGenes_Kang<-KangModules_v1[which(KangModules_v1$symbol %in% GOI_Pheno$Symbol),] #GOI which are in KANG

AllGenes_Kang_Uniq<-AllGenes_Kang[duplicated(AllGenes_Kang$symbol)==FALSE,]

AllGenes<-data.frame(AllGenes_Kang_Uniq$symbol)
colnames(AllGenes)<-"V1"

Module=df$Module[df$Module=="Module 1"]
printMOIgraph=function(Module){
  
  GenesInMod<-KangModules_v1[which(gsub("M", "Module ",KangModules_v1$Module) %in% Module),][,2] #Respective Entrez gene ids in the enriched module
  GenesInModEntrez=KangModules_v1[which(gsub("M", "Module ",KangModules_v1$Module) %in% Module),][,10] ###Respective gene Symbols in the enriched module
  num_nodes<- min(c(50, length(GenesInMod)))
  adjMat1 = bicor(datExprPlot[,which(colnames(datExprPlot) %in% as.character(GenesInMod))])
  diag(adjMat1)=NA
  index=order(colSums(adjMat1, na.rm = T), decreasing = T)[1:num_nodes]
  reducedTOM=adjMat1[index, index]
  
  
  ### define layouts
  n=ceiling(num_nodes/6)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[1:n,1:n]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[(n+1):(3*n),(n+1):(3*n)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[((3*n)+1):ncol(reducedTOM),((3*n)+1):ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatc <- layout.circle(g0)
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*1,layoutMatb*2.5, layoutMatc*4); ##Plot 3 circles
  
  
  ### define the colors
  tmpcol<-c()  
  
  for(i in 1:length(V(g1)$name)){
    Gene<-V(g1)$name[i]
    Faccol=ifelse(Gene %in% AllGenes$V1,"darkorange3","grey")
    tmpcol[i]<-Faccol
  }
  
  
  colormap = colorRampPalette(c("red","white","blue"))(201);
  
  E(g1)$weight[is.na(E(g1)$weight)]=0
  roundedweight = round(E(g1)$weight*100,0)
  E(g1)$color = paste0(colormap[roundedweight+100], as.hexmode(floor(abs(roundedweight/100*255))))
  
  plot(g1,cex.lab=10,
       vertex.color=tmpcol, 
       vertex.size=5+6*(V(g1)$name %in% as.character(AllGenes$V1))*1,
       vertex.label=as.character(colnames(reducedTOM)),
       #vertex.label.cex=5+5*(V(g1)$name %in% AllGenes$V1),
       vertex.label.cex=5,
       vertex.label.dist=0.4+0.3*(V(g1)$name %in% AllGenes$V1)*1,
       vertex.label.degree=90,
       vertex.label.color="black",
       vertex.label.font=2, 
       vertex.label.family="sans",
       edge.width=2,
       layout= layoutMat, 
       main=Module)
  DFinMod=KangUnivers[(KangUnivers$Symbol %in% GenesInMod) &
                        (KangUnivers$Symbol %in% AllGenes_Kang_Uniq$symbol),]
  write.table(DFinMod, file=paste0(pheno,"_Genes_in_",Module,".txt"), sep="\t", row.names=F)
  #        }
}
printMOIgraph(Module)