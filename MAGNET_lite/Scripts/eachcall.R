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
