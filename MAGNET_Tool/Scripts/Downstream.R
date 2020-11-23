cat("Checking if Rversion installed is R>=3.5")

RVer<-R.Version()
x<-RVer[13][1]
if(grepl("3.5|3.6|4.0",x[1])==FALSE) {
  print("Please use R-version>=3.5")
} else {
  print(paste0(x,"is installed",collapse=" "))
}


cat("Checking if all required packages are installed")
source("RefData/mergewithorder.R")


.libPaths("./tmpRlib")

ifelse(!dir.exists(file.path("./tmpRlib")), dir.create(file.path("./tmpRlib")), FALSE)


options(download.file.method="libcurl",useHTTPS=FALSE,BioC_mirror = "http://bioconductor.org",INSTALL_opts = c('--no-lock'),repos = c(CRAN = "http://cran.r-project.org", CRANextra = "http://cran.uk.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE))
	BiocManager::install(version = "3.12")
	#pkgs <- rownames(installed.packages())
	#BiocManager::install(pkgs, type = "source", checkBuilt = TRUE,ask=FALSE,lib="./tmpRlib")

package<-c("codetools","igraph","gprofiler2")

if(length(setdiff(package, rownames(installed.packages(lib.loc="./tmpRlib")))) > 0)	{
install.packages(setdiff(package, rownames(installed.packages(lib.loc="./tmpRlib"))),INSTALL_opts = c('--no-lock'),lib="./tmpRlib",ask=FALSE)  
                                                              		       }
require("BiocManager",lib.loc="./tmpRlib")


packages2_1<-c("annotate","AnnotationDbi","base64","Biobase", "bit64","cluster","crosstalk","devtools","digest","foreign")
if (length(setdiff(packages2_1, rownames(installed.packages(lib.loc="./tmpRlib")))) > 0) {
.libPaths("./tmpRlib")
BiocManager::install(setdiff(packages2_1, rownames(installed.packages(lib.loc="./tmpRlib"))),site_repository="http://bioconductor.org/packages/3.9/bioc",lib="./tmpRlib",ask=FALSE)
			                                                                  } 


packages2_2<-c("gtools","gplots","gtable","GO.db","Hmisc","illuminaio","impute","lattice","limma","Matrix","mime","mvtnorm","nlme","openssl",
"org.Hs.eg.db","plyr","preprocessCore","rmarkdown","rpart","rrcov","RSQLite","scales","shiny","survival","tables","viridisLite","WGCNA","XML")

if (length(setdiff(packages2_2, rownames(installed.packages(lib.loc="./tmpRlib")))) > 0) {
	.libPaths("./tmpRlib")
	BiocManager::install(setdiff(packages2_2, rownames(installed.packages(lib.loc="./tmpRlib"))),dependencies=TRUE,site_repository="http://bioconductor.org/packages/3.9/bioc",lib="./tmpRlib",ask=FALSE)
                                                                                          } 






##### install and load/load libraries


require("devtools",lib.loc="./tmpRlib")
require("gtable",lib.loc="./tmpRlib")
require("gprofiler2",lib.loc="./tmpRlib")
require("igraph",lib.loc="./tmpRlib")
install_github("ethanbahl/cerebroViz")
require("WGCNA",lib.loc="./tmpRlib")
require("annotate",lib.loc="./tmpRlib")
require("org.Hs.eg.db",lib.loc="./tmpRlib")
require("limma",lib.loc="./tmpRlib")
require("cerebroViz",lib.loc="./tmpRlib")
require("gplots",lib.loc="./tmpRlib")
require("plyr",lib.loc="./tmpRlib")
require("scales",lib.loc="./tmpRlib")
require("cerebroViz",lib.loc="./tmpRlib")
require("gtools",lib.loc="./tmpRlib")
require(org.Hs.eg.db,lib.loc="./tmpRlib")
require("viridisLite",lib.loc="./tmpRlib")
# #############################################
 
options(stringsAsFactors = F)
 
args=commandArgs(trailingOnly = T)
 
if(!length(args)==6){
  cat("downstream analysis did not have enough arguments supplied setting to standard\n")
  Pathwaysfolder="./Output"
  outputfolder="OutputDownstream"
  pheno="Phenotype"
  MAGNETHome="../MAGNET"
  Genelist=c("Genelist")
  UniversalGenes=c("UniversalGenes")
                    } else{
  Pathwaysfolder=args[1]
  Outputfolder=args[2]
  pheno=args[3]
  MAGNETHome=args[4]
  Genelist=args[5]
  UniversalGenes=args[6]
                          }

 
Home=getwd()

#Generating Biological pathway Plots
 
#GO plots

Genelist<-read.table(args[5],sep="\t",header = TRUE)
Genes<-Genelist[,1]

setwd(Home)

UniversalGenes<-read.table(args[6],sep = "\t",header=TRUE)

UG<-UniversalGenes[,1]


setwd(Pathwaysfolder)

#1) Gprofiler

gostres <- gost(query = as.vector(Genes), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = strtoi(UG), 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

pdf(paste0(pheno,"_GProfiler.pdf"),height=4,width=5);
p<-gostplot(gostres, capped = FALSE, interactive = FALSE)
publish_gostplot(p,width = NA, height = NA, filename = NULL )
dev.off()

#Table<-publish_gosttable(gostres)
#write.table(Table,paste0("Table_",pheno,".txt"),sep="\t",row.names=FALSE,quote=FALSE)

#2) GO Plots
GOoutput = gostres$result
if(identical(GOoutput$term_id,logical(0))==FALSE) {
  pdf(paste0(pheno,"_GO_Plot.pdf"),height=4,width=5);
  par(oma=c(3,6,3,6));
  par(mgp=c(1, 0, 0))
  GO = GOoutput[(GOoutput$source=="GO:BP"|GOoutput$source=="GO:BP"|GOoutput$source=="GO:CC"),];
  if(length(GO$term_id)>0){
  bp = barplot(GO[order(GO$p_value,decreasing = TRUE),][,3],
               main=paste("GO Ontology"),
               horiz=TRUE,
               yaxt='n',col="#D95F02",
               xlab="p_value",
               cex.main=0.5,cex.axis=0.3,cex.lab=0.5,tck=-0.05);
  axis(2,at=bp,labels=GO[order(GO$p_value,decreasing = TRUE),][,11],
       tick=FALSE,las=2,cex.axis=0.3);
  abline(v=0.05,col="red",lwd=2,lty=1);
  dev.off()
  }else {
    cat("Not enough enriched terms for GO plot")
  }
  
#3) KEGG plot   
  pdf(paste0(pheno,"_KEGGPlot.pdf"),height=4,width=5);
  par(oma=c(3,6,3,6));
  par(mgp=c(1, 0, 0))
  kegg = GOoutput[(GOoutput$source=="KEGG"),];
  if(length(kegg$term_id)>0){
  KG = barplot(kegg[order(kegg$p_value,decreasing = TRUE),][,3],
               main=paste("KEGG pathways"),
               horiz=TRUE,
               yaxt='n',col="#D95F02",
               xlab="p_value",
               cex.main=0.5,cex.axis=0.3,cex.lab=0.5,tck=-0.05);
  axis(2,at=KG,labels=kegg[order(kegg$p_value,decreasing = TRUE),][,11],
       tick=FALSE,las=2,cex.axis=0.3);
  abline(v=0.05,col="red",lwd=2,lty=1);
  dev.off()
  } else {
    cat("Not enough enriched terms for KEGG plot")
  }
}

 
# ##############################################################################
# ########################## TRANSCRIPTOME ANALYSIS ############################
# ##############################################################################
# 
setwd(Home)
 
#Enrichment in Kang Dataset
 
setwd(MAGNETHome)
# 
#1)All Kang data set genes
KangUnivers<- read.table("RefData/KangUnivers.txt", sep="\t", header=T) 
colnames(KangUnivers)<-c("EntrezId","Symbol") #Name the columns
 
#2)All Kang genes that are in Modules
Kang_genes<-read.table("RefData/Kang_dataset_genesMod_version2.txt",sep="\t",header=TRUE)
 
#3)Generate Gene universe to be used for single gene lists
tmp=merge(KangUnivers,Kang_genes,by.y="EntrezGene",by.x="EntrezId",all=TRUE) #18826
KangUni_Final<-tmp[duplicated(tmp$EntrezId)==FALSE,] #18675
 
#4)Read the Gene list of interest
#magmagenes_significant<-read.table(Genelist,sep="\t",header=TRUE)
magmagenes_significant<-Genelist


#hs <- org.Hs.eg.db
#my.symbols <- read.table("Genelist.txt")
#select(hs, 
#       keys = as.character(my.symbols$V1),
#       columns = c("ENTREZID", "SYMBOL"),
#       keytype = "SYMBOL")
 
#5)Convert the gene list to gene symbols
 
All_genes_Entrez<-data.frame(magmagenes_significant[,1])
colnames(All_genes_Entrez)<-"EntrezId"


Annotation_list<-as.data.frame(matrix(NA,nrow=length(All_genes_Entrez$EntrezId),ncol=2))#Create an empty matrix
 
colnames(Annotation_list)<-c("EntrezId","Symbol")
 
Annotation_list$EntrezId<-All_genes_Entrez$EntrezId
 
#.libPaths("./tmpRlib")
for(i in 1:nrow(Annotation_list)){
	.libPaths("./tmpRlib")
	Annotation_list$Symbol[i] <- unlist(mget(as.character(All_genes_Entrez$EntrezId[i]),org.Hs.egSYMBOL,ifnotfound=NA))#Extract gene symbols of the entrez ids
				}


setwd(Home)

setwd(Outputfolder)


Annotation_list$Module = Kang_genes$Module[match(Annotation_list$EntrezId,Kang_genes$EntrezGene)] #Attach module information
 
GOI_Pheno<-Annotation_list[duplicated(Annotation_list$EntrezId)==FALSE,]

Sys.sleep(1)
 
#Create list of genes of interest
 
GOI1<-list(GOI_Pheno$Symbol)

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
 
Sys.sleep(2)
setwd(Home)

setwd(Outputfolder)
 
pdf(paste0(pheno,"Enrichmentplot_Kang.pdf"))
betas=df$beta
SEs=df$SE
idx=is.infinite(betas)
betas[idx]=NaN
SEs[idx]=NaN
CIUpper=betas+1.96*SEs
CILower=betas-1.96*SEs
pval=convertpval(df$P)
pval[df$Padj<0.05]="adj.p**"
pval[is.na(betas)]="n.a."
xlim=range(c(CIUpper, CILower), na.rm=T)*1.2
par(mar=c(5,12,5,3))
plot(x=betas, y=1:length(betas),
	type="n", panel.first = grid(ny=NA),
	yaxt = "n", ylab="",
	xlim=xlim,
	xlab=expression(paste('log(OR)'%+-%95,"%CI")),
	main=paste(pheno, "associated Genes"))
	abline(v=0,col="black", lty=3)
	axis(2, at=1:30,
	labels=c("All brain genes", paste0("Module ",formatC(29:1, width=2, flag=0))),
	las=1)
	arrows(x0=CILower, x1=CIUpper, y0=30:1, y1=30:1, col=rainbow(30), length=0, lwd=2,code = 3)
	points(y=30:1, x=betas, pch=18, col="black")
	betas[is.na(betas)]=0
	text(y=(30:1)+0.5, x=betas, labels=pval, cex=0.7)
	dev.off()
Sys.sleep(2)
 
# ##############################################################################
# #
# #		Networks
# ##############################################################################
Sys.sleep(1)
 
setwd(Home)
setwd(MAGNETHome)
load("RefData/DataPreprocessing.RData") #Load the Kang expression data of all genes 
datExprPlot=matriz #Expression data of Kang loaded as Rdata object DataPreprocessing.RData
setwd(Home)
setwd(Outputfolder)
 
Sys.sleep(1)
 
KangModules_v1<-Kang_genes[is.na(Kang_genes$Module)==FALSE,] #Selecting dataset without NA modules
 
AllGenes_Kang<-KangModules_v1[which(KangModules_v1$symbol %in% GOI_Pheno$Symbol),] #GOI which are in KANG
 
AllGenes_Kang_Uniq<-AllGenes_Kang[duplicated(AllGenes_Kang$symbol)==FALSE,]
 
AllGenes<-data.frame(AllGenes_Kang_Uniq$symbol)
colnames(AllGenes)<-"V1"

 
Sys.sleep(1)
.libPaths("./tmpRlib")
printMOIgraph=function(Module){
  pdf("Graphs_of_Enriched_Modules.pdf") 
	for(Module in MOI){
	GenesInMod<-KangModules_v1[which(gsub("M", "Module ",KangModules_v1$Module) %in% Module),][,2] #Respective Entrez gene ids in the enriched module
	GenesInModEntrez=KangModules_v1[which(gsub("M", "Module ",KangModules_v1$Module) %in% Module),][,10] ###Respective gene Symbols in the enriched module
  num_nodes<- min(c(50, length(GenesInMod)))
  adjMat1 = bicor(datExprPlot[,which(colnames(datExprPlot) %in% as.character(GenesInMod))])
  diag(adjMat1)=NA
  index=order(colSums(adjMat1, na.rm = T), decreasing = T)[1:num_nodes]
  reducedTOM=adjMat1[index, index]

  ## define layouts
  n=ceiling(num_nodes/6)
  g0 <- graph.adjacency(as.matrix(reducedTOM[1:n,1:n]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)
  g0 <- graph.adjacency(as.matrix(reducedTOM[(n+1):(3*n),(n+1):(3*n)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)
  g0 <- graph.adjacency(as.matrix(reducedTOM[((3*n)+1):ncol(reducedTOM),((3*n)+1):ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatc <- layout.circle(g0)
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*1,layoutMatb*2.5, layoutMatc*4); ##Plot 3 circles
     
  ## define the colors
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
  plot(g1, 
  vertex.color=tmpcol, 
	vertex.size=5+6*(V(g1)$name %in% as.character(AllGenes$V1))*1,
	vertex.label=as.character(colnames(reducedTOM)),
	vertex.label.cex=0.7+0.3*(V(g1)$name %in% AllGenes$V1),
	vertex.label.dist=0.4+0.3*(V(g1)$name %in% AllGenes$V1)*1,
	vertex.label.degree=90,
	vertex.label.color="black",
	vertex.label.font=1, 
	vertex.label.family="sans",
	edge.width=2,
	layout= layoutMat, 
	main=Module)
	DFinMod=KangUnivers[(KangUnivers$Symbol %in% GenesInMod) &
	(KangUnivers$Symbol %in% AllGenes_Kang_Uniq$symbol),]
	write.table(DFinMod, file=paste0(pheno,"_Genes_in_",Module,".txt"), sep="\t", row.names=F)
			   }
  dev.off()
                                }
 
 
 
MOI=df$Module[df$P<=0.05 & df$Module!="Module all"]
 
MS<-as.numeric(gsub("\\D", "", MOI)) #list of significant modules
 
if(length(MOI)<1){
  cat(" No Significant Enrichment in Brain Expressed Modules detected \n Program exits")
  #quit("no")	} else
 	              } else{
   			                for(Module in MOI){
    			              printMOIgraph(MOI)
  					                              }
 			                  }
 
 
# #################################################################End of Enrichment plots############################################################################################
# 
cerebroScale2<-
  function (x, clamp, divData, center_zero)
  {
    xmed <- median(x, na.rm = TRUE)
    xmad <- mad(x, constant = 1, na.rm = TRUE)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    avoidClamp <- max(abs(xmed - xmin), abs(xmed - xmax))/xmad
    fillMatrix <- x
    if (is.null(clamp)) {
      clamp <- avoidClamp + 1
    }
    outlrs <- clamp * xmad
    if (clamp <= 0)
      stop("clamp must be >0")
    pctOL <- round(length(which(x[!is.na(x)] <= (xmed - (outlrs)) |
                                  x[!is.na(x)] >= (xmed + (outlrs))))/length(x[!is.na(x)]) *
                     100, 2)
    if (pctOL > 0) {
      warning(paste("The clamp value of ", clamp, " will clamp ",
                    pctOL, "% of input values (outliers) to the min or max of the scaled range.",
                    sep = ""))
    }
    if (divData == TRUE & center_zero == FALSE) {
      abvMed <- x[x >= xmed & x <= (xmed + outlrs) & !is.na(x)]
      belMed <- x[x <= xmed & x >= (xmed - outlrs) & !is.na(x)]
      if (length(which(!is.na(x)))%%2 == 0) {
        rightsc <- rescale(c(xmed, abvMed), c(0.5, 1))[-1]
        fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(c(xmed, belMed), c(-1, 0))[-1]
        fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
      if ((length(which(!is.na(x))))%%2 == 1) {
        rightsc <- rescale(abvMed, c(-1, 1))
        fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(belMed, c(-1, 0))
        fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
    }
    if (divData == TRUE & center_zero == TRUE) {
      abvMed <- x[x >= -1 & x <= (-1 + outlrs) & !is.na(x)]
      belMed <- x[x <= -1 & x >= (-1 - outlrs) & !is.na(x)]
      if (length(which(!is.na(x)))%%2 == 0) {
        rightsc <- rescale(c(-1, abvMed), c(-1, 1))[-1]
        fillMatrix[x >= -1 & x <= (-1 + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(c(-1, belMed), c(-1, 0))[-1]
        fillMatrix[x <= -1 & x >= (-1 - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (-1 - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (-1 + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
      if ((length(which(!is.na(x))))%%2 == 1) {
        rightsc <- rescale(abvMed, c(-1, 1))
        fillMatrix[x >= -1 & x <= (-1 + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(belMed, c(-1, 1))
        fillMatrix[x <= -1 & x >= (-1 - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (-1 - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (-1 + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
    }
    if (divData == FALSE) {
      nonoutlrs <- x[x >= (xmed - outlrs) & x <= (xmed + outlrs) &
                       !is.na(x)]
      xsc <- rescale(nonoutlrs, c(-1, 1))
      fillMatrix[x >= (xmed - outlrs) & x <= (xmed + outlrs) &
                   !is.na(x)] <- xsc
      fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
      fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
      xScaled <- fillMatrix
    }
    return(xScaled)
  }

##########################################################################################################################################################################################################

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


  pdf(paste0("Module",j,"Heatmap.pdf"), paper="a4r")
  toplot<-toplot[,1:998]
  date<-c(1:998)
  dateY<-paste0(round(date/365,2),"_Years")
  names(toplot)<-dateY
  cerebroScale2(as.matrix(toplot),clamp=NULL,divData = FALSE,center_zero=FALSE) -> ex1_scaled
  heatmap.2(as.matrix(toplot), col = cols,
            trace = "none",
            na.color = "grey",
            Colv = F, Rowv = F,labCol = labvec,
            breaks = seq(-0.1,0.1, length.out=101),
            symkey =F,
            key.title = "",dendrogram = "none",
            key.xlab = "eigengene",
            density.info = "none",
            main=paste("Module",j),
            srtCol=90,tracecol = "none",
            add.expr=eval.parent(abline(h=0, v=282),axis(1,at=c(1,77,282,392,640,803,996),labels =FALSE)),cexCol = 1)
  dev.off()
  dir.create(paste0("BrainView_Module_",j))
  setwd(paste0("BrainView_Module_",j))
  cerebroViz(ex1_scaled,palette = cols,customNames = cnm,figLabel = TRUE,timePoint = c(1,77,282,392,640,803,996),filePrefix = paste0("Module_Brain",j),divData = FALSE)
  setwd("../")

}



