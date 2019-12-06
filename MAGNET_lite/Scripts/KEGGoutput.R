pheno="Phenotype"
setwd(Home)
Outputfolder=setwd("Output/goelite_output")
FilesKEGG=list.files("../goelite_input/")
ZscoresKEGG=FilesKEGG[grep("GOElite.input", FilesKEGG)]
dir="../goelite_input/"
KEGGoutput = read.table(paste0(dir,ZscoresKEGG),sep="\t",header=T,fill=T, quote="")
#KEGG
dir="../goelite_denom/"
setwd(Home)
setwd("Output/goelite_denom")
Universal_Set<-read.table(paste0(dir,"GOElite.universe.txt"),sep="\t",header=TRUE)
Pheno<-kegga(KEGGoutput$SourceIdentifier,species = "Hs",universe=Universal_Set$SourceIdentifier)
topKEGG_Pheno<-topKEGG(Pheno,number=20L)

ColNames<-topKEGG_Pheno$Pathway[10:1]
par(mar=c(10,35,7,7))
bp=barplot(topKEGG_Pheno$P.DE[10:1],
           main=paste("GO Ontology", pheno),
           horiz=TRUE,
           yaxt='n',col="#fec44f",names.arg=ColNames,axisnames=TRUE,
           xlab="P-value",
           cex.main=2,cex.axis=1.5,cex.lab=1.5,tck=-0.01)
axis(2,at=bp,labels=topKEGG_Pheno$Pathway[10:1], tick=FALSE,las=1.58,cex.axis=1.5)
abline(v=0.05,col="red",lwd=2,lty=1)
dev.off()