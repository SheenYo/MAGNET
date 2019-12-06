GOI1 <- GOI_1()
setwd(Home)
setwd("Outputfolder")
GL_all=list()
for(i in 1:29){
  GL_all[[i]]=Kang_genes$symbol[Kang_genes$Module==paste0("M",i)]
}
setwd(Home)

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
setwd("Outputfolder")

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

par(mar=c(5.1, 10, 4.1, 2))
bp=plot(x=betas, y=1:length(betas),
        type="n", panel.first = grid(ny=NA),
        yaxt = "n", ylab="",
        xlim=xlim,
        xlab=expression(paste(log(OR)%+-%95,"%CI")),
        main=paste(pheno, "associated Genes"))
abline(v=0,col="black", lty=3);
axis(2, at=1:30,cex.axis=0.7,
     labels=c("All Brain Genes", paste0("Module ",formatC(29:1, width=2, flag=0))),
     las=1);
arrows(x0=CILower, x1=CIUpper, y0=30:1, y1=30:1, col=rainbow(30), length=0, lwd=2,code = 3)
points(y=30:1, x=betas, pch=18, col="black")
betas[is.na(betas)]=0
text(y=(30:1)+0.5, x=betas, labels=pval, cex=1)
dev.off()
