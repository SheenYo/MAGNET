options(stringsAsFactors = F)

args=commandArgs(trailingOnly = T)

options(scipen = 999) #Turning off scientific notation only for this block

args=c("QC1_report.lmiss","QC1_report.imiss","PreQC_hardy.hwe","PreQC_AlleleFreq.frq","QC2_Sexcheck.sexcheck","PreQC_Inbreeding.het", "PostQC_report.lmiss","PostQC_report.imiss","PostQC_Hardyreport.hwe","PostQC_AlleleFreq.frq","Final_Sexcheck.sexcheck","PostQC_Inbreeding.het",
       "SampleInds.txt","../../RefData/Hapmap_siteinfo.txt","HM3mds.mds")

cat("...Plotting the data before QC...\n")

cat("...Plotting missingness per SNP...\n")

if(file.exists(args[1])){
lmiss<-read.table(args[1],header=T) #Read in the lmiss file


pdf("PreQC_MissignessReport_SNP.pdf",width = 10,height=10)
hist(lmiss$F_MISS,col="blue",cex.lab=1.3,breaks = seq(min(lmiss$F_MISS), max(lmiss$F_MISS), length.out = 100),ylab="SNPs",xlab="Fraction of missing data",main="Pre-QC_SNP")
abline(v=0.05,lty=2)
dev.off()

} else {
cat("Missingness Report file per SNP was not generated")
}



cat("...Plotting missingness per individual...")

if(file.exists(args[2])){
imiss=read.table(args[2],head=T) #Read in the imiss file
pdf("PreQC_MissignessReport_Individual.pdf",width = 10,height=10)
hist(imiss$F_MISS,col="green",cex.lab=1.3,breaks = seq(min(imiss$F_MISS), max(imiss$F_MISS), length.out = 100),ylab="Samples",xlab="Fraction of genotype call missing",main="Pre-QC_Individuals")
abline(v=0.05,lty=2)
dev.off()
} else {
cat("Missingness Report file per individual was not generated")
}


cat("...Plotting Hardy Weinberg distribution...")


if(file.exists(args[3])){
hardy=read.table(args[3],head=T)
pdf("PreQC_HWE.pdf",width = 10,height=10)
hist(hardy$P,xlab="P-value from HWE Test",main="Pre-QC_HWE",col="red",cex.lab=1.3)
abline(v=0.00000001,lty=2)
dev.off()
} else {
cat("HWE file was not generated")
}

cat("...Plotting minor allele frequency...")

if(file.exists(args[4])){
freq=read.table(args[4],header=T)
pdf("PreMAF.pdf")
hist(freq$MAF,n=100,xlab="MAF",main="Pre-QC_MAF",col="yellow")
abline(v=0.02,lty=2)
dev.off()
} else {
cat("minor allele frequency file was not generated")
}

cat("...Infer Sex of Individuals from SNP Genotypes...")

if(file.exists(args[5])){
sex=read.table(args[5],header=TRUE,stringsAsFactors=FALSE)
a<-0
b<-1.48
c<-1.6
pdf("PreQC_sexcheck.pdf")
plot(sex[,6],col=as.numeric(sex[,4])+1,pch=19,
     ylab="Estimated Heterozygosity",xlab="Individual",main="Gender check")
abline(h=0.8,lty=2)
text(mean(sex[,6])*2.5,(b*max(sex[,6])/2+a),"FEMALES",srt=0.2,pos=3,cex=0.5,col="blue")
text(mean(sex[,6])*3,(c*max(sex[,6])/2+a),"MALES",srt=0.2,pos=3,cex=0.5,col="blue")
dev.off()
} else {
cat("PreQC_sexcheck file was not generated")
}

cat("...Plotting Pre QC Inbreeding output...")

if(file.exists(args[6])){
PreInbred<-read.table(args[6],header=T)
pdf("PreQC_Inbreeding.pdf")
hist(PreInbred$F,xlab="Inbreeding Coefficient",col="orange",main="Pre-QC Histogram Inbreeding Co-efficent")
abline(v=c(-0.15,0.25),lty=2)
dev.off()
} else {
cat("Pre QC Inbreeding file was not generated")
}


pdf("Main Pre-QC plot.pdf")
par(mfrow=c(3,2))
if(exists("lmiss")==TRUE){
hist(lmiss$F_MISS,col="blue",cex.lab=1.3,breaks = seq(min(lmiss$F_MISS), max(lmiss$F_MISS), length.out = 100),ylab="SNPs",xlab="Fraction of missing data",main="Pre-QC_SNP")
abline(v=0.05,lty=2)
} else if (exists("imiss")==TRUE) {
hist(imiss$F_MISS,col="green",cex.lab=1.3,breaks = seq(min(imiss$F_MISS), max(imiss$F_MISS), length.out = 100),ylab="Samples",xlab="Fraction of genotype call missing",main="Pre-QC_Individuals")
abline(v=0.05,lty=2)
} else if (exists("hardy")==TRUE) {
hist(hardy$P,xlab="P-value from HWE Test",main="Pre-QC_HWE",col="red",cex.lab=1.3)
abline(v=0.00000001,lty=2)
} else if (exists("freq")==TRUE) {
hist(freq$MAF,n=100,xlab="MAF",main="Pre-QC_MAF",col="yellow")
abline(v=0.02,lty=2)
} else if (exists("PreInbred")==TRUE) {
hist(PreInbred$F,xlab="Inbreeding Coefficient",col="orange",main="Pre-QC Inbreeding Co-efficent")
abline(v=c(-0.15,0.25),lty=2)
dev.off()
} else {print("Pre QC files do not exists")}






#########################################################################################################################################################################################
cat("...Plotting the data after QC...")
#########################################################################################################################################################################################


cat("...Plotting missingness per SNP...")

PostQClmiss<-read.table(args[7],header=T) #Read in the lmiss file
pdf("PostQC_MissignessReport_SNP.pdf",width = 10,height=10)
hist(PostQClmiss$F_MISS,col="blue",cex.lab=1.3,breaks = seq(min(PostQClmiss$F_MISS), max(PostQClmiss$F_MISS), length.out = 100),ylab="SNPs",xlab="Fraction of missing data",main="PostQC_SNP")
abline(v=0.05,lty=2)
dev.off()


cat("...Plotting missingness per individual...")

PostQCimiss=read.table(args[8],head=T) #Read in the imiss file
pdf("PostQC_MissignessReport_Individual.pdf",width = 10,height=10)
hist(PostQCimiss$F_MISS,col="green",cex.lab=1.3,breaks = seq(min(PostQCimiss$F_MISS), max(PostQCimiss$F_MISS), length.out = 100),ylab="Samples",xlab="Fraction of genotype call missing",main="PostQC_Individuals")
abline(v=0.05,lty=2)
dev.off()



cat("...Plotting Hardy Weinberg distribution...")

Posthardy=read.table(args[9],head=T)
png("PostQC_HWE.png")
#hist(hardy$P,xlab="P-value from HWE Test",main="PreQC_HWE",col="red",cex.lab=1.3)
plot(-log10(Posthardy[,8]), -log10(Posthardy[,7]), asp=1,xlab="-log10 Expected Heterozygosity",ylab="-log10 Observed Heterozygosity",main="PostQC_hardy",col="blue")
abline(0,1)
abline(v=0.00000001,lty=2)
dev.off()


cat("...Plotting minor allele frequency...")


Postfreq=read.table(args[10],header=T)
pdf("PostQC_MAF.pdf")
hist(Postfreq$MAF,xlab="MAF",main="Pre-QC_MAF",col="yellow")
dev.off()


cat("...Infer Sex of Individuals from SNP Genotypes...")

sex=read.table(args[11],header=TRUE,stringsAsFactors=FALSE)
a<-0
b<-1.48
c<-1.6
pdf("PreQC_sexcheck.pdf")
plot(sex[,6],col=as.numeric(sex[,4])+1,pch=19,
     ylab="Estimated Heterozygosity",xlab="Individual")
abline(h=0.8,lty=2)
text(mean(sex[,6])*2.5,(b*max(sex[,6])/2+a),"FEMALES",srt=0.2,pos=3,cex=0.5,col="blue")
text(mean(sex[,6])*3,(c*max(sex[,6])/2+a),"MALES",srt=0.2,pos=3,cex=0.5,col="blue")
dev.off()



cat("....Plotting post QC Inbreeding output...")

PostInbred<-read.table(args[12],header=T)
pdf("PostQC_Inbreeding.pdf")
hist(PostInbred$F,xlab="Inbreeding Coefficient",col="orange",main="Histogram Inbreeding Co-efficent")
abline(v=c(-0.15,0.25),lty=2)
dev.off()




cat("...Filter individuals by heterozygosity and missingness...")

missf=as.numeric(PostQCimiss[,6])
hetf=1-as.numeric(PostInbred[,3])/as.numeric(PostInbred[,5])

pdf("Inds by heterozygosity and missingness")
par(mfrow=c(1,2))
hist(missf,xlab="Missingness",col="blue");hist(hetf,xlab="Heterozygosity",col="orange")
plot(missf,hetf,xlab="Missingness",ylab="Heterozygosity")
abline(v=.007,col=2);abline(h=c(.365,.28),col=2)
dev.off()




pdf("Main Post-QC plot.pdf")
par(mfrow=c(3,2))

hist(PostQClmiss$F_MISS,col="blue",cex.lab=1.3,breaks = seq(min(PostQClmiss$F_MISS), max(PostQClmiss$F_MISS), length.out = 100),ylab="SNPs",xlab="Fraction of missing data",main="PostQC_SNP")
abline(v=0.05,lty=2)

hist(PostQCimiss$F_MISS,col="green",cex.lab=1.3,breaks = seq(min(PostQCimiss$F_MISS), max(PostQCimiss$F_MISS), length.out = 100),ylab="Samples",xlab="Fraction of genotype call missing",main="PostQC_Individuals")
abline(v=0.05,lty=2)

plot(-log10(Posthardy[,8]), -log10(Posthardy[,7]), asp=1,xlab="-log10 Expected Heterozygosity",ylab="-log10 Observed Heterozygosity",main="PostQC_hardy",col="blue")
abline(0,1)
abline(v=0.00000001,lty=2)

hist(Postfreq$MAF,xlab="MAF",main="Pre-QC_MAF",col="yellow")

hist(PostInbred$F,xlab="Inbreeding Coefficient",col="orange",main="Histogram Inbreeding Co-efficent")
abline(v=0.05,lty=2)

dev.off()




                                                                        ###########################################
                                                                        #               MDS-plot                  #
                                                                        ###########################################

StudyFID<-read.table(args[13],sep="\t")
colnames(StudyFID)=c("FID","IID")
StudyFID$Population<-"NA"
  

HMFID<-read.table(args[14],sep=" ",header=TRUE)

DE_HM<-rbind(StudyFID,HMFID)
#DE_HM_PopInfo<-merge(DE_HM,pop_info,by.x="V1",by.y="FID")

mds.cluster<-read.table(args[15],header=TRUE)

FINAL=merge(DE_HM,mds.cluster,by="FID") 

pdf("Study_HM_samples_Group.pdf",height=8,width=8)
plot(x=FINAL$C1,y=FINAL$C2,xlab="MDS Componenet1",ylab="MDS Component2",type ="n")
points(x=FINAL$C1[FINAL$Population=="NA"],y=FINAL$C2[FINAL$Population=="NA"],col="blue",pch=16)
points(x=FINAL$C1[FINAL$Population=="CEU"],y=FINAL$C2[FINAL$Population=="CEU"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="YRI"],y=FINAL$C2[FINAL$Population=="YRI"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="JPT"],y=FINAL$C2[FINAL$Population=="JPT"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="CHB"],y=FINAL$C2[FINAL$Population=="CHB"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="ASW"],y=FINAL$C2[FINAL$Population=="ASW"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="CHD"],y=FINAL$C2[FINAL$Population=="CHD"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="GIH"],y=FINAL$C2[FINAL$Population=="GIH"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="LWK"],y=FINAL$C2[FINAL$Population=="LWK"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="MEX"],y=FINAL$C2[FINAL$Population=="MEX"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="MKK"],y=FINAL$C2[FINAL$Population=="MKK"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="TSI"],y=FINAL$C2[FINAL$Population=="TSI"],col="red",pch=17)
points(x=FINAL$C1[FINAL$Population=="YRI"],y=FINAL$C2[FINAL$Population=="YRI"],col="red",pch=17)
legend("bottomright",legend=c("StudyData","Hapmap"),
       col=c("blue",'red')
       ,pch=c(16,17))
dev.off()  

cat("...MDS plots based on HM population and Study data...")

pdf("Study_HM_samples_Population.pdf",height = 8,width = 8)
plot(x=FINAL$C1,y=FINAL$C2,xlab="MDS Component1",ylab="MDS Component2",type="n")
points(x=FINAL$C1[FINAL$Population=="NA"],y=FINAL$C2[FINAL$Population=="NA"],col="blue",pch=17)
points(x=FINAL$C1[FINAL$Population=="CEU"],y=FINAL$C2[FINAL$Population=="CEU"],col="burlywood4",pch=16)
points(x=FINAL$C1[FINAL$Population=="YRI"],y=FINAL$C2[FINAL$Population=="YRI"],col="red",pch=16)
points(x=FINAL$C1[FINAL$Population=="JPT"],y=FINAL$C2[FINAL$Population=="JPT"],col="purple",pch=16)
points(x=FINAL$C1[FINAL$Population=="CHB"],y=FINAL$C2[FINAL$Population=="CHB"],col="green",pch=16)
points(x=FINAL$C1[FINAL$Population=="ASW"],y=FINAL$C2[FINAL$Population=="ASW"],col="violet",pch=16)
points(x=FINAL$C1[FINAL$Population=="CHD"],y=FINAL$C2[FINAL$Population=="CHD"],col="aquamarine1",pch=16)
points(x=FINAL$C1[FINAL$Population=="GIH"],y=FINAL$C2[FINAL$Population=="GIH"],col="bisque",pch=16)
points(x=FINAL$C1[FINAL$Population=="LWK"],y=FINAL$C2[FINAL$Population=="LWK"],col="orange",pch=16)
points(x=FINAL$C1[FINAL$Population=="MEX"],y=FINAL$C2[FINAL$Population=="MEX"],col="coral1",pch=16)
points(x=FINAL$C1[FINAL$Population=="MKK"],y=FINAL$C2[FINAL$Population=="MKK"],col="gray",pch=16)
points(x=FINAL$C1[FINAL$Population=="TSI"],y=FINAL$C2[FINAL$Population=="TSI"],col="deeppink",pch=16)
points(x=FINAL$C1[FINAL$Population=="YRI"],y=FINAL$C2[FINAL$Population=="YRI"],col="forestgreen",pch=16)
legend("bottomright",legend=c("Study Data","CEU","YRI","JPT","CHB","ASW","CHD","GIH","LWK","MEX","MKK","TSI","YRI"),
       col=c("blue","burlywood4","red","purple","green","violet","aquamarine1","bisque","orange","coral1","gray","deeppink","forestgreen"),pch=c(17,16,16,16,16,16,16,16,16,16,16,16,16))
dev.off()   

cat("...Study Data only...")

pdf("Study_Data Only.pdf",height = 8,width = 8)
plot(x=FINAL$C1,y=FINAL$C2,xlab="MDS Componenet1",ylab="MDS Component2",type ="n")
points(x=FINAL$C1[FINAL$Population=="NA"],y=FINAL$C2[FINAL$Population=="NA"],col="blue",pch=16)
legend("bottomright",legend="Study Data",col="blue",pch=16)
dev.off()
