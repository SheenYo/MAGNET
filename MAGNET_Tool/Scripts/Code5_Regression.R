
##### install and laod or just load libraries

packagesBioC=c()
packagesCRAN=c("lme4", "readxl")

pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    { dir.create("./tmpRlib", showWarnings=F)
      install.packages(x,dep=TRUE, lib="./tmpRlib", repos="http://cloud.r-project.org")
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

pkgTestBioC <- function(x)
  {
    if (!require(x,character.only = TRUE))
    { dir.create("./tmpRlib", showWarnings=F)
        source("https://bioconductor.org/biocLite.R")
        biocLite(x, lib="./tmpRlib", dependencies=T, update="n")
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }


if(length(packagesCRAN)>0) { sapply(packagesCRAN, pkgTest) }
if(length(packagesBioC)>0) { sapply(packagesBioC, pkgTestBioC) }



options(echo=TRUE, stringsAsFactors = F) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)

print (args)

if (length(args)!=6){
  cat("not all options are provided \n")
  cat("arg 1= Rawdatafile 2=Phenofile 3=ColNamePhenoToTest 4=Covariates 5=fixedEffects[none] 6=Outputname\n")
  stop()
  
}
source("../../RefData/mergewithorder.R")

cat("processing", args[1], "\n\n")

##########################################################
##Read the file from PLINK after recode AD
###########################################################

gc(TRUE)

rawdata<-read.table(args[1],sep=" ", header=T)

dim(rawdata)
names(rawdata)
Phenofile<-read.table(args[2],sep="\t",header=TRUE) #Regression gathered file
dim(Phenofile)
head(Phenofile)

Geno_Pheno<-merge(Phenofile,rawdata,by.x ="ID_Genetik", by.y="IID")
dim(Geno_Pheno)


Geno_Pheno_Phenotype<-subset(Geno_Pheno,Geno_Pheno$PHENOTYPE==-9) 

dim(Geno_Pheno_Phenotype)

#Select only one individual per family

Complete_PhenoGenoData<-Geno_Pheno_Phenotype[duplicated(Geno_Pheno_Phenotype$FID.x)==FALSE,] 


dim(Complete_PhenoGenoData)
snps<-Complete_PhenoGenoData[,colnames(rawdata)[-c(1:6)]]

colnames(snps)<- gsub('_[A-z ]*', '' , colnames(snps))
colnames(snps)<- gsub('_[0-9 ]*', '' , colnames(snps))

colnames(snps)<-sub("X", "", colnames(snps))

#args[3]="JA"

y<-Complete_PhenoGenoData[,args[3]]

#####################################
##Extract SNPs B-coefficients and fits
#####################################

pvars <- c("C1","C2",
           "C3","C4")
datsc <- Complete_PhenoGenoData
datsc[pvars] <- lapply(datsc[pvars],scale)

options(lmerControl=list(check.nobs.vs.rankZ = "warningSmall",
check.nobs.vs.nlev = "warning",
check.nobs.vs.nRE = "warning",
check.nlev.gtreq.5 = "ignore",
check.nlev.gtr.1 = "warning"))



########################################################
### run regression choose y (your phenotype of interest)
########################################################
gc(TRUE)
## Regress over SNPS in apply function

#args[4]="Age,Sex"
#args[5]="Site"


covars=unlist(strsplit(args[4], ","))

for(cov in covars ) {
  if(  length(unique(datsc[,cov]))<4 | is.character(datsc[1,cov]) ) {
    cat(cov, "is treated as factor")
    datsc[,cov]=as.factor(datsc[,cov])
  }
  
}



if(args[5]!="none"){
  cat(" doing the fixed effect models\n")
  FixEff=args[5]
  datsc[,FixEff]=as.factor(datsc[,FixEff])
  cat("Converted Fixed effect to factor")
  formula=as.formula(paste0("y ~ x+", paste(covars, collapse="+"), "+C1+C2+C3+C4+ (1|",FixEff,")"))
  snps=rbind(as.data.frame(snps), 1:ncol(snps))
  Res=apply(snps[,], 2, function(x)	{
    cat(x[length(x)],"\r")
    datsc$x=x[-length(x)]
    fit<- try(lmer(formula ,data=datsc, REML=FALSE,na.action=na.exclude,control=lmerControl(optimizer="Nelder_Mead")))
    fit
  }
  )
  
} else {
  cat(" doing the linear models")
  formula=as.formula(paste0("y ~ x+", paste(covars, collapse="+"), "+C1+C2+C3+C4"))
  Res=apply(snps[,], 2, function(x)
  {
    datsc$x=x
    fit<- try(lm(formula ,data=datsc, na.action=na.exclude))
    fit
  }
  )
}


## extract coefficients 
Comp_Res<-list()
m <-as.data.frame(matrix(NA,nrow=ncol(snps),ncol=4))
names(m)<-c("snp","beta", "se", "t")
m[[1]]<- colnames(snps)	  
#get a summary and the coefficients out of the model 
for (col in 1:length(Res)) {
  smry <- summary(Res[[col]])
  Comp_Res[[length(Comp_Res)+1]] = smry
  coefs <- smry$coefficients
  # Populate the results data frame with the appropriate values
  diditrun <- (nrow(coefs)>=2 & row.names(coefs)[2]=="x")
  m$beta[col] <- ifelse(diditrun, round(coefs["x",1],4), NA)
  m$se[col] <- ifelse(diditrun, round(coefs["x",2],4), NA)
  m$t[col] <- ifelse(diditrun, round(coefs["x",3],4), NA)
}


## output 
p <- pnorm(abs(m$t), lower.tail = FALSE) * 2
ctable <- cbind(m, "p_value" = p)
ctable$adj.pvalues<-p.adjust(ctable$p,method="fdr",n=length(p))
## output 
names(Comp_Res)<-colnames(snps)
#dir.create("Output",showWarnings = F)
#setwd("Output")
capture.output(Comp_Res, file = args[6])
write.table(ctable,paste0("Results",args[6]),sep="\t",row.names=FALSE)

cat("processed", args[1], "\n\n")
