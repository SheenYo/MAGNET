# Generating Manhattan plots for each Chromosome

#Read the linear regression file

# needs pheno and color from inputline

args = commandArgs(trailingOnly=TRUE)


if(!length(args)==2){
    cat("arguments for Manhattanplots not set, setting to standard\n")
    args=c("Phenotype", "#FF0000")}

##########
pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    { dir.create("./tmpRlib", showWarnings=F)
      install.packages(x,dep=TRUE, lib="./tmpRlib", repos="http://cloud.r-project.org")
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }




pkgTest("data.table")
pkgTest("qqman")
pkgTest("lattice")

library("data.table")
library("lattice")
library("qqman")

options(stringsAsFactors=FALSE)

LR_output<-fread(paste0("Complete_SNP_info_",args[1]), header=T, data.table=F)
colnames(LR_output)<-c("SNP","beta","se","t","P","adj.pvalues","CHR","CHR2","GP","BP","A1","A2")

GWsnps_pheno<-LR_output[which(LR_output$P<0.00000005),] 
Nomsnps_pheno<-LR_output[which(LR_output$P<0.05),]
write.table(GWsnps_pheno[,c(1:7,10:12)],paste0("GWsnps",args[1],".txt"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(Nomsnps_pheno[,c(1:7,10:12)],paste0("Nomsnps",args[1],".txt"),sep="\t",row.names=FALSE,quote=FALSE)


LR_output<-LR_output[which(LR_output$CHR!="NA"),]





snpsOfInterest<-c("")


png(paste0(args[1], "_Manhattan_Final.png"), width=1250,height=702)
par(mar=c(5, 8,3,2))
manhattan(x=LR_output, col=c("grey50", "#8ac51b"), cex.axis=1.5, cex.lab=1.5, cex.main=2,main=args[1],genomewideline=FALSE,suggestiveline=FALSE,ylim=c(0,10),highlight=snpsOfInterest)
abline(h=2, lwd=3, col="blue") # 
abline(h=7.30103, lwd=3, col="red") 
dev.off()




#### function for nice QQ plot

qqunif.plot<-function(pvalues, 
	should.thin=T, thin.obs.places=2, thin.exp.places=2, 
	xlab=list(expression(paste("Expected ",-log[10], "(",italic("p"),")")),cex=1.5),
	ylab=list(expression(paste("Observed ",-log[10], "(",italic("p"),")")),cex=1.5), 
	draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
	already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
	par.settings=list(superpose.symbol=list(pch=pch)), ...) {
	
	
	#error checking
	if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
	if(!(class(pvalues)=="numeric" || 
		(class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
		stop("pvalue vector is not numeric, can't draw plot")
	if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
	if (already.transformed==FALSE) {
		if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
	} else {
		if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
	}
	
	
	grp<-NULL
	n<-1
	exp.x<-c()
	if(is.list(pvalues)) {
		nn<-sapply(pvalues, length)
		rs<-cumsum(nn)
		re<-rs-nn+1
		n<-min(nn)
		if (!is.null(names(pvalues))) {
			grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
			names(pvalues)<-NULL
		} else {
			grp=factor(rep(1:length(pvalues), nn))
		}
		pvo<-pvalues
		pvalues<-numeric(sum(nn))
		exp.x<-numeric(sum(nn))
		for(i in 1:length(pvo)) {
			if (!already.transformed) {
				pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
				exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
			} else {
				pvalues[rs[i]:re[i]] <- pvo[[i]]
				exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
			}
		}
	} else {
		n <- length(pvalues)+1
		if (!already.transformed) {
			exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
			pvalues <- -log10(pvalues)
		} else {
			exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
		}
	}


	#this is a helper function to draw the confidence interval
	panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
		require(grid)
		conf.points = min(conf.points, n-1);
		mpts<-matrix(nrow=conf.points*2, ncol=2)
        	for(i in seq(from=1, to=conf.points)) {
            		mpts[i,1]<- -log10((i-.5)/n)
            		mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
            		mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
            		mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        	}
        	grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
    	}

	#reduce number of points to plot
	if (should.thin==T) {
		if (!is.null(grp)) {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places),
				grp=grp))
			grp = thin$grp
		} else {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places)))
		}
		pvalues <- thin$pvalues
		exp.x <- thin$exp.x
	}
	gc()
	
	prepanel.qqunif= function(x,y,...) {
		A = list()
		A$xlim = range(x, y)*1.02
		A$xlim[1]=0
		A$ylim = A$xlim
		return(A)
	}

	#draw the plot
	xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
               prepanel=prepanel,
               scales=list(axs="i",cex=1.5),
               pch=pch,
               panel = function(x, y, ...) {
			if (draw.conf) {
				panel.qqconf(n, conf.points=conf.points, 
					conf.col=conf.col, conf.alpha=conf.alpha)
			};
			panel.xyplot(x,y, ...);
			panel.abline(0,1);
		}, par.settings=par.settings, ...
	)
}




png(paste0(args[1], "_QQ_Final.png"), width=702,height=702)
par(mar=c(5, 8,3,2))
qqunif.plot(LR_output$P, main=list(label=args[1], cex=2))
dev.off()



