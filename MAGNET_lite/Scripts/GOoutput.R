GOoutput = GOoutput[(GOoutput$Ontology.Type=="biological_process" |
                       GOoutput$Ontology.Type=="molecular_function"),]
ColNames<-GOoutput$Ontology.Name[10:1]
par(mgp=c(3,0,0),mar=c(10,17,10,5))
#par(mar = rep(2, 4))
bp=barplot(GOoutput$Z.Score[10:1],
           main=paste("GO Ontology", pheno),
           horiz=TRUE,
           yaxt='n',col="#D95F02",
           xlab="Z-score",names.arg=ColNames,axisnames=TRUE,
           cex.main=1,cex.axis=0.7,cex.lab=0.7,tck=-0.01)
axis(2,at=bp,labels=GOoutput$Ontology.Name[10:1],
     tick=FALSE,las=2,cex.axis=0.7);
abline(v=1.96,col="red",lwd=2,lty=1)
dev.off()