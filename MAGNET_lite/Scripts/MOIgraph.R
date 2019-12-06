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
  
  plot(g1, 
       vertex.color=tmpcol, 
       vertex.size=5+6*(V(g1)$name %in% as.character(AllGenes$V1))*1,
       vertex.label=as.character(colnames(reducedTOM)),
       #vertex.label.cex=1.5+0.7*(V(g1)$name %in% AllGenes$V1),
       vertex.label.cex=0.8+0.05*(V(g1)$name %in% AllGenes$V1),
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