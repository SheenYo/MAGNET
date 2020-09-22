cat("Checking if Rversion installed is R>=3.5")

RVer<-R.Version()
x<-RVer[13][1]
if(grepl("3.5|3.6",x[1])==FALSE) {# make an or with 3.4, 3.3, 3.2 and 3.1
  print("Please use R-version>=3.5")
} else {
  print(paste0(x,"is installed",collapse=" "))
}

Home=getwd()
cat("Checking if all required packages are installed")

#SourceRmerge<-paste(Home,"/RefData/mergewithorder.R",sep="")

source(paste(Home,'/RefData/mergewithorder.R',sep=""), local = TRUE)

.libPaths(.libPaths())


if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
     #pkgs <- rownames(installed.packages())
     #BiocManager::install(pkgs, type = "source", checkBuilt = TRUE,ask=FALSE)
 package<-c("igraph","DT")
     if(length(setdiff(package, rownames(installed.packages()))) > 0)	{
     install.packages(setdiff(package, rownames(installed.packages())),ask=FALSE)  
                                                                      }
 
packages2_2<-c("shiny","methods","WGCNA","annotate","org.Hs.eg.db","limma","scales","gplots")
if (length(setdiff(packages2_2, rownames(installed.packages()))) > 0) {
   BiocManager::install(setdiff(packages2_2, rownames(installed.packages())),site_repository="http://bioconductor.org/packages/3.9/bioc",ask=FALSE,dependencies=TRUE)
                                                                       } 
 


# Call the required libraries
require("shiny")
require("DT")
require("igraph")
require("WGCNA")
require("annotate")
require("org.Hs.eg.db")
require("limma")
require("scales")
require("gplots")

Home=getwd()

dir.create("Output", showWarnings = FALSE)
setwd("Output")

dir.create("goelite_input", showWarnings = FALSE)
dir.create("goelite_denom", showWarnings = FALSE)
dir.create("goelite_output", showWarnings = FALSE)

setwd(Home)

# identify the folders
current.folder <- paste(Home,"/Input",sep="")
new.folder <- paste(Home,"/Output/goelite_input",sep="")

setwd(current.folder)
# find the files that you want
list.of.files <- list.files(current.folder, "\\.txt$")

# copy the files to the new folder
file.copy(list.of.files,new.folder,recursive = FALSE)

setwd(Home)
setwd("Output/goelite_input")
file.rename(list.files(pattern="*.txt"), "GOElite.input.txt")
#file.append(list.files(pattern="*.txt"), "GOElite.input.txt")
setwd(Home)




ui <- fluidPage(
  
  # App title ----
  titlePanel("MAGNET"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select the random distribution type ----
      fileInput('file', 'Upload File',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )),
      
      # Taken from: http://shiny.rstudio.com/gallery/file-upload.html
      tags$hr(),
      checkboxInput('header', 'Header', FALSE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   '\t'),
      
      ################################################################
      
      actionButton("choice", "upload")
      
      ,width=4,height=2),
    
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      selectInput("columns", "Select Columns", choices = NULL), # no choices before uploading
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("MyGenelist",DT::dataTableOutput("table_display")),
                  #tabPanel("GoElitePlot", fluidRow(column(4,offset = 5,align="center",plotOutput(width = "1000px", height="500px","plot"), downloadButton(outputId = "down", label = "Download the plot")))),
                  tabPanel("GoElitePlot",plotOutput("plot"), downloadButton(outputId = "Ontology", label = "Download")),
                  tabPanel("KEGG",plotOutput("plot2"),downloadButton(outputId = "KEGGenrichment", label = "Download")),
                  tabPanel("Enrichment", plotOutput("plot3"), downloadButton(outputId = "ModuleEnrichment", label = "Download")),
                  tabPanel("Brainnetwork",tabsetPanel(type = "tabs",
                                                      tabPanel("MOI_1", plotOutput("MOI1"),downloadButton(outputId = "Module1", label = "Download")),
                                                      tabPanel("MOI_2", plotOutput("MOI2"),downloadButton(outputId = "Module2", label = "Download")),
                                                      tabPanel("MOI_3", plotOutput("MOI3"),downloadButton(outputId = "Module3", label = "Download")),
                                                      tabPanel("MOI_4", plotOutput("MOI4"),downloadButton(outputId = "Module4", label = "Download")),
                                                      tabPanel("MOI_5", plotOutput("MOI5"),downloadButton(outputId = "Module5", label = "Download")),
                                                      tabPanel("MOI_6", plotOutput("MOI6"),downloadButton(outputId = "Module6", label = "Download")),
                                                      tabPanel("MOI_7", plotOutput("MOI7"),downloadButton(outputId = "Module7", label = "Download")),
                                                      tabPanel("MOI_8", plotOutput("MOI8"),downloadButton(outputId = "Module8", label = "Download")),
                                                      tabPanel("MOI_9", plotOutput("MOI9"),downloadButton(outputId = "Module9", label = "Download")),
                                                      tabPanel("MOI_10", plotOutput("MOI10"),downloadButton(outputId = "Module10", label = "Download")),
                                                      tabPanel("MOI_11", plotOutput("MOI11"),downloadButton(outputId = "Module11", label = "Download")),
                                                      tabPanel("MOI_12", plotOutput("MOI12"),downloadButton(outputId = "Module12", label = "Download")),
                                                      tabPanel("MOI_13", plotOutput("MOI13"),downloadButton(outputId = "Module13", label = "Download")),
                                                      tabPanel("MOI_14", plotOutput("MOI14"),downloadButton(outputId = "Module14", label = "Download")),
                                                      tabPanel("MOI_15", plotOutput("MOI15"),downloadButton(outputId = "Module15", label = "Download")),
                                                      tabPanel("MOI_16", plotOutput("MOI16"),downloadButton(outputId = "Module16", label = "Download")),
                                                      tabPanel("MOI_17", plotOutput("MOI17"),downloadButton(outputId = "Module17", label = "Download")),
                                                      tabPanel("MOI_18", plotOutput("MOI18"),downloadButton(outputId = "Module18", label = "Download")),
                                                      tabPanel("MOI_19", plotOutput("MOI19"),downloadButton(outputId = "Module19", label = "Download")),
                                                      tabPanel("MOI_20", plotOutput("MOI20"),downloadButton(outputId = "Module20", label = "Download")),
                                                      tabPanel("MOI_21", plotOutput("MOI21"),downloadButton(outputId = "Module21", label = "Download")),
                                                      tabPanel("MOI_22", plotOutput("MOI22"),downloadButton(outputId = "Module22", label = "Download")),
                                                      tabPanel("MOI_23", plotOutput("MOI23"),downloadButton(outputId = "Module23", label = "Download")),
                                                      tabPanel("MOI_24", plotOutput("MOI24"),downloadButton(outputId = "Module24", label = "Download")),
                                                      tabPanel("MOI_25", plotOutput("MOI25"),downloadButton(outputId = "Module25", label = "Download")),
                                                      tabPanel("MOI_26", plotOutput("MOI26"),downloadButton(outputId = "Module26", label = "Download")),
                                                      tabPanel("MOI_27", plotOutput("MOI27"),downloadButton(outputId = "Module27", label = "Download")),
                                                      tabPanel("MOI_28", plotOutput("MOI28"),downloadButton(outputId = "Module28", label = "Download")),
                                                      tabPanel("MOI_29", plotOutput("MOI29"),downloadButton(outputId = "Module29", label = "Download")))),
                  tabPanel("Heatmaps",tabsetPanel(type = "tabs",
                                                  tabPanel("Heatmap_1", plotOutput("Heatmap1"), downloadButton(outputId = "Module1Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_2", plotOutput("Heatmap2"), downloadButton(outputId = "Module2Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_3", plotOutput("Heatmap3"), downloadButton(outputId = "Module3Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_4", plotOutput("Heatmap4"), downloadButton(outputId = "Module4Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_5", plotOutput("Heatmap5"), downloadButton(outputId = "Module5Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_6", plotOutput("Heatmap6"), downloadButton(outputId = "Module6Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_7", plotOutput("Heatmap7"), downloadButton(outputId = "Module7Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_8", plotOutput("Heatmap8"), downloadButton(outputId = "Module8Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_9", plotOutput("Heatmap9"), downloadButton(outputId = "Module9Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_10", plotOutput("Heatmap10"), downloadButton(outputId = "Module10Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_11", plotOutput("Heatmap11"), downloadButton(outputId = "Module11Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_12", plotOutput("Heatmap12"), downloadButton(outputId = "Module12Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_13", plotOutput("Heatmap13"), downloadButton(outputId = "Module13Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_14", plotOutput("Heatmap14"), downloadButton(outputId = "Module14Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_15", plotOutput("Heatmap15"), downloadButton(outputId = "Module15Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_16", plotOutput("Heatmap16"), downloadButton(outputId = "Module16Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_17", plotOutput("Heatmap17"), downloadButton(outputId = "Module17Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_18", plotOutput("Heatmap18"), downloadButton(outputId = "Module18Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_19", plotOutput("Heatmap19"), downloadButton(outputId = "Module19Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_20", plotOutput("Heatmap20"), downloadButton(outputId = "Module20Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_21", plotOutput("Heatmap21"), downloadButton(outputId = "Module21Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_22", plotOutput("Heatmap22"), downloadButton(outputId = "Module22Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_23", plotOutput("Heatmap23"), downloadButton(outputId = "Module23Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_24", plotOutput("Heatmap24"), downloadButton(outputId = "Module24Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_25", plotOutput("Heatmap25"), downloadButton(outputId = "Module25Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_26", plotOutput("Heatmap26"), downloadButton(outputId = "Module26Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_27", plotOutput("Heatmap27"), downloadButton(outputId = "Module27Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_28", plotOutput("Heatmap28"), downloadButton(outputId = "Module28Heatmap", label = "Download")),
                                                  tabPanel("Heatmap_29", plotOutput("Heatmap29"), downloadButton(outputId = "Module29Heatmap", label = "Download"))))
      )
    )
  ))
server <- function(input, output, session) { # added session for updateSelectInput
  

  
  info <- eventReactive(input$choice, {
    inFile <- input$file
    # Instead # if (is.null(inFile)) ... use "req"
    req(inFile)
    
    # Changes in read.table 
    f <- read.table(inFile$datapath, header = input$header, sep = input$sep, quote = input$quote)
    vars <- names(f)
    # Update select input immediately after clicking on the action button. 
    updateSelectInput(session, "columns","Select Columns", choices = vars)
    f
  })
  
  output$table_display <- DT::renderDataTable({
    f <- info()
    f <- subset(f, select = input$columns) #subsetting takes place here
    g <- rep("L", length(f$V1))
    GOElite.input <- data.frame(f, g)
    colnames(GOElite.input)<- c("SourceIdentifier", "SourceCode")
    ## dir.create("goelite_input", showWarnings = F)
    setwd(Home)
    setwd("Output/goelite_input")
    write.table(GOElite.input, "GOElite.input.txt",quote=F, sep="\t",row.names = F)
    colnames(f) <- "Gen"
    DT::datatable(f, options = list(pageLength = 15, info = FALSE,
                                    lengthMenu = list(c(15, -1), c("15", "All"))))
  })
  
  
  ################################################################    Create Universal geneset    #################################################################      
  if(file.exists("Output/goelite_input/GOElite.input.txt")==TRUE & file.exists("Output/goelite_denom/GOElite.universe.txt")==FALSE) {
  setwd(Home)
  setwd("RefData")
  EG <- read.table("KangUnivers.txt", sep="\t", header=T)
  tmp2 <- rep("L", length(EG$Entrez_id))
  UniversalGene <- data.frame(EG$Entrez_id, tmp2)
  colnames(UniversalGene)<- c("SourceIdentifier", "Sourcecode")
  setwd(Home)
  ## dir.create("goelite_denom", showWarnings = F)
  setwd("Output/goelite_denom")
  write.table(UniversalGene, "GOElite.universe.txt",quote=F, sep="\t",row.names = F)
                                                                  }
  

  
  
  ## Sys.sleep(3)
  setwd(Home)
  if(file.exists("Output/goelite_output/GO-Elite_results/CompleteResults/ORA_pruned/GOElite.input-GO_z-score_elite.txt")==FALSE){  
    setwd(Home)
    if(file.exists("Output/goelite_input/GOElite.input.txt")==TRUE) {
      command="python"
      setwd(Home)
      Home= getwd()
      url<-"https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/go-elite/GO-Elite_v.1.2.5-Py.zip"
      download.file(url, destfile = basename(url))
      unzip("GO-Elite_v.1.2.5-Py.zip", exdir=Home)
      OutputF=paste0(Home,"/Output/goelite_output/")
      InputF=paste0(Home,"/Output/goelite_input/")
      DenomF=paste0(Home,"/Output/goelite_denom/")
      Species="Hs"
      args1=c(paste("--update Official",collapse=" "),paste("--species Hs",collapse=" "),paste("--version EnsMart62Plus",collapse=" "))
      args=c(paste("--input",InputF,collapse=" "),paste("--denom",DenomF,collapse=" "),paste("--output",OutputF,collpase=" "),paste("--species",Species,collapse=" "))
      #allArgs1=c("/home/afsheenyousaf/GO-Elite_v.1.2.5-Py/GO_Elite.py",args1)
      #allArgs=c("/home/afsheenyousaf/GO-Elite_v.1.2.5-Py/GO_Elite.py",args)
      Path=paste0(Home,"/GO-Elite_v.1.2.5-Py/GO_Elite.py")
      allArgs1=c(Path,args1)
      allArgs=c(Path,args)
      
      
      withProgress(message = 'Updating Species in GOelite', value = 0, {
        # Number of times we'll go through the loop
        n <- 3
        
        for (i in 1:n) {
          # Each time through the loop, add another row of data. This is
          # a stand-in for a long-running computation.
          RunPy1=system2(command,args = allArgs1, stdout=TRUE)
          
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("Doing part", i))
          
          # Pause for 0.1 seconds to simulate a long computation.
          #Sys.sleep(0.1)
        }
      })
      
      withProgress(message = 'Performing GO-term analysis', value = 0, {
        # Number of times we'll go through the loop
        n <- 3
        
        for (i in 1:n) {
          # Each time through the loop, add another row of data. This is
          # a stand-in for a long-running computation.
          RunPy2= system2(command, args=allArgs, stdout=TRUE)
          
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("Doing part", i))
          
          # Pause for 0.1 seconds to simulate a long computation.
          #Sys.sleep(0.1)
        }
      })
      
      
    } else {print("Your input file does not exists")}
    
  }
  if(file.exists("Output/goelite_output/GO-Elite_results/CompleteResults/ORA_pruned/GOElite.input-GO_z-score_elite.txt")==TRUE){  
  #############################################################    Generate the GO-elite Plot     #################################################################    
  if(file.exists("Output/goelite_output/GO-Elite_results/CompleteResults/ORA_pruned/GOElite.input-GO_z-score_elite.txt")==TRUE){
  setwd(Home)
  setwd("Output/goelite_output/GO-Elite_results/CompleteResults")
  Files=list.files("ORA_pruned")
  ZscoresGO=Files[grep("GO_z-score", Files)]
  GOoutput = read.table(paste0("ORA_pruned/",ZscoresGO),
                        sep="\t",header=T,fill=T, quote="")
  
  if(identical(GOoutput$Ontology.ID,logical(0))==FALSE) {
    
    setwd(Home)
    dir.create("Outputfolder", showWarnings = F)
    
    GoEliteplot <- reactive({
      pheno="Phenotype"
      ## dir.create("Outputfolder", showWarnings = F)
      setwd(Home)
      #setwd("Outputfolder")
      GOoutput = GOoutput[(GOoutput$Ontology.Type=="biological_process" |
                             GOoutput$Ontology.Type=="molecular_function"),]
      ColNames<-GOoutput$Ontology.Name[10:1]
      par(mar=c(10,35,7,7))
      bp= barplot(GOoutput$Z.Score[10:1],
                  main=paste("GO Ontology", pheno),
                  horiz=TRUE,
                  yaxt='n',col="#D95F02",
                  xlab="Z-score",names.arg=ColNames,axisnames=TRUE,
                  cex.main=2,cex.axis=1.5,cex.lab=1.5,tck=-0.01)
      axis(2,at=bp,labels=GOoutput$Ontology.Name[10:1],
           tick=FALSE,las=2,cex.axis=1.5);
      
      abline(v=1.96,col="red",lwd=2,lty=1);
    })}
    
    output$plot <- renderPlot({
      print(GoEliteplot())
    })
  }
    #############################################################    Download the GO-elite Plot     #################################################################  
    output$Ontology<- downloadHandler(
      filename = function(){
        paste("pheno","GO_ElitePlot.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source("Scripts/GOoutput.R",local=TRUE)
      } )
    
    #################################################################################################################################################################   
    ################################################################################################################################################################# 
    
    
    #################################################################    Generate the KEGG Plot     #################################################################    
    KEGG<- reactive({
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
      write.table(topKEGG_Pheno,paste0("topKEGG_",pheno,".txt"),sep="\t",row.names=FALSE,quote=FALSE)
      
      
      ColNames<-topKEGG_Pheno$Pathway[10:1]
      par(mar=c(10,35,7,7))
      bp=barplot(topKEGG_Pheno$P.DE[10:1],
                 main=paste("KEGG pathways", pheno),
                 horiz=TRUE,
                 yaxt='n',col="#fec44f",names.arg=ColNames,axisnames=TRUE,
                 xlab="P-value",
                 cex.main=2,cex.axis=1.5,cex.lab=1.5,tck=-0.01)
      axis(2,at=bp,labels=topKEGG_Pheno$Pathway[10:1], tick=FALSE,las=1.58,cex.axis=1.5)
      abline(v=0.05,col="red",lwd=2,lty=1)
    })
    
    output$plot2 <- renderPlot({
      print(KEGG())
    },height = 800, width = 1500)
    
    #############################################################    Download the KEGG Plot     #################################################################      
    output$KEGGenrichment<- downloadHandler(
      filename = function(){
        paste("pheno","keggplot.pdf", sep = '')},
      content = function(file){
        pdf(file,width = 20, height = 20)
        source("../../Scripts/KEGGoutput.R",local=TRUE)
      } )
    
    ####################################################################################################################################################################################
    ############################################################################## TRANSCRIPTOME ANALYSIS ############################################################################## 
    ####################################################################################################################################################################################
    
    
    setwd(Home)
    
    geneList <- reactive({
      geneList <- info()
      g <- subset(geneList, select = input$columns)
      All_genes_Entrez<- data.frame(g$V1)
      colnames(All_genes_Entrez)<-"EntrezId"
      All_genes_Entrez
    })
    
    setwd(Home)
    KangUnivers<- read.table("RefData/KangUnivers.txt", sep="\t", header=T) 
    colnames(KangUnivers)<-c("EntrezId","Symbol") #Name the columns
    
    #2)All Kang genes that are in Modules
    Kang_genes<-read.table("RefData/Kang_dataset_genesMod_version2.txt",sep="\t",header=TRUE)
    
    #3)Generate Gene universe to be used for single gene lists
    tmp=merge(KangUnivers,Kang_genes,by.y="EntrezGene",by.x="EntrezId",all=TRUE) #18826
    KangUni_Final<-tmp[duplicated(tmp$EntrezId)==FALSE,] #18675
    
    geneList2 <- reactive({
      All_genes_Entrez <- geneList()
      Annotation_list<-as.data.frame(matrix(NA,nrow=length(All_genes_Entrez$EntrezId),ncol=2))#Create an empty matrix
      colnames(Annotation_list)<-c("EntrezId","Symbol")
      Annotation_list$EntrezId<-All_genes_Entrez$EntrezId
      for(i in 1:nrow(Annotation_list)){
        Annotation_list$Symbol[i] <- unlist(mget(as.character(All_genes_Entrez$EntrezId[i]),org.Hs.egSYMBOL,ifnotfound=NA))#Extract gene symbols of the entrez ids
      }
      Annotation_list })
    
    
    GOI_P <- reactive({
      Annotation_list <- geneList2()
      Annotation_list$Module = Kang_genes$Module[match(Annotation_list$EntrezId,Kang_genes$EntrezGene)] #Attach module information
      GOI_Pheno<-Annotation_list[duplicated(Annotation_list$EntrezId)==FALSE,] 
      GOI_Pheno })
    
    GOI_1 <- reactive({
      GOI_Pheno <- GOI_P()
      #Create list of genes of interest
      GOI1<-list(GOI_Pheno$Symbol)
      GOI1
    })
    
    
    pheno="Phenotype"
    
    
    
    #################################################################    Generate the Enrichment Plot     #################################################################        
    Enrichment<- reactive({
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
      
      par(mar=c(10,17,10,10))
      bp=plot(x=betas, y=1:length(betas),
              type="n", panel.first = grid(ny=NA),
              yaxt = "n", ylab="",
              xlim=xlim,
              xlab=expression(paste(log(OR)%+-%95,"%CI")),
              main=paste(pheno, "associated Genes"))
      
      abline(v=0,col="black", lty=3);
      axis(2, at=1:30,cex.axis=0.7,
           labels=c("All Brain Genes", paste0("Module",formatC(29:1, width=2, flag=0))),
           las=1);
      arrows(x0=CILower, x1=CIUpper, y0=30:1, y1=30:1, col=rainbow(30), length=0, lwd=2,code = 3)
      points(y=30:1, x=betas, pch=18, col="black")
      betas[is.na(betas)]=0
      text(y=(30:1)+0.5, x=betas, labels=pval, cex=1)
    })
    
    output$plot3 <- renderPlot({
      print(Enrichment())
    },height = 800, width = 900)
    
    
    #############################################################    Download the Enrichment Plot     #################################################################      
    output$ModuleEnrichment<- downloadHandler(
      filename = function(){
        paste("pheno","Enrichmentplot.pdf", sep = '')},
      content = function(file){
        pdf(file)
        #source("../Scripts/EnrichmentModules.R",local=TRUE)
        source(paste(Home,'/Scripts/EnrichmentModules.R',sep=""), local = TRUE)
      })
    
 
    ###############################################################    Generate the Enrichment Plot  MOI 1   #################################################################            
    output$MOI1<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 1"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module1<- downloadHandler(
      filename = function(){
        paste("pheno","Module1.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 1"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    
    ###############################################################    Generate the Enrichment Plot  MOI 2   #################################################################            
    output$MOI2<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 2"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module2<- downloadHandler(
      filename = function(){
        paste("pheno","Module2.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 2"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 3   #################################################################            
    output$MOI3<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 3"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module3<- downloadHandler(
      filename = function(){
        paste("pheno","Module3.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 3"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 4   #################################################################            
    output$MOI4<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 4"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module4<- downloadHandler(
      filename = function(){
        paste("pheno","Module4.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 4"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 5   #################################################################            
    output$MOI5<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 5"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module5<- downloadHandler(
      filename = function(){
        paste("pheno","Module5.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 5"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 6   #################################################################            
    output$MOI6<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 6"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module6<- downloadHandler(
      filename = function(){
        paste("pheno","Module6.pdf", sep = '')},
      content = function(file){
        pdf(file)
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 6"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 7   #################################################################            
    output$MOI7<- renderPlot ({
      source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 7"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module7<- downloadHandler(
      filename = function(){
        paste("pheno","Module7.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 7"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 8   #################################################################            
    output$MOI8<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 8"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module8<- downloadHandler(
      filename = function(){
        paste("pheno","Module8.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 8"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 9   #################################################################            
    output$MOI9<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 9"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module9<- downloadHandler(
      filename = function(){
        paste("pheno","Module9.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 9"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 10   #################################################################            
    output$MOI10<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 10"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module10<- downloadHandler(
      filename = function(){
        paste("pheno","Module10.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 10"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 11   #################################################################            
    output$MOI11<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 11"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module11<- downloadHandler(
      filename = function(){
        paste("pheno","Module11.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 11"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 12   #################################################################            
    output$MOI12<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 12"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module12<- downloadHandler(
      filename = function(){
        paste("pheno","Module12.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 12"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 13   #################################################################            
    output$MOI13<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 13"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module13<- downloadHandler(
      filename = function(){
        paste("pheno","Module13.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 13"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 14   #################################################################            
    output$MOI14<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 14"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module14<- downloadHandler(
      filename = function(){
        paste("pheno","Module14.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 14"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 15   #################################################################            
    output$MOI15<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 15"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module15<- downloadHandler(
      filename = function(){
        paste("pheno","Module15.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 15"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 16   #################################################################            
    output$MOI16<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 16"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module15<- downloadHandler(
      filename = function(){
        paste("pheno","Module16.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 16"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 17   #################################################################            
    output$MOI17<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 17"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module17<- downloadHandler(
      filename = function(){
        paste("pheno","Module17.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 17"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 18   #################################################################            
    output$MOI18<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 18"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module18<- downloadHandler(
      filename = function(){
        paste("pheno","Module18.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 18"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 19   #################################################################            
    output$MOI19<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 19"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module19<- downloadHandler(
      filename = function(){
        paste("pheno","Module19.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 19"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 20   #################################################################            
    output$MOI20<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 20"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module20<- downloadHandler(
      filename = function(){
        paste("pheno","Module20.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 20"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 21   #################################################################            
    output$MOI21<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 21"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module21<- downloadHandler(
      filename = function(){
        paste("pheno","Module21.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 21"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 22   #################################################################            
    output$MOI22<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 22"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module22<- downloadHandler(
      filename = function(){
        paste("pheno","Module22.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 22"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 23   #################################################################            
    output$MOI23<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 23"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module23<- downloadHandler(
      filename = function(){
        paste("pheno","Module23.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 23"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 24   #################################################################            
    output$MOI24<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 24"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module24<- downloadHandler(
      filename = function(){
        paste("pheno","Module24.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 24"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 25   #################################################################            
    output$MOI25<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 25"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module25<- downloadHandler(
      filename = function(){
        paste("pheno","Module25.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 25"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 26   #################################################################            
    output$MOI26<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 26"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module26<- downloadHandler(
      filename = function(){
        paste("pheno","Module26.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 26"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 27   #################################################################            
    output$MOI27<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 27"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module27<- downloadHandler(
      filename = function(){
        paste("pheno","Module27.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 27"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 28   #################################################################            
    output$MOI28<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 28"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module28<- downloadHandler(
      filename = function(){
        paste("pheno","Module28.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 28"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    ###############################################################    Generate the Enrichment Plot  MOI 29   #################################################################            
    output$MOI29<- renderPlot ({
        source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
      Module=df$Module[df$Module=="Module 29"]
        source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
    },height = 1000, width = 1000) 
    
    
    output$Module29<- downloadHandler(
      filename = function(){
        paste("pheno","Module29.pdf", sep = '')},
      content = function(file){
        pdf(file)
          source(paste(Home,'/Scripts/eachcall.R',sep=""), local = TRUE)
        Module=df$Module[df$Module=="Module 29"]
          source(paste(Home,'/Scripts/MOIgraph.R',sep=""), local = TRUE)
        #source("ModulePlot.R",local=TRUE)
        dev.off()
      })
    
    ########################################################################### Heatmaps ###########################################################################################################
    
    
    ##########################################################################################################################################################################################################
    Heatmaps<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 1"]
      source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
      source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap1<- renderPlot({
      print(Heatmaps())
    },height = 1000, width = 1000)
    
    output$Module1Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module1Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 1"]
        source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
        source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps2<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 2"]
      source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
      source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap2<- renderPlot({
      print(Heatmaps2())
    },height = 1000, width = 1000)
    
    output$Module2Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module2Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 2"]
        source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
        source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################
    Heatmaps3<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 3"]
      source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
      source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap3<- renderPlot({
      print(Heatmaps3())
    },height = 1000, width = 1000)
    
    output$Module3Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module3Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 3"]
        source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
        source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps4<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 4"]
      source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
      source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap4<- renderPlot({
      print(Heatmaps4())
    },height = 1000, width = 1000)
    
    output$Module4Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module4Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 4"]
        source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
        source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps5<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 5"]
      source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
      source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap5<- renderPlot({
      print(Heatmaps5())
    },height = 1000, width = 1000)
    
    output$Module5Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module5Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 5"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps6<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 6"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap6<- renderPlot({
      print(Heatmaps6())
    },height = 1000, width = 1000)
    
    output$Module6Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module6Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 6"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps7<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 7"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap7<- renderPlot({
      print(Heatmaps7())
    },height = 1000, width = 1000)
    
    output$Module7Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module7Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 7"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps8<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 8"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap8<- renderPlot({
      print(Heatmaps8())
    },height = 1000, width = 1000)
    
    output$Module8Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module8Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 8"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps9<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 9"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap9<- renderPlot({
      print(Heatmaps9())
    },height = 1000, width = 1000)
    
    output$Module9Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module9Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 9"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps10<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 10"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap10<- renderPlot({
      print(Heatmaps10())
    },height = 1000, width = 1000)
    
    output$Module10Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module10Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 10"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps11<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 11"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap11<- renderPlot({
      print(Heatmaps11())
    },height = 1000, width = 1000)
    
    output$Module11Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module11Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 11"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps12<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 12"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap12<- renderPlot({
      print(Heatmaps12())
    },height = 1000, width = 1000)
    
    output$Module12Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module12Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 12"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps13<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 13"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap13<- renderPlot({
      print(Heatmaps13())
    },height = 1000, width = 1000)
    
    output$Module13Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module13Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 13"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps14<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 14"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap14<- renderPlot({
      print(Heatmaps14())
    },height = 1000, width = 1000)
    
    output$Module14Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module14Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 14"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps15<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 15"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap15<- renderPlot({
      print(Heatmaps15())
    },height = 1000, width = 1000)
    
    output$Module15Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module15Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 15"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps16<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 16"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap16<- renderPlot({
      print(Heatmaps16())
    },height = 1000, width = 1000)
    
    output$Module16Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module16Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 16"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps17<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 17"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap17<- renderPlot({
      print(Heatmaps17())
    },height = 1000, width = 1000)
    
    output$Module17Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module17Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 17"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps18<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 18"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap18<- renderPlot({
      print(Heatmaps18())
    },height = 1000, width = 1000)
    
    output$Module18Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module18Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 18"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps19<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 19"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap19<- renderPlot({
      print(Heatmaps19())
    },height = 1000, width = 1000)
    
    output$Module19Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module19Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 19"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps20<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 20"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap20<- renderPlot({
      print(Heatmaps20())
    },height = 1000, width = 1000)
    
    output$Module20Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module20Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 20"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps21<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 21"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap21<- renderPlot({
      print(Heatmaps21())
    },height = 1000, width = 1000)
    
    output$Module21Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module21Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 21"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps22<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 22"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap22<- renderPlot({
      print(Heatmaps22())
    },height = 1000, width = 1000)
    
    output$Module22Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module22Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 22"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps23<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 23"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap23<- renderPlot({
      print(Heatmaps23())
    },height = 1000, width = 1000)
    
    output$Module23Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module23Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 23"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps24<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 24"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap24<- renderPlot({
      print(Heatmaps24())
    },height = 1000, width = 1000)
    
    output$Module24Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module24Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 24"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps25<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 25"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap25<- renderPlot({
      print(Heatmaps25())
    },height = 1000, width = 1000)
    
    output$Module25Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module25Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 25"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps26<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 26"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap26<- renderPlot({
      print(Heatmaps26())
    },height = 1000, width = 1000)
    
    output$Module26Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module26Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 26"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps27<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 27"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap27<- renderPlot({
      print(Heatmaps27())
    },height = 1000, width = 1000)
    
    output$Module27Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module27Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 27"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps28<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 28"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap28<- renderPlot({
      print(Heatmaps28())
    },height = 1000, width = 1000)
    
    output$Module28Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module28Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 28"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    ############################################ Module ##############################################################  
    Heatmaps29<- reactive({
      source("../Scripts/eachcall.R",local=TRUE)
      Module=df$Module[df$Module=="Module 29"]
    
     source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
     source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
    })
    
    output$Heatmap29<- renderPlot({
      print(Heatmaps29())
    },height = 1000, width = 1000)
    
    output$Module29Heatmap<- downloadHandler(
      filename = function(){
        paste("Brain","Module29Heatmap.pdf", sep = '_')},
      content = function(file){
        pdf(file)
        source("../Scripts/eachcall.R",local=TRUE)
        Module=df$Module[df$Module=="Module 29"]
      
       source(paste(Home,'/Scripts/cerebroScale.R',sep=""), local = TRUE)
       source(paste(Home,'/Scripts/BrainModule.R',sep=""), local = TRUE)
        dev.off()
      })
    #################################################################################################################################################################################################
  
  }
}
shinyApp(ui, server)
