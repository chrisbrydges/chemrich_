library(shiny)
library(shinythemes)
library(shinyjs)
library(devtools)
library(RCurl)
library(BiocManager)
options(repos = BiocManager::repositories())
library(GGally)
library(DT)
library(RJSONIO)
library(ape)
library(devEMF)
library(dynamicTreeCut)
library(extrafont)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(htmlwidgets)
library(igraph)
#library(magrittr)
library(network)
library(officer)
library(openxlsx)
library(phytools)
library(plotly)
library(plotrix)
library(rcdk)
library(rvg)
library(sna)
library(visNetwork)

options(shiny.maxRequestSize=50*1024^2)

# Functions that need to run to get MeSH ontologies from SMILES codes:
getCNames <- function(x) {
  if(length(which(treenames.df$MESHTREE==x))>0){
    as.character(treenames.df$ClassName[which(treenames.df$MESHTREE==x)])
  } else {
    x
  }
}

makeSmiles.clean <- function (smi) {
  smi <- gsub("@","",smi)
  smi <- gsub("O-","O",smi)
  smi <- gsub("H[1-4][+]","",smi)
  smi <- gsub("N[+]","N",smi)
  smi <- gsub("/|[\\]","",smi)
  smi <- gsub("[[]C(@{1,}|)H[]]","C",smi)
  smi <- gsub("[[]CH2[]]","C",smi)
  smi
}

load("mesh_bit_loc_list.RData")
load("cidmesh_smiles_fpbit.RData")
load("treenames.df.RData")

df.mega.mesh$CompoundName <- tolower(df.mega.mesh$CompoundName)
treenames.df$ClassName <- tolower(treenames.df$ClassName)


# Define UI for application
ui <- shinyUI(fluidPage(theme = shinytheme("cerulean"),
                        useShinyjs(),
                        # Application title
                        titlePanel("ChemRICH for Multiple Comparisons"),
                        p("A Shiny app to conduct ChemRICH analysis on multiple comparisons"),
                        
                        # Sidebar with column selection input
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("file1", "Upload a csv file.",
                                      multiple = FALSE, accept = c("text/csv", ".csv", 
                                                                   "text/comma-separated-values,text/plain")),
                            selectInput("compound_names_column", "Select compound names column", c("Upload data, then select column")),
                            selectInput("smiles_column", "Select SMILES column", c("Upload data, then select column")),
                            selectInput("cid_column", "Select PubChem CID column", c("Upload data, then select column")),
                            hr(),
                            actionButton("go", "Start"),
                            downloadButton("download_data", "Download")
                            ),
                          
                          # Main panel will have a completion message or error message
                          mainPanel(h3("Instructions:"),
                                    br(),
                                    p("Upload a csv file with columns of compound names, PubChem CIDs, SMILES, p values, and effect sizes."),
                                    br(),
                                    p("You will specify which column contains compound names, PubChem CIDs, and SMILES codes, so they can be named anything."),
                                    br(),
                                    p("The columns of p values must end with '_pvalue' and the effect size column must come immediately after it's corresponding p value column."),
                                    br(),
                                    helpText(a("Click here for an example data file.", href = "example_chemrich_file.csv")),
                                    hr(),
                                    textOutput("text")
                                    )
                          )
                        )
              )




# Define server logic
server <- shinyServer(function(input, output, session) {
  
  # Disable the two buttons
  disable("download_data")
  disable("go")
  
  
  # Once a csv file is uploaded, populate the column names selections and enable the "go" button
  observe({
    x <- input$file1
    if(is.null(x)){
    }else{
      req(input$file1)
      enable("go")
      file_location <- input$file1$datapath
      dat1 <- data.frame(read.csv(file_location))
      
      updateSelectInput(session, "compound_names_column", 
                        label = "Select compound names column", 
                        choices = colnames(dat1))
      
      updateSelectInput(session, "smiles_column", 
                        label = "Select SMILES column", 
                        choices = colnames(dat1))
      
      updateSelectInput(session, "cid_column", 
                        label = "Select PubChem CID column", 
                        choices = colnames(dat1))
            
    }
  })
  
  # Reload the data in again, and select the already-specified columns
  dat2 <- eventReactive(input$go,{
    showNotification("Reading dataset", duration = 5)
    
    req(input$file1)
    
    file_location <- input$file1$datapath
    ndf <- data.frame(read.csv(file_location))
    compound_names_column_name <- input$compound_names_column
    smiles_column_name <- input$smiles_column
    cid_column_name <- input$cid_column
    
    #########################
    ## Fatty acid labels ####
    #########################
    
    smi.all.fas <- as.character(sapply(ndf[[smiles_column_name]], makeSmiles.clean))
    falabelvec <- sapply(smi.all.fas, function(x) {
      elecount <- table(strsplit(gsub("[0-9]|[)]|[(]|=","",x),"")[[1]])
      falabel <-  ""
      if (length(table(c("c","o")%in%tolower(names(elecount))  )  )==1) {
        if(length(grep("n",x,ignore.case = T))==0) {
          if (elecount['C']>7 & length(grep("CCCC",x))==1 & length(grep("C2",x))!=1  ) { # long carbon but not aromatic or cyclic.
            if (elecount['O']==2) {
              dlen <- length(strsplit(x,"=")[[1]])-2
              falabel <- paste(c("FA",elecount['C'],dlen), collapse="_")
            }
            if(elecount['O']>=3) { ## Put Rules here. How many O and then how many carbon chain. That will make the class.
              if( length(grep("C1",x))==1) {
                if (length(strsplit(x,"C1")[[1]]) ==3 ) {
                  dlen <- length(strsplit(x,"=")[[1]])-2
                  #falabel <- paste(c("Prostaglandin",elecount['C']), collapse="_")
                } else {
                  dlen <- length(strsplit(x,"=")[[1]])-2
                  falabel <- paste(c("Epoxy FA",elecount['C']), collapse="_")
                }
              } else {
                if (length(strsplit(x,"=O|CO|OC")[[1]])-2==0){
                  dlen <- length(strsplit(x,"=")[[1]])-2
                  falabel <- paste(c("OH-FA",elecount['C'],dlen,(elecount['O']-2)), collapse="_")
                } else {
                  if (length(strsplit(x,"OC|CO")[[1]]) <3 ) {
                    dlen <- length(strsplit(x,"=")[[1]])-2
                    falabel <- paste(c("O=FA",elecount['C'],dlen), collapse="_")
                  }
                }
              }
            }
          }
        }
      }
      falabel
    })
    
    falabelvec[which(falabelvec=="OH-FA_20_3_2")] <- "DiHETrE"
    falabelvec[which(falabelvec=="OH-FA_20_4_2")] <- "DiHETE"
    falabelvec[which(falabelvec=="O=FA_18_3")] <- "oxo-ODE"
    falabelvec[which(falabelvec=="O=FA_20_5")] <- "oxo-ETE"
    falabelvec[which(falabelvec=="OH-FA_18_1_2")] <- "DiHOME"
    falabelvec[which(falabelvec=="OH-FA_18_1_3")] <- "TriHOME"
    falabelvec[which(falabelvec=="OH-FA_18_2_1")] <- "HODE"
    falabelvec[which(falabelvec=="OH-FA_18_2_2")] <- "DiHODE"
    falabelvec[which(falabelvec=="OH-FA_18_3_1")] <- "HOTrE"
    falabelvec[which(falabelvec=="OH-FA_20_3_1")] <- "HETrE"
    falabelvec[which(falabelvec=="OH-FA_20_4_1")] <- "HETE"
    falabelvec[which(falabelvec=="OH-FA_20_5_1")] <- "HEPE"
    falabelvec[which(falabelvec=="OH-FA_22_5_2")] <- "DiHDPE"
    falabelvec[which(falabelvec=="Epoxy FA_22")] <- "EpDPE"
    falabelvec[which(falabelvec=="Epoxy FA_18")] <- "EpETrE"
    falabelvec[which(falabelvec=="Epoxy FA_20")] <- "EpODE"
    falabelvec[grep("^FA_[0-9]{1,2}_0$", falabelvec)] <- "Saturated FA"
    falabelvec[grep("^FA_[0-9]{1,2}_[1-9]$", falabelvec)] <- "UnSaturated FA"
    
    showNotification("Computing sub-structure fingerprint (this can take a while)", duration = 30)
    
    smiles_list <- sapply(1:nrow(ndf),function(x){
      rcdk::parse.smiles(ndf[[smiles_column_name]][x])[[1]]
    })
    
    fps <- t(sapply(1:nrow(ndf), function(x) {
      gc()
      library(rcdk)
      xy <- 0
      xy <- as.character(rcdk::get.fingerprint(smiles_list[[x]],type="pubchem"))
      xy
    }))
    
    cid.mesh.df <- data.frame(CID=df.mega.mesh$CID[df.mega.mesh$CID%in%ndf[[cid_column_name]]],MESHTREE= df.mega.mesh$MESSTREE[df.mega.mesh$CID%in%ndf[[cid_column_name]]], stringsAsFactors = F)
    
    ## Covered by MeSH ###
    
    inmesh.vec <- rep("No",nrow(ndf))
    inmesh.vec[ndf[[cid_column_name]]%in%cid.mesh.df$CID] <- "Yes"
    
    ## chemical similarity matrix.
    df1.bitmat <- do.call(rbind,lapply(fps,function(x) as.integer(strsplit(x,"")[[1]][1:881])))
    df1.bitmat.location <- lapply(1:nrow(df1.bitmat), function(x) { which(df1.bitmat[x,]==1) })
    only_b <- sapply(1:length(bitloclist), function(x) {  length(bitloclist[[x]]) })
    bitmeans <- sapply(1:length(bitloclist), function(x) {  median(bitloclist[[x]]) })
    fpsmeans <- sapply(df1.bitmat.location, function(x){median(x)})
    
    showNotification("Obtaining MESH class annotation", duration = 30)
    
    
    directlabels <- sapply(tolower(ndf$compound_name), function(x) {
      clabel="Not Found"
      findind <- which(df.mega.mesh$CompoundName==x)
      if(length( findind)>0) {
        classvec <- as.character(df.mega.mesh$MESSTREE[findind])
        classvec <- strsplit(classvec[1],";")[[1]]
        if( length(grep("^D01[.]",classvec)) > 0 ) {
          classvec <- classvec[-grep("^D01[.]|^D03[.]",classvec)]
        }
        clabel <- names(which.max(sapply(classvec,nchar)))
      }
      clabel
    })
    
    labelvec <- sapply(1:nrow(ndf), function(i) {
      clabel <- "Not Found"
      if(falabelvec[i]=="" & inmesh.vec[i]=="No" & directlabels[i]=="Not Found") {
        #if(falabelvec[i]=="") {
        #print(i)
        meanindex <- which(bitmeans < (fpsmeans[i]+5) & bitmeans > (fpsmeans[i]-5))
        bitloclist.sb <- bitloclist[meanindex]
        only_b.sb <- only_b[meanindex]
        overlapvec <- sapply(1:length(bitloclist.sb), function(x) {  length(which(bitloclist.sb[[x]]%in%df1.bitmat.location[[i]]==TRUE)) })
        tmvec <- overlapvec/((length(df1.bitmat.location[[i]])+only_b.sb)-overlapvec)
        if(length(which(tmvec>0.90))>0) {
          if(length(which(tmvec>0.98))>0){
            cidindex <- meanindex[which(tmvec>0.98)]
            if(length(cidindex)==1) {
              clabel <- df.mega.mesh$MESSTREE[cidindex]
            } else {
              clabel.table <- sort(table(unlist(sapply(unique(df.mega.mesh[cidindex[order(tmvec[which(tmvec>0.98)])],"MeSHUID"])[1:10], function(x) {if(!is.na(x)) { strsplit(df.mega.mesh$MESSTREE[which(df.mega.mesh$MeSHUID==x)][1],";")  }}))),decreasing = T)
              clabel.table <- which(clabel.table==clabel.table[1])
              clabel.table.sort <- sort(sapply(names(clabel.table),nchar),decreasing = T)
              clabel.table.sort.max <- which.max(clabel.table.sort)
              if(length(clabel.table.sort.max==1)) {
                clabel <- names(clabel.table.sort.max)
              } else {
                clabel <- sort(names(clabel.table.sort.max))[1]
              }
            }
          } else {
            ## CID_MESH
            cidindex <- meanindex[which(tmvec>0.90)]
            if(length(cidindex)==1) {
              clabel <- df.mega.mesh$MESSTREE[cidindex]
            } else {
              clabel.table <- sort(table(unlist(sapply(unique(df.mega.mesh[cidindex[order(tmvec[which(tmvec>0.90)])],"MeSHUID"])[1:10], function(x) {if(!is.na(x)) { strsplit(df.mega.mesh$MESSTREE[which(df.mega.mesh$MeSHUID==x)][1],";")  }}))),decreasing = T)
              clabel.table <- which(clabel.table==clabel.table[1])
              clabel.table.sort <- sort(sapply(names(clabel.table),nchar),decreasing = T)
              clabel.table.sort.max <- which.max(clabel.table.sort)
              if(length(clabel.table.sort.max==1)) {
                clabel <- names(clabel.table.sort.max)
              } else {
                clabel <- sort(names(clabel.table.sort.max))[1]
              }
            }
          }
        }
      }
      clabel
    })
    
    finalMesh.df <- cid.mesh.df
    finalMesh.df <- finalMesh.df[which(finalMesh.df$CID%in%ndf[[cid_column_name]][which(falabelvec!="")]==FALSE),]
    finalMesh.df <- finalMesh.df[which(finalMesh.df$CID%in%ndf[[cid_column_name]][which(directlabels!="Not Found")]==FALSE),]
    finalMesh.df <- rbind(finalMesh.df,data.frame(CID=ndf[[cid_column_name]], MESHTREE=labelvec ))
    finalMesh.df$NewMesh <- finalMesh.df$MESHTREE
    finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh!="Not Found"),]
    finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh!=""),]
    finalMesh.df <- finalMesh.df[!duplicated(finalMesh.df),]
    
    # Calculate the chemical similarity and clusters.
    
    showNotification("Computing chemical similarity", duration = 30)
    
    m <- df1.bitmat
    mat <- m%*%t(m)
    len <- length(m[,1])
    s <- mat.or.vec(len,len)
    for (i in 1:len) {
      for (j in 1:len){
        s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
      }
    }
    diag(s) <- 0
    hc <- hclust(as.dist(1-s), method="ward.D2") # ward method provide better clusters.
    clust1 <- cutreeDynamic(hc,distM = as.matrix(1-s),deepSplit =4, minClusterSize = 3)  # can be used to merge cluster, but it is better to keep them as it is. # Clusters are detected using the average linkage hclust.
    ndf$ClusterNumber <- clust1
    ndf$xlogp <- as.numeric(sapply(ndf[[smiles_column_name]], function(x)  { rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]]) }))
    
    ## creating the final label data frame.
    
    finalterm.df <- data.frame(CID=ndf[[cid_column_name]], Clabel=falabelvec,stringsAsFactors = F) # first add the fatty acid labels.
    directlabindex <- as.integer(which(directlabels!="Not Found"))[which(as.integer(which(directlabels!="Not Found"))%in%which(finalterm.df$Clabel=="")==TRUE)] ## then add the direct labels found by names matching
    finalterm.df$Clabel[directlabindex] <- as.character(directlabels[directlabindex])
    
    for (i in 1:nrow(finalterm.df)) {
      if(finalterm.df$Clabel[i]=="" & length(which(finalMesh.df$CID==ndf[[cid_column_name]][i]))>0 ) {
        finalterm.df$Clabel[i] <- names(which.max(sapply(unlist(strsplit(finalMesh.df$NewMesh[which(finalMesh.df$CID==ndf[[cid_column_name]][i])],";")),nchar)))
      }
    }
    
    ##########################################
    ####  Detect for new compound clusters ###
    ##########################################
    showNotification("Detecting new compound classes", duration = 30)
    
    finalterm.df.2 <- finalterm.df
    newClustVec <- names(which(table(ndf$ClusterNumber[which(finalterm.df.2$Clabel=="")])>3))
    clustMeanvec <- sapply(newClustVec, function(x) {  mean(s[which(ndf$ClusterNumber==x),which(ndf$ClusterNumber==x)])  }  )
    newClustVec <- newClustVec[which(clustMeanvec>0.70)]
    if(length(newClustVec)>0) {
      for(i in which(finalterm.df.2$Clabel=="")) {
        if(ndf$ClusterNumber[i]%in%newClustVec){
          finalterm.df$Clabel[i] <- paste0("NewCluster_",ndf$ClusterNumber[i])
        }
      }
    }
    
    ##### Map the compounds that have at least 0.75 similarity to others. Only for compounds that do not have any labels.
    
    for ( i in which(finalterm.df$Clabel=="")) { ## if there is a metabolite that has score higher than 0.80 then we get the class using that compound.
      print(i)
      if(max(s[i,])>0.75) {
        simorder <- order(s[i,],decreasing = T)[which(s[i,][order(s[i,],decreasing = T)]>0.75)]
        simorder.class <- sapply(simorder, function(x) { finalterm.df$Clabel[x]})
        simorder.class <- simorder.class[!is.na(simorder.class)]
        if(!is.na(simorder.class[1])){
          if(simorder.class[1]!=""){
            finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
          } else if(length(simorder.class)>1) {
            finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
          }
        }
      }
    }
    
    showNotification("Detecting non-overlapping class definitions", duration = 30)
    
    
    finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
    finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))
    
    exclusionVec <- c("D02","D03.383","D03.633.100","D03.633.300","D03.633.400","D03.633","D03.605","D02.241.081") ## we will have some static list of terms that need to be excluded.
    exclusionVec <- c(exclusionVec, unique(falabelvec)[-1]) ## if we see Fatty acid label, we dont touch them.
    
    for ( i in which(finalterm.df$gCount<3)) { ## Drop the compound to the neareast one.
      qpat <- gsub("[.][0-9]{2,3}$","",finalterm.df$Clabel[i])
      if(length(grep(qpat,finalterm.df$Clabel))>2 & !qpat%in%exclusionVec){
        finalterm.df$Clabel[i] <- qpat
      }
    }
    
    finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
    finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))
    
    for ( i in which(finalterm.df$gCount<3)) { ## Map to the closest ones.
      print(i)
      if(max(s[i,])>0.85) {
        simorder <- order(s[i,],decreasing = T)[which(s[i,][order(s[i,],decreasing = T)]>0.85)]
        simorder.class <- sapply(simorder, function(x) { finalterm.df$Clabel[x]})
        simorder.class <- simorder.class[!is.na(simorder.class)]
        if(simorder.class[1]!=""){
          finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
        } else if(length(simorder.class)>1) {
          finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
        }
      }
    }
    
    finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
    finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))
    
    # Repeat it one more time.
    for ( i in which(finalterm.df$gCount<3)) { ## Drop the compound to the neareast one.
      qpat <- gsub("[.][0-9]{2,3}$","",finalterm.df$Clabel[i])
      if(length(grep(qpat,finalterm.df$Clabel))>2 & !qpat%in%exclusionVec){
        finalterm.df$Clabel[i] <- qpat
      }
    }
    finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
    finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))
    
    for ( i in which(finalterm.df$gCount<3)) { ## Map to the closest ones.
      if(max(s[i,])>0.85) {
        simorder <- order(s[i,],decreasing = T)[which(s[i,][order(s[i,],decreasing = T)]>0.85)]
        simorder.class <- sapply(simorder, function(x) { finalterm.df$Clabel[x]})
        simorder.class <- simorder.class[!is.na(simorder.class)]
        if(simorder.class[1]!=""){
          finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
        } else if(length(simorder.class)>1) {
          finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
        }
      }
    }
    
    finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
    finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))
    finalterm.df$Clabel[which(finalterm.df$Count<3)] <- finalterm.df.2$Clabel[which(finalterm.df$Count<3)] ### We reverse back the original labels as this did not create any higher labels.
    
    finallabelvec <- finalterm.df$Clabel
    
    HasSaturatedFats <- names(which(table(finallabelvec[grep("D10|D09.400.410",finallabelvec)[which(sapply(grep("D10|D09.400.410",finallabelvec), function(x)  { length(grep("C=C",ndf[[smiles_column_name]][x]))  })==0)]])>2)) ### we are only selecting lipid classes that has atleast 3 saturated lipids.
    
    for (i in 1:nrow(finalterm.df)) {
      if(finallabelvec[i]%in%HasSaturatedFats){
        if(length(grep("C=C",ndf[[smiles_column_name]][i]))==0) {
          finallabelvec[i] <- paste0("Saturated_",getCNames(finallabelvec[i]))
        } else {
          finallabelvec[i] <- paste0("Unsaturated_",getCNames(finallabelvec[i]))
        }
      }
    }
    clusterids <- sapply(as.character(finallabelvec),getCNames)
    ndf$MeSH_Class <- clusterids
    
    classVariable <- "MeSH_Class"
    compoundName <- "compound_name"
    
    data_dict <- ndf # Data Dictionary
    pvalvec <- grep("pvalue",names(data_dict))
    classVec <- names(which(table(data_dict[[classVariable]])>2))
    clusterids <- classVec[which(classVec!="")]
    
    l <- list("ChemRICH_Input" = data_dict)
    
    my_pptx <- read_pptx()
    
    for(k in pvalvec) {
      df1 <- data.frame(Compound = data_dict[[compoundName]], Class = data_dict[[classVariable]], xlogp = data_dict$xlogp, SMILES = data_dict[[smiles_column_name]], pvalue= as.numeric(data_dict[[k]]), foldchange= as.numeric(data_dict[[k+1]]), stringsAsFactors = F)
      
      ## Get CHEMRICH computation
      cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
        cl.member <- which(df1$Class==x)
        if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
          pval.cl.member <- df1$pvalue[cl.member]
          p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
          p.test.results$p.value
        } else {
          1
        }
      })
      
      cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
      clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)
      
      clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
        dfx <- df1[which(df1$Class==x),]
        dfx$Compound[which.min(dfx$pvalue)]
      })
      
      altrat <- sapply(clusterdf$name, function (k) {
        length(which(df1$Class==k & df1$pvalue<0.05))/length(which(df1$Class==k))
      })
      
      uprat <-sapply(clusterdf$name, function (k) {
        length(which(df1$Class==k & df1$pvalue<0.05 & df1$foldchange > 1.00000000))/length(which(df1$Class==k & df1$pvalue<0.05))
      })
      
      clust_s_vec <- sapply(clusterdf$name, function (k) {
        length(which(df1$Class==k))
      })
      
      clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(df1$Class==k & df1$pvalue<0.05))})
      clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(df1$Class==k & df1$pvalue<0.05 & df1$foldchange > 1.00000000))})
      clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(df1$Class==k & df1$pvalue<0.05 & df1$foldchange < 1.00000000))})
      clusterdf$upratio <- uprat
      clusterdf$altratio <- altrat
      clusterdf$csize <- clust_s_vec
      clusterdf <- clusterdf[which(clusterdf$csize>2),]
      clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
      clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(df1$xlogp[which(df1$Class==x)]) }))
      #clusterdf
      clusterdf$Compounds <- sapply(clusterdf$name, function(x) {
        dfx <- df1[which(df1$Class==x),]
        paste(dfx$Compound,collapse="<br>")
      }) ## this one is the label on the tooltip of the ggplotly plot.
      clustdf <- clusterdf[which(clusterdf$pvalues!=1),]
      
      #################################################
      ########## Impact Visualization Graph ###########
      #################################################
      
      clustdf.alt.impact <- clustdf[which(clustdf$pvalues<0.05 & clustdf$csize>1 & clustdf$alteredMetabolites>1) ,]
      clustdf.alt.impact <- clustdf.alt.impact[order(clustdf.alt.impact$xlogp),]
      clustdf.alt.impact$order <- 1:nrow(clustdf.alt.impact) ### Order is decided by the hclust algorithm.
      clustdf.alt.impact$logPval <- -log(clustdf.alt.impact$pvalues)
      
      p2 <- ggplot(clustdf.alt.impact,aes(x=xlogp,y=-log(pvalues)))
      p2 <- p2 + geom_point(aes(size=csize, color=upratio)) +
        #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
        scale_color_gradient(low = "blue", high = "red", limits=c(0,1))+
        scale_size(range = c(5, 30)) +
        scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+4  )) +
        scale_x_continuous(" Lipophilicity ", limit=c(min(df1$xlogp)-2,max(df1$xlogp)+2)) +
        theme_bw() +
        labs(title = paste0("ChemRICH cluster impact plot for ", gsub("_pvalue","",names(data_dict)[k]))) +
        geom_label_repel(aes(label = name), color = "gray20",family="Arial",data=subset(clustdf.alt.impact, csize>2),force = 5)+
        theme(text=element_text(family="Arial Black"))+
        theme(
          plot.title = element_text(face="bold", size=30,hjust = 0.5),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20, angle=90),
          panel.grid.major = element_blank(), # switch off major gridlines
          panel.grid.minor = element_blank(), # switch off minor gridlines
          legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
          legend.title = element_blank(), # switch off the legend title
          legend.text = element_text(size=12),
          legend.key.size = unit(1.5, "lines"),
          legend.key = element_blank(), # switch off the rectangle around symbols in the legend
          legend.spacing = unit(.05, "cm"),
          axis.text.x = element_text(size=10,angle = 0, hjust = 1),
          axis.text.y = element_text(size=15,angle = 0, hjust = 1)
        )
      my_pptx <- add_slide(my_pptx, layout = "Title and Content", master = "Office Theme") |>
        ph_with(dml(ggobj = p2), location = ph_location(type = "body",width=10, height=8,left = 0, top = 0))
      
      clustdf.e <- clusterdf[order(clusterdf$pvalues),]
      clustdf.e$pvalues <- signif(clustdf.e$pvalues, digits = 2)
      clustdf.e$adjustedpvalue <- signif(clustdf.e$adjustedpvalue, digits = 2)
      clustdf.e$upratio <- signif(clustdf.e$upratio, digits = 1)
      clustdf.e$altratio <- signif(clustdf.e$altratio, digits = 1)
      clustdf.e <- clustdf.e[,c("name","csize","pvalues","adjustedpvalue","keycpdname","alteredMetabolites","upcount","downcount","upratio","altratio")]
      names(clustdf.e) <- c("Cluster name","Cluster size","p-values","FDR","Key compound","Altered metabolites","Increased","Decreased","Increased ratio","Altered Ratio")
      l[[gsub("_pvalue","",names(data_dict)[k])]] <- clustdf.e
    }
    
    
    print(my_pptx, target = "ChemRICH_bubble_plots.pptx") |>
      invisible()
    openxlsx::write.xlsx(l, file = "ChemRICH_results.xlsx", asTable = TRUE)
    showNotification("ChemRICH_results.xlsx has been saved", duration = 5)
    
    
    
    
    
    
    
    library(dplyr)
    library(purrr)
    library(tidyr)
    
    
    sheets <- names(l)[-1]
    df <- l[-1]
    
    # Throw out the non-significant compound clusters for each pairwise comparison
    for (j in 1:length(df)) {
      df[[j]] <- dplyr::filter(df[[j]], `p-values` < 0.05)
      df[[j]] <- dplyr::select(df[[j]], `Cluster name`, `Increased ratio`)
    }
    
    # Turn the list of data frames into one merged data frame
    df <- reduce(df, dplyr::full_join, by = "Cluster name")
    colnames(df)[2:ncol(df)] <- sheets
    df <- tidyr::pivot_longer(df, cols = -`Cluster name`)
    df$`Cluster name` <- factor(df$`Cluster name`)
    
    # Make the heat map
    p <- ggplot2::ggplot(df, aes(x = name, y = `Cluster name`, fill = value)) + 
      geom_tile() +
      scale_fill_gradient(name = "Effect Direction", 
                          low = "blue", 
                          high = "red", 
                          limits = c(0, 1), 
                          na.value = "white", 
                          breaks = c(0, 0.5, 1), 
                          labels = c("All negative", "Mixed", "All positive")) +
      theme_classic() +
      scale_x_discrete(name = "Comparison") +
      scale_y_discrete(limits = rev(levels(df$`Cluster name`))) +
      theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90))
    
    ggplot2::ggsave("ChemRICH_heatmap.png", p, height = 8, width = 8, dpi = 300)
    showNotification("ChemRICH_heatmap.png has been created")
    
    
    
    # Enable the download button for file download
    enable("download_data")
    
    return("ChemRICH is complete. Please click the 'Download' button to download the zip file.")
    
  })
  
  # Put the text in line 616 on the main panel of the app
  output$text <- renderText({dat2()})
  
    
    output$download_data <- downloadHandler(
      filename = function() {"ChemRICH_output.zip"},
      content = function(fname) {zip(fname, c("ChemRICH_bubble_plots.pptx", "ChemRICH_results.xlsx", "ChemRICH_heatmap.png"))},
      contentType = "application/zip")
    
})

# Run the application 
shinyApp(ui = ui, server = server)
