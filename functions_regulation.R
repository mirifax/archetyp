library(dplyr)

## some of the functions used for the analysis scripts of novosparc output data 

#want to check the high/liw weights
feat_from_arche <- function(nmf, archetype=2, weights="high", n=25){
  
  #nmf: is the nmf object imported with readnmf()
  #archetye: the number of archetype u want to inspect
  #weights: looking at high or low weights
  #n: number of features zou want to get 
  
  data = nmf$data
  if(weights=="high"){decreasing<-T}
  if(weights=="low"){decreasing<-F}
  
  w <- nmf$weights[,archetype] #subset to archetype specific weight
  features <- nmf$features
  n_features <- length(features)
  w_order <- order(w, decreasing = decreasing)
  
  plot(c(1:n_features),w[w_order],main="ranked features")
  points(1:n,w[w_order][1:n],col = "red")
  
  features_out <- features[w_order][1:n]
  features_out_embryo <- data[features_out,]
  return(features_out_embryo)
}


 #get for a specific feature the according weights for the archetype
wei_from_feat<- function(nmf, feature,locations,type="RNA", lf=25){
 #nmf is the nmf read in with readnmf()
  #feature is the gene or ATAC peak/region u want to inspect
  #locations are the coordinates for the embryo plot 
 
  if(!feature %in% nmf$features){print(paste("the feature:",feature,"is not present!"))}
  weights <- nmf$weights
  data <- nmf$data
  
  barplot(as.numeric(weights[feature,]), main=paste("weight on each LF from:",feature),xlab="Latent Factor", ylab="weight", names.arg = c(1:lf))
  pl_embryos(data[feature,],locations,type=type,psize=3)
}




# load in nmf data
readnmf <- function(folder){
  nobject <- list(data= read.table(file.path(folder,"input_data.txt")),
                features= as.character(read.table(file.path(folder,"features.txt"))$x),
                latfact = read.table(file.path(folder,"LatentFactors.txt")),
                weights = read.table(file.path(folder,"weights.txt")))
  return(nobject)
  }



#get the maximum score of each granges in the bigwig file 
max_per_bigwig_peak <- function(bigwigfile1, bigwigfile2, granges){
  #bigwigfile should be the path to the bigwig file
  #granges the granges object of peaks we want to get the scores for
  vec <- c(length(granges)) # store the scores in
  for (g in c(1:length(granges))){
    bw_1 <- import(bigwigfile1,format="bw", which=granges[g])
    bw_2 <- import(bigwigfile2,format="bw", which=granges[g])
    vec[g] <- (max(bw_1$score)+max(bw_2$score))/2
  }
  return(vec)
}


make_binary <- function(x) {
  #make eyerything above quantile 0.75 to 1 , rest to 0
  quant <- quantile(x, 0.8)
  #med <- median(x)
  #quant <- 4.066944e-05/2
  x[x>quant]<-1
  x[x!=1] <-0
  return(x)
}


perc_covered <- function(GEX, CRM) {
  #function for retreiving percent of overlap of accessibile areas within the GEX (not the other way around!)
  # GEX has to be a binarised digital gene expression dim: (1,3039)
  # CR< has to be a binarised digital acccessibility profil  dim: (1,3039)
  a <- length(which(CRM==1)) #how many accessible locations 
  b <- length(which((CRM==1)&(GEX==1))) #how many accessible areas overlap with expression
  c <- b/a #percent of accessibility covered
  return(c)
}


ReScale <- function(x,first,last){(last-first)/(max(x)-min(x))*(x-min(x))+first}




################     plotting  multiple spatial patterns #########


pl_embryos_dorsal <- function(patterns, coordinates, title="", type="RNA", psize=1){
  
  ##plot from dorsal view
  stopifnot(is.character(title))
  if (class(patterns)!="numeric"){
    if (dim(patterns)[2]==3039){
      patterns <- t(patterns)
    }
  }
  #double the patterns and reduce to no negative z vvalues
  locations2 <- coordinates
  locations2$y <- locations2$y * -1
  locations3 <- rbind(coordinates,locations2)
  
  test<- append(as.numeric(patterns),as.numeric(patterns))
  test2 <- cbind(test,locations3)
  test_dorsal<- test2[test2$z>2.25,][1:3039,]
  
  
  if(type=="RNA"){option='viridis';begin=0;end=1; legendname="vGEX"}
  if(type=="ATAC"){option='magma';begin=0;end=1; legendname="vCA"}
  
  forplot <- test_dorsal %>% reshape2::melt( id.vars= c("x","y","z"))#plot the gene
  plot <- ggplot(forplot, aes(x = x, y = y, color = value)) +
    ggtitle(title)+ labs( col=legendname) +
    scale_colour_viridis_c(option=option,direction=-1,begin=begin, end=end)+
    geom_point(size = psize)+
    facet_wrap(~variable) + theme_classic()+#theme(aspect.ratio = 0.5)+
    theme(strip.text = element_text(margin= margin(1,0),size = 10, color="white"), 
          strip.background = element_rect(fill = "white", color="white"),
          aspect.ratio = 0.45,
          axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    ) 
  return(plot)
  
}

####################

pl_fig_embryos <- function(patterns, coordinates, title="", type="RNA", psize=1,nrow=1){
  #patterns are x by 3039 matrices 
  #are the coordinates for the loactions (always x,y,z)
  stopifnot(is.character(title))
  if (class(patterns)!="numeric"){
    if (dim(patterns)[2]==3039){
      patterns <- t(patterns)
    }
  }
  
  if(type=="RNA"){option='viridis';begin=0;end=1; legendname="vGEX"}
  if(type=="ATAC"){option='magma';begin=0;end=1; legendname="vCA"}
  if(type=="LF"){option='rocket';begin=0;end=1}
  if(type=="resp_TF"){option='magma';begin=0.2;end=0.9}
  if(type=="qc"){option='cividis';begin=0;end=1;legendname="sum of probabilities"}
  
  forplot <- cbind(coordinates, patterns) %>% reshape2::melt( id.vars= c("x","y","z"))#plot the gene
  plot <- ggplot(forplot, aes(x = x, y = z, color = value)) +
    ggtitle(title)+ labs( col=legendname) +
    scale_colour_viridis_c(option=option,direction=-1,begin=begin, end=end)+
    geom_point(size = psize)+
    facet_wrap(~variable,nrow=nrow) + theme_classic()+#theme(aspect.ratio = 0.5)+
    theme(strip.text = element_text(margin= margin(1,0),size = 10, color="white"), 
          strip.background = element_rect(fill = "white", color="white"),
          aspect.ratio = 0.45,
          axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5),
          #legend.position = "none"
          legend.text = element_blank(),
          legend.title = element_blank()
    ) 
  return(plot)
}
  

pl_embryos <- function(patterns, coordinates, title="", type="RNA", psize=1,nrow=1){
  stopifnot(is.character(title))
  if (class(patterns)!="numeric"){     #patterns are x by 3039 matrices # are the coordinates for the loactions (always x,y,z)
    if (dim(patterns)[2]==3039){
      patterns <- t(patterns)
    }
  }
  
  if(type=="RNA"){option='viridis';begin=0;end=1; legendname="vGEX"}
  if(type=="ATAC"){option='magma';begin=0;end=1; legendname="vCA"}
  if(type=="LF_ATAC"){option='magma';begin=0;end=0.9}
  if(type=="resp_TF"){option='magma';begin=0.2;end=0.9}
  if(type=="qc"){option='cividis';begin=0;end=1;legendname="coverage"}

  forplot <- cbind(coordinates, patterns) %>% reshape2::melt( id.vars= c("x","y","z"))#plot the gene
  plot <- ggplot(forplot, aes(x = x, y = z, color = value)) +
    ggtitle(title)+ labs( col=legendname) +
    scale_colour_viridis_c(option=option,begin=begin,end=end,direction=-1)+
    geom_point(size = psize)+
    facet_wrap(~variable, nrow = nrow) + theme_classic()+
    theme(strip.text = element_text(margin= margin(1,0),size = 10, color="black"),
          strip.background = element_rect(fill = "light grey", color="white"),
          aspect.ratio = 0.45,
          axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5)
          #legend.position = "none"
    )
  return(plot)
}

pl_LF_embryos <- function(patterns, coordinates, title="",type="RNA",psize=1,nrow=4){
  stopifnot(is.character(title))
  if (class(patterns)!="numeric"){
    if (dim(patterns)[2]==3039){
      patterns <- t(patterns)
    }
  }
  
  if(type=="RNA"){option='viridis';begin=0;end=1; legendname="values"}
  if(type=="ATAC"){option='magma';begin=0;end=1; legendname="values"}
  #if (type=="RNA"){low = "yellow";high= "dark blue"}
  #if (type=="ATAC"){low="yellow";high="dark blue"}
  #if (type=="grey"){low="light gray";high='black'} 
  
  forplot <- cbind(coordinates, patterns) %>% reshape2::melt( id.vars= c("x","y","z"))#plot the gene
  plot <- ggplot(forplot, aes(x = x, y = z, color = value)) +
    ggtitle(title)+ labs( col=legendname) +
    geom_point(size = psize)+
    scale_colour_viridis_c(option=option,begin=begin,end=end,direction=-1)+
    facet_wrap(~variable,nrow=nrow) + theme_classic()+
    theme(strip.text = element_text(margin= margin(1,0),size = 12, color="black",face="bold"), 
          strip.background = element_rect(fill = "light grey", color="light grey"),
          aspect.ratio = 0.5,
          axis.line.x=element_line(color = "light grey"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "light grey"),
          panel.spacing = unit(0,"lines")
          
          #legend.position = "none"

    )
  return(plot)
}



##############  translate fbgn gene names both directions #######
translateG <- function(genes, dict){
  col_fbgn <- grep("FBgn",dict[1,])
  col_n <-1
  while (is.character(dict[1,col_n]) && grepl("FBgn",dict[1,col_n])){
    col_n <- col_n+1
  }
  
  
  if (grepl("FBgn",as.vector(genes)[1])){
    trans_genes <- dict[[col_n]][match(genes,dict[[col_fbgn]])]
  } else {
    trans_genes <- dict[[col_fbgn]][match(genes, dict[[col_n]])]
  }
  return(trans_genes)
}
  
  
  
#translate signs from glmnet to division and multiplication 

translateSigns<- function(x){
  if(x=="POS"){
    '*'
  } else if (x=="NEG"){
    '/'
  }
}
