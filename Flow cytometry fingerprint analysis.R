# Author: Niklas Gador, Kristianstad University Sweden

library(flowCore)
library(lattice)
library(flowViz)
library(png)
library(MASS)
library(matlib)
dir.create("images")

myColors<-gray.colors(255)
myColors<-rev(myColors)
myColRamp<-colorRampPalette(myColors)

gate_x<-c(1,1,5.5,2.9)
gate_y<-c(3.1,6,6,3.1)

binStorlek <- 64
filval<-"TB1"

normalisera_inom_gate<-"F"

volume_normalize<-"T"

datafolder<-"/All"

yeardata<-matrix(numeric(36),nrow = 6,ncol = 6,
                 dimnames = list(c("cluster1_Meanconf_int_x_min",
                                   "cluster1_Meanconf_int_x_max",
                                   "cluster2_Meanconf_int_x_min",
                                   "cluster2_Meanconf_int_x_max",
                                   "cluster_separation",
                                   "cluster_mean_separation P-value"
                 ),
                 c("18","19","20","21","22","23"))) # add years

cluster1<-c("SP1","SP2","SP3","SP4")  # add sampling points here, adjust the code below !
cluster2<-c("SP5","SP6","SP7","SP8")
number_of_viz_boot_elipses<-3
N<-100 #number of bootstraps for calculations

# Start 1. ## testplot an fcs file 2D histogram ####################################################################################
setwd(paste(".",datafolder,sep=""))
frame<-read.FCS(paste(".","A01 SP1 191111 TB1 250 m.fcs",sep="/"),alter.names=TRUE,transformation=FALSE) #test plot of one sample, check if gate fits

one<-truncateTransform("truncate at 1",a=1)
frame2<-{transform(frame,transformList(c('FL3.A','FL1.A'), one) )}  #look name specififations in FCM file
frame3<-{transform(frame2,transformList(c('FL3.A','FL1.A'), log10) )}
pg1<-{polygonGate(cbind("FL3.A"=gate_x,"FL1.A"=gate_y))}
xyplot(FL3.A~FL1.A,frame3, smooth = FALSE, xbin = binStorlek , filter=pg1, colramp =  myColRamp, par.settings = list(panel.background=list(col = "white")))
setwd("..")
# End 1. ###########################################################################################################################



# Start 2. ############################ Make png files #############################################################################
setwd(paste(".",datafolder,sep=""))
namevector<-list.files(".", full=F, pattern = "*.fcs") ## inspect all file names in Console !

## namevector volumes, search for 4th space from left,select number of rest to the right of 4th space #####################
namevector2<-substr(namevector,unlist(gregexpr(pattern="\\ ",namevector))[4]+1,nchar(namevector))
x <- unlist(regmatches(namevector2, gregexpr('\\(?[0-9]+', namevector2))) #only select integer
x <- as.numeric(gsub('\\(', '-', gsub(',', '', x)))
namevector_volumes<- x
namevector_volumes_scale<-namevector_volumes/min(namevector_volumes)


# hitta minst antal rader inom gaten av alla fcs filer ###############
fcs_fil<-read.FCS(namevector[1], alter.names = TRUE)
one<-truncateTransform("truncate at 1",a=1)
fcs_fil2<-transform(fcs_fil,transformList(c('FL3.A','FL1.A'), one) )
fcs_fil3<-transform(fcs_fil2,transformList(c('FL3.A','FL1.A'), log10) )
pg1<-{polygonGate(cbind("FL3.A"=gate_x,"FL1.A"=gate_y))}
fcs_fil_filtered<-{Subset(fcs_fil3,pg1)}
rader<-nrow(fcs_fil_filtered)
for (j in 2:length(namevector)){
  # j<-2
  ntl<-read.FCS(namevector[j], alter.names = TRUE)
  one<-truncateTransform("truncate at 1",a=1)
  ntl2<-transform(ntl,transformList(c('FL3.A','FL1.A'), one) )
  ntl3<-transform(ntl2,transformList(c('FL3.A','FL1.A'), log10) )
  pg1<-{polygonGate(cbind("FL3.A"=gate_x,"FL1.A"=gate_y))}
  ntl_filtered<-{Subset(ntl3,pg1)}
  
  if (nrow(ntl_filtered)<rader){
    rader<-nrow(ntl_filtered)
  }
}
######################################################################
#length(namevector)
for(k in 1:length(namevector)){
  #k<-2
  
  if (unlist(gregexpr(pattern=filval,namevector[k]))>-1){
    
    frame<-read.FCS(namevector[k],alter.names=TRUE,transformation=FALSE)
    one<-truncateTransform("truncate at 1",a=1)
    frame2<-transform(frame,transformList(c('FL3.A','FL1.A'), one) )
    frame3<-transform(frame2,transformList(c('FL3.A','FL1.A'), log10) )
    pg1<-{polygonGate(cbind("FL3.A"=gate_x,"FL1.A"=gate_y))}
    framefiltered<-{Subset(frame3,pg1)}
    
    namn4<-namevector[k]
    
    ######### normalizing within gate #######################
    if(normalisera_inom_gate=="T"){
      sf <- sampleFilter(filterId="mySampleFilter", size=rader) ## randomly select 'rader' number of events without replacement
      result<-filter(framefiltered,sf)
      framefiltered<-Subset(framefiltered,result)
    }
    #########################################################
    
    
    ######### normalizing volumes #######################
    if(volume_normalize=="T"){
      
      rader<-as.integer(nrow(framefiltered)/namevector_volumes_scale[k])    
      
      sf <- sampleFilter(filterId="mySampleFilter", size=rader) ## randomly select 'rader' number of events without replacement
      result<-filter(framefiltered,sf)
      framefiltered<-Subset(framefiltered,result)
    }
    ###########################################
    
    
    
    
    setwd("..")
    setwd(paste(".","/images",sep=""))
    
    name<-paste(namn4,".png", sep="")
    png(name, width=300, height=300, res=30, units="px")
    print(xyplot(FL3.A~FL1.A,framefiltered , smooth = FALSE, xbin = binStorlek , colramp = myColRamp, par.settings = list(panel.background=list(col = "white")), scales=list(x=list(draw=FALSE),y=list(draw=FALSE)), xlab="", ylab=""))
    dev.off()
    
    setwd("..")
    setwd(paste(".",datafolder,sep=""))
    rm(frame,one,frame2,frame3,framefiltered)
    gc()
  }
}

setwd("..") 

# End 2. ###########################################################################################################################

### Paus until all png files are made in folder 'images' ##################

# Start 3. #######  make a data.frame from windowed png histogram images, with and without kategori factor column  #################

fonstersize<-20 

setwd(paste(".","/images",sep=""))
subsets <-list.files(path=paste(".",sep="/"),pattern="*.png")

for (m in 1:length(subsets)){
  
  kategori<-substr(subsets[m],unlist(gregexpr(pattern="\\ ",subsets[m]))[1]+1,unlist(gregexpr(pattern="\\ ",subsets[m]))[2]-1)
  
  bild1<-readPNG(subsets[m])
  bild1<-bild1[,,1]
  #bild1low<-matrix(1, nrow = 300, ncol = 300)
  bild1low<-matrix(1, nrow = 300/fonstersize, ncol = 300/fonstersize)
  for (i in 1:(300/fonstersize)){
    for (j in 1:(300/fonstersize)){
      bigpix<-0
      for (k in 0:(fonstersize-1)){
        for (l in 0:(fonstersize-1)){
          bigpix<-bigpix+bild1[fonstersize*i-fonstersize+1+k,fonstersize*j-fonstersize+1+l]
        }
      }
      bild1low[i,j]<-bigpix
    }
  }
  enrad<- as.vector(bild1low)
  bild1df<-as.data.frame(t(enrad))
  bild1df<-cbind(bild1df, "kategori"=kategori)
  
  if (m==1) {
    bilddf<-bild1df
  } else {
    bilddf<-rbind(bilddf,bild1df)
  }
}
#str(bilddf)
bilddf_nokategori<- subset (bilddf, select = -kategori)

pngbild<-readPNG(subsets[1])
pngbild<-pngbild[,,1]
pngbild<-t(pngbild)[,dim(pngbild)[1]:1]
image(pngbild)
grid(nx=300/fonstersize,ny=300/fonstersize)
kolumnnamn<-names(bilddf)
for (i in 1:(300/fonstersize)){
  for (j in 1:(300/fonstersize)){
    
    text(x=(fonstersize/300*j)-(fonstersize/600),y=1-(fonstersize/300*i)+(fonstersize/600),kolumnnamn[j+(300/fonstersize)*(i-1)],col=0,cex=0.3*(fonstersize/300)*15)
    
  }
}
#View(bilddf)
# End 3. ############################################################################################################################

# Start 4. ################# pca on bilddf_nokategori ###############################################################################

former<-subset (bilddf, select = kategori)
former3<-integer(dim(former)[1])+1
for (m in 1:dim(former)[1]){
  # m<-2
  if (former[m,1]=="SP1"){
    former3[m]<-15
  }else if(former[m,1]=="SP2" ){
    former3[m]<-17
  }else if(former[m,1]=="SP3" ){
    former3[m]<-19
  }else if(former[m,1]=="SP4" ){
    former3[m]<-18
  }else if(former[m,1]=="SP5" ){
    former3[m]<-8
  }else if(former[m,1]=="SP6" ){
    former3[m]<-0
  }else if(former[m,1]=="SP7" ){
    former3[m]<-2
  }else if(former[m,1]=="SP8" ){
    former3[m]<-1
  }
}

lab=subsets
for (i in 1:length(lab)){
  #i<-1  
  lab[i]<- paste(substr(lab[i],unlist(gregexpr(pattern="\\ ",subsets[m]))[1]+1,unlist(gregexpr(pattern="\\ ",subsets[m]))[2]-1) ,
                 substr(lab[i],unlist(gregexpr(pattern="\\ ",subsets[m]))[2]+1,unlist(gregexpr(pattern="\\ ",subsets[m]))[3]-1),
                 sep="")
}

datum1<-subsets
for (i in 1:length(datum1)){
  # i<-1  
  datum1[i]<-substr(datum1[i],((unlist(gregexpr(pattern="\\ ",datum1[i]))))[2]+1,((unlist(gregexpr(pattern="\\ ",datum1[i]))))[2]+2)
}
datum1<-as.numeric(datum1)
if(max(datum1)-min(datum1)>0){
  datum2<-1*(datum1-min(datum1))/((max(datum1))-min(datum1))} else{datum2<-datum1/max(datum1)}
datum3<-gray(datum2)

######## PCA calculation ###################
pca<-prcomp(bilddf_nokategori,scale=FALSE)

####### year lightup ################################
min(datum1)
max(datum1)

for(i in min(datum1):max(datum1)){
  #i=18
  datum4<-ifelse(datum1==i,4,8)
  
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  plot(pca$x[,1],pca$x[,2],pch=former3,cex=1,col=datum4,main=paste("All",as.character(i)), xlab=paste("PC1 - ",pca.var.per[1],"%",sep=""),ylab=paste("PC2 - ",pca.var.per[2],"%",sep=""))
  #plot(pca$x[,1],pca$x[,2],xlim=c(-100,100),ylim=c(-100,100),pch=former3,cex=1,col=datum3,main="All", xlab="PC1_80%",ylab="PC2_10%")
  legend("bottomright",c("SP1","SP2","SP3","SP4","SP5","SP6","SP7","SP8"),pch=c(15,17,19,18,8,0,2,1))
  #text(pca$x,labels=lab, cex=0.3,pos=1,offset =0.2 )  ## ads datapoint labels
}
########################################################


########################## PCA plot ######################
par(bg ='white')
plot(pca$x[,1],pca$x[,2],pch=former3,cex=1,col=datum3,main="All", xlab=paste("PC1 - ",pca.var.per[1],"%",sep=""),ylab=paste("PC2 - ",pca.var.per[2],"%",sep=""))
#plot(pca$x[,1],pca$x[,2],xlim=c(-100,100),ylim=c(-100,100),pch=former3,cex=1,col=datum3,main="All", xlab="PC1_80%",ylab="PC2_10%")
legend("bottomright",c("SP1","SP2","SP3","SP4","SP5","SP6","SP7","SP8"),pch=c(15,17,19,18,8,0,2,1))
#text(pca$x,labels=lab, cex=0.3,pos=1,offset =0.2 )  ## ads datapoint labels


##########################################################

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#pca.var.per
sum(pca.var.per[1:2])
#barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

loading_scores <- pca$rotation[,1] #for PC1#
pixelwindow_scores <- abs(loading_scores) ## get the magnitudes
pixelwindow_score_ranked <- sort(pixelwindow_scores, decreasing=TRUE);pixelwindow_score_ranked[1:10]
top_10_pixelwindows <- names(pixelwindow_score_ranked[1:10])
top_10_pixelwindows ## show the names of the top 10 pixelwindows
write.csv(top_10_pixelwindows, file="top_10_pixelwindows_PC1.csv")

loading_scores <- pca$rotation[,2] #for PC2#
pixelwindow_scores <- abs(loading_scores) ## get the magnitudes
pixelwindow_score_ranked <- sort(pixelwindow_scores, decreasing=TRUE)
top_10_pixelwindows <- names(pixelwindow_score_ranked[1:10])
top_10_pixelwindows ## show the names of the top 10 pixelwindows


setwd("..")
# End 4. ############################################################################################################################

# Start 5. ##########################################################################################################################

pcadatan<-pca$x[,1:2]
### plats och ?r array ####
plats_ar<-subsets
for (i in 1:length(plats_ar)){
  # i<-1  
  plats_ar[i]<-substr(plats_ar[i],((unlist(gregexpr(pattern="\\ ",plats_ar[i]))))[1]+1,((unlist(gregexpr(pattern="\\ ",plats_ar[i]))))[2]+2)
}
#plats_ar
# filenames syncronized with pcadatan is 'subsets' #


for(year in as.integer(dimnames(yeardata)[[2]])){
#year=18
  
  mu_orig<-cbind(c(0,0),c(0,0))
  s2_orig<-cbind(c(0),c(0))
  dimnames(mu_orig)<-list(c("mu_x","mu_y" ),c("cluster1","cluster2" ))
  mues<-cbind(c(0:(N-1)),c(0:(N-1)),c(0:(N-1)),c(0:(N-1)))
  
  ####### lightup ?r year ################################
  
  datum4<-ifelse(datum1==year,4,8)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  plot(pca$x[,1],pca$x[,2],pch=former3,cex=1,col=datum4,main=paste("All",as.character(year)), xlab=paste("PC1 - ",pca.var.per[1],"%",sep=""),ylab=paste("PC2 - ",pca.var.per[2],"%",sep=""))
  
  ########################################################
  
  ### plocka ut platserna f?r ?r year, cluster1  ###
  pcadata_cluster1<-plats_ar==0  
  for(i in 1:length(cluster1)){
    # i<-2 
    pcadata_onepoint<-plats_ar==paste(cluster1[i],as.character(year)) 
    pcadata_cluster1<-pcadata_cluster1|pcadata_onepoint
  }  
  pcadatan_cluster1<-pcadatan[pcadata_cluster1,]
  N_cluster1<-dim(pcadatan_cluster1)[1]
  
  ### plocka ut platserna f?r ?r year, cluster2  ###
  pcadata_cluster2<-plats_ar==0  
  for(i in 1:length(cluster2)){
    # i<-2 
    pcadata_onepoint<-plats_ar==paste(cluster2[i],as.character(year)) 
    pcadata_cluster2<-pcadata_cluster2|pcadata_onepoint
  }  
  pcadatan_cluster2<-pcadatan[pcadata_cluster2,]
  N_cluster2<-dim(pcadatan_cluster2)[1]
  
  
  ## 2D Gaussian fit, distribution and 95% probabillity contour for both clusters  ##############################################################################
  for (cluster in 1:2){
    #cluster<-1
    
    Nmatrix<-100
    ifelse(cluster==1,data<-pcadatan_cluster1,data<-pcadatan_cluster2)
    #plot(data[,1],data[,2],xlim=c(-100,0),ylim=c(-20,20))
    ########################################################
    
    # sample mean #
    mu<-cbind(c(0,0))
    mu[1,1]<-(1/dim(data)[1])*sum(data[,1])
    mu[2,1]<-(1/dim(data)[1])*sum(data[,2])
    
    ## mean values of original cluster data ##########################
    ifelse(cluster==1,mu_orig[,"cluster1"]<-mu,mu_orig[,"cluster2"]<-mu)
    
    ## marginal x standard deviation square of original cluster data##
    s2_orig[cluster]<-((sum((data[,1]-mu[1,1])^2)/(dim(data)[1]-1))^0.5)^2 
    #################################  
    
    # sample standard deviations square#
    sigma<-cbind(c(0,0),c(0,0))       
    x_minus_mu<-data
    x_minus_mu[,1]<-data[,1]-matrix(mu[1,1],nrow=dim(data)[1],ncol=1)
    x_minus_mu[,2]<-data[,2]-matrix(mu[2,1],nrow=dim(data)[1],ncol=1)
    sigma<-(1/dim(data)[1])*t(x_minus_mu)%*%x_minus_mu # maximum likelihood sigma^2?
    
    range_x<-cbind(c(0),c(0))
    range_x[1,1]<-mu[1,1]-4*sigma[1,1]^0.5
    range_x[1,2]<-mu[1,1]+4*sigma[1,1]^0.5
    range_y<-cbind(c(0),c(0))
    range_y[1,1]<-mu[2,1]-4*sigma[2,2]^0.5
    range_y[1,2]<-mu[2,1]+4*sigma[2,2]^0.5
    ###################################################################################
    x <- seq(range_x[1,1], range_x[1,2], length= Nmatrix)
    y<-seq(range_y[1,1], range_y[1,2], length= Nmatrix)
    z<-matrix(1, nrow = Nmatrix, ncol = Nmatrix)
    for (i in 1:Nmatrix){
      for (j in 1:Nmatrix){
        xmat<-cbind(c(x[i],y[j])) 
        z[i,j]<- 1/(2*pi*sqrt(det(sigma)))*exp(as.numeric(-0.5*t(xmat-mu)%*%inv(sigma)%*%(xmat-mu)))
      } 
    }
    
    
    ######### centervalue of probabillity distribution ###########################
    xmat<-mu
    z0<- 1/(2*pi*sqrt(det(sigma)))*exp(as.numeric(-0.5*t(xmat-mu)%*%inv(sigma)%*%(xmat-mu)))
    ##############################################################################
    
    ##### finding roughly 95% probabillity contour      check total probabillity !!!                          #########
    ##### ( a new sample would end up inside ellipse with 95% probabillity ) #########
    
    k<-z0
    while (sum(z[z>k])*((range_x[1,2]-range_x[1,1])/Nmatrix)*((range_y[1,2]-range_y[1,1])/Nmatrix)<0.95)
    {
      k<- k-(z0/10000)
    }
    contour(x, y, z,xlab="x1",ylab="x2",nlevels=1,levels = c(k),drawlabels = FALSE,add = TRUE)
    
    
    
    
    
    ######################################################################################
    
    ## bootstrapping number_of_viz_boot_elipses ellipses ################################################################################################################### 
    
    
    for (m in 2:number_of_viz_boot_elipses)
    {
      data_ny<-data
      for (i in 1:dim(data)[1])
      {
        random<-floor(runif(1,min=1,dim(data)[1]+1))
        data_ny[i,]<-data[random,]
      }
      #plot(data_ny[,1],data_ny[,2],xlim=c(0,10),ylim=c(0,10))
      
      # sample mean #
      mu<-cbind(c(0,0))
      mu[1,1]<-(1/dim(data_ny)[1])*sum(data_ny[,1])
      mu[2,1]<-(1/dim(data_ny)[1])*sum(data_ny[,2])
      
      # sample standard deviations square#
      sigma<-cbind(c(0,0),c(0,0))       
      x_minus_mu<-data_ny
      x_minus_mu[,1]<-data_ny[,1]-matrix(mu[1,1],nrow=dim(data_ny)[1],ncol=1)
      x_minus_mu[,2]<-data_ny[,2]-matrix(mu[2,1],nrow=dim(data_ny)[1],ncol=1)
      
      sigma<-(1/dim(data_ny)[1])*t(x_minus_mu)%*%x_minus_mu # maximum likelihood sigma^2?
      
      ###################################################################################
      x <- seq(range_x[1,1], range_x[1,2], length= Nmatrix)
      y<-seq(range_y[1,1], range_y[1,2], length= Nmatrix)
      
      z<-matrix(1, nrow = Nmatrix, ncol = Nmatrix)
      
      for (i in 1:Nmatrix){
        for (j in 1:Nmatrix){
          xmat<-cbind(c(x[i],y[j])) 
          z[i,j]<- 1/(2*pi*sqrt(det(sigma)))*exp(as.numeric(-0.5*t(xmat-mu)%*%inv(sigma)%*%(xmat-mu)))
        } 
      }
      
      
      ######### centervalue of probabillity distribution ###########################
      xmat<-mu
      z0<- 1/(2*pi*sqrt(det(sigma)))*exp(as.numeric(-0.5*t(xmat-mu)%*%inv(sigma)%*%(xmat-mu)))
      ##############################################################################
      
      ##### finding roughly 95% probabillity contour                              #########
      k<-z0
      while (sum(z[z>k])*((range_x[1,2]-range_x[1,1])/Nmatrix)*((range_y[1,2]-range_y[1,1])/Nmatrix)<0.95)
      {
        k<- k-(z0/10000)
      }
      contour(x, y, z,xlab="x1",ylab="x2",nlevels=1,levels = c(k),drawlabels = FALSE,col=m,add = TRUE)
      
      
    }
    ## end bootstrapping number_of_viz_boot_elipses ellipses ######################### 
    
    
    ####  bootstrapped mean uncertainty #########
    # bootstrap N=100 times #
    
    for (n in 1:N)
      # n<-1  
    {
      
      data_ny<-data
      for (i in 1:dim(data)[1])
      {
        random<-floor(runif(1,min=1,dim(data)[1]+1))
        data_ny[i,]<-data[random,]
      }
      
      # sample mean #
      mu<-cbind(c(0),c(0))
      mu[1,1]<-(1/dim(data_ny)[1])*sum(data_ny[,1])
      mu[1,2]<-(1/dim(data_ny)[1])*sum(data_ny[,2])
      
      mues[n,(cluster*2-1):(cluster*2)]<-mu
      
    }  
    
    # error of means #
    
    mean_bar_x<-sum(mues[,(cluster*2-1)])/N
    mean_bar_y<-sum(mues[,(cluster*2)])/N
    
    meanstand_x<-(sum((mues[,(cluster*2-1)]-mean_bar_x)^2/(N-1)))^0.5
    meanstand_y<-(sum((mues[,(cluster*2)]-mean_bar_y)^2/(N-1)))^0.5
    
    error_mean_x<-2*meanstand_x ## about 95% certainty #
    error_mean_y<-2*meanstand_y
    
    Meanconf_int_x<-cbind(c(0),c(0))
    Meanconf_int_x[1,1]<-mu_orig[1,cluster]-error_mean_x
    Meanconf_int_x[1,2]<-mu_orig[1,cluster]+error_mean_x
    Meanconf_int_x
    Meanconf_int_y<-cbind(c(0),c(0))
    Meanconf_int_y[1,1]<-mu_orig[2,cluster]-error_mean_y
    Meanconf_int_y[1,2]<-mu_orig[2,cluster]+error_mean_y
    Meanconf_int_y
    
    yeardata[paste("cluster",as.character(cluster),"_Meanconf_int_x_min",sep=""),as.character(year)]<-Meanconf_int_x[1,1] # only in x (PC1)
    yeardata[paste("cluster",as.character(cluster),"_Meanconf_int_x_max",sep=""),as.character(year)]<-Meanconf_int_x[1,2]
    
    ####  end bootstrapped mean uncertainty #########
    
  }
  
  ###### Seperation in PC1 ,smallest mu difference from meanconf_int_x #####################################################################
  
  
  if (yeardata[2,as.character(year)]< yeardata[3,as.character(year)] | yeardata[4,as.character(year)]< yeardata[1,as.character(year)]   ){
    
    if(yeardata[2,as.character(year)]< yeardata[3,as.character(year)] ){
      mu1<-yeardata[2,as.character(year)]
      mu2<-yeardata[3,as.character(year)]
      separation<-(mu1-mu2)^2/(s2_orig[1]+s2_orig[2])
    }else{
      mu1<-yeardata[4,as.character(year)]
      mu2<-yeardata[1,as.character(year)]
      separation<-(mu1-mu2)^2/(s2_orig[1]+s2_orig[2])
    }
  }else{separation<-0}
  
  
  yeardata[5,as.character(year)]<-separation
  
  
  ###  t test for difference of means ###
  
  N_cluster1
  N_cluster2
  mu_orig[1,1] # in x #
  s2_orig[1] ## s square in x##
  mu_orig[1,2] # in x #
  s2_orig[2] ## s square in x##
  t<- abs((mu_orig[1,1]-mu_orig[1,2])/((s2_orig[1]/N_cluster1)+(s2_orig[2]/N_cluster2))^0.5)
  P<-2*pt(q=t, df=36, lower.tail = F) 
  
  yeardata[6,as.character(year)]<-P
  
}


# end 5. ###############################################################################################################
View(yeardata)
  
write.csv(yeardata, file="yeardata.csv")
par(bg ='white')
plot(as.integer(dimnames(yeardata)[[2]]),yeardata["cluster_separation",],main="separation",xlab="year",ylab="seperation")



