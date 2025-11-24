##### multivariate #####
options(java.parameters = "-Xmx32g")

rm(list=ls()) ; graphics.off() ; gc()
require(MBC)
library(loadeR)
library(visualizeR)
library(downscaleR)
library(climate4R.value)
library(loadeR.2nc)
library(climate4R.indices)
require(lattice)
require(abind)
# install.packages("MBC")
# install.packages("https://cran.r-project.org/src/contrib/Archive/gsl/gsl_2.1-6.tar.gz", repos = NULL, type = "source")
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico"

setwd(path_1)

ratio.seq <- c(FALSE, FALSE)
trace <- c(Inf,0)
pp.type=6
method='mbcr'
method='mbcn'
files2=list.files(pattern=paste0('Datos.*.','hus'))
files=list.files(pattern=paste0('Datos.*.','tasmax'))
files=files[c(2,3,1,4:9)]
files2=files2[c(2,4,3,5:10)]
name=c("ETA_HadGEM2","ETA_MIROC5", "ETA_CanESM2",           
       "RegCM_HadGEM2",  "RegCM_MPI-M-MPI", "RegCM_NorESM1",   
       "REMO_HadGEM2",   "REMO_MPI-M-MPI",  "REMO_NorESM1") 

years=1981:2004

experiment="Train"
years=1981:2004

obs.huss=get(load("Datos_CV_BC_huss_MSWX.rda"))
obs.huss=subsetGrid(obs.huss,years=1981:2004,season=c(12,1,2))#2166 6 bisciestos

obs.Tx=get(load("Datos_CV_BC_Tmax_MSWX.rda"))
obs.Tx=subsetGrid(obs.Tx,years=1981:2004,season=c(12,1,2))#2166 6 bisciestos

obs_1.huss=getTemporalIntersection('obs'=obs.huss,get(load(files2[4])),which.return = 'obs')
obs_1.Tx=getTemporalIntersection('obs'=obs.Tx,get(load(files[4])),which.return = 'obs')

#### Futuro ###
####ETA#####
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/ETA/"
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/"

setwd(path_2)

files2_fut=list.files(path=path_2,pattern  = paste0('Datos_huss_ETA_RCP8.5_.*.rda'))
files2_fut=files2_fut[c(2,3,1)]

files_fut=list.files(path=path_2,pattern  = paste0('Datos_tasmax_ETA_RCP8.5_.*.rda'))
files_fut=files_fut[c(2,3,1)]

##opc A ajusto huss
var='huss'
y.obs=list('Tx'=obs_1.Tx,'hus'=obs_1.huss)
y.obs[["Tx"]][["Variable"]][["varName"]]='tasmax'
y.obs[["hus"]][["Variable"]][["varName"]]='huss'


for (z in 1:length(files_fut)){##opc A
  
  Tx <-get(load(paste0(path_1,files[z])))
  hus <-get(load(paste0(path_1,files2[z])))
  
  Tx=getTemporalIntersection(obs =y.obs$Tx,prd=Tx,which.return = 'prd')
  hus=getTemporalIntersection(obs =y.obs$Tx,prd=hus,which.return = 'prd')
  x=list(Tx,hus)
  rm(Tx,hus);gc()
  y=y.obs; x=x
  
  Tx_fut <-get(load(paste0(path_2,files_fut[z])))
  hus_fut <-get(load(paste0(path_2,files2_fut[z])))
  hus_fut=getTemporalIntersection(obs =Tx_fut,prd=hus_fut,which.return = 'prd')
  Tx_fut=getTemporalIntersection(obs =Tx_fut,prd=hus_fut,which.return = 'obs')
  
  x_fut=list(Tx_fut,hus_fut)
  rm(Tx_fut,hus_fut);gc()
  
########################## BC ######################33
    oo <- redim(makeMultiGrid(y[[1]],y[[2]]), drop=TRUE) 
    pp <- redim(makeMultiGrid(x[[1]],x[[2]]), drop=TRUE) 
    ss<- redim(makeMultiGrid(x_fut[[1]],x_fut[[2]]), drop=TRUE)
    
##apply MBC correction
    newdata.mbc <- array(NA,dim=dim(ss$Data))
    
    for(i in 1:length(ss[["xyCoords"]][["y"]])){
      for(j in 1:length(ss[["xyCoords"]][["x"]])){
        message("BC point ", i, ", ", j)
        if(all(is.na(oo$Data[,,i,j]))){
          next
        } 
        
        o <- t(oo$Data[,,i,j])
        p <- t(pp$Data[,,i,j])
        s <- t(ss$Data[,,i,j])
        if(method== "mbcr"){
        fit.mbc <- MBCr(o.c=o, m.c=p,
                        m.p=s, ratio.seq=ratio.seq, trace=trace,pp.type=pp.type)}
        if (method== "mbcn"){
          fit.mbc <- MBCn(o.c=o, m.c=p,
                          m.p=s, ratio.seq=ratio.seq, trace=trace,pp.type=pp.type)
          
        }
        newdata.mbc[,,i,j]<-t(fit.mbc$mhat.p)
      }
    }
    
    fut.mbc <- ss
    fut.mbc$Data <- newdata.mbc
    attr(fut.mbc$Data,"dimensions")<- c("var","time" ,"lat", "lon")
  
  save(fut.mbc,file=paste0('BC/',var,"_",name[z],'_',method,"_RCP8.5_SESA.rda"))
  
}  


####REMO ####
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/"

setwd(path_2)

files2_fut=list.files(path=path_2,pattern  = paste0('Datos_huss_REMO_RCP8.5_.*.rda'))

files_fut=list.files(path=path_2,pattern  = paste0('Datos_tasmax_REMO_RCP8.5_.*.rda'))

##opc A ajusto huss
var='huss'
y.obs=list('Tx'=obs_1.Tx,'hus'=obs_1.huss)
y.obs[["Tx"]][["Variable"]][["varName"]]='tasmax'
y.obs[["hus"]][["Variable"]][["varName"]]='huss'


for (z in 1:length(files_fut)){##opc A
  
  Tx <-get(load(paste0(path_1,files[z+6])))
  hus <-get(load(paste0(path_1,files2[z+6])))
  
  Tx=getTemporalIntersection(obs =y.obs$Tx,prd=Tx,which.return = 'prd')
  hus=getTemporalIntersection(obs =y.obs$Tx,prd=hus,which.return = 'prd')
  x=list(Tx,hus)
  rm(Tx,hus);gc()
  y=y.obs; x=x
  
  Tx_fut <-get(load(paste0(path_2,files_fut[z])))
  hus_fut <-get(load(paste0(path_2,files2_fut[z])))
  hus_fut=getTemporalIntersection(obs =Tx_fut,prd=hus_fut,which.return = 'prd')
  Tx_fut=getTemporalIntersection(obs =Tx_fut,prd=hus_fut,which.return = 'obs')
  
  x_fut=list(Tx_fut,hus_fut)
  rm(Tx_fut,hus_fut);gc()
  
########################## BC ######################
  oo <- redim(makeMultiGrid(y[[1]],y[[2]]), drop=TRUE) 
  pp <- redim(makeMultiGrid(x[[1]],x[[2]]), drop=TRUE) 
  ss<- redim(makeMultiGrid(x_fut[[1]],x_fut[[2]]), drop=TRUE)
  
  #apply MBC correction
  newdata.mbc <- array(NA,dim=dim(ss$Data))
  
  for(i in 1:length(ss[["xyCoords"]][["y"]])){
    for(j in 1:length(ss[["xyCoords"]][["x"]])){
      message("BC point ", i, ", ", j)
      if(all(is.na(oo$Data[,,i,j]))){
        next
      } 
      
      o <- t(oo$Data[,,i,j])
      p <- t(pp$Data[,,i,j])
      s <- t(ss$Data[,,i,j])
      if(method== "mbcr"){
        fit.mbc <- MBCr(o.c=o, m.c=p,
                        m.p=s, ratio.seq=ratio.seq, trace=trace,pp.type=pp.type)}
      if (method== "mbcn"){
        fit.mbc <- MBCn(o.c=o, m.c=p,
                        m.p=s, ratio.seq=ratio.seq, trace=trace,pp.type=pp.type)
        
      }
      newdata.mbc[,,i,j]<-t(fit.mbc$mhat.p)
    }
  }
  
  fut.mbc <- ss
  fut.mbc$Data <- newdata.mbc
  attr(fut.mbc$Data,"dimensions")<- c("var","time" ,"lat", "lon")
  
  save(fut.mbc,file=paste0('REMO/BC_',var,"_",name[z+6],'_',method,"_RCP8.5_SESA.rda"))
  
}  
########################### RegCM ####################################

files2_fut=list.files(path=paste0(path_2,'RegCM'),pattern  = paste0('Datos_huss_RegCM_RCP8.5_.*.rda'))

files_fut=list.files(path=paste0(path_2,'RegCM'),pattern  = paste0('Datos_tasmax_RegCM_RCP8.5_.*.rda'))

##opc A ajusto huss
var='huss'
y.obs=list('Tx'=obs_1.Tx,'hus'=obs_1.huss)
y.obs[["Tx"]][["Variable"]][["varName"]]='tasmax'
y.obs[["hus"]][["Variable"]][["varName"]]='huss'


for (z in 1:length(files_fut)){##opc A
  
  Tx <-get(load(paste0(path_1,files[z+3])))
  hus <-get(load(paste0(path_1,files2[z+3])))
  
  Tx=getTemporalIntersection(obs =y.obs$Tx,prd=Tx,which.return = 'prd')
  hus=getTemporalIntersection(obs =y.obs$Tx,prd=hus,which.return = 'prd')
  x=list(Tx,hus)
  rm(Tx,hus);gc()
  y=y.obs; x=x
  
  Tx_fut <-get(load(paste0(path_2,'RegCM/',files_fut[z])))
  hus_fut <-get(load(paste0(path_2,'RegCM/',files2_fut[z])))
  hus_fut=getTemporalIntersection(obs =Tx_fut,prd=hus_fut,which.return = 'prd')
  Tx_fut=getTemporalIntersection(obs =Tx_fut,prd=hus_fut,which.return = 'obs')
  
  x_fut=list(Tx_fut,hus_fut)
  rm(Tx_fut,hus_fut);gc()
  
  ########################## BC ######################
  oo <- redim(makeMultiGrid(y[[1]],y[[2]]), drop=TRUE) 
  pp <- redim(makeMultiGrid(x[[1]],x[[2]]), drop=TRUE) 
  ss<- redim(makeMultiGrid(x_fut[[1]],x_fut[[2]]), drop=TRUE)
  
  #apply MBC correction
  newdata.mbc <- array(NA,dim=dim(ss$Data))
  
  for(i in 1:length(ss[["xyCoords"]][["y"]])){
    for(j in 1:length(ss[["xyCoords"]][["x"]])){
      message("BC point ", i, ", ", j)
      if(all(is.na(oo$Data[,,i,j]))){
        next
      } 
      
      o <- t(oo$Data[,,i,j])
      p <- t(pp$Data[,,i,j])
      s <- t(ss$Data[,,i,j])
      if(method== "mbcr"){
        fit.mbc <- MBCr(o.c=o, m.c=p,
                        m.p=s, ratio.seq=ratio.seq, trace=trace,pp.type=pp.type)}
      if (method== "mbcn"){
        fit.mbc <- MBCn(o.c=o, m.c=p,
                        m.p=s, ratio.seq=ratio.seq, trace=trace,pp.type=pp.type)
        
      }
      newdata.mbc[,,i,j]<-t(fit.mbc$mhat.p)
    }
  }
  
  fut.mbc <- ss
  fut.mbc$Data <- newdata.mbc
  attr(fut.mbc$Data,"dimensions")<- c("var","time" ,"lat", "lon")
  
  save(fut.mbc,file=paste0('RegCM/BC_',var,"_",name[z+3],'_',method,"_RCP8.5_SESA.rda"))
  
}  
