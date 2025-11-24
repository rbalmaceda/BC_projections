options(java.parameters = "-Xmx32g")

rm(list=ls()) ; graphics.off() ; gc()
library(loadeR)
library(visualizeR)
library(downscaleR)
library(climate4R.value)
library(loadeR.2nc)
library(climate4R.indices)
require(lattice)
require(abind)


###### pruebas BC con otros periodos de estudio ######
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/"

setwd(path_1)
vars=c('tasmax','hus')
GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
rcms=c("REMO","RegCM",'ETA')
# season=c(12,1,2)
# set.seed(123)

var=vars[1]

experiment="Train"
years=1981:2004

if(var=='hus'){
  
  obs=get(load("Datos_CV_BC_huss_MSWX.rda"))
  
}

if(var=='tasmax'){
  
  obs=get(load("Datos_CV_BC_Tmax_MSWX.rda"))
  
}

obs=subsetGrid(obs,years=years,season=c(12,1,2))#2166 6 bisciestos
source("source_BC_update_dqm.R")
methods=c('eqm','dqm','qdm')

########## ETA########################
GCM=c("MIROC5","CanESM2","HadGEM")
setwd(path_2)

for (j in 1:3){
GCM_u=GCM[j]
message(paste0('GCMs.. ',GCM_u))
files_Datos=list.files(path=paste0(path_2,'ETA/'),pattern  = paste0('Datos_',var,"s_ETA_RCP8.5_",GCM_u,".rda"))
SESA_fut=get(load(paste0(path_2,'ETA/',files_Datos)))
#rm(data_SESA_fut)

files_Datos=list.files(path=path_1,pattern  = paste0('Datos.*.',var,".*.ETA.*.",GCM_u))

x=get(load(paste0(path_1,files_Datos)))
x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')
rm(data_SESA)  
gc()
  
for (s in c(1,2,3)){
method=methods[s]
message(paste0('method...',method))

ETA.BC <- biasCorrection_update(y=obs_1, x=x,newdata = SESA_fut, precipitation = FALSE,
                                   method=method, extrapolation = "constant", 
                                   n.quantiles=100)
  gc()
  save(ETA.BC,file=paste0("ETA/BC/",var,"_ETA_",method,"_",GCM_u,"_RCP8.5_SESA.rda"))
  rm(ETA.BC);gc()
}
  }


# ########## REMO ######################
GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
setwd(path_2)

for (j in 1:3){
  GCM_u=GCM[j]
  message(paste0('GCMs.. ',GCM_u))
  
##futuro
# files_Datos=list.files(path=path_2,pattern  = paste0('Datos_',var,"s_REMO_RCP8.5_",GCM_u,".rda"))
files_Datos=list.files(path=path_2,pattern  = paste0('Datos_',var,"_REMO_RCP8.5_",GCM_u,".rda"))

SESA_fut=get(load(paste0(path_2,files_Datos)))
  
  
###historico   
files_Datos=list.files(path=path_1,pattern  = paste0('Datos.*.',var,".*.REMO.*.",GCM_u))
x=get(load(paste0(path_1,files_Datos)))
x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')

for (k in c(1:3)){
  
  method=methods[k]
  message(paste0('method ',method))
  
  REMO.BC <- biasCorrection_update(y=obs_1, x=x,newdata = SESA_fut, precipitation = FALSE,
                                  method=method, extrapolation = "constant", 
                                  n.quantiles=100)
  gc()
  save(REMO.BC,file=paste0("REMO/BC_",var,"_REMO_",method,"_",GCM_u,"_RCP8.5_SESA.rda"))
  rm(REMO.BC);gc()  
}
rm(SESA_fut,data_SESA_fut)
}

####### RegCM #######################

GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
setwd(path_2)
k=2
methods=c('eqm','dqm','qdm')

for (j in 1:3){
  GCM_u=GCM[j]
  message(paste0('GCMs.. ',GCM_u))
  
  ##futuro
  # files_Datos=list.files(path=path_2,pattern  = paste0('Datos_',var,"s_RegCM_RCP8.5_",GCM_u,".rda"))
  files_Datos=list.files(path=paste0(path_2,rcms[k]),pattern  = paste0('Datos_',var,".*._RegCM_RCP8.5_",GCM_u,".rda"))
  
  SESA_fut=get(load(paste0(path_2,rcms[k],'/',files_Datos)))
  
  
  ###historico   
  files_Datos=list.files(path=path_1,pattern  = paste0('Datos.*.',var,".*.RegCM.*.",GCM_u))
  x=get(load(paste0(path_1,files_Datos)))
  x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
  obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')
  
  for (s in 1:3){
    
    method=methods[s]
    message(paste0('method ',method))
    
    RegCM.BC <- biasCorrection_update(y=obs_1, x=x,newdata = SESA_fut, precipitation = FALSE,
                                     method=method, extrapolation = "constant", 
                                     n.quantiles=100)
    gc()
    save(RegCM.BC,file=paste0("RegCM/BC_",var,"_RegCM_",method,"_",GCM_u,"_RCP8.5_SESA.rda"))
    rm(RegCM.BC);gc()  
  }
  rm(SESA_fut,data_SESA_fut)
}
