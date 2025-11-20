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
################ Datos #############
setwd("/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico")
# load("Datos_CV_BC_Tmax_RegCM_SESA_NorESM1.rda");load("Datos_CV_BC_Tmax_REMO_SESA_NorESM1.rda")
vars=c('tasmax','hus')
GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
rcms=c("REMO","RegCM",'ETA')
years=1981:2004
season=c(12,1,2)
set.seed(123)
pp=sample(years,size=24,replace = F)
# [1] 1995 1999 1994 1983 1990 1998 1991 1985 1984 2002 1986 1989 2000 2003 1997 2001 1993 1981 1996 1987
# [21] 1992 2004 1988 1982

experiment="CrosVal"
folds=list(pp[1:8],pp[9:16],pp[17:24])

print(folds)

var=vars[2]

if(var=='hus'){
  
  obs=get(load("Datos_CV_BC_huss_MSWX.rda"))
  
}

if(var=='tasmax'){
  
  obs=get(load("Datos_CV_BC_Tmax_MSWX.rda"))
  
}

obs=subsetGrid(obs,years=1981:2004,season=c(12,1,2))#2166 6 bisciestos
source("source_BC_update_dqm.R")
methods=c('eqm','dqm','qdm')
k=2
nn=100 #Al cambiar la funcion de BC esto vale para todos los metodos
# nn=100
for (j in 1:3){
  GCM_u=GCM[j]
  message(paste0('GCMs ',GCM_u))
  
################################################################################

##########################
# *** Bias correction ****
##########################

### Hago una prueba para ver que esta haciendo el BC, aplico y entreno en el mismo periodo ####

# files_Datos=list.files(pattern  = paste0('Datos.*.',var,".*.REMO.*.",GCM_u))
# x=get(load(files_Datos))
# x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
# obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')

# REMO.eqm.training <- biasCorrection(y=obs_1, x=x,newdata = x, precipitation = FALSE,
#                            method="eqm", extrapolation = "constant", n.quantiles=99 ,
#                            consecutive = F)
# 
# bias.P98 <- valueMeasure(x=REMO.eqm.training, y=obs_1, measure.code =  "bias",index.code = "P98",return.NApercentage = F)
# 
# spatialPlot(bias.P98, backdrop.theme = "countries",
#             rev.colors = T, main= "bias eqm P98 tasmax",  colorkey = list(space = "left", height = 0.8,width=0.5),
#             ylim=c(-45,-16),#at= seq(-1.8,1.8,0.1),
#             as.table=TRUE, xlim=c(-70,-40),strip = strip.custom(bg='grey90'),cex=0.3)#, scales = list(draw = TRUE)
# 
# # calcular el sesgo de ambos en el número de días por encima de 35 grados, año a año
# obs.P90 <- indexGrid(tx = obs_1,  index.code = "P",value=30)
# 
# bc.eqm.P90 <- indexGrid(tx = REMO.eqm.training,  index.code = "P")
# 
# bias.eqm.P90 <- gridArithmetics(climatology(bc.eqm.P90), climatology(obs.P90), operator="-")

#########scaling########
# RegCM.scaling <- biasCorrection(y=tasmax_MSWX.masked, x=Tmax_RegCM_SESA, newdata= Tmax_RegCM_SESA, precipitation = FALSE, 
#                              method="scaling", scaling.type="additive",cross.val = "kfold",
#                              folds = folds) # correción solo de la media (aditiva para temperatura)
# 
# save(RegCM.scaling,file=paste0(var,"_RegCM_",GCM_u,experiment,"_SESA.rda"))

# REMO.scaling <- biasCorrection(y=tasmax_MSWX.masked, x=Tmax_REMO_SESA, precipitation = FALSE, 
#                                 method="scaling", scaling.type="additive",cross.val = "kfold",
#                                 folds = folds,consecutive = F) # correción solo de la media (aditiva para temperatura)
# 
# 
# 
# 
# save(REMO.scaling,file=paste0(var,"_REMO_",GCM_u,experiment,"_SESA.rda"))

#Aca empieza la Cross Validacion 
##################### EQM ####################################################3
### En el caso de los eqm,y qdm y dqm n.quantile es la forma en que parten los percentiles, por eso uso
  # n.quatile 100 para que queden particiones de percentiles 0.01,0.2,...98
  
# ########## RegCM ######################

files_Datos=list.files(pattern  = paste0('Datos.*.',var,".*.RegCM.*.",GCM_u))
x=get(load(files_Datos))
#Esto ya lo hace el biascorrection function tambien
x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')

for (k in c(1,3)){
  
  method=methods[k]
  message(paste0('method ',method))
  
RegCM.eqm <- biasCorrection_update(y=obs_1, x=x, precipitation = FALSE,
                         method=method, extrapolation = "constant", n.quantiles=nn,cross.val = "kfold",
                         folds = folds,consecutive = F)

save(RegCM.eqm,file=paste0(var,"_RegCM_",method,"_",GCM_u,experiment,"_SESA.rda"))

}
# ########## REMO ######################
files_Datos=list.files(pattern  = paste0('Datos.*.',var,".*.REMO.*.",GCM_u))
x=get(load(files_Datos))
x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')

for (k in c(1,3)){
  
  method=methods[k]
  message(paste0('method ',method))

REMO.eqm <- biasCorrection_update(y=obs_1, x=x, precipitation = FALSE,
                         method=method, extrapolation = "constant", n.quantiles=nn ,cross.val = "kfold",
                         folds = folds,consecutive = F)

save(REMO.eqm,file=paste0(var,"_REMO_",method,"_",GCM_u,experiment,"_SESA.rda"))
}
 }

################## ETA ###########################3
#Lo hago por separado porque tiene otros GCMs

GCMs=c("MIROC5","CanESM2","HadGEM")

for (j in 1:3){
GCM_u=GCMs[j]

files_Datos=list.files(pattern  = paste0('Datos.*.',var,".*.ETA.*.",GCM_u))
x=get(load(files_Datos))
x=getTemporalIntersection(obs =obs,prd=x,which.return = 'prd')
obs_1=getTemporalIntersection(obs =obs,prd=x,which.return = 'obs')

for (k in c(1,3)){
  
  method=methods[k]
  message(paste0('method ',method))

  ETA.eqm <- biasCorrection_update(y=obs_1, x=x, precipitation = FALSE,
                           method=method, extrapolation = "constant", n.quantiles=nn ,cross.val = "kfold",
                           folds = folds,consecutive = F)

save(ETA.eqm,file=paste0(var,"_ETA_",method,"_",GCM_u,experiment,"_SESA.rda"))
}
}
