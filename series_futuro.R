options(java.parameters = "-Xmx32g")

library(loadeR)
library(visualizeR)
library(downscaleR)
library(climate4R.value)
library(loadeR.2nc)
library(climate4R.indices)
require(lattice)
require(abind)
require(dplyr)
library(purrr)

########################## Series temporales ##############################
path_0="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/Pruebas_CV/Train/"
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/"
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"

vars=c('tasmax','hus')
GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
rcms=c("REMO","RegCM",'ETA')
methods=c('eqm','dqm','qdm',"mbcr","mbcn")
name=c("ETA_HadGEM2","ETA_MIROC5", "ETA_CanESM2",           
       "RegCM_HadGEM2",  "RegCM_MPI-M-MPI", "RegCM_NorESM1",   
       "REMO_HadGEM2",   "REMO_MPI-M-MPI",  "REMO_NorESM1") 
index='mean'
####observacion #####
setwd(path_1)

obs=get(load(paste0(path_1,"/Datos_CV_BC_Tmax_MSWX.rda")))%>% subsetGrid(years=1981:2004,season=c(12,1,2))
  
obs=aggregateGrid(obs,aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                season = c(12, 1, 2)))
  
tasmax_obs_esc=scaleGrid(grid=obs,base = subsetGrid(obs,years = 1984:2003))
  
####### RAW #########
for (k in 1:3) {
message("Processing RCM: ", RCMs[k])

files_train_raw=list.files(path=path_1,pattern=paste0('Datos.*.','tasmax'))
files_train_raw=files_train_raw[c(2,3,1,4:9)]

tasmax_raw_train <- files_train_raw %>% 
  map(~ get(load(.x)) %>% subsetGrid(years=1981:2004,season=c(12,1,2)))
      
tasmax_raw_train <-lapply(1:length(files_train_raw), function(i){

pp=aggregateGrid(tasmax_raw_train[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                    season = c(12, 1, 2)))
}) %>% bindGrid(dimension = "member")

# tasmax_raw_train <- subsetDimension(tasmax_raw_train,dimension = 'member',indices=1:3)
# tasmax_raw_train$Members=name[1:3]
# tasmax_raw_train <- subsetDimension(tasmax_raw_train,dimension = 'member',indices=7:9)
# tasmax_raw_train$Members=name[7:9]
tasmax_raw_train <- subsetDimension(tasmax_raw_train,dimension = 'member',indices=4:6)
tasmax_raw_train$Members=name[4:6]

tasmax_raw_train_esc=scaleGrid(grid=tasmax_raw_train, base = subsetGrid(tasmax_raw_train,years = 1984:2003),by.member = T)

##### futuro ########

files_fut_raw=list.files(path=paste0(path_2,rcms[k],'/'),pattern  = paste0("Datos_tasmax_.*.RCP8.5_.*.rda"))
# files_fut_raw=files_fut_raw[c(2,3,1)]

setwd(paste0(path_2,rcms[k],'/'))

tasmax_raw_fut= files_fut_raw %>% 
  map(~ get(load(.x))) 

tasmax_raw_fut <-lapply(1:length(files_fut_raw), function(i){
  
  pp=aggregateGrid(tasmax_raw_fut[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                       season = c(12, 1, 2))) %>% 
    subsetGrid(years = 2006:2098)
    
}) %>% bindGrid(dimension = "member")

# tasmax_raw_fut$Members=name[7:9]
# tasmax_raw_fut$Members=name[1:3]
 tasmax_raw_fut$Members=name[4:6]

tasmax_raw_fut_esc=scaleGrid(grid=tasmax_raw_fut, base = subsetGrid(tasmax_raw_train,years = 1984:2003),by.member = T)

raw_serie_list=list(tasmax_raw_train,tasmax_raw_fut,tasmax_raw_fut_esc,tasmax_raw_train_esc)
names(raw_serie_list)=c("tasmax_raw_train","tasmax_raw_fut","tasmax_raw_fut_esc","tasmax_raw_train_esc")
assign(raw_serie_list,x=paste0("raw_serie_list_",rcms[k]))
save(raw_serie_list,file=paste0("tasmax_raw_serie_list_",rcms[k],".rda"))

####################### BC #####################################
methods=c('eqm','dqm','qdm','mbcn','mbcr')

for (j in 1:3){ #length(methods)
method=methods[j]
message(paste0('method..',method))

###historico ########
if(j %in% c(4,5)){
  
files_train=list.files(path=path_0,pattern=paste0("huss_.*._",method,"_Train_SESA.rda"))
# files_train=files_train[c(2,3,1)]
# files_train=files_train[c(7,8,9)]
files_train=files_train[c(4,5,6)]

setwd(path_0)
tasmax_BC_train <- files_train  %>%
map(~ get(load(.x)))


tasmax_BC_train <-lapply(1:length(files_train), function(i){
  
  pp=subsetGrid(tasmax_BC_train[[i]] ,var = 'tasmax')
  
  pp=aggregateGrid(pp,aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                      season = c(12, 1, 2)))
}) %>% bindGrid(dimension = "member")


} else {
files_train=list.files(path=path_0,pattern=paste0("tasmax_.*._",method,"_.*.Train","_SESA.rda"))
# files_train=files_train[c(2,3,1)]
files_train=files_train[c(7,8,9)]
# files_train=files_train[c(4,5,6)]

setwd(path_0)
tasmax_BC_train <- files_train %>%
  map(~ get(load(.x)))

tasmax_BC_train <-lapply(1:length(files_train), function(i){
  
  pp=aggregateGrid(tasmax_BC_train[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                     season = c(12, 1, 2)))
}) %>% bindGrid(dimension = "member")

}
# tasmax_BC_train$Members=name[1:3]
tasmax_BC_train$Members=name[7:9]
# tasmax_BC_train$Members=name[4:6]

#### Climatology 1986-2005 

tasmax_BC_train_esc=scaleGrid(grid=tasmax_BC_train,base = subsetGrid(tasmax_BC_train,years = 1984:2003),by.member = T)

###futuro ###
if (j %in% c(4,5)){
  
  files_fut=list.files(path=paste0(path_2,rcms[k],"/BC/"),paste0("hus.*.",method,"_RCP8.5_SESA.rda"))
  # files_fut=files_fut[c(2,3,1)]
  
  setwd(paste0(path_2,rcms[k],"/BC/"))
  
  tasmax_BC_fut <- files_fut %>%
    map(~ get(load(.x)))
  
  tasmax_BC_fut <-lapply(1:length(files_fut), function(i){
    pp=subsetGrid(tasmax_BC_fut[[i]] ,var = 'tasmax')
    pp=aggregateGrid(pp,aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                      season = c(12, 1, 2)))%>% 
      subsetGrid(years = 2006:2098)
  }) %>% bindGrid(dimension = "member")
  
  
  
}else {

files_fut=list.files(path=paste0(path_2,rcms[k],"/BC/"),paste0("tasmax_.*.",method,"_.*._RCP8.5_SESA.rda"))
# files_fut=files_fut[c(2,3,1)]

setwd(paste0(path_2,rcms[k],"/BC/"))

tasmax_BC_fut <- files_fut %>%
  map(~ get(load(.x)))

tasmax_BC_fut <-lapply(1:length(files_fut), function(i){
  
  pp=aggregateGrid(tasmax_BC_fut[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                      season = c(12, 1, 2)))%>% 
    subsetGrid(years = 2006:2098)
}) %>% bindGrid(dimension = "member")

}
# tasmax_BC_fut$Members=name[1:3]
tasmax_BC_fut$Members=name[7:9]
# tasmax_BC_fut$Members=name[4:6]

gc()

tasmax_BC_fut_esc=scaleGrid(grid=tasmax_BC_fut,base = subsetGrid(tasmax_BC_train,years = 1984:2003),by.member = T)

setwd(path_2)
BC_serie_list=list(tasmax_BC_train,tasmax_BC_train_esc,tasmax_BC_fut,tasmax_BC_fut_esc)
names(BC_serie_list)=c("tasmax_BC_train","tasmax_BC_train_esc","tasmax_BC_fut","tasmax_BC_fut_esc")
assign(BC_serie_list,x=paste0("BC_serie_list_",method))


rm(tasmax_BC_train,tasmax_BC_fut_esc,tasmax_BC_fut,tasmax_BC_train_esc)
rm(files_fut,files_train)
gc()

}

save(list = ls(pattern = "^BC_serie_list_"), file = paste0(rcms[k],"_tasmax_BC_serie_list.RData"))

}


####### grafico #########
#solo EQM
temporalPlot("OBS"=tasmax_obs_esc,#"ETA_BC_RCP8.5"=BC_serie_list_eqm$tasmax_BC_fut_esc,"ETA_hist"=raw_serie_list$tasmax_raw_train_esc,ETA_BC_hist=BC_serie_list_eqm$tasmax_BC_train_esc,
             ETA_RAW_RCP8.5=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))

temporalPlot("OBS"=tasmax_obs_esc,"REMO_BC_RCP8.5"=BC_serie_list_eqm$tasmax_BC_fut_esc,"REMO_hist"=raw_serie_list$tasmax_raw_train_esc,
             REMO_BC_hist=BC_serie_list_eqm$tasmax_BC_train_esc,
             REMO_RAW_RCP8.5=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))

temporalPlot("OBS"=tasmax_obs_esc,"RegCM_BC_RCP8.5"=BC_serie_list_eqm$tasmax_BC_fut_esc,"RegCM_hist"=raw_serie_list$tasmax_raw_train_esc,
             RegCM_BC_hist=BC_serie_list_eqm$tasmax_BC_train_esc,
             RegCM_RAW_RCP8.5=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))


temporalPlot("OBS"=tasmax_obs_esc,"ETA_BC_RCP8.5"=BC_serie_list_dqm$tasmax_BC_fut_esc,"ETA_hist"=raw_serie_list$tasmax_raw_train_esc,ETA_BC_hist=BC_serie_list_dqm$tasmax_BC_train_esc,
             ETA_RAW_RCP8.5=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))


temporalPlot("OBS"=tasmax_obs_esc,"ETA_QDM_RCP8.5"=BC_serie_list_qdm$tasmax_BC_fut_esc,
             "ETA_DQM_RCP8.5"=BC_serie_list_dqm$tasmax_BC_fut_esc,
             "ETA_EQM_RCP8.5"=BC_serie_list_eqm$tasmax_BC_fut_esc,
             "ETA_RAW_RCP8.5"=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))

temporalPlot("OBS"=tasmax_obs_esc,"REMO_QDM_RCP8.5"=BC_serie_list_qdm$tasmax_BC_fut_esc,
             "REMO_DQM_RCP8.5"=BC_serie_list_dqm$tasmax_BC_fut_esc,
             "REMO_EQM_RCP8.5"=BC_serie_list_eqm$tasmax_BC_fut_esc,
             "REMO_RAW_RCP8.5"=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))


temporalPlot("OBS"=tasmax_obs_esc,"RegCM_QDM_RCP8.5"=BC_serie_list_qdm$tasmax_BC_fut_esc,
             "RegCM_DQM_RCP8.5"=BC_serie_list_dqm$tasmax_BC_fut_esc,
             "RegCM_EQM_RCP8.5"=BC_serie_list_eqm$tasmax_BC_fut_esc,
             "RegCM_RAW_RCP8.5"=raw_serie_list$tasmax_raw_fut_esc,cols = c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))

####una run en particular 

temporalPlot("OBS"=tasmax_obs_esc,
             "ETA_BC_eqm"=subsetDimension(BC_serie_list_eqm$tasmax_BC_fut_esc,indices = 3,dimension = 'member'),
             "ETA_BC_dqm"=subsetDimension(BC_serie_list_dqm$tasmax_BC_fut_esc,indices = 3,dimension = 'member'),
             "ETA_BC_qdm"=subsetDimension(BC_serie_list_qdm$tasmax_BC_fut_esc,indices = 3,dimension = 'member'),
             "ETA_RAW"=subsetDimension(tasmax_raw_fut_esc,indices = 3,dimension = 'member'),
             cols=c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))


temporalPlot("OBS"=tasmax_obs_esc,
             "REMO_BC_eqm"=BC_serie_list_eqm$tasmax_BC_fut_esc,
             "REMO_BC_dqm"=BC_serie_list_dqm$tasmax_BC_fut_esc,
             "REMO_BC_qdm"=BC_serie_list_qdm$tasmax_BC_fut_esc,
             "REMO_RAW"=raw_serie_list$tasmax_raw_fut_esc,
             cols=c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))

temporalPlot("OBS"=tasmax_obs_esc,
             "RegCM_BC_eqm"=BC_serie_list_eqm$tasmax_BC_fut_esc,
             "RegCM_BC_dqm"=BC_serie_list_dqm$tasmax_BC_fut_esc,
             "RegCM_BC_qdm"=BC_serie_list_qdm$tasmax_BC_fut_esc,
             "RegCM_RAW"=raw_serie_list$tasmax_raw_fut_esc,
             cols=c('black','red','blue','green','orange'),xyplot.custom = list(ylim=c(-2,10)))

##tomo el punto 
coordenadas=c(110,89)
tasmax_obs_esc$xyCoords$x[89]
tasmax_obs_esc$xyCoords$y[110]

temporalPlot("OBS"=subsetDimension(subsetDimension(tasmax_obs_esc,indices = c(25),dimension = 'lon'),indice=21,dimension = 'lat'),
             "ETA_fut"=subsetDimension(subsetDimension(tasmax_raw_fut_esc,indices = c(25),dimension = 'lon'),indice=21,dimension = 'lat'),
             "ETA_EQM_fut"=subsetDimension(subsetDimension(BC_serie_list_eqm$tasmax_BC_fut_esc,indices = c(25),dimension = 'lon'),indice=21,dimension = 'lat'),
             "ETA_QDM_fut"=subsetDimension(subsetDimension(BC_serie_list_qdm$tasmax_BC_fut_esc,indices = c(25),dimension = 'lon'),indice=21,dimension = 'lat'),
            "ETA_mbcn_fut"=subsetDimension(subsetDimension(BC_serie_list_mbcn$tasmax_BC_fut_esc,indices = c(25),dimension = 'lon'),indice=21,dimension = 'lat'),
            cols=c('black','red','blue','green','orange'))


#########################

setwd(path_2)
list=ls(pattern = "^BC_serie_list_")
list_years=list("HadGEM2"=years_HadGEM2,"MIROC5"=years_MIROC5,"CanESM2"=years_CanESM2)

