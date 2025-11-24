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

####observacion ########################################################33
setwd(path_1)

obs=get(load(paste0(path_1,"/Datos_CV_BC_huss_MSWX.rda")))%>% subsetGrid(years=1981:2004,season=c(12,1,2))%>% 
  aggregateGrid(aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                season = c(12, 1, 2)))
  
hus_obs_esc=scaleGrid(grid=obs,base = subsetGrid(obs,years = 1984:2003))
  
####### RAW ######################################################################
k=3#ETA
k=1#REMO
k=2#RegCM
files_train_raw=list.files(path=path_1,pattern=paste0('Datos.*.','hus'))
files_train_raw=files_train_raw[c(2,4,3,5:10)]
files_train_raw=files_train_raw[1:3]
# files_train_raw=files_train_raw[7:9]
# files_train_raw=files_train_raw[4:6]

hus_raw_train <- files_train_raw %>% 
  map(~ get(load(.x)) %>% subsetGrid(years=1981:2004,season=c(12,1,2)))
      
hus_raw_train <-lapply(1:length(files_train_raw), function(i){

pp=aggregateGrid(hus_raw_train[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                    season = c(12, 1, 2)))
}) %>% bindGrid(dimension = "member")

# hus_raw_train <- subsetDimension(hus_raw_train,dimension = 'member',indices=1:3)
hus_raw_train$Members=name[1:3]

# hus_raw_train <- subsetDimension(hus_raw_train,dimension = 'member',indices=7:9)
# hus_raw_train$Members=name[7:9]

# hus_raw_train$Members=name[4:6]

hus_raw_train_esc=scaleGrid(grid=hus_raw_train, base = subsetGrid(hus_raw_train,years = 1984:2003),by.member = T,type='ratio')

##### futuro #####
files_fut_raw=list.files(path=paste0(path_2,rcms[k],'/'),pattern  = paste0("Datos_hus.*.RCP8.5_.*.rda"))
if(k==3){files_fut_raw=files_fut_raw[c(2,3,1)]}

setwd(paste0(path_2,rcms[k],'/'))


hus_raw_fut= files_fut_raw %>% 
  map(~ get(load(.x))) 

hus_raw_fut <-lapply(1:length(files_fut_raw), function(i){
  
  pp=aggregateGrid(hus_raw_fut[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                       season = c(12, 1, 2)))%>% 
    subsetGrid(years = 2006:2098)
}) %>% bindGrid(dimension = "member")

hus_raw_fut$Members=name[1:3]
# hus_raw_fut$Members=name[7:9]
# hus_raw_fut$Members=name[4:6]

hus_raw_fut_esc=scaleGrid(grid=hus_raw_fut, base = subsetGrid(hus_raw_train,years = 1984:2003),by.member = T)

raw_serie_list=list(hus_raw_train,hus_raw_fut,hus_raw_train_esc,hus_raw_fut_esc)
names(raw_serie_list)=c("hus_raw_train","hus_raw_fut","hus_raw_train_esc","hus_raw_fut_esc")
assign(raw_serie_list,x=paste0("raw_serie_list_",rcms[k]))
save(raw_serie_list,file=paste0("huss_raw_serie_list_",rcms[k],".rda"))

####################### BC #####################################
methods=c('eqm','dqm','qdm','mbcn','mbcr')
k=3#ETA
k=1#REMO
k=2#RegCM
scale='ratio'

for (j in 1:length(methods)){ #
method=methods[j]
message(paste0('method..',method))

###historico ########
if(j %in% c(4,5)){
  
files_train=list.files(path=path_0,pattern=paste0("huss_.*._",method,"_Train_SESA.rda"))
# files_train=files_train[c(2,3,1)]
# files_train=files_train[c(7,8,9)]
files_train=files_train[c(4,5,6)]

setwd(path_0)
hus_BC_train <- files_train  %>%
map(~ get(load(.x)))


hus_BC_train <-lapply(1:length(files_train), function(i){
  
  pp=subsetGrid(hus_BC_train[[i]] ,var = hus_BC_train[[i]]$Variable$varName[2])
  
  pp=aggregateGrid(pp,aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                      season = c(12, 1, 2)))
}) %>% bindGrid(dimension = "member")


} else {
files_train=list.files(path=path_0,pattern=paste0("hus_.*._",method,"_.*.Train","_SESA.rda"))
files_train=files_train[c(2,3,1)]
# files_train=files_train[c(7,8,9)]
# files_train=files_train[c(4,5,6)]

# pp=aggregateGrid(get(load(paste0(path_0,files_train[z]))),aggr.s=list(FUN = list("mean", na.rm = TRUE),season=c(12,1,2)))
# tasmax_BC_train=get(load(files_train[z]))%>% aggregateGrid(aggr.s=list(FUN = list("mean", na.rm = TRUE),season=c(12,1,2))
#                                            %>% bindGrid(dimension = 'time')

setwd(path_0)
hus_BC_train <- files_train %>%
  map(~ get(load(.x)))

hus_BC_train <-lapply(1:length(files_train), function(i){
  
  pp=aggregateGrid(hus_BC_train[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                     season = c(12, 1, 2)))
}) %>% bindGrid(dimension = "member")

}
# hus_BC_train$Members=name[1:3]
# hus_BC_train$Members=name[7:9]
hus_BC_train$Members=name[4:6]

#### Climatology 1986-2005 ####

hus_BC_train_esc=scaleGrid(grid=hus_BC_train,base = subsetGrid(hus_BC_train,years = 1984:2003),by.member = T)

###futuro #######

if (j %in% c(4,5)){
  
  files_fut=list.files(path=paste0(path_2,rcms[k],"/BC/"),paste0("BC_hus.*.",method,"_RCP8.5_SESA.rda"))
  # files_fut=files_fut[c(2,3,1,4:9)]
  # files_fut=files_fut[c(2,3,1)]
  
  setwd(paste0(path_2,rcms[k],"/BC/"))
  
  hus_BC_fut <- files_fut %>%
    map(~ get(load(.x)))
  
  hus_BC_fut <-lapply(1:length(files_fut), function(i){
    pp=subsetGrid(hus_BC_fut[[i]] ,var = hus_BC_fut[[i]]$Variable$varName[2])
    pp=aggregateGrid(pp,aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                      season = c(12, 1, 2)))%>% 
      subsetGrid(years = 2006:2098)
  }) %>% bindGrid(dimension = "member")
  
  
  
}else {

files_fut=list.files(path=paste0(path_2,rcms[k],"/BC/"),paste0("hus_.*.",method,"_.*._RCP8.5_SESA.rda"))
files_fut=files_fut[c(2,3,1)]

setwd(paste0(path_2,rcms[k],"/BC/"))

hus_BC_fut <- files_fut %>%
  map(~ get(load(.x)))

hus_BC_fut <-lapply(1:length(files_fut), function(i){
  
  pp=aggregateGrid(hus_BC_fut[[i]],aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                                      season = c(12, 1, 2)))%>% 
    subsetGrid(years = 2006:2098)
}) %>% bindGrid(dimension = "member")

}
# hus_BC_fut$Members=name[1:3]
# hus_BC_fut$Members=name[7:9]
hus_BC_fut$Members=name[4:6]

gc()

hus_BC_fut_esc=scaleGrid(grid=hus_BC_fut,base = subsetGrid(hus_BC_train,years = 1984:2003),by.member = T)

if(scale=='ratio'){
  hus_BC_fut_esc=scaleGrid(grid=hus_BC_fut_esc, base = subsetGrid(hus_BC_train,years = 1984:2003),by.member = T,type='ratio')
}


BC_serie_list=list(hus_BC_train,hus_BC_train_esc,hus_BC_fut,hus_BC_fut_esc)
names(BC_serie_list)=c("hus_BC_train","hus_BC_train_esc","hus_BC_fut","hus_BC_fut_esc")
assign(BC_serie_list,x=paste0("BC_serie_list_",method))


rm(hus_BC_train,hus_BC_fut_esc,hus_BC_fut,hus_BC_train_esc)
rm(files_fut,files_train)
gc()

}

setwd(path_2)
save(list = ls(pattern = "^BC_serie_list_"), file = paste0(rcms[k],"_hus_BC_serie_list.RData"))
save(list = ls(pattern = "^BC_serie_list_"), file = paste0(rcms[k],"_hus_BC_serie_list_ratio.RData"))

####### grafico #########
load(paste0(path_2,"REMO_hus_BC_serie_list.RData"))
########################## Series temporales ##############################

#solo EQM
hus_obs_esc_1=hus_obs_esc
hus_obs_esc_1$Data=hus_obs_esc$Data*1000

png("temporal_plot.png", width = 1200, height = 600, res = 150)
temporalPlot("OBS"=gridArithmetics(hus_obs_esc,1000,operator = '*'),
             "ETA_RAW_RCP8.5"=gridArithmetics(raw_serie_list$hus_raw_fut_esc,1000,operator = '*'), 
             "ETA_BC_RCP8.5_eqm"=gridArithmetics(BC_serie_list_eqm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_hist"=gridArithmetics(raw_serie_list$hus_raw_train_esc,1000,operator = '*'),
             # "ETA_BC_hist"=gridArithmetics(BC_serie_list_eqm$hus_BC_train_esc,1000,operator = '*'),
             "ETA_BC_RCP8.5_qdm"=gridArithmetics(BC_serie_list_qdm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_BC_RCP8.5_dqm"=gridArithmetics(BC_serie_list_dqm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_BC_RCP8.5_mcbr"=gridArithmetics(BC_serie_list_mbcr$hus_BC_fut_esc,1000,operator = '*'),
             "ETA_BC_RCP8.5_mcbn"=gridArithmetics(BC_serie_list_mbcn$hus_BC_fut_esc,1000,operator = '*'),
             xyplot.custom = list(ylim=c(-2,5)),cols=c('black','red','blue','green','orange'))
             

dev.off()

REMO_RAW_RCP8.5=gridArithmetics(raw_serie_list$hus_raw_fut_esc,1000,operator = '*')
REMO_BC_RCP8.5_eqm=gridArithmetics(BC_serie_list_eqm$hus_BC_fut_esc,1000,operator = '*')
ETA_RAW_RCP8.5=gridArithmetics(raw_serie_list$hus_raw_fut_esc,1000,operator = '*')
ETA_BC_RCP8.5_eqm=gridArithmetics(BC_serie_list_eqm$hus_BC_fut_esc,1000,operator = '*')

temporalPlot("OBS"=gridArithmetics(hus_obs_esc,1000,operator = '*'),
             "REMO_RAW_RCP8.5"=REMO_RAW_RCP8.5, 
             "REMO_BC_RCP8.5_eqm"=REMO_BC_RCP8.5_eqm,
             "ETA_BC_RCP8.5_eqm"=ETA_BC_RCP8.5_eqm,
             "ETA_RAW_RCP8.5"=ETA_RAW_RCP8.5,
             xyplot.custom = list(ylim=c(-2,5)),cols=c('black','red','blue','green','orange'))

temporalPlot("OBS"=gridArithmetics(hus_obs_esc,1000,operator = '*'),
             "REMO_RAW_RCP8.5"=gridArithmetics(raw_serie_list$hus_raw_fut_esc,1000,operator = '*'), 
             "REMO_BC_RCP8.5_eqm"=gridArithmetics(BC_serie_list_eqm$hus_BC_fut_esc,1000,operator = '*'),
             "REMO_hist"=gridArithmetics(raw_serie_list$hus_raw_train_esc,1000,operator = '*'),
             "REMO_BC_hist"=gridArithmetics(BC_serie_list_eqm$hus_BC_train_esc,1000,operator = '*'),
             # "REMO_BC_RCP8.5_qdm"=gridArithmetics(BC_serie_list_qdm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_BC_RCP8.5_dqm"=gridArithmetics(BC_serie_list_dqm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_BC_RCP8.5_mcbr"=gridArithmetics(BC_serie_list_mbcr$hus_BC_fut_esc,1000,operator = '*'),
             # "REMO_BC_RCP8.5_mcbn"=gridArithmetics(BC_serie_list_mbcn$hus_BC_fut_esc,1000,operator = '*'),
             xyplot.custom = list(ylim=c(-2,5)),cols=c('black','red','blue','green','orange'))


temporalPlot("OBS"=gridArithmetics(hus_obs_esc,1000,operator = '*'),
             "RegCM_RAW_RCP8.5"=gridArithmetics(raw_serie_list$hus_raw_fut_esc,1000,operator = '*'), 
             "RegCM_BC_RCP8.5_eqm"=gridArithmetics(BC_serie_list_eqm$hus_BC_fut_esc,1000,operator = '*'),
             "RegCM_hist"=gridArithmetics(raw_serie_list$hus_raw_train_esc,1000,operator = '*'),
             "RegCM_BC_hist"=gridArithmetics(BC_serie_list_eqm$hus_BC_train_esc,1000,operator = '*'),
             # "RegCM_BC_RCP8.5_qdm"=gridArithmetics(BC_serie_list_qdm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_BC_RCP8.5_dqm"=gridArithmetics(BC_serie_list_dqm$hus_BC_fut_esc,1000,operator = '*'),
             # "ETA_BC_RCP8.5_mcbr"=gridArithmetics(BC_serie_list_mbcr$hus_BC_fut_esc,1000,operator = '*'),
             # "RegCM_BC_RCP8.5_mcbn"=gridArithmetics(BC_serie_list_mbcn$hus_BC_fut_esc,1000,operator = '*'),
             xyplot.custom = list(ylim=c(-2,5)),cols=c('black','red','blue','green','orange'))




# 
temporalPlot("OBS"=gridArithmetics(hus_obs_esc,1000,operator = '*'),xyplot.custom = list(ylim=c(-1,1)))
####################################### mapas ########################################
setwd(path_2)
list=ls(pattern = "^BC_serie_list_")
list_years=list("HadGEM2"=years_HadGEM2,"MIROC5"=years_MIROC5,"CanESM2"=years_CanESM2)



temporalPlot("OBS"=gridArithmetics(hus_obs_esc,1000,operator = '*'),
             "RegCM_RAW_RCP8.5"=gridArithmetics(raw_serie_list$hus_raw_fut_esc,1000,operator = '*'), 
             "RegCM_BC_RCP8.5_eqm"=gridArithmetics(BC_serie_list_eqm$hus_BC_fut_esc,1000,operator = '*'),
             # "RegCM_hist"=gridArithmetics(raw_serie_list$hus_raw_train_esc,1000,operator = '*'),
             # "RegCM_BC_hist"=gridArithmetics(BC_serie_list_eqm$hus_BC_train_esc,1000,operator = '*'),
             "RegCM_BC_RCP8.5_qdm"=gridArithmetics(BC_serie_list_qdm$hus_BC_fut_esc,1000,operator = '*'),
             # "RegCM_BC_RCP8.5_dqm"=gridArithmetics(BC_serie_list_dqm$hus_BC_fut_esc,1000,operator = '*'),
             # "RegCM_BC_RCP8.5_mcbr"=gridArithmetics(BC_serie_list_mbcr$hus_BC_fut_esc,1000,operator = '*'),
             "RegCM_BC_RCP8.5_mcbn"=gridArithmetics(BC_serie_list_mbcn$hus_BC_fut_esc,1000,operator = '*'),
             xyplot.custom = list(ylim=c(-2,5)),cols=c('black','red','blue','green','orange'))
