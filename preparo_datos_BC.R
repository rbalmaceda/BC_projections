
####Data cordex 
df <- read.csv("https://data.meteo.unican.es/inventory.csv")
subset <- subset(df, activity == 'CORDEX' & domain=='SAM-22')
unique(subset$model)

subset=subset(subset, model %in% unique(subset$model)[2:5])
subset=subset(subset, variable %in% unique(subset$variable)[c(1:3,7)])

###### datos #######

#####Voy a comenzar probando con la variable Tmax ########
unique(subset$variable)
GCM=as.character(unique(subset$model))
var=as.character(unique(subset$variable))

####REMO2015######

Tmax_REMO=loadGridData("//media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/Proyecciones_Regionales/Tmax/crudos/BiasCorr/tasmax_SAM-22_NCC-NorESM1-M_historical_r1i1p1_GERICS-REMO2015_v1_day_1970_2005ll.nc",
                       year=1981:2005,season=c(12,1,2),var="tasmax")
gc()
grid_REMO=getGrid(Tmax_REMO)

Tmax_REMO_SESA<- subsetGrid(data.masked,lonLim = c(-64,-45),latLim = c(-40,-18))
#enmascarar y luego recortar
rm(Tmax_REMO)
Tmax_REMO_SESA <- gridArithmetics(Tmax_REMO_SESA, 273.15, operator="-")
gc()

spatialPlot(climatology(Tmax_REMO_SESA) , backdrop.theme = "countries", 
            rev.colors = TRUE, main= "tasmax-REMO", as.table=TRUE, scales = list(draw = TRUE),ylim=c(-45,-16),
            xlim=c(-70,-40),at= seq(15,36)) # pinta ambos en las latitudes que tocan
save(Tmax_REMO_SESA,file="Datos_CV_BC_Tmax_REMO_SESA_NorESM1.rda")

######################## RegCM######################################################

Tmax_RegCM=loadGridData("/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/Proyecciones_Regionales/Tmax/crudos/BiasCorr/tasmax_SAM-22_NCC-NorESM1-M_historical_r1i1p1_ICTP-RegCM4-7_v0_day_1970_2005_pp.nc",
                        year=1981:2005,season=c(12,1,2),var="tasmax")

#Enmascarar el oceano levantando el archivo de rcm.mask
Tmax_RegCM_SESA=subsetGrid(data.masked,lonLim = c(-64,-45),latLim = c(-40,-18))
Tmax_RegCM_SESA <- gridArithmetics(Tmax_RegCM_SESA, 273.15, operator="-")

spatialPlot(climatology(Tmax_RegCM_SESA) , backdrop.theme = "countries", 
            rev.colors = TRUE, main= "tasmax-RegCM", as.table=TRUE, scales = list(draw = TRUE),ylim=c(-45,-16),
            xlim=c(-70,-40),at= seq(15,36)) # pinta ambos en las latitudes que tocan
save(Tmax_RegCM_SESA,file="Datos_CV_BC_Tmax_RegCM_SESA_NorESM1.rda")
######################## Observacion ######################################################

# huss_MSWX=loadGridData("/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/specific_humidity/specific_humidity_1979_2022_up.nc",
#                        year=1981:2005,season=c(12,1,2),var="specific_humidity",lonLim = c(-64,-45),latLim = c(-40,-18))

tasmax_MSWX=loadGridData("/media/usuario/Seagate Expansion Drive/MSWEX/Tmax/Tmax_1980_2022_DJF.nc",
                       year=1981,season=c(12,1,2),var="air_temperature")

grid_MSWX=getGrid(tasmax_MSWX)
save(grid_MSWX,file="grid_MSWX.rda")
#Enmascaro

######Veo que pasa con la laguna de mar chiquita #######

spatialPlot(climatology(tasmax_MSWX.masked),backdrop.theme = 'countries',scales = list(draw = TRUE),
            main= "tasmax-MSWX-up",rev.colors = T,at= seq(15,36),xlim=c(-70,-40),ylim=c(-45,-16))

save(tasmax_MSWX.masked,file="Datos_CV_BC_Tmax_MSWX.rda")
gc()


######### UNIDADES ####
#kg/kg y Eta es g/kg
huss_MSWX=loadGridData("/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/specific_humidity/specific_humidity_1979_2022_up.nc",
                         year=1981:2005,season=c(12,1,2),var="specific_humidity",lonLim = c(-64,-45),latLim = c(-40,-18))

#Enmascaro
spatialPlot(climatology(huss_MSWX.masked),backdrop.theme = 'countries',scales = list(draw = TRUE),
            main= "hus-MSWX-up",rev.colors = T,xlim=c(-70,-40),ylim=c(-45,-16))

gc()
