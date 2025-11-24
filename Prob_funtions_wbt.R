
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
require(HeatStress)

path_0="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/Pruebas_CV/Train/"
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/"
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"

var='wbt'
GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
rcms=c("REMO","RegCM",'ETA')
methods=c('eqm','dqm','qdm',"mbcr","mbcn")
name=c("ETA_HadGEM2","ETA_MIROC5", "ETA_CanESM2",           
       "RegCM_HadGEM2",  "RegCM_MPI-M-MPI", "RegCM_NorESM1",   
       "REMO_HadGEM2",   "REMO_MPI-M-MPI",  "REMO_NorESM1") 

####################### funciones #######################33
load_rda_list <- function(files, path) {
  full_paths <- file.path(path, files)
  map(full_paths, ~ get(load(.x)))
}

calc_wbt <- function(tas_grid, hus_grid) {
  HR <- huss2hurs(huss = as.vector(hus_grid$Data),
                  tas  = as.vector(tas_grid$Data), ps = 1013)
  sw <- wbt.Stull(as.vector(tas_grid$Data), hurs = HR)
  hus_grid$Data <- array(sw, dim = dim(hus_grid$Data))
  attributes(hus_grid$Data)=attributes(tas_grid$Data)
  hus_grid
}
# PAFeff=prctile(BC_TEST, pct)- prctile(RAW_TEST, pct);
# PAF=prctile(OBS, pct)- prctile(RAW_TRAIN, pct);
####observacion ########################################################
setwd(path_1)

obs.huss=get(load(paste0(path_1,"/Datos_CV_BC_huss_MSWX.rda")))%>% subsetGrid(years=1981:2004,season=c(12,1,2))

obs.Tx=get(load(paste0(path_1,"/Datos_CV_BC_Tmax_MSWX.rda")))%>% subsetGrid(years=1981:2004,season=c(12,1,2))#2166 6 bisciestos

HR_aprox_Tmax.obs=huss2hurs(huss=as.vector(obs.huss$Data),tas=as.vector(obs.Tx$Data),ps=1013)

wbt.obs=wbt.Stull(as.vector(obs.Tx$Data), hurs=HR_aprox_Tmax.obs)

wbt.obs_data=obs.huss
wbt.obs_data$Data=array(wbt.obs,dim=dim(obs.huss$Data))
attributes(wbt.obs_data$Data)=attributes(obs.huss$Data)
rm(wbt.obs,HR_aprox_Tmax.obs)
gc()

obs=wbt.obs_data

####### RAW ######################################################################
k=3#ETA
k=1#REMO
k=2#RegCM

#Ejemplo coordenadas al azar

files_Tx  <- list.files(path_1,pattern=paste0('Datos.*.','tasmax','.*.',rcms[k]))
files_Hus <- list.files(path_1,pattern=paste0('Datos.*.','hus','.*.',rcms[k]))

if(k==3){
  files_Tx=files_Tx[c(2,3,1)];files_Hus=files_Hus[c(1,3,2)]}

tasmax_hist <- load_rda_list(files_Tx, path_1)
hus_hist    <- load_rda_list(files_Hus, path_1)


coordenadas=c(50,21)
obs$xyCoords$y[coordenadas[1]]
obs$xyCoords$x[coordenadas[2]]

coordenadas=c(100,71)
obs$xyCoords$y[coordenadas[1]]
obs$xyCoords$x[coordenadas[2]]

hus_raw_train <- map2(tasmax_hist, hus_hist, function(tx, hu) {
  
  # Calcular WBT y subset DJF
  raw <- calc_wbt(tx, hu) %>% 
    subsetGrid(years = 1981:2004, season = c(12,1,2))
  
  # Alinear temporalmente con OBS
  raw_aligned <- getTemporalIntersection(obs = obs, prd = raw, which.return = 'prd')
  obs_aligned <- getTemporalIntersection(obs = obs, prd = raw, which.return = 'obs')
  
  # Extraer series temporales
  val_raw <- raw_aligned$Data[, coordenadas[1], coordenadas[2]]
  val_obs <- obs_aligned$Data[, coordenadas[1], coordenadas[2]]
  
  # Percentiles
  p <- seq(0, 100, 1)
  q_raw <- quantile(val_raw, probs = p/100, na.rm=TRUE)
  q_obs <- quantile(val_obs, probs = p/100, na.rm=TRUE)
  
  data.frame(
    percentil = p,
    dif = q_obs - q_raw,
    Lon = obs$xyCoords$x[coordenadas[2]],
    Lat = obs$xyCoords$y[coordenadas[1]],
    Period = "Train"
  )
})

hus_raw_df <-rbind(bind_rows(hus_raw_train[[1]]),bind_rows(hus_raw_train[[2]]),bind_rows(hus_raw_train[[3]]))
row.names(hus_raw_df)=NULL
hus_raw_df$modelo = rep(name[ list_rcms[[k]] ], each = 101)
hus_raw_df$method = "train"


plot(hus_raw_train[[1]]$percentil, hus_raw_train[[1]]$dif, type="l",ylim=c(-1,1.5),
     xlab="Percentil", ylab="Obs - RAW")
abline(h=0, lty=2)
lines(hus_raw_train[[3]]$dif,col='red')
lines(hus_raw_train[[2]]$dif,col='green')

################### futuro ####################

k <- 1  # 1=REMO, 2=RegCM, 3=ETA


# ---  RAW futuro ---
files_fut_raw_hus <- list.files(
  path = paste0(path_2, rcms[k]),
  pattern = paste0("Datos_hus.*RCP8.5_.*.rda"),
  full.names = TRUE
)

files_fut_raw_Tx <- list.files(
  path = paste0(path_2, rcms[k]),
  pattern = paste0("Datos_tasmax.*RCP8.5_.*.rda"),
  full.names = TRUE
)

if (k == 3) {
  files_fut_raw_Tx <- files_fut_raw_Tx[c(2,3,1)]
files_fut_raw_hus <- files_fut_raw_hus[c(2,3,1)]
                                       }

# --- Crear data.frame vacío donde guardar resultados ---
pp_all <- data.frame()

for (f in seq_along(files_fut_raw_Tx)) {
  
  message(paste0("Procesando RAW: ", basename(files_fut_raw_Tx[f])))
  
  tas_fut <- get(load(files_fut_raw_Tx[f]))
  hus_fut <- get(load(files_fut_raw_hus[f]))
  
  raw <- calc_wbt(tas_fut, hus_fut) %>%
    subsetGrid(years = 2007:2099, season = c(12,1,2))
  
  gc()
  # --- Diferentes  BC ---
  for (j in seq_along(methods)) {
    method <- methods[j]
    message(paste0("   method.. ", method))
    special_case <- method %in% c("mbcr", "mbcn")
    
    # Archivos BC futuros 
    if (!special_case) {
    hus_files_fut_method <- list.files(
      path = paste0(path_2, rcms[k], "/BC"),
      pattern = paste0("hus_.*_", method, "_.*_RCP8.5_SESA.rda"),
      full.names = TRUE
    )
    
    # Archivos BC futuros de tasmax
    Tx_files_fut_method <- list.files(
      path = paste0(path_2, rcms[k], "/BC"),
      pattern = paste0("tasmax_.*_", method, "_.*_RCP8.5_SESA.rda"),
      full.names = TRUE
    )
    # Si ETA: ordenar
    if (k == 3) {
      hus_files_fut_method <- hus_files_fut_method[c(2,3,1)]
      Tx_files_fut_method  <- Tx_files_fut_method[c(2,3,1)]
    }
    
    # Cargar BC hus y tasmax para el modelo f
    hus_BC <- get(load(hus_files_fut_method[f]))
    Tx_BC  <- get(load(Tx_files_fut_method[f]))
    }
    # Para mbcr/mbcn: seleccionar la variable correcta del objeto multi-variable
    if (special_case) {
      hus_files_fut_method <- list.files(
        path = paste0(path_2, rcms[k], "/BC"),
        pattern = paste0("hus.*_", method, "_RCP8.5_SESA.rda"),
        full.names = TRUE
      )
      
      
      if (k == 3) {
        hus_files_fut_method <- hus_files_fut_method[c(2,3,1)]}
        
        hus_BC <- get(load(hus_files_fut_method[f]))
        Tx_BC  <- subsetGrid(hus_BC,  var = hus_BC$Variable$varName[1])
        hus_BC <- subsetGrid(hus_BC, var = hus_BC$Variable$varName[2])
    }
    
    # Calcular WBT_BC
    wbt_BC <- calc_wbt(Tx_BC, hus_BC) %>%
      subsetGrid(years = 2007:2099, season = c(12,1,2))
    
    # Alinear temporalmente WBT_BC con WBT_RAW pre-computado
    raw_aligned <- getTemporalIntersection(obs = wbt_BC, prd = raw, which.return = 'prd')
    BC_aligned  <- getTemporalIntersection(obs = wbt_BC, prd = raw, which.return = 'obs')
    
    # Extraer punto
    val_raw <- raw_aligned$Data[, coordenadas[1], coordenadas[2]]
    val_BC  <- BC_aligned$Data[, coordenadas[1], coordenadas[2]]
    
    # Percentiles y diferencia
    p <- seq(0, 100, 1)
    q_raw <- quantile(val_raw, probs = p/100, na.rm = TRUE)
    q_BC  <- quantile(val_BC,  probs = p/100, na.rm = TRUE)
    q_diff <- q_BC - q_raw
    
    # Guardar info
    pp <- data.frame(
      percentil = p,
      dif = q_diff,
      Lon = raw$xyCoords$x[coordenadas[2]],
      Lat = raw$xyCoords$y[coordenadas[1]],
      Period = "Fut",
      method = method,
      modelo = name[list_rcms[[k]][f]]
    )
    
    pp_all <- rbind(pp_all, pp)
    gc()
  }
  
}

# head(pp_all)

pp_all_all=rbind(pp_all,hus_raw_df)
# pp_all_all_p1=pp_all_all
library(ggplot2)
my_theme <- theme_bw() + theme(legend.position = "bottom",legend.title = element_text(size=8,colour = "black",face = "bold"),
                               plot.title = element_text(hjust = 0.5,size=7,colour = "black",face = "bold"),axis.text = element_text(size=7,colour = "black"),
                               axis.title = element_text(size=6,face="bold"),legend.key.size = unit(0.1,"line"),legend.text = element_text(size = 8),legend.key.height = unit(0.1,"inch"),legend.key.width = unit(0.1,"inch"),plot.margin = margin(0.15,0.15,0.15,0,"cm")
                               ,axis.text.x = element_text(angle = 360,hjust = 1,size = 6,color = "black",face="bold"),
                               axis.text.y = element_text(hjust = 1,size = 7,color = "black",face="bold"),strip.text.x = element_text(size = 8, colour = "black",face = "bold"),
                               strip.text.y = element_text(size = 6, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank(),panel.grid.major = element_blank(),
                               )

p=ggplot(pp_all_all, aes(x = percentil, y = dif, color = method)) + scale_color_manual(values=colores_metodos)+
  geom_line() +my_theme+
  facet_wrap(.~ modelo) +scale_y_continuous(breaks = seq(-1.2,2.5,by=0.5),limits=c(-1.2,2.5))+
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Percentil", y="Diferencia (BC - RAW)") 
p
ggsave(p,filename=paste0("P_fut_punto_wbt_",rcms[k],"_p1.png"),device="png",dpi=300,width=6,height = 4,scale=0.7)

  colores_metodos <- c(
    eqm  = "#1b9e77",  # verde
    dqm  = "#d95f02",  # naranja
    qdm  = "#7570b3",  # violeta
    mbcr = "#e7298a",  # magenta
    mbcn = "#66a61e" ,  # verde claro
    train = "black"   # verde claro
    
  )
  