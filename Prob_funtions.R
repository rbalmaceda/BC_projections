
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

var='hus'
GCM=c("HadGEM2","MPI-M-MPI", "NorESM1")
rcms=c("REMO","RegCM",'ETA')
methods=c('eqm','dqm','qdm',"mbcr","mbcn")
name=c("ETA_HadGEM2","ETA_MIROC5", "ETA_CanESM2",           
       "RegCM_HadGEM2",  "RegCM_MPI-M-MPI", "RegCM_NorESM1",   
       "REMO_HadGEM2",   "REMO_MPI-M-MPI",  "REMO_NorESM1") 


# PAFeff=prctile(BC_TEST, pct)- prctile(RAW_TEST, pct);
# PAF=prctile(OBS, pct)- prctile(RAW_TRAIN, pct);
####observacion ########################################################33
setwd(path_1)

obs=get(load(paste0(path_1,"/Datos_CV_BC_huss_MSWX.rda")))%>% subsetGrid(years=1981:2004,season=c(12,1,2))

####### RAW ######################################################################
k=3#ETA
k=1#REMO
k=2#RegCM
#Ejemplo coordenadas al azar

setwd(path_1)

files_train_raw=list.files(path=path_1,pattern=paste0('Datos.*.','hus'))
files_train_raw=files_train_raw[-1]
# files_train_raw=files_train_raw[c(1,3,2)]
files_train_raw=files_train_raw[7:9]
# files_train_raw=files_train_raw[4:6]

coordenadas=c(50,21)
obs$xyCoords$y[coordenadas[1]]
obs$xyCoords$x[coordenadas[2]]

coordenadas=c(100,71)
obs$xyCoords$y[coordenadas[1]]
obs$xyCoords$x[coordenadas[2]]


hus_raw_train <- map(files_train_raw, function(f) {
  raw <- get(load(f)) %>% subsetGrid(years=1981:2004, season=c(12,1,2))
  
  # Alinear temporalmente
  raw_aligned <- getTemporalIntersection(obs=obs, prd=raw, which.return='prd')
  obs_aligned <- getTemporalIntersection(obs=obs, prd=raw, which.return='obs')
  
  # Extraer serie temporal del punto
  val_raw <- raw_aligned$Data[, coordenadas[1], coordenadas[2]]
  val_obs <- obs_aligned$Data[, coordenadas[1], coordenadas[2]]
  
  # Calcular percentiles 0-100
  p <- seq(0, 100, by = 1)[-c(1,101)]
  q_raw <- quantile(val_raw, probs = p/100, na.rm=TRUE)
  q_obs <- quantile(val_obs, probs = p/100, na.rm=TRUE)
  
  # Diferencia por percentil
  q_diff <- q_obs - q_raw
  
  data.frame(percentil=p, dif=q_diff, Lon  = rep(obs$xyCoords$x[coordenadas[2]],length(p)),
             Lat  = rep(obs$xyCoords$y[coordenadas[1]], length(p)),
             Period = 'Train')
})

# names(hus_raw_train) <- name[1:3]
# names(hus_raw_train) <- name[7:9]
hus_raw_df <-rbind(bind_rows(hus_raw_train[[1]]),bind_rows(hus_raw_train[[2]]),bind_rows(hus_raw_train[[3]]))
row.names(hus_raw_df)=NULL
hus_raw_df$modelo=rep(name[7:9],each=99)
hus_raw_df$method=rep("train",297)


plot(hus_raw_train[[2]]$percentil, hus_raw_train[[2]]$dif*1000, type="l",ylim=c(-5,5),
     xlab="Percentil", ylab="Obs - RAW")
abline(h=0, lty=2)
lines(hus_raw_train[[1]]$dif*1000,col='red')
lines(hus_raw_train[[3]]$dif*1000,col='green')

################### futuro ####################
list_rcms=list(7:9,4:6,1:3)

k <- 1  # 1=REMO, 2=RegCM, 3=ETA

# ---  RAW futuro ---
files_fut_raw <- list.files(
  path = paste0(path_2, rcms[k]),
  pattern = paste0("Datos_hus.*RCP8.5_.*.rda"),
  full.names = TRUE
)

if (k == 3) files_fut_raw <- files_fut_raw[c(2,3,1)]

# --- Crear data.frame vacío donde guardar resultados ---
pp_all <- data.frame()

for (f in seq_along(files_fut_raw)) {
  
  message(paste0("Procesando RAW: ", basename(files_fut_raw[f])))
  raw <- get(load(files_fut_raw[f])) %>% subsetGrid(years = 2007:2099, season = c(12,1,2))
  # raw <- get(load(files_fut_raw[f])) %>% subsetGrid(years = 2040:2059, season = c(12,1,2))
  
    gc()
  # --- Diferentes  BC ---
  for (j in seq_along(methods)) {
    method <- methods[j]
    message(paste0("   method.. ", method))
    special_case <- method %in% c("mbcr", "mbcn")
   
     files_fut_method <- list.files(
      path = paste0(path_2, rcms[k], "/BC"),
      pattern = paste0("hus.*.", method, ".*._RCP8.5_SESA.rda"), 
      full.names = TRUE
    )
    if(length(files_fut_method)==0) {
      files_fut_method <- list.files(
        path = paste0(path_2, rcms[k], "/BC"),
        pattern = paste0("hus.*.", method, "_RCP8.5_SESA.rda"), 
        full.names = TRUE
      )
      
    }
    if (k == 3) files_fut_method <- files_fut_method[c(2,3,1)]
    
    BC <- get(load(files_fut_method[f])) 
    if (special_case) {
      BC <- subsetGrid(BC, var = BC[["Variable"]][["varName"]][2])
    }
    # --- Alinear temporalmente ---
    raw_aligned <- getTemporalIntersection(obs = BC, prd = raw, which.return = 'prd')
    BC_aligned  <- getTemporalIntersection(obs = BC, prd = raw, which.return = 'obs')
    
    # --- Extraer serie del punto ---
    val_raw <- raw_aligned$Data[, coordenadas[1], coordenadas[2]]
    val_BC  <- BC_aligned$Data[, coordenadas[1], coordenadas[2]]
    
    # --- Calcular percentiles y diferencias ---
    p <- seq(0, 100, by = 1)[-c(1,101)]
    q_raw <- quantile(val_raw, probs = p/100, na.rm = TRUE)
    q_BC  <- quantile(val_BC,  probs = p/100, na.rm = TRUE)
    q_diff <- q_BC - q_raw
    
    # --- Crear dataframe con etiquetas ---
    pp <- data.frame(
      percentil = p,
      dif = q_diff,
      Lon = rep(raw$xyCoords$x[coordenadas[2]], length(p)),
      Lat = rep(raw$xyCoords$y[coordenadas[1]], length(p)),
      Period = "Fut",
      method = rep(method, length(p)),
      modelo = rep(name[list_rcms[[k]][f]], length(p))
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

# pp=subset(pp_all_all,method %in% methods[c(4)])
p=ggplot(pp_all_all, aes(x = percentil, y = dif*1000, color = method)) + scale_color_manual(values=colores_metodos)+
  geom_line() +my_theme+
  facet_wrap(.~ modelo) +scale_y_continuous(breaks = seq(-1.2,2.5,by=0.5))+
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Percentil", y="Diferencia (BC - RAW)") 
ggsave(p,filename=paste0("P_fut_punto_",rcms[k],"_p1.png"),device="png",dpi=300,width=6,height = 4,scale=0.7)

  colores_metodos <- c(
    eqm  = "#1b9e77",  # verde
    dqm  = "#d95f02",  # naranja
    qdm  = "#7570b3",  # violeta
    mbcr = "#e7298a",  # magenta
    mbcn = "#66a61e" ,  # verde claro
    train = "black"   # verde claro
    
  )
  