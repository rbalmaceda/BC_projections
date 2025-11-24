################## serie wbt ##########################
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
library(transformeR)  
library(HeatStress)
library(climate4R.indices)

# --- Paths y configuraciones generales ---
base_path <- "/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC"
path_train <- file.path(base_path, "Datos_historico/Pruebas_CV/Train")
path_hist  <- file.path(base_path, "Datos_historico")
path_fut   <- file.path(base_path, "Datos_rcp8.5")
list_rcms=list(7:9,4:6,1:3)

vars <- c("tasmax", "hus")
GCMs <- c("HadGEM2", "MPI-M-MPI", "NorESM1")
RCMs <- c("REMO", "RegCM", "ETA")
methods <- c("eqm", "dqm", "qdm", "mbcr", "mbcn")
members_name <- c("ETA_HadGEM2", "ETA_MIROC5", "ETA_CanESM2",
                  "RegCM_HadGEM2", "RegCM_MPI-M-MPI", "RegCM_NorESM1",
                  "REMO_HadGEM2", "REMO_MPI-M-MPI", "REMO_NorESM1")
var='wbt'
# index='Ind27'
index='max'

# --- Funciones  ---

## Cargar una lista de archivos .rda
load_rda_list <- function(files, path) {
  full_paths <- file.path(path, files)
  map(full_paths, ~ get(load(.x)))
}

## Calcular wbt a partir de tasmax y hus/huss
calc_wbt <- function(tas_grid, hus_grid) {
  HR <- huss2hurs(huss = as.vector(hus_grid$Data),
                  tas  = as.vector(tas_grid$Data), ps = 1013)
  sw <- wbt.Stull(as.vector(tas_grid$Data), hurs = HR)
  hus_grid$Data <- array(sw, dim = dim(hus_grid$Data))
  attributes(hus_grid$Data)=attributes(tas_grid$Data)
  hus_grid
}

# ## Promediar en DJF
aggregate_djf <- function(grid) {
  aggregateGrid(grid, aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                    season = c(12, 1, 2)))
}


my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

aggregate_djf_max <- function(grid) {
  aggregateGrid(grid, aggr.s = list(FUN = list('my.max'),
                                    season = c(12, 1, 2)))}

# fecha_Train=aggregate_djf(tasmax_train[[2]])%>% subsetGrid(years = 1980:2003)


## Indice dias por encima de 30 
# aggregate_djf_Ind27 <- function(grid) {
#   indexGrid(tx = grid, time.resolution = "year", index.code = "TXth", th=30)%>% 
#     gridArithmetics(operator = '/',90) %>% gridArithmetics(operator = '*',100) 
# }

aggregate_djf_Ind27 <- function(grid) {
  ind <- indexGrid(tx = grid, time.resolution = "year", index.code = "TXth", th = 27)
  
  ind <- gridArithmetics(ind, operator = '/', e2 = 90) #Dividir solo donde el valor es > 0
  ind <- gridArithmetics(ind, operator = '*', e2 = 100) #Porcentual
  
  # Reemplazar infinitos o NaN por 0
  ind$Data[!is.finite(ind$Data)] <- 0
  
  return(ind)
}

#for (k in 1:3) {
  k=2
  message("Processing RCM: ", RCMs[k])
  
# --- BC ---################
for (method in methods) {
  message("Processing method: ", method)
  
  special_case <- method %in% c("mbcr", "mbcn")
  
  # --- TRAIN/HISTORICO ---
  if (special_case) {
    files_train <- list.files(path_train, pattern = paste0("huss_.*_", method, "_Train_SESA.rda"))[list_rcms[[k]]]
    
    if(k==3){files_train=files_train[c(2,3,1)]}
    grids_train <- load_rda_list(files_train, path_train)
    
    wbt_train <- map(1:length(grids_train), function(g) {
      Tx  <- subsetGrid(grids_train[[g]], var = grids_train[[g]][["Variable"]][["varName"]][1])
      Hus <- subsetGrid(grids_train[[g]], var = grids_train[[g]][["Variable"]][["varName"]][2])
      if(index=='Ind27'){
       out <- calc_wbt(Tx, Hus) %>% aggregate_djf_Ind27()
       out$Dates <- aggregate_djf(Tx)$Dates  }
      if(index=='max'){
      out <- calc_wbt(Tx, Hus) %>% aggregate_djf_max()}
        return(out)
      }) %>%  bindGrid(dimension = "member")
      
  
  } else {
    files_Tx  <- list.files(path_train, pattern = paste0("tasmax_.*_", method, "_.*Train_SESA.rda"))[list_rcms[[k]]]
    files_Hus <- list.files(path_train, pattern = paste0("hus_.*_", method, "_.*Train_SESA.rda"))[list_rcms[[k]]]
    
    if(k==3){files_Tx=files_Tx[c(2,3,1)]; files_Hus=files_Hus[c(2,3,1)]}
    
    tasmax_train <- load_rda_list(files_Tx, path_train)
    hus_train    <- load_rda_list(files_Hus, path_train)
    
    wbt_train <- map2(tasmax_train, hus_train, function(tx, hu) {
      if(index=='Ind27'){
      out <- calc_wbt(tx, hu) %>% aggregate_djf_Ind27()
      out$Dates <- aggregate_djf(tx)$Dates
      }
      if(index=='max'){
      out <- calc_wbt(tx, hu) %>% aggregate_djf_max()
      }
      return(out)
    }) %>% bindGrid(dimension = "member")
    
  } 
  wbt_train$Members <- members_name[list_rcms[[k]]]
  wbt_train[["Variable"]][["varName"]] <- "wbt"
  
  wbt_train_esc <- scaleGrid(
    grid = wbt_train,
    base = subsetGrid(wbt_train, years = 1984:2003),
    by.member = TRUE
  )
  
  # --- FUTURE ---

  if (special_case) {
    files_fut <- list.files(file.path(path_fut, RCMs[k], "BC"),
                            pattern = paste0("hus.*_", method, "_RCP8.5_SESA.rda"))
    if(k==3){files_fut=files_fut[c(2,3,1)]}
    
    grids_fut <- load_rda_list(files_fut, file.path(path_fut, RCMs[k], "BC"))
    
    wbt_fut <- map(1:length(grids_fut), function(g) {
      Tx  <- subsetGrid(grids_fut[[g]], var = grids_fut[[g]][["Variable"]][["varName"]][1])
      Hus <- subsetGrid(grids_fut[[g]], var = grids_fut[[g]][["Variable"]][["varName"]][2])
      if(index=='Ind27'){
      out <-aggregate_djf_Ind27(calc_wbt(Tx, Hus))
      out$Dates <- aggregate_djf(Tx)$Dates }
      if(index=='max'){
        out <-aggregate_djf_max(calc_wbt(Tx, Hus))}
      out=subsetGrid(out,years = 2006:2098)
      return(out)
    }) %>% bindGrid(dimension = "member")
    
  } else {
    files_Tx <- list.files(file.path(path_fut, RCMs[k], "BC"),
                           pattern = paste0("tasmax_.*_", method, "_.*_RCP8.5_SESA.rda"))
    files_Hus <- list.files(file.path(path_fut, RCMs[k], "BC"),
                            pattern = paste0("hus_.*_", method, "_.*_RCP8.5_SESA.rda"))
    if(k==3){files_Tx=files_Tx[c(2,3,1)]; files_Hus=files_Hus[c(2,3,1)]}
    
    tasmax_fut <- load_rda_list(files_Tx, file.path(path_fut, RCMs[k], "BC"))
    hus_fut    <- load_rda_list(files_Hus, file.path(path_fut, RCMs[k], "BC"))
    

  wbt_fut <- map2(tasmax_fut, hus_fut, function(tx, hu) {
    if(index=='Ind27'){
    out <- calc_wbt(tx, hu) %>% 
      aggregate_djf_Ind27()
    out$Dates <- aggregate_djf(tx)$Dates  }
    if(index=='max'){
      out <- calc_wbt(tx, hu) %>% 
        aggregate_djf_max()
      
    }
    out=subsetGrid(out,years = 2006:2098)
    return(out)
  }) %>% bindGrid(dimension = "member")
  
  }
  wbt_fut$Members <- members_name[list_rcms[[k]]]
  
  wbt_fut_esc <- scaleGrid(
    grid = wbt_fut,
    base = subsetGrid(wbt_train, years = 1984:2003),
    by.member = TRUE
  )
  
  # --- Guardo la lista ---
  BC_serie_list <- list(
    wbt_BC_train     = wbt_train,
    wbt_BC_train_esc = wbt_train_esc,
    wbt_BC_fut       = wbt_fut,
    wbt_BC_fut_esc   = wbt_fut_esc
  )
  
  assign(paste0("wbt_BC_serie_list_", method), BC_serie_list)
}

save(list = ls(pattern = "^wbt_BC_serie_list_"), file = paste0("wbt_",index,"_BC_serie_list_",RCMs[k],".rda"))

###### RAW #####
##################### historico ####################################
files_Tx  <- list.files(path_hist,pattern=paste0('Datos.*.','tasmax','.*.',RCMs[k]))
files_Hus <- list.files(path_hist,pattern=paste0('Datos.*.','hus','.*.',RCMs[k]))

if(k==3){
  files_Tx=files_Tx[c(2,3,1)];files_Hus=files_Hus[c(1,3,2)]}
    
tasmax_hist <- load_rda_list(files_Tx, path_hist)
hus_hist    <- load_rda_list(files_Hus, path_hist)

wbt_raw_hist <- map2(tasmax_hist, hus_hist, function(tx, hu) {

  if(index=='max'){
    out=aggregate_djf_max(calc_wbt(tx, hu)) %>% 
      subsetGrid(years = 1981:2004)
  }
  if (index=='Ind27'){
  out=aggregate_djf_Ind27(calc_wbt(tx,hu)) %>%
    subsetGrid(years = 1981:2005)}
  
  # if (index=='Ind27'){
  #   out=calc_wbt(tx,hu) 
  #   out=indexGrid(tx=out,time.resolution = "year", index.code = "TXth", th=30)%>%
  #    climatology()}

    return(out)
}) %>% bindGrid(dimension = "member")

# pp$Dates <- grid$Dates

wbt_raw_hist$Members <- members_name[list_rcms[[k]]]
wbt_raw_hist[["Variable"]][["varName"]] <- "wbt"

if(index=='max'){
  
  wbt_hist_esc <- scaleGrid(
    grid = wbt_raw_hist,
    base = subsetGrid(wbt_raw_hist, years = 1984:2003),
    # base = subsetGrid(wbt_raw_hist, years = 1986:2005),
    by.member = TRUE )
}

if(index=='Ind27'){
wbt_hist_esc <- scaleGrid(
  grid = wbt_raw_hist,
  base = subsetGrid(wbt_raw_hist, years = 1985:2004),
  # base = subsetGrid(wbt_raw_hist, years = 1986:2005),
  
  by.member = TRUE
)
}
############ futuro ###############
files_Tx <- list.files(file.path(path_fut, RCMs[k]),
                       pattern = paste0("Datos_tasmax_.*.RCP8.5_.*.rda"))
files_Hus <- list.files(file.path(path_fut, RCMs[k]),
                        pattern  = paste0("Datos_hus.*.RCP8.5_.*.rda"))

if(k==3){files_Tx=files_Tx[c(2,3,1)]; files_Hus=files_Hus[c(2,3,1)]}

tasmax_fut <- load_rda_list(files_Tx, file.path(path_fut, RCMs[k]))
hus_fut    <- load_rda_list(files_Hus, file.path(path_fut, RCMs[k]))

wbt_raw_fut <- map2(tasmax_fut, hus_fut, function(tx, hu) {
  if(index=='max'){
   out= aggregate_djf_max(calc_wbt(tx, hu))%>% 
      subsetGrid(years = 2006:2098)
  }
  if(index=='Ind27'){
 out= aggregate_djf_Ind27(calc_wbt(tx, hu))%>% 
    subsetGrid(years = 2007:2099)}
  return(out)
  

}) %>% bindGrid(dimension = "member")

# pp=indexGrid(tx = hus_fut[[2]], time.resolution = "year", index.code = "TXth", th=35)
# pp=subsetGrid(pp,years = 2007:2099)

wbt_raw_fut$Members <-members_name[list_rcms[[k]]]

if(index=='Ind27'){
wbt_raw_fut_esc <- scaleGrid(
  grid = wbt_raw_fut,
  base = subsetGrid(wbt_raw_hist, years = 1985:2004),
  # base = subsetGrid(wbt_raw_hist, years = 1986:2005),
  by.member = TRUE
)
}

if(index=='max'){
  wbt_raw_fut_esc <- scaleGrid(
    grid = wbt_raw_fut,
    base = subsetGrid(wbt_raw_hist, years = 1984:2003),
    # base = subsetGrid(wbt_raw_hist, years = 1986:2005),
    by.member = TRUE
  )
}

raw_serie_list <- list(
  wbt_raw_train     = wbt_raw_hist,
  wbt_raw_train_esc = wbt_hist_esc,
  wbt_raw_fut       = wbt_raw_fut,
  wbt_raw_fut_esc   = wbt_raw_fut_esc
)

assign("wbt_indice_raw_serie_list", raw_serie_list)
save(list = ls(pattern = "^wbt_indice_raw_serie_list"), file = paste0("wbt_",index,"_raw_serie_list_",RCMs[k],".rda"))
