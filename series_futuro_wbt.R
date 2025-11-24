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

## Promediar en DJF
aggregate_djf <- function(grid) {
  aggregateGrid(grid, aggr.s = list(FUN = list("mean", na.rm = TRUE),
                                    season = c(12, 1, 2)))
}

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
      aggregate_djf(calc_wbt(Tx, Hus))
    }) %>% bindGrid(dimension = "member")
    
  } else {
    files_Tx  <- list.files(path_train, pattern = paste0("tasmax_.*_", method, "_.*Train_SESA.rda"))[list_rcms[[k]]]
    files_Hus <- list.files(path_train, pattern = paste0("hus_.*_", method, "_.*Train_SESA.rda"))[list_rcms[[k]]]
    
    if(k==3){files_Tx=files_Tx[c(2,3,1)]; files_Hus=files_Hus[c(2,3,1)]}
    
    tasmax_train <- load_rda_list(files_Tx, path_train)
    hus_train    <- load_rda_list(files_Hus, path_train)
    
    wbt_train <- map2(tasmax_train, hus_train, ~ {
      aggregate_djf(calc_wbt(.x, .y))
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
      aggregate_djf(calc_wbt(Tx, Hus))%>% 
        subsetGrid(years = 2006:2098)
    }) %>% bindGrid(dimension = "member")
    
  } else {
    files_Tx <- list.files(file.path(path_fut, RCMs[k], "BC"),
                           pattern = paste0("tasmax_.*_", method, "_.*_RCP8.5_SESA.rda"))
    files_Hus <- list.files(file.path(path_fut, RCMs[k], "BC"),
                            pattern = paste0("hus_.*_", method, "_.*_RCP8.5_SESA.rda"))
    if(k==3){files_Tx=files_Tx[c(2,3,1)]; files_Hus=files_Hus[c(2,3,1)]}
    
    tasmax_fut <- load_rda_list(files_Tx, file.path(path_fut, RCMs[k], "BC"))
    hus_fut    <- load_rda_list(files_Hus, file.path(path_fut, RCMs[k], "BC"))
    
    wbt_fut <- map2(tasmax_fut, hus_fut, ~ {
      aggregate_djf(calc_wbt(.x, .y))%>% 
        subsetGrid(years = 2006:2098)
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

save(list = ls(pattern = "^wbt_BC_serie_list_"), file = paste0("wbt_BC_serie_list_",RCMs[k],".rda"))

###### RAW #####

files_Tx  <- list.files(path_hist,pattern=paste0('Datos.*.','tasmax','.*.',RCMs[k]))
files_Hus <- list.files(path_hist,pattern=paste0('Datos.*.','hus','.*.',RCMs[k]))

if(k==3){
  files_Tx=files_Tx[c(2,3,1)];files_Hus=files_Hus[c(1,3,2)]}
    
tasmax_hist <- load_rda_list(files_Tx, path_hist)
hus_hist    <- load_rda_list(files_Hus, path_hist)

wbt_raw_hist <- map2(tasmax_hist, hus_hist, ~ {
  aggregate_djf(calc_wbt(.x, .y))
}) %>% bindGrid(dimension = "member")

wbt_raw_hist$Members <- members_name[list_rcms[[k]]]
wbt_raw_hist[["Variable"]][["varName"]] <- "wbt"

wbt_hist_esc <- scaleGrid(
  grid = wbt_raw_hist,
  base = subsetGrid(wbt_raw_hist, years = 1984:2003),
  by.member = TRUE
)

############ futuro ###############
files_Tx <- list.files(file.path(path_fut, RCMs[k]),
                       pattern = paste0("Datos_tasmax_.*.RCP8.5_.*.rda"))
files_Hus <- list.files(file.path(path_fut, RCMs[k]),
                        pattern  = paste0("Datos_hus.*.RCP8.5_.*.rda"))

if(k==3){files_Tx=files_Tx[c(2,3,1)]; files_Hus=files_Hus[c(2,3,1)]}

tasmax_fut <- load_rda_list(files_Tx, file.path(path_fut, RCMs[k]))
hus_fut    <- load_rda_list(files_Hus, file.path(path_fut, RCMs[k]))

wbt_raw_fut <- map2(tasmax_fut, hus_fut, ~ {
  aggregate_djf(calc_wbt(.x, .y))%>% 
    subsetGrid(years = 2006:2098)
}) %>% bindGrid(dimension = "member")


wbt_raw_fut$Members <-members_name[list_rcms[[k]]]

wbt_raw_fut_esc <- scaleGrid(
  grid = wbt_raw_fut,
  base = subsetGrid(wbt_raw_hist, years = 1984:2003),
  by.member = TRUE
)

raw_serie_list <- list(
  wbt_raw_train     = wbt_raw_hist,
  wbt_raw_train_esc = wbt_hist_esc,
  wbt_raw_fut       = wbt_raw_fut,
  wbt_raw_fut_esc   = wbt_raw_fut_esc
)

assign("wbt_raw_serie_list", raw_serie_list)
save(list = ls(pattern = "^wbt_raw_serie_list"), file = paste0("wbt_raw_serie_list_",RCMs[k],".rda"))
