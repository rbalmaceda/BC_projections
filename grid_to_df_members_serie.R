require(lubridate)
grid_to_df_members_serie <- function(grid) {
  
  
  # Extraer componentes
  data_array <- grid$Data
  lon <- grid$xyCoords$x
  lat <- grid$xyCoords$y
  time <- grid$Dates$start
  members <- grid$Members
  
  # Dimensiones esperadas: [member, lat, lon, time]
  stopifnot(length(dim(data_array)) == 4)
  
  n_mem <- dim(data_array)[1]
  
  df_list <- vector("list", n_mem)
  
  for (m in 1:n_mem) {
    df_list[[m]] <- data.frame(
      Lon   = rep(lon, each = length(lat)),
      Lat   = rep(lat, length(lon)),
      Data  = as.vector(data_array[m,,,]),
      Member = rep(members[m],length(lon)),
      Time=year(time)
    )
  }
  
  df_final <- do.call(rbind, df_list)
  return(df_final)
}



grid_to_df_members_serie_means <- function(grid) {
  
  grid=aggregateGrid(grid,aggr.spatial = list(FUN='mean',na.rm=T))
  # Extraer componentes
  data_array <- grid$Data
  time <- grid$Dates$start
  members <- grid$Members
  
  # Dimensiones esperadas: [member, time]
  stopifnot(length(dim(data_array)) == 2)
  
  n_mem <- dim(data_array)[1]
  
  df_list <- vector("list", n_mem)
  
  for (m in 1:n_mem) {
    df_list[[m]] <- data.frame(
      Data  = as.vector(data_array[m,]),
      Member = rep(members[m],length(time)),
      Time=year(time)
    )
  }
  
  df_final <- do.call(rbind, df_list)
  return(df_final)
}
