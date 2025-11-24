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
library(RColorBrewer)

####################################### mapas ########################################
# list_years=list("HadGEM2"=years_HadGEM2,"MIROC5"=years_MIROC5,"CanESM2"=years_CanESM2)
list_years=list("cercano"=2040:2059,"lejano"=2080:2098)
# methods=c('dqm','eqm','mbcn',"mbcr","qdm")
rcms=c("REMO","RegCM",'ETA')
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"
rcms=c("REMO","RegCM",'ETA')
var='hus'
index='mean'
porcentual=TRUE #revisar con el resto de los scrips 

for (k in 1:3) {
  message("Processing RCM: ", rcms[k])
  setwd(paste0(path_2,rcms[k]))
  
  load(paste0(var,"_BC_serie_list_",rcms[k],".rda"))
  list=ls(pattern = "BC_serie_list_")
  
  if( k==1){methods=substr(x = list,stop=nchar(list),start = 15)}
  
  load(paste0(var,"_raw_serie_list_",rcms[k],".rda"))

df_clim_mediano <- map_dfr(seq_along(methods), function(i) {
  
  pp_all <- get(list[[i]])[[paste0(var,"_BC_fut_esc")]]

    reference <-get(list[[i]])[[paste0(var,"_BC_train")]]%>%
    subsetGrid(years=1984:2003)%>%
      climatology()%>%redim(drop = T)

    # Sepaaro los  miembros y combinamos
  map_dfr(1:3, function(j) {
    reference2 <- reference%>%
      subsetDimension(dimension = "member", indice = j)
    
    pp <- pp_all %>%
      subsetDimension(dimension = "member", indice = j) %>%
      # subsetGrid(years = list_years[[j]]) %>%
      subsetGrid(years = list_years[[1]]) %>%
      climatology() 
    
    pp <-gridArithmetics(pp,reference2,operator='/')
    
    data.frame(
      Lon  = rep(pp$xyCoords$x, each = length(pp$xyCoords$y)),
      Lat  = rep(pp$xyCoords$y,  length(pp$xyCoords$x)),
      Data = as.vector(pp$Data[1,,]),
      RCM  = rep(pp[["Members"]], length(as.vector(pp$Data[1,,]))),
      Method = methods[i]
    )
  })
})


df_clim_lejano <- map_dfr(seq_along(methods), function(i) {
  
  pp_all <- get(list[[i]])[[paste0(var,"_BC_fut_esc")]]
 
   reference <-get(list[[i]])[[paste0(var,"_BC_train")]]%>%
    subsetGrid(years=1984:2003)%>%
    climatology()%>%redim(drop = T)
  
  # Sepaaro los  miembros y combinamos
  map_dfr(1:3, function(j) {
    reference2 <- reference%>%
      subsetDimension(dimension = "member", indice = j)
    pp <- pp_all %>%
      subsetDimension(dimension = "member", indice = j) %>%
      # subsetGrid(years = list_years[[j]]) %>%
      subsetGrid(years = list_years[[2]]) %>%
      climatology()
    pp <-gridArithmetics(pp,reference2,operator='/')
    
    data.frame(
      Lon  = rep(pp$xyCoords$x, each = length(pp$xyCoords$y)),
      Lat  = rep(pp$xyCoords$y,  length(pp$xyCoords$x)),
      Data = as.vector(pp$Data[1,,]),
      RCM  = rep(pp[["Members"]], length(as.vector(pp$Data[1,,]))),
      Method = methods[i]
    )
  })
})

####Agrego el raw##################

pp_all <-raw_serie_list[[paste0(var,"_raw_fut_esc")]]

reference <-raw_serie_list[[paste0(var,"_raw_train")]]%>%
  subsetGrid(years=1984:2003)%>%
  climatology()%>%redim(drop = T)

df_clim_raw_mediano <- map_dfr(1:3, function(j) {
 
   reference2 <- reference%>%
    subsetDimension(dimension = "member", indice = j)
  
  pp <- pp_all %>%
    subsetDimension(dimension = "member", indice = j) %>%
    subsetGrid(years = list_years[[1]]) %>%
    climatology()
  pp <-gridArithmetics(pp,reference2,operator='/')
  
  data.frame(
    Lon  = rep(pp$xyCoords$x, each = length(pp$xyCoords$y)),
    Lat  = rep(pp$xyCoords$y,  length(pp$xyCoords$x)),
    Data = as.vector(pp$Data[1,,]),
    RCM  = rep(pp[["Members"]], length(as.vector(pp$Data[1,,]))),
    Method = 'raw'
  )
})

df_clim_raw_lejano <- map_dfr(1:3, function(j) {
  reference2 <- reference%>%
    subsetDimension(dimension = "member", indice = j)
  
  pp <- pp_all %>%
    subsetDimension(dimension = "member", indice = j) %>%
    subsetGrid(years = list_years[[2]]) %>%
    climatology()
  pp <-gridArithmetics(pp,reference2,operator='/')
  
  data.frame(
    Lon  = rep(pp$xyCoords$x, each = length(pp$xyCoords$y)),
    Lat  = rep(pp$xyCoords$y,  length(pp$xyCoords$x)),
    Data = as.vector(pp$Data[1,,]),
    RCM  = rep(pp[["Members"]], length(as.vector(pp$Data[1,,]))),
    Method = 'raw'
  )
})

#######junto ############################
df_clim_all_mediano=rbind(df_clim_raw_mediano,df_clim_mediano)
df_clim_all_lejano=rbind(df_clim_raw_lejano,df_clim_lejano)

df_clim_all_mediano$family=rep(rcms[k],length(df_clim_all_mediano$RCM))
df_clim_all_lejano$family=rep(rcms[k],length(df_clim_all_lejano$RCM))

assign(paste0(var,'_',index,"_serie_list_lejano_", rcms[k]), df_clim_all_lejano)
assign(paste0(var,'_',index,"_serie_list_mediano_", rcms[k]), df_clim_all_mediano)
rm(df_clim_all_lejano,df_clim_all_mediano,list)
rm(list=ls(pattern="BC_serie_list_"))
rm(list=ls(pattern="raw_serie_list"))
rm(list=ls(pattern="df_clim"))
gc()
}


hus_df_clim_all_lejano_mean=rbind(hus_mean_serie_list_lejano_ETA,hus_mean_serie_list_lejano_RegCM,hus_mean_serie_list_lejano_REMO)
hus_df_clim_all_mediano_mean=rbind(hus_mean_serie_list_mediano_ETA,hus_mean_serie_list_mediano_RegCM,hus_mean_serie_list_mediano_REMO)
save(hus_df_clim_all_lejano_mean,file="hus_df_clim_all_lejano_mean.rda")
save(hus_df_clim_all_mediano_mean,file="hus_df_clim_all_mediano_mean.rda")


################################## grafico ########################
require(ggplot2)
require(metR)
# source(paste0(path_2,'mapas_graficos.R'))

my_theme <- theme_bw() + theme(legend.position = "bottom",legend.title = element_text(size=8,colour = "black",face = "bold"),
                               plot.title = element_text(hjust = 0.5,size=7,colour = "black",face = "bold"),axis.text = element_text(size=7,colour = "black"),
                               axis.title = element_text(size=6,face="bold"),legend.key.size = unit(0.1,"line"),legend.text = element_text(size = 7),legend.key.height = unit(0.1,"inch"),legend.key.width = unit(0.1,"inch"),plot.margin = margin(0.15,0.15,0.15,0,"cm")
                               ,axis.text.x = element_text(angle = 360,hjust = 1,size = 6,color = "black",face="bold"),
                               axis.text.y = element_text(hjust = 1,size = 7,color = "black",face="bold"),strip.text.x = element_text(size = 8, colour = "black",face = "bold"),
                               strip.text.y = element_text(size = 6, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())

scale <- scale_fill_gradientn(
  colours = c("#fff7bc","#b8e186", "#4dac26", "#2b83ba", "#053061"),
  limits = c(0, 20),
  # limits = c(0, 4),
  oob = scales::squish,
  breaks = pretty_breaks(n = 10),
  # name = expression(Delta~(g/kg)),na.value = "white",
  name = expression(Delta~("%")),na.value = "white",
  guide = guide_colorsteps(barwidth = 10, barheight = 0.6)
)

df_clim_all_mediano$Method=factor(df_clim_all_mediano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

# gg=ggplot(df_clim_all) +geom_contour_fill(aes(x=Lon,y=Lat,z=Data*1000),bins = 50)+
gg= ggplot(df_clim_all_mediano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs( x = "Longitude",y = "Latitude")

ggsave(gg,filename=paste0(var,"_fut_mediano_",rcms[k],"_porcentual_.png"),device="png",dpi=300,width=7,height = 6,scale=0.7)


##### Lejano #######
df_clim_all_lejano$Method=factor(df_clim_all_lejano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

scale <- scale_fill_gradientn(
  colours = c("#fff7bc","#b8e186", "#4dac26", "#2b83ba", "#053061"),
  limits = c(0, 40),
  # limits = c(5, 25),
  # limits = c(0, 4),
  oob = scales::squish,
  breaks = pretty_breaks(n = 10),
  # name = expression(Delta~(g/kg)),na.value = "white",
  name = expression(Delta~("%")),na.value = "white",
  guide = guide_colorsteps(barwidth = 10, barheight = 0.6)
)


gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Hus', x = "Longitude",y = "Latitude")

ggsave(gg,filename=paste0(var,"_fut_lejano_",rcms[k],"_porcentual_.png"),device="png",dpi=300,width=7,height = 6,scale=0.7)

df_clim_all_lejano$horizonte=rep('lejano',length(df_clim_all_lejano$Lat))
df_clim_all_mediano$horizonte=rep('intermedio',length(df_clim_all_mediano$Lat))

df_clim_all_hor=rbind(df_clim_all_mediano,df_clim_all_lejano)

save(df_clim_all_hor,file=paste0(var,"_fut_df_",rcms[k],".rda"))

# scale=scale_fill_gradientn(
#   # colours = c("#543005", "#8c510a", "#bf812d", "#dfc27d",
#   #             "#f6e8c3", "#f7f7f7", "#c7eae5", "#80cdc1",
#   #             "#35978f", "#01665e", "#003c30"),
#   # limits = c(min(df_clim_all$Data, na.rm=TRUE)*1000, max(df_clim_all$Data, na.rm=TRUE)*1000),
#   
#   colours = c(  "#fdd0a2", "#f16913", "#a63603", "#7f2704","#f7f7f7",
#                  "#c6dbef","#6baed6","#2171b5","#08306b"),
#   # limits = c(-max(df_clim_all_mediano$Data, na.rm=TRUE)*1000, max(df_clim_all_mediano$Data, na.rm=TRUE)*1000),
#   limits = c(-3, 3),
#   oob = scales::squish,breaks=pretty_breaks(n=10),
#   name = expression(Delta~(g/kg~x~1000)),
#   guide = guide_colorsteps(barwidth = 10, barheight = 0.6)
# )
# colfunc<-colorRampPalette(c("#EEDD82","green","#104E8B",'Midnight Blue'))
# colfunc <- colorRampPalette(c("#EEDD82", "#66C2A5", "green", "#104E8B", "MidnightBlue"))
# colfunc <- colorRampPalette(c("#b8e186", "#4dac26", "#2b83ba", "#053061"))

scale <- scale_fill_gradientn(
  colours = c("#fff7bc","#b8e186", "#4dac26", "#2b83ba", "#053061"),
  limits = c(3, 20),
  # limits = c(0, 4),
  oob = scales::squish,
  breaks = pretty_breaks(n = 10),
  # name = expression(Delta~(g/kg)),na.value = "white",
  name = expression(Delta~("%")),na.value = "white",
  guide = guide_colorsteps(barwidth = 10, barheight = 0.6)
)

# 
# scale <- scale_fill_gradientn(
#   colours =color_hus(50),
#   # colours = c(
#   #   "#A1D99B",
#   #   "#eafaea",  # verde muy claro
#   #   # "#a1d99b",  # verde medio
#   #   "#31a354",  # verde oscuro
#   #   # "#c7e9b4",  # verde-agua
#   #   "#7fcdbb",  # turquesa
#   #   "#41b6c4",  # celeste medio
#   #   "#1d91c0",  # celeste fuerte
#   #   "#225ea8",  # azul medio
#   #   "#0c2c84"   # azul profundo
#   # ),
#   limits = c(0, 3),
#   oob = scales::squish,
#   breaks = pretty_breaks(n = 10),
#   name = expression(Delta~(g/kg~x~1000)),
#   guide = guide_colorsteps(barwidth = 10, barheight = 0.6)
# )



############################# diferencia en la senal ##########################

df_diff <-lapply(1:length(methods), function(j){
  pp=subset(df_clim_all_mediano,Method==methods[j])
  pp$Dif=pp$Data-df_clim_raw_mediano$Data
  return(pp)
})%>% do.call(rbind, .)
df_diff_mediano=df_diff

scale=scale_fill_gradient2(low = "#8C510A", mid = "white", high = "#01665E", midpoint = 0,
                           na.value = "white",breaks=pretty_breaks(n=9),
                           limits=c(-10,5),oob = scales::squish,
                           guide = guide_colorsteps(barwidth =0.5, barheight=5))

df_diff$Method=factor(df_diff$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='hus',x = "Longitude",y = "Latitude") +
 guides(fill=guide_legend(reverse = F,title ="%",ncol =12,nrow=1,label.position = "bottom"))

ggsave(gg,filename=paste0(var,"_dif_fut_mediano_",rcms[k],"_porcentual.png"),device="png",dpi=300,width=7,height = 6,scale=0.6)

# scale = scale_fill_distiller(palette = "RdBu",direction = -1,
#                              limits = c(-2.5,2.5),
#                              oob = scales::squish,
#                              breaks = scales::pretty_breaks(n = 10),
#                              name = expression(Delta~(degC)),
#                              guide = guide_colorsteps(barwidth = 10, barheight = 0.6))
# scale=scale_fill_gradient2(
#   low = "blue",   
#   mid = "white",   
#   high = "red", limits=c(-1,2.5),
#   midpoint = 0, oob = scales::squish,guide = guide_colorsteps(barwidth =8, barheight=0.5),  
#   name = expression(Delta~(degC)))
# 


###########################
df_diff <-lapply(1:length(methods), function(j){
  pp=subset(df_clim_all_lejano,Method==methods[j])
  pp$Dif=pp$Data-df_clim_raw_lejano$Data
  return(pp)
})%>% do.call(rbind, .)

df_diff_lejano=df_diff

# scale=scale_fill_gradient2(low = "#8C510A", mid = "white", high = "#01665E", midpoint = 0,
#                            na.value = "white",breaks=seq(-10,5,by=2),
#                            limits=c(-10,5),oob = scales::squish,
#                            guide = guide_colorsteps(barwidth =0.5, barheight=5))

df_diff$Method=factor(df_diff$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='hus',x = "Longitude",y = "Latitude") +
  guides(fill=guide_legend(reverse = F,title ="%",ncol =15,nrow=1,label.position = "bottom"))

ggsave(gg,filename=paste0(var,"_dif_fut_lejano_",rcms[k],"_porcentual.png"),device="png",dpi=300,width=7,height = 6,scale=0.6)



df_diff_lejano$horizonte=rep('lejano',length(df_diff_lejano$Lat))
df_diff_mediano$horizonte=rep('intermedio',length(df_diff_mediano$Lat))

diff_clim_all_hor=rbind(df_diff_mediano,df_diff_lejano)

save(diff_clim_all_hor,file=paste0(var,"_fut_diff_",rcms[k],".rda"))
