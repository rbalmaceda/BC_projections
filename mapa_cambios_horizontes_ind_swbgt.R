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
setwd(paste0(path_2,rcms[k]))
var="swbgt"
# index="max"
index="Ind30"
porcentual=FALSE

for (k in 1:3) {
message("Processing RCM: ", rcms[k])
setwd(paste0(path_2,rcms[k]))
  
load(paste0(var,"_",index,"_BC_serie_list_",rcms[k],".rda"))
list=ls(pattern = "_BC_serie_list_")

if( k==1){methods=substr(x = list,stop=nchar(list),start = 21)}

load(paste0(var,"_",index,"_raw_serie_list_",rcms[k],".rda"))

df_clim_mediano <- map_dfr(seq_along(methods), function(i) {

  pp_all <- get(list[[i]])[[paste0(var,"_BC_fut_esc")]]

  if (porcentual){
    reference <-get(list[[i]])[[paste0(var,"_BC_train")]]%>%
    subsetGrid(years=1984:2003)%>%
      climatology()%>%redim(drop = T)}

    # Sepaaro los  miembros y combinamos
  map_dfr(1:3, function(j) {
    if (porcentual){
    reference2 <- reference%>%
      subsetDimension(dimension = "member", indice = j)}
    # 
    pp <- pp_all %>%
      subsetDimension(dimension = "member", indice = j) %>%
      # subsetGrid(years = list_years[[j]]) %>%
      subsetGrid(years = list_years[[1]]) %>%
      climatology() 
   if(porcentual){ 
    pp <-gridArithmetics(pp,reference2,operator='/')}
    
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
 
  if(porcentual){
   reference <-get(list[[i]])[[paste0(var,"_BC_train")]]%>%
    subsetGrid(years=1984:2003)%>%
    climatology()%>%redim(drop = T)
  }
  # Sepaaro los  miembros y combinamos
  map_dfr(1:3, function(j) {
    if(porcentual){
      reference2 <- reference%>%
      subsetDimension(dimension = "member", indice = j)}
    
    pp <- pp_all %>%
      subsetDimension(dimension = "member", indice = j) %>%
      # subsetGrid(years = list_years[[j]]) %>%
      subsetGrid(years = list_years[[2]]) %>%
      climatology()
    
    if(porcentual){
      pp <-gridArithmetics(pp,reference2,operator='/')}
    
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

# raw_serie_list=swbgt_raw_serie_list
# names(swbgt_indice_raw_serie_list)=c("swbgt_raw_train","swbgt_raw_train_esc", "swbgt_raw_fut", "swbgt_raw_fut_esc")
pp_all <-swbgt_indice_raw_serie_list[[paste0(var,"_raw_fut_esc")]]

if(porcentual){
  
reference <-swbgt_indice_raw_serie_list[[paste0(var,"_raw_train")]]%>%
  subsetGrid(years=1984:2003)%>%
  climatology()%>%redim(drop = T)

}
df_clim_raw_mediano <- map_dfr(1:3, function(j) {
 
  if(porcentual){
   reference2 <- reference%>%
    subsetDimension(dimension = "member", indice = j)}

  pp <- pp_all %>%
    subsetDimension(dimension = "member", indice = j) %>%
    subsetGrid(years = list_years[[1]]) %>%
    climatology()
  if(porcentual){
  pp <-gridArithmetics(pp,reference2,operator='/')
  }
  data.frame(
    Lon  = rep(pp$xyCoords$x, each = length(pp$xyCoords$y)),
    Lat  = rep(pp$xyCoords$y,  length(pp$xyCoords$x)),
    Data = as.vector(pp$Data[1,,]),
    RCM  = rep(pp[["Members"]], length(as.vector(pp$Data[1,,]))),
    Method = 'raw'
  )
})

df_clim_raw_lejano <- map_dfr(1:3, function(j) {
  if(porcentual){
    
  reference2 <- reference%>%
    subsetDimension(dimension = "member", indice = j)
  }
  pp <- pp_all %>%
    subsetDimension(dimension = "member", indice = j) %>%
    subsetGrid(years = list_years[[2]]) %>%
    climatology()
  
  if(porcentual){
    
  pp <-gridArithmetics(pp,reference2,operator='/')
  }
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

####### Uno para calcular los ensambles #######
# df_clim_all_lejano=rbind(swbgt_max_serie_list_lejano_ETA,swbgt_max_serie_list_lejano_RegCM,swbgt_max_serie_list_lejano_REMO)
# df_clim_all_mediano=rbind(swbgt_max_serie_list_mediano_ETA,swbgt_max_serie_list_mediano_RegCM,swbgt_max_serie_list_mediano_REMO)

df_clim_all_lejano=rbind(swbgt_Ind30_serie_list_lejano_ETA,swbgt_Ind30_serie_list_lejano_RegCM,swbgt_Ind30_serie_list_lejano_REMO)
df_clim_all_mediano=rbind(swbgt_Ind30_serie_list_mediano_ETA,swbgt_Ind30_serie_list_mediano_RegCM,swbgt_Ind30_serie_list_mediano_REMO)

df_clim_all_lejano$Data[which(df_clim_all_lejano$Data==0)] <- NaN
df_clim_all_mediano$Data[which(df_clim_all_mediano$Data==0)] <- NaN

################################## grafico ########################
require(sf)
require(ggplot2)
require(scales)
require(shadowtext)
require(hrbrthemes)
require(maps)
library("rnaturalearth")
library("rnaturalearthdata")
require(RColorBrewer)
require(metR)

regions <- sf::read_sf("https://raw.githubusercontent.com/IPCC-WG1/Atlas/devel/reference-regions/IPCC-WGI-reference-regions-v4.geojson")
regions$area <- as.numeric(sf::st_area(regions))
regions_piola=subset(regions,Continent=="SOUTH-AMERICA")  
world <- ne_countries(scale = "medium", returnclass = "sf")

my_theme <- theme_bw() + theme(legend.position = "bottom",legend.title = element_text(size=8,colour = "black",face = "bold"),
                               plot.title = element_text(hjust = 0.5,size=7,colour = "black",face = "bold"),axis.text = element_text(size=7,colour = "black"),
                               axis.title = element_text(size=6,face="bold"),legend.key.size = unit(0.1,"line"),legend.text = element_text(size = 7),legend.key.height = unit(0.1,"inch"),legend.key.width = unit(0.1,"inch"),plot.margin = margin(0.15,0.15,0.15,0,"cm")
                               ,axis.text.x = element_text(angle = 360,hjust = 1,size = 6,color = "black",face="bold"),
                               axis.text.y = element_text(hjust = 1,size = 7,color = "black",face="bold"),strip.text.x = element_text(size = 8, colour = "black",face = "bold"),
                               strip.text.y = element_text(size = 6, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())

####swbgt##############

# colfunc<-colorRampPalette(c('green',"springgreen","yellow",'orange','red'))
# colfunc <- colorRampPalette(c("#66c2a5","#b8e186", "#fee08b", "#f46d43", "#d73027"))
# colfunc <- colorRampPalette(c("#66c2a5", "#fee08b", "#f46d43", "#d73027"))
# colfunc <- colorRampPalette(c(
#   "#006837",  # verde oscuro adicional
#   "#1a9850",  # verde medio
#   "#66c2a5",  # verde agua
#   "#b8e186",  # verde amarillento
#   "#fee08b",  # amarillo suave
#   "#fdae61",  # naranja claro
#   "#f46d43",  # naranja fuerte
#   "#d73027"   # rojo oscuro
# ))
# 
# colfunc <- colorRampPalette(c(
#   "#006837",  # verde oscuro
#   "#1a9850",  # verde medio
#   "#66c2a5",  # verde agua
#   "#b8e186",  # verde amarillento
#   "#fee08b",  # amarillo suave
#   "#fdae61",  # naranja claro
#   "#f46d43",  # naranja fuerte
#   "#d73027",  # rojo oscuro
#   "#762a83"   # violeta intenso para valores extremos
# ))

# colfunc <- colorRampPalette(c(
#   "#a6d96a",  # verde claro
#   "#d9ef8b",  # verde-amarillento
#   "#ffffbf",  # amarillo suave
#   "#fee08b",  # amarillo-anaranjado
#   "#fdae61",  # naranja medio
#   "#f46d43",  # naranja rojizo
#   "#d73027",  # rojo fuerte
#   "#b2182b",  # rojo oscuro
#   "#762a83"   # violeta profundo
# ))
###### escala de colores elegida ####
colfunc <- colorRampPalette(c(
  "#a6d96a",  # verde claro
  "#d9ef8b",  # verde-amarillento
  "#ffffbf",  # amarillo suave
  "#fee08b",  # amarillo-anaranjado
  "#fdae61",  # naranja medio
  "#f46d43",  # naranja rojizo
  "#d73027",  # rojo fuerte
  "#e082c1",  # violeta claro
  "#762a83"   # violeta oscuro
))

################################### OBS############################
path_1="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_historico/"
require(HeatStress)
swbgt30.obs <-aggregate_djf_Ind30(calc_swbgt(tas=obs.Tx,hus_grid = obs.huss))
df_obs=grid_to_df_members(climatology(bindGrid(swbgt30.obs,swbgt30.obs,dimension='member')))
df_obs$Member=as.character(df_obs$Member)
df_obs$Member[which(df_obs$Member=='Member_1')]='OBS'
df_obs$Data[which(df_obs$Data==0)]=NA
quantile(df_obs$Data,na.rm=T,probs=c(0.025,0.975))

my_scale <- scale_fill_gradientn(limits=c(0,35),colours=colfunc(7),#(brewer.pal(n=9,"YlOrRd")),
                                 na.value = "white",aesthetics = c("fill"),#,breaks=seq(0,60,by=10)
                                 name = expression('% days in summer'),breaks=seq(0,35,by=5),
                                 oob = scales::squish,
                                 guide = guide_colorsteps(barwidth = 10, barheight = 0.6,show.limits = T,ticks = 50))

df$Member='OBS'
gg= ggplot(subset(df_obs,Member=='OBS')) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
  facet_grid(.~ Member) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  scale_x_continuous(breaks = c(-65,-58,-50)) +
  scale_y_continuous(breaks = c(-35,-30,-25,-20)) +my_scale +my_theme+ 
  labs(main=var, x = "Longitude",y = "Latitude")+
  geom_point(data=df,aes(x=Lon,y=Lat),color ='black')+ geom_text(data=df,aes(x=Lon,y=Lat,label=label),nudge_x = 0.9,color ='black')
gg
ggsave(gg,filename=paste0(var,"_",index,"_OBS.png"),device="png",dpi=300,width=5,height = 6,scale=0.7)
###################################
# df_clim_all_mediano$Data[!is.finite(df_clim_all_mediano$Data)] <- NaN
df_clim_all_mediano$Method=factor(df_clim_all_mediano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

quantile(df_clim_all_mediano$Data,na.rm=T,probs=c(0.025,0.975))
quantile(df_clim_all_lejano$Data,na.rm=T,probs=c(0.025,0.975))

my_scale <- scale_colour_gradientn(limits=c(0,60),colours=colfunc(9),#(brewer.pal(n=9,"YlOrRd")),
                                   na.value = "white",breaks=seq(0,60,by=10),aesthetics = c("fill"),
                                   name = expression(Delta~('% days in summer')),
                                   oob = scales::squish,
                                   guide = guide_colorsteps(barwidth = 10, barheight = 0.6,show.limits = T))


# gg=ggplot(df_clim_all_mediano) +geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
gg= ggplot(df_clim_all_mediano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  scale_x_continuous(breaks = c(-65,-58,-50)) +
  scale_y_continuous(breaks = c(-35,-30,-25,-20)) +my_scale +my_theme+ 
  theme(legend.key.width = unit(0.2,"inch"),legend.margin = margin(0, 0, 0, 0),legend.position = 'right',
                       )+
  labs(main='Tx',x = "Longitude",y = "Latitude") +
  guides(fill=guide_legend(reverse = F,title ="% days in summer",ncol =1,nrow=10,label.position = "right"))

gg
ggsave(gg,filename=paste0(var,"_",index,"_fut_mediano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.7)
ggsave(gg,filename=paste0(var,"_",index,"_fut_mediano_all.png"),device="png",dpi=300,width=6,height = 8,scale=0.9)


##### Lejano #######
# my_scale = scale_fill_gradientn(colors= colfunc(50),#palette = "YlOrRd",direction = 1,
#                              limits = c(0, 5),
#                              oob = scales::squish,
#                              breaks = seq(0,5,by=0.5),
#                              name = expression(Delta~(degC)),
#                              guide = guide_colorsteps(barwidth = 10, barheight = 0.6,show.limits = T))
df_clim_all_lejano$Data[which(df_clim_all_lejano$Data==0)] <- NaN
quantile(df_clim_all_lejano$Data,na.rm=T,probs=c(0.025,0.975))

# my_scale = scale_fill_gradientn(colors= colfunc(50),#palette = "YlOrRd",direction = 1,
#                                 limits = c(0, 6),
#                                 oob = scales::squish,
#                                 breaks = seq(0,5,by=0.5),
#                                 name = expression(Delta~(degC)),
#                                 guide = guide_colorsteps(barwidth = 10, barheight = 0.6,show.limits = T))


df_clim_all_lejano$Method=factor(df_clim_all_lejano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
# gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +my_scale +my_theme+ 
  labs( x = "Longitude",y = "Latitude")
gg
ggsave(gg,filename=paste0(var,"_",index,"indice_fut_lejano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.7)

# df_clim_all_lejano$horizonte=rep('lejano',length(df_clim_all_lejano$Lat))
# df_clim_all_mediano$horizonte=rep('intermedio',length(df_clim_all_mediano$Lat))
# 
# df_clim_all_hor=rbind(df_clim_all_mediano,df_clim_all_lejano)
# 
# save(df_clim_all_hor,file=paste0(var,"_fut_df_",rcms[k],".rda"))

############################# diferencia en la senal mediano ##########################

df_clim_all_mediano=rbind(df_clim_raw_mediano,df_clim_mediano)

df_clim_all_lejano=rbind(df_clim_raw_lejano,df_clim_lejano)

df_diff <-lapply(1:length(methods), function(j){
  pp=subset(df_clim_all_mediano,Method==methods[j])
  pp$Dif=pp$Data-df_clim_raw_mediano$Data
  return(pp)
})%>% do.call(rbind, .)

df_diff$Method=factor(df_diff$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

quantile(df_diff$Dif,na.rm=T,probs=c(0.025,0.975))


# scale=scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
#                            na.value = "white",breaks=pretty_breaks(n=9),
#                            limits=c(-20,9),oob = scales::squish,
#                            guide = guide_colorsteps(barwidth =0.5, barheight=5))

scale=scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           na.value = "white",breaks=pretty_breaks(n=9),
                           limits=c(-30,10),oob = scales::squish,#values = scales::rescale(c(-40, 10, 5)),
                           guide = guide_colorsteps(barwidth =0.5, barheight=5))


# gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Tx',x = "Longitude",y = "Latitude") +
 guides(fill=guide_legend(reverse = F,title ="days %",ncol =15,nrow=1,label.position = "bottom"))
gg
ggsave(gg,filename=paste0(var,"_",index,"_dif_fut_mediano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.6)

# scale = scale_fill_distiller(palette = "RdBu",direction = -1,
#                              limits = c(-2.5,2.5),
#                              oob = scales::squish,
#                              breaks = scales::pretty_breaks(n = 10),
#                              name = expression(Delta~(degC)),
#                              guide = guide_colorsteps(barwidth = 10, barheight = 0.6))


###########################
df_diff <-lapply(1:length(methods), function(j){
  pp=subset(df_clim_all_lejano,Method==methods[j])
  pp$Dif=pp$Data-df_clim_raw_lejano$Data
  return(pp)
})%>% do.call(rbind, .)

df_diff_lejano=df_diff

df_diff$Method=factor(df_diff$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))
df_diff$Data[which(df_diff$Data==0)] <- NaN

quantile(df_diff$Dif,na.rm=T,probs=c(0.025,0.975))

# gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Tx',x = "Longitude",y = "Latitude") +
  guides(fill=guide_legend(reverse = F,title ="% days in summer",ncol =15,nrow=1,label.position = "bottom"))
gg
ggsave(gg,filename=paste0(var,"_",index,"_ind_dif_fut_lejano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.6)



df_diff_lejano$horizonte=rep('lejano',length(df_diff_lejano$Lat))
df_diff_mediano$horizonte=rep('intermedio',length(df_diff_mediano$Lat))

diff_clim_all_hor=rbind(df_diff_mediano,df_diff_lejano)

save(diff_clim_all_hor,file=paste0(var,"_fut_diff_",rcms[k],".rda"))
