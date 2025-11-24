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
k=3
setwd(paste0(path_2,rcms[k]))
# load('HR_BC_serie_list_REMO.rda');load('HR_raw_serie_list_REMO.rda')
list=ls(pattern = "_BC_serie_list_")
methods=substr(x = list,stop=nchar(list),start = 18)
var="HR"

df_clim_mediano <- map_dfr(seq_along(methods), function(i) {
  
  pp_all <- get(list[[i]])[[paste0(var,"_BC_fut_esc")]]

    # reference <-get(list[[i]])[[paste0(var,"_BC_train")]]%>%
    # subsetGrid(years=1984:2003)%>%
    #   climatology()%>%redim(drop = T)

    # Sepaaro los  miembros y combinamos
  map_dfr(1:3, function(j) {
    # reference2 <- reference%>%
    #   subsetDimension(dimension = "member", indice = j)
    # 
    pp <- pp_all %>%
      subsetDimension(dimension = "member", indice = j) %>%
      # subsetGrid(years = list_years[[j]]) %>%
      subsetGrid(years = list_years[[1]]) %>%
      climatology() 
    
    # pp <-gridArithmetics(pp,reference2,operator='/')
    
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
 
   # reference <-get(list[[i]])[[paste0(var,"_BC_train")]]%>%
   #  subsetGrid(years=1984:2003)%>%
   #  climatology()%>%redim(drop = T)
  
  # Sepaaro los  miembros y combinamos
  map_dfr(1:3, function(j) {
    # reference2 <- reference%>%
    #   subsetDimension(dimension = "member", indice = j)
    pp <- pp_all %>%
      subsetDimension(dimension = "member", indice = j) %>%
      # subsetGrid(years = list_years[[j]]) %>%
      subsetGrid(years = list_years[[2]]) %>%
      climatology()
    # pp <-gridArithmetics(pp,reference2,operator='/')
    
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
# raw_serie_list=HR_raw_serie_list
# names(raw_serie_list)=c("HR_raw_train","HR_raw_train_esc", "HR_raw_fut", "HR_raw_fut_esc")  
pp_all <-HR_raw_serie_list[[paste0(var,"_raw_fut_esc")]]
# reference <-raw_serie_list[[paste0(var,"_raw_train")]]%>%
#   subsetGrid(years=1984:2003)%>%
#   climatology()%>%redim(drop = T)

df_clim_raw_mediano <- map_dfr(1:3, function(j) {
 
   # reference2 <- reference%>%
   #  subsetDimension(dimension = "member", indice = j)

  pp <- pp_all %>%
    subsetDimension(dimension = "member", indice = j) %>%
    subsetGrid(years = list_years[[1]]) %>%
    climatology()
  # pp <-gridArithmetics(pp,reference2,operator='/')
  
  data.frame(
    Lon  = rep(pp$xyCoords$x, each = length(pp$xyCoords$y)),
    Lat  = rep(pp$xyCoords$y,  length(pp$xyCoords$x)),
    Data = as.vector(pp$Data[1,,]),
    RCM  = rep(pp[["Members"]], length(as.vector(pp$Data[1,,]))),
    Method = 'raw'
  )
})

df_clim_raw_lejano <- map_dfr(1:3, function(j) {
  # reference2 <- reference%>%
  #   subsetDimension(dimension = "member", indice = j)
  
  pp <- pp_all %>%
    subsetDimension(dimension = "member", indice = j) %>%
    subsetGrid(years = list_years[[2]]) %>%
    climatology()
  # pp <-gridArithmetics(pp,reference2,operator='/')
  
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

################################## grafico ########################
require(ggplot2)
require(metR)
my_theme <- theme_bw() + theme(legend.position = "bottom",legend.title = element_text(size=8,colour = "black",face = "bold"),
                               plot.title = element_text(hjust = 0.5,size=7,colour = "black",face = "bold"),axis.text = element_text(size=7,colour = "black"),
                               axis.title = element_text(size=6,face="bold"),legend.key.size = unit(0.1,"line"),legend.text = element_text(size = 7),legend.key.height = unit(0.1,"inch"),legend.key.width = unit(0.1,"inch"),plot.margin = margin(0.15,0.15,0.15,0,"cm")
                               ,axis.text.x = element_text(angle = 360,hjust = 1,size = 6,color = "black",face="bold"),
                               axis.text.y = element_text(hjust = 1,size = 7,color = "black",face="bold"),strip.text.x = element_text(size = 8, colour = "black",face = "bold"),
                               strip.text.y = element_text(size = 6, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())

####HR##############
# colfunc<-colorRampPalette(c("springgreen","yellow",'orange','red'))
colfunc <- colorRampPalette(c("#b8e186", "#fee08b", "#f46d43", "#d73027"))
# colfunc <- colorRampPalette(c("#66c2a5", "#fee08b", "#f46d43", "#d73027"))
df_clim_all_mediano$Method=factor(df_clim_all_mediano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

scale = scale_fill_distiller(palette = "PiYG",direction = 1,
                             limits = c(-5, 10),
                             values = scales::rescale(c(-5, 0, 10)),
                             # limits = c(-20, 10),
                             # values = scales::rescale(c(-20, 0, 10)),
                             oob = scales::squish,
                             breaks = scales::pretty_breaks(n = 10),
                             name = expression(Delta~("%")),
                             guide = guide_colorsteps(barwidth = 10, barheight = 0.6))


# gg=ggplot(df_clim_all_mediano) +geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
gg= ggplot(df_clim_all_mediano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main=var, x = "Longitude",y = "Latitude")
gg
ggsave(gg,filename=paste0(var,"_fut_mediano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.7)


##### Lejano #######
# scale = scale_fill_gradientn(colors= colfunc(50),#palette = "YlOrRd",direction = 1,
#                              limits = c(0, 6),
#                              oob = scales::squish,
#                              breaks = scales::pretty_breaks(n = 10),
#                              name = expression(Delta~(degC)),
#                              guide = guide_colorsteps(barwidth = 10, barheight = 0.6))

df_clim_all_lejano$Method=factor(df_clim_all_lejano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
# gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs( x = "Longitude",y = "Latitude")
gg
ggsave(gg,filename=paste0(var,"_fut_lejano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.7)

df_clim_all_lejano$horizonte=rep('lejano',length(df_clim_all_lejano$Lat))
df_clim_all_mediano$horizonte=rep('intermedio',length(df_clim_all_mediano$Lat))

df_clim_all_hor=rbind(df_clim_all_mediano,df_clim_all_lejano)

save(df_clim_all_hor,file=paste0(var,"_fut_df_",rcms[k],".rda"))

############################# diferencia en la senal ##########################33

df_diff <-lapply(1:length(methods), function(j){
  pp=subset(df_clim_all_mediano,Method==methods[j])
  pp$Dif=pp$Data-df_clim_raw_mediano$Data
  return(pp)
})%>% do.call(rbind, .)

df_diff_mediano=df_diff
df_diff$Method=factor(df_diff$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

scale=scale_fill_gradient2(low = "#8C510A", mid = "white", high = "#01665E", midpoint = 0,
                           na.value = "white",breaks=pretty_breaks(n=9),#breaks=seq(-10,5,by=2),
                           limits=c(-10,5),oob = scales::squish,
                           guide = guide_colorsteps(barwidth =0.5, barheight=5))



# gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Tx',x = "Longitude",y = "Latitude") +
 guides(fill=guide_legend(reverse = F,title ="%",ncol =12,nrow=1,label.position = "bottom"))
gg
ggsave(gg,filename=paste0(var,"_dif_fut_mediano_",rcms[k],"_porcentual.png"),device="png",dpi=300,width=7,height = 6,scale=0.6)

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

# gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Tx',x = "Longitude",y = "Latitude") +
  guides(fill=guide_legend(reverse = F,title ="%",ncol =15,nrow=1,label.position = "bottom"))
gg
ggsave(gg,filename=paste0(var,"_dif_fut_lejano_",rcms[k],"_porcentual.png"),device="png",dpi=300,width=7,height = 6,scale=0.6)



df_diff_lejano$horizonte=rep('lejano',length(df_diff_lejano$Lat))
df_diff_mediano$horizonte=rep('intermedio',length(df_diff_mediano$Lat))

diff_clim_all_hor=rbind(df_diff_mediano,df_diff_lejano)

save(diff_clim_all_hor,file=paste0(var,"_fut_diff_",rcms[k],".rda"))
