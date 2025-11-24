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
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"
rcms=c("REMO","RegCM",'ETA')
var="swbgt"
index='mean'
for (k in 1:3) {
  message("Processing RCM: ", rcms[k])
  setwd(paste0(path_2,rcms[k]))
  
  load(paste0(var,"_BC_serie_list_",rcms[k],".rda"))
  list=ls(pattern = "_BC_serie_list_")
  
  if( k==1){methods=substr(x = list,stop=nchar(list),start = 21)}
  
  load(paste0(var,"_raw_serie_list_",rcms[k],".rda"))
  


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
# raw_serie_list=swbgt_raw_serie_list
names(swbgt_raw_serie_list)=c("swbgt_raw_train","swbgt_raw_train_esc", "swbgt_raw_fut", "swbgt_raw_fut_esc")
pp_all <-swbgt_raw_serie_list[[paste0(var,"_raw_fut_esc")]]
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

swbgt_df_clim_all_lejano_mean=rbind(swbgt_mean_serie_list_lejano_ETA,swbgt_mean_serie_list_lejano_RegCM,swbgt_mean_serie_list_lejano_REMO)
swbgt_df_clim_all_mediano_mean=rbind(swbgt_mean_serie_list_mediano_ETA,swbgt_mean_serie_list_mediano_RegCM,swbgt_mean_serie_list_mediano_REMO)
save(swbgt_df_clim_all_lejano_mean,file="swbgt_df_clim_all_lejano_mean.rda")
save(swbgt_df_clim_all_mediano_mean,file="swbgt_df_clim_all_mediano_mean.rda")


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

####swbgt##############
df_clim_all_mediano$Method=factor(df_clim_all_mediano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

colfunc<-colorRampPalette(c('green',"springgreen","yellow",'orange','red'))
colfunc <- colorRampPalette(c("#66c2a5","#b8e186", "#fee08b", "#f46d43", "#d73027"))
# colfunc <- colorRampPalette(c("#66c2a5", "#fee08b", "#f46d43", "#d73027"))
colfunc <- colorRampPalette(c(
  "#006837",  # verde oscuro adicional
  "#1a9850",  # verde medio
  "#66c2a5",  # verde agua
  "#b8e186",  # verde amarillento
  "#fee08b",  # amarillo suave
  "#fdae61",  # naranja claro
  "#f46d43",  # naranja fuerte
  "#d73027"   # rojo oscuro
))

colfunc <- colorRampPalette(c(
  "#006837",  # verde oscuro
  "#1a9850",  # verde medio
  "#66c2a5",  # verde agua
  "#b8e186",  # verde amarillento
  "#fee08b",  # amarillo suave
  "#fdae61",  # naranja claro
  "#f46d43",  # naranja fuerte
  "#d73027",  # rojo oscuro
  "#762a83"   # violeta intenso para valores extremos
))

colfunc <- colorRampPalette(c(
  "#a6d96a",  # verde claro
  "#d9ef8b",  # verde-amarillento
  "#ffffbf",  # amarillo suave
  "#fee08b",  # amarillo-anaranjado
  "#fdae61",  # naranja medio
  "#f46d43",  # naranja rojizo
  "#d73027",  # rojo fuerte
  "#b2182b",  # rojo oscuro
  "#762a83"   # violeta profundo
))

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




my_scale <- scale_colour_gradientn(limits=c(0,3),colours=colfunc(50),#(brewer.pal(n=9,"YlOrRd")),
                                   na.value = "white",breaks=seq(0,3,by=0.5),aesthetics = c("fill"),
                                   name = expression(Delta~(degC)),oob = scales::squish,values = scales::rescale(c(-1, 0, 3)),
                                   guide = guide_colorsteps(barwidth = 10, barheight = 0.6,show.limits = T))

# gg=ggplot(df_clim_all_mediano) +geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
gg= ggplot(df_clim_all_mediano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +my_scale +my_theme+ 
  labs(main=var, x = "Longitude",y = "Latitude")
gg
ggsave(gg,filename=paste0(var,"_fut_mediano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.7)


##### Lejano #######
my_scale = scale_fill_gradientn(colors= colfunc(50),#palette = "YlOrRd",direction = 1,
                             limits = c(0, 5),
                             oob = scales::squish,
                             breaks = seq(0,5,by=0.5),
                             name = expression(Delta~(degC)),
                             guide = guide_colorsteps(barwidth = 10, barheight = 0.6,show.limits = T))

df_clim_all_lejano$Method=factor(df_clim_all_lejano$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data),bins = 50)+
# gg= ggplot(df_clim_all_lejano) + geom_contour_fill(aes(x=Lon,y=Lat,z=Data*100),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +my_scale +my_theme+ 
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

quantile(df_diff_mediano$Dif,na.rm=T,probs=c(0.025,0.975))


scale=scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           na.value = "white",breaks=pretty_breaks(n=9),
                           limits=c(-0.75,0.5),oob = scales::squish,
                           guide = guide_colorsteps(barwidth =0.5, barheight=5))

# gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Tx',x = "Longitude",y = "Latitude") +
 guides(fill=guide_legend(reverse = F,title ="degC",ncol =15,nrow=1,label.position = "bottom"))
gg
ggsave(gg,filename=paste0(var,"_dif_fut_mediano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.6)

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
gg<- ggplot(df_diff) + geom_contour_fill(aes(x=Lon,y=Lat,z=Dif),bins = 50)+
  facet_grid(RCM~ Method) +geom_sf(data = world,fill=NA,lwd=0.2)+
  coord_sf(xlim=c(-65,-44),ylim=c(-39,-17),expand = FALSE)+ 
  xlab("longitude") + ylab("latitude")+scale_x_continuous(breaks = seq(-65, -44, by = 15))+
  scale_y_continuous(breaks = seq(-40,-18, by = 10)) +scale +my_theme+ 
  labs(main='Tx',x = "Longitude",y = "Latitude") +
  guides(fill=guide_legend(reverse = F,title ="degC",ncol =15,nrow=1,label.position = "bottom"))
gg
ggsave(gg,filename=paste0(var,"_dif_fut_lejano_",rcms[k],".png"),device="png",dpi=300,width=7,height = 6,scale=0.6)



df_diff_lejano$horizonte=rep('lejano',length(df_diff_lejano$Lat))
df_diff_mediano$horizonte=rep('intermedio',length(df_diff_mediano$Lat))

diff_clim_all_hor=rbind(df_diff_mediano,df_diff_lejano)

save(diff_clim_all_hor,file=paste0(var,"_fut_diff_",rcms[k],".rda"))


################################## boxplot ###########################

p <- ggplot(data=data,aes(x=Simulation,y=Index*100,fill=Method,group=interaction(Method,Simulation)))+ stat_summary(fun.data = f, geom="boxplot",na.rm=T,show.legend = T,lwd=0.25,position="dodge")
p <-p +  theme_bw() + theme(text=element_text(face="bold",size=8),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face="bold",color = "black"),
                            strip.text.x = element_text(size = 8, colour = "black",face = "bold"),axis.text.y = element_text(size = 6,color = "black"),axis.title.x = element_blank(),
                            strip.text.y = element_text(size = 9, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank())
p <- p    + geom_hline(aes(yintercept=0),color='grey30',lwd=0.2)+ ylab('Delta (degC)') 
p <- p + scale_fill_manual(aesthetics = c("fill"),values = values=c("grey95","#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C")) #+ facet_grid_sc(Metric~., scales = list(y = scales_y))

swbgt_df_clim_all_mediano_mean$Method=factor(swbgt_df_clim_all_mediano_mean$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))

ggplot(swbgt_df_clim_all_mediano_mean, aes(x = RCM, y = Data, fill = Method,color=Method,group=interaction(Method,RCM))) +
  geom_boxplot(outlier.size = 0.1, alpha = 0.3,outlier.shape = 20,color = "grey15",lwd=0.25,width=0.9) +
  theme_bw() + scale_fill_manual(aesthetics = c("fill"),values =c("grey50","#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C")) +
   labs( x = "RCMs", y = "Delta(degC)") +scale_y_continuous(breaks=seq(0,4,1),limits=c(0,4)) +
   theme(text=element_text(face="bold",size=8),axis.text.x = element_text(angle = 45, hjust=1,face="bold"),
          strip.text.x = element_text(size = 8, colour = "black",face = "bold"),axis.text.y = element_text(size = 6,color = "black"),axis.title.x = element_blank(),
          strip.text.y = element_text(size = 9, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank(),
         text = element_text(size = 8, family = "AvantGarde"))



p <- ggplot(data=data,aes(x=Simulation,y=Index*100,fill=Method,group=interaction(Method,Simulation)))+ stat_summary(fun.data = f, geom="boxplot",na.rm=T,show.legend = T,lwd=0.25,position="dodge")
p <-p +  theme_bw() + theme(text=element_text(face="bold",size=8),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face="bold",color = "black"),
                            strip.text.x = element_text(size = 8, colour = "black",face = "bold"),axis.text.y = element_text(size = 6,color = "black"),axis.title.x = element_blank(),
                            strip.text.y = element_text(size = 9, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank())
p <- p    + geom_hline(aes(yintercept=0),color='grey30',lwd=0.2)+ ylab('biasRel(%)') 
p <- p + scale_fill_manual(aesthetics = c("fill"),values = paleta) + facet_grid_sc(Metric~., scales = list(y = scales_y))
# p <- p + theme_bw() + theme(strip.background = element_rect(fill="grey92"),axis.text.y = element_text(color = "black",size = 8,face = "bold"),legend.position = "bottom",legend.title = element_text(size=6),axis.text.x = element_text(angle = 0,hjust = 0.5,size =9 ,color = "black",face ="bold" ),axis.title.x = element_blank(),
#                             strip.text.x = element_text(face="bold",size=9),strip.text.y = element_blank(),legend.text = element_text(face = "bold"))
