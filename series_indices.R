######################### Series temporales ################################3
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
library(tidyr)
library(ggplot2)
library(ggnewscale)
####################################### mapas ########################################
rcms=c("REMO","RegCM",'ETA')
path_2="/media/usuario/9304a642-45f1-44f5-a662-fc30d7f0d1a2/home/usuario/PostDOC/BC/Datos_rcp8.5/"
var="swbgt"
# index="max"
index="Ind30"
source(paste0(path_2,"../grid_to_df_members_serie.R"))
for (k in 1:3) {
  message("Processing RCM: ", rcms[k])
  setwd(paste0(path_2,rcms[k]))
  
  load(paste0(var,"_",index,"_BC_serie_list_",rcms[k],".rda"))
  list=ls(pattern = "_BC_serie_list_")
  
  if( k==1){methods=substr(x = list,stop=nchar(list),start = 21)}
  
  load(paste0(var,"_",index,"_raw_serie_list_",rcms[k],".rda"))

serie_BC=lapply (1:length(list), function(i) {

  pp=grid_to_df_members_serie_means(get(list[[i]])$swbgt_BC_fut_esc)
  pp1=grid_to_df_members_serie_means(get(list[[i]])$swbgt_BC_train_esc)
  pp=rbind(pp,pp1)
  pp$Method=rep(methods[i],length(pp$Time))
  pp$family=rep(rcms[k],length(pp$Time))
  return(pp)
}) %>% rbind()     

serie_BC=bind_rows(serie_BC)

###### raw######
pp=grid_to_df_members_serie_means(swbgt_indice_raw_serie_list$swbgt_raw_train_esc)
pp1=grid_to_df_members_serie_means(swbgt_indice_raw_serie_list$swbgt_raw_fut_esc)
pp=rbind(pp,pp1)
pp$Method=rep('raw',length(pp$Time))
pp$family=rep(rcms[k],length(pp$Time))
serie_raw=pp
rm(pp,pp1)

Serie=rbind(serie_BC,serie_raw)
assign(paste0("Serie_",index,"_",rcms[k]),Serie)
}

###### obs ######
if(index=='Ind30'){
swbgt30.obs <-aggregate_djf_Ind30(calc_swbgt(tas=obs.Tx,hus_grid = obs.huss))%>%
  subsetGrid(years = 1981:2004)

  swbgt30.obs_esc <- scaleGrid(
  grid = swbgt30.obs,
  base = subsetGrid(swbgt30.obs, years = 1985:2004),
  by.member = TRUE
)

}

if(index=='max'){
  swbgt30.obs <-aggregate_djf_max(calc_swbgt(tas=obs.Tx,hus_grid = obs.huss))

  swbgt30.obs_esc <- scaleGrid(
  grid = swbgt30.obs,
  base = subsetGrid(swbgt30.obs, years = 1984:2003),
  by.member = TRUE
)

}
pp=grid_to_df_members_serie_means(bindGrid(swbgt30.obs_esc,swbgt30.obs_esc,dimension='member'))
pp$Method=rep('OBS',length(pp$Time))
pp$Member=as.character(pp$Member)
pp$Member[which(pp$Member=='Member_1')]='OBS'
pp$Data[which(pp$Data==0)]=NA
pp$Member=rep('OBS')
serie_obs=subset(pp,Member=='OBS')
rm(pp)

names(serie_obs)
########

Serie_Ind30=rbind(Serie_Ind30_ETA,Serie_Ind30_RegCM,Serie_Ind30_REMO)
Serie_Ind30_BC=subset(Serie_Ind30,Method %in% methods)
rm(Serie_Ind30_ETA,Serie_Ind30_RegCM,Serie_Ind30_REMO)
gc()
Serie_Ind30$Data[which(Serie_Ind30$Data==0)]=NA


Serie_max=rbind(Serie_max_ETA,Serie_max_RegCM,Serie_max_REMO)
Serie_max_BC=subset(Serie_max,Method %in% methods)
rm(Serie_max_ETA,Serie_max_RegCM,Serie_max_REMO)
gc()

################ Ensamble de modelos #######################3

##Mediana
my.median <- function(x) ifelse( !all(is.na(x)), median(x, na.rm=T), NA)
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)),min(x, na.rm=T), NA)
my.mean <- function(x) ifelse( !all(is.na(x)), mean(x, na.rm=T), NA)


Serie_Ind30_BC_median=aggregate(x=Serie_Ind30_BC$Data,by=list(Serie_Ind30_BC$Member,Serie_Ind30_BC$Time),FUN=my.median)
colnames(Serie_Ind30_BC_median)=c('Model','Model','Data')
gc()

Serie_Ind30=Serie_max
######
Serie_Ind30_median=aggregate(x=Serie_Ind30$Data,by=list(Serie_Ind30$Method,Serie_Ind30$Time,Serie_Ind30$family),FUN=my.median)
colnames(Serie_Ind30_median)=c('Method','Time','RCM','Data')
Serie_Ind30_median$Type='median'
gc()

######
Serie_Ind30_max=aggregate(x=Serie_Ind30$Data,by=list(Serie_Ind30$Method,Serie_Ind30$Time,Serie_Ind30$family),FUN=my.max)
colnames(Serie_Ind30_max)=c('Method','Time','RCM','Data')
Serie_Ind30_max$Type='max'
gc()

######
Serie_Ind30_min=aggregate(x=Serie_Ind30$Data,by=list(Serie_Ind30$Method,Serie_Ind30$Time,Serie_Ind30$family),FUN=my.min)
colnames(Serie_Ind30_min)=c('Method','Time','RCM','Data')
gc()

######
df_all <- Serie_Ind30_median
df_all$min=Serie_Ind30_min$Data
df_all$max=Serie_Ind30_max$Data


grafico=df_all

grafico$Method=factor(grafico$Method,levels=c('OBS','raw','eqm','dqm','qdm',"mbcr","mbcn"))

names(Serie_Ind30)[5]="RCM"

grafico=subset(grafico,Method %in% c('raw','eqm','dqm','mbcn'))

################# suavizado para el grafico #######
library(dplyr)
library(zoo)

grafico_suav <- grafico %>%
  arrange(Method, RCM, Time) %>%   
  group_by(Method, RCM) %>%        
  mutate(Data_suav = rollapply(
    Data,
    width = 5,              # ventana
    FUN = mean,            
    align = "center",       # centrado
    fill = NA               
  )) %>%
  ungroup()

names(serie_obs)[2]='RCM'

obs_suav <- serie_obs %>%
  arrange(Method, RCM, Time) %>%   
  group_by(Method, RCM) %>%        
  mutate(Data_suav = rollapply(
    Data,
    width = 5,              # ventana
    FUN = mean,            
    align = "center",       # centrado
    fill = NA               
  )) %>%
  ungroup()

obs_suav=rbind(obs_suav,obs_suav,obs_suav)
obs_suav$RCM=rep(rcms,each=length(serie_obs$Time))
####################################### escala colores #####
my_scale_fill=scale_color_manual(values = c(
  "raw" = "black",
  "eqm" = "#1f78b4",
  # "qdm" = "#e31a1c",
  # "mbcr" = "#6a3d9a",
  "mbcn" = "#ff7f00",
  "dqm" = "#33a02c"
))

my_scale_color=scale_fill_manual(values = c(
  "raw" = alpha("grey85",0.3),
  "eqm" = alpha("#1f78b480",0.2),
  # "qdm" = alpha("#e31a1c40",0.1),
  # "mbcr" = "#6a3d9a40",
  "mbcn" = alpha("#ff7f0040",0.2),
  "dqm" = alpha("#33a02c40",0.2)
))

scale=c(' (%)',' (degC)')

######## grafico#####
pp=ggplot(grafico_suav, aes(x = Time, color = Method)) +
  
  geom_ribbon(grafico_suav,mapping=aes(ymin = min, ymax = max, fill = Method),
             color = NA) +
  geom_line(aes(y = Data), size = 0.2) +
  my_scale_color+my_scale_fill+
  # scale_color_manual(values=c("black","#1F78B4","#B2DF8A" ,"#33A02C" ,"#FB9A99",'violet'))+
  # scale_fill_manual(values=c("grey80","#1F78B4","#B2DF8A" ,"#33A02C" ,"#FB9A99",'violet'))+
  
   labs(x = "Years", y = paste0("Delta swbgt",index,scale[2])) +
  # geom_line(subset(Serie_Ind30,Member==),aes(y = Data), size = 0.3,color='magenta')
  # geom_line(data=subset(Serie_Ind30,Member=="ETA_HadGEM2" & Method=='raw'), mapping=aes(y = Data,x = Time),size = 0.1,color="orange")+
  # geom_line(data=subset(Serie_Ind30,Member=="RegCM_HadGEM2" & Method=='raw'), mapping=aes(y = Data,x = Time),size = 0.1,color="orange")+
  # geom_line(data=subset(Serie_Ind30,Member=="REMO_HadGEM2" & Method=='raw'), mapping=aes(y = Data,x = Time),size = 0.1,color="orange")+
  geom_line(obs_suav, mapping=aes(y = Data,x = Time), size = 0.2,color="#E31A1C")+
  
  theme_bw() +facet_wrap(.~ RCM)+
  theme(axis.text.x = element_text(size=6),axis.text.y = element_text(size=6),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )+ guides(size = guide_legend(override.aes = list(size = 1))) 
ggsave(pp,filename = paste0("serie_swbgt_",index,".png"),device = "png",dpi=300,width=8,height = 4,scale =0.7)

###########################

list_years=list("cercano"=2040:2059,"lejano"=2080:2098)

Serie_Ind30_l <- Serie_Ind30 %>% 
  filter(Time %in% list_years[[2]]) %>% 
  group_by(Member, Method) %>% 
  summarise(Data = mean(Data, na.rm = TRUE)) %>% 
  ungroup()

Serie_Ind30_c <- Serie_Ind30 %>% 
  filter(Time %in% list_years[[1]]) %>% 
  group_by(Member, Method) %>% 
  summarise(Data = mean(Data, na.rm = TRUE)) %>% 
  ungroup()

names(Serie_Ind30_c)[1]='RCM'
Serie_Ind30_c$Period='2041-2060'
Serie_Ind30_c2=subset(Serie_Ind30_c, Method %in% c('raw','eqm','dqm','mbcn'))

names(Serie_Ind30_l)[1]='RCM'
Serie_Ind30_l$Period='2081-2100'
Serie_Ind30_l2=subset(Serie_Ind30_l, Method %in% c('raw','eqm','dqm','mbcn'))

Serie_Ind30_all=rbind(Serie_Ind30_l2,Serie_Ind30_c2)
########## 
# Serie_Ind30_c$Method=factor(Serie_Ind30_c$Method,levels=c('raw','eqm','dqm','qdm',"mbcr","mbcn"))
# Serie_Ind30_c2$Method=factor(Serie_Ind30_c2$Method,levels=c('raw','eqm','dqm','mbcn'))
Serie_Ind30_all$Method=factor(Serie_Ind30_all$Method,levels=c('raw','eqm','dqm','mbcn'))
#################3
#### escala

my_scale_fill=scale_fill_manual(values = c(
  "raw" = "grey85",
  "eqm" = "#1f78b4",
  # "qdm" = "#e31a1c",
  # "mbcr" = "#6a3d9a",
  "mbcn" = "#ff7f00",
  "dqm" = "#33a02c"
))

pp=ggplot(Serie_Ind30_all, aes(x = RCM, y = Data, fill = Method,group=interaction(Method,RCM),shape=Period)) +
    geom_point(Serie_Ind30_all,mapping=aes(fill = Method),size=2,color='black',stroke=0.3) +
    geom_hline(aes(yintercept = 0),color="grey50",lwd=0.2,lty=2)+ my_scale_fill+
    geom_vline(
      xintercept = seq(1.5, length(unique(Serie_Ind30_all$RCM)) - 0.5, 1),
      color = "grey70",
      lwd = 0.2)+ 
      # new_scale_fill()+
      theme_bw() +scale_shape_manual(values=c(21,23))+ new_scale_fill()+
     # scale_fill_manual(aesthetics = c("fill"),values =c("grey50","#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99", "#E31A1C")) +
    labs( x = "RCMs",y='%') +#scale_y_continuous(breaks=breaks[[j]],limits=scales[[j]]) + 
    theme(text=element_text(face="bold",size=7),axis.text.x = element_text(angle = 45, hjust=1,face="bold"),
          strip.text.x = element_text(size = 6, colour = "black",face = "bold"),axis.text.y = element_text(size = 6,color = "black",face="bold"),axis.title.x = element_blank(),
          strip.text.y = element_text(size = 9, colour = "black",face = "bold"),strip.background.x = element_blank(),strip.background.y = element_blank(),
          legend.margin = margin(0,0,0,0,"cm"),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank())#+ facet_grid(Period~.,scales = 'free_y')
  
ggsave(pp,filename = paste0("swbgt_",index,".png"),device = "png",dpi=300,width=8,height = 4,scale =0.6)
