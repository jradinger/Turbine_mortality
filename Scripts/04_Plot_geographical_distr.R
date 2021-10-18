##############################################
### Mortality of fish passing turbines
### Map plot of investigated turbine mortality studies
### Author: Johannes Radinger
### October 2021
##############################################


library(ggplot2)
library(scales)
library(ggsci) #https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
library(xlsx)
library(dplyr)
library(plyr)
library(see)
library(rgdal)
library(gtable)
library(gridExtra)
library(sp)
library(ggspatial)
library(patchwork)
#devtools::install_github("ropenscilabs/rnaturalearth")
devtools::install_github("burchill/zplyr")
library(zplyr)
library("rnaturalearth")
library("rnaturalearthdata")

################################################################
# Load data
turb_mort_df <- read.csv("./Data/turb_mort_df.csv",stringsAsFactors=T)
turb_df <- read.csv("./Data/turb_df.csv",stringsAsFactors=T)


################################################################
summary_ls <- list()
count=1
for(i in 1:length(unique(turb_df$location))){
  message(paste("This is location ",unique(turb_df$location)[i],sep=""))
  turb_df_i <- turb_df[turb_df$location==unique(turb_df$location)[i],]
  turb_df_i <- droplevels(turb_df_i)
  turb_n <- levels(as.factor(turb_df_i$turbine_type))
  for(j in 1:length(turb_n)){
    summary_ls[[count]] <- data.frame(location = unique(turb_df$location)[i],
                                  Lat = mean(turb_df_i$Lat,na.rm=T),
                                  Long = mean(turb_df_i$Long,na.rm=T),
                                  capacity = mean(turb_df_i$total_capacity_kW,na.rm=T),
                                  turbine_type = turb_n[j],
                                  capacity_class = turb_df_i$capacity_class[1])
    count <- count+1
  }
}

summary_df <- do.call(rbind,summary_ls)
summary_df$capacity_class_numeric <- as.numeric(as.character(mapvalues(summary_df$capacity_class,
                                             from=c("vSHP","SHP","LHP"), to = c(1,2,3))))
summary_df$turbine_type <- ordered(summary_df$turbine_type, levels = c("Kaplan","Francis","VLH","Screw","Water wheel","Cross-flow","Other"))

#Export data as shape to make maps in GIS
summary_df_GIS <- summary_df[!is.na(summary_df$Lat),]
summary_df_GIS$turbine_type <- factor(summary_df_GIS$turbine_type,ordered=FALSE)
summary(summary_df_GIS)
coordinates(summary_df_GIS)=~Long+Lat
proj4string(summary_df_GIS)<- CRS("+init=epsg:4326")

writeOGR(summary_df_GIS,
         "./Results_Figs_Tabs/Turbine_locations_SHP/",
         layer="Turbine_locations",
         driver="ESRI Shapefile")

### Load country maps
europe <- ne_countries(scale = "medium", returnclass = "sf",
                      continent = "europe")
world <- ne_countries(scale = "medium", returnclass = "sf")


### Plot maps for continents
Europe <- ggplot(data = world) +
  geom_sf(fill="grey70",color="grey90",size=0.3) +
  geom_point(data=summary_df, aes(x=Long,y=Lat,
                                  size=capacity_class_numeric,
                                  fill=capacity_class),
             shape=21,
             #fill=material_colors("deep purple"),
             color="white",
             alpha=0.7)+
  scale_fill_manual(values=c("#ffcd38","#cb9904","#997304"))+
  scale_radius(name="Hydropower scale", range = c(3, 9),breaks=c(1,2,3))+
  coord_sf(xlim = c(-12, 33), ylim = c(40, 65), expand = T,datum = NA)+
  theme_minimal()+
  annotation_scale(location = "br", width_hint = 0.5,
                   bar_cols = c("grey40", "grey90"),
                   line_width = 0,
                   height=unit(0.15,"cm")) +
  geom_abs_text(aes(xpos=0.05,ypos=0.95,label="b"),size=5.5,data = data.frame())+
  theme(axis.title=element_blank(),
        legend.position = "none",
        legend.justification = c(0, 0),
        legend.box.background=element_rect(fill="white",colour="grey50"),
        panel.border = element_rect(colour = "grey50", fill=NA),
        plot.background = element_rect(colour = "grey50", fill=NA),
        plot.margin = unit(c(0, 0, 0, 0), "null"))+
  guides(shape = guide_legend(ncol = 1, override.aes = list(size = 3, color = material_colors("deep purple"))),
         size = guide_legend(override.aes = list(fill=c("#cb9904","#ffcd38","#997304"))),
         fill = "none")

NorthAmerica <- ggplot(data = world) +
  geom_sf(fill="grey70",color="grey90",size=0.3) +
  geom_point(data=summary_df, aes(x=Long,y=Lat,
                                  size=capacity_class_numeric,
                                  fill=capacity_class),
             shape=21,
             #fill=material_colors("deep purple"),
             color="white",
             alpha=0.7)+
  #scale_fill_manual(values=c("#ffcd38","#997304","#cb9904"))+
  scale_fill_manual(values=c("#ffcd38","#cb9904","#997304"))+
  scale_radius(name="Hydropower scale", range = c(3, 9),breaks=c(1,2,3),
               labels=c(expression(""<1~'MW'),
                        expression(1-10~'MW'),
                        expression("">=10~'MW')))+
  #annotate(geom = 'text', label = 'a', x = -Inf, y = Inf, hjust = 0.1, vjust = 0.95)+
  geom_abs_text(aes(xpos=0.05,ypos=0.95,label="a"),size=5.5,data = data.frame())+
  coord_sf(xlim = c(-127, -50), ylim = c(17, 66.5), expand = T, datum = NA,clip="on")+
  theme_minimal()+
  annotation_scale(location = "br", width_hint = 0.5,
                   bar_cols = c("grey40", "grey90"),
                   line_width = 0,
                   height=unit(0.15,"cm")) +
  theme(axis.title=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.background = element_rect(fill=alpha('white',0.75),colour="grey50"),
        panel.border = element_rect(colour = "grey50", fill=NA),
        plot.background = element_rect(colour = "grey50", fill=NA),
        plot.margin = unit(c(0, 0, 0, 0), "null"))+
  guides(shape = guide_legend(ncol = 1, override.aes = list(size = 3, color = material_colors("deep purple"))),
         size = guide_legend(override.aes = list(fill=c("#cb9904","#ffcd38","#997304"))),
         fill = "none")


NorthAmerica + Europe



cairo_pdf("./Results_Figs_Tabs/Map_turbines.pdf",width=9,height=4.5)
NorthAmerica + Europe
dev.off()
