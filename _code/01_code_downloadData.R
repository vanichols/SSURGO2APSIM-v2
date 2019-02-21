#############################
##
## Date: March 28 2018
## Purpose: Use lat/long from field site to get SSURGO data in APSIM format
## Author: Rafa, modified by Gina on Nov 13 2018
## Inputs: user defined txt file from _dataIn in tools_SSURGO folder
##         hydrogroup.rds from _dataIn
## Outputs: soilData_raw.rds and soilData_raw.csv in _dataOut
##          sites$abv.jpeg in _figs
##
##############################

rm(list=ls()) #Clears environment, not always necessary

# Load packages
#library(sp)
library(tidyverse)
library(ggmap)
library(lubridate)
library(FedData)
library(raster)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(maptools)

# Set working directory to wherever active document is kept ######################################################
#
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

# Read data for locations ######################################################
#getwd()
setwd("../_dataOut")

# in illinois 1 degree lat is 110,000 m and long is 90,000 m  
# so a 0.001 represents an area of 1 ha.

w <- 0.001

read.csv("../_dataIn/FACTS_Stations2.txt") %>%
  mutate(north = lat + w,
         south = lat - w,
         east = long - w,
         west = long + w,
         abv = as.character(abv)) %>%
  dplyr::select(name, abv, north:west) -> sites

#data.frame(name = "Ames", abv = "AMES", 
#           north = 42.021697,
#           south = 42.020702,
#           east = -93.775228,
#           west =-93.775996) -> sites

#sites <- sites[1,]

out <- data.frame()

# Download ssurgo data with fed data package and extract needed info ###########

for(i in 1:length(sites$abv)){
  
  a <- polygon_from_extent(extent(sites$east[i], sites$west[i],
                                  sites$south[i],sites$north[i]),
                           proj4string="+proj=longlat")
  
  x <- get_ssurgo(template=a, label=sites$abv[i])
  
  print(as.character(sites$abv[i]))
  
  ### Calculate % of area #################################################
  area <- x$spatial@data
  area$id <- rownames(x$spatial@data)
  area$area <- NA
  
  for (j in 1:length(area$id)) area$area[j] <- x$spatial@polygons[[j]]@area
  
  area %>% 
    group_by(MUSYM,MUKEY) %>% 
    summarise(area=sum(area)) %>%
    group_by() %>%
    mutate(MUKEY = as.numeric(as.character(MUKEY)),
           area = area/sum(area)) %>%
    as.data.frame() %>%
    `names<-`(c("musym","mukey","area")) -> area
  
  ### Create map for site ##################################################
  
  get_map(location = c(lon = sites$east[i] + w,
                       lat = sites$south[i] + w ),
          zoom = 18, maptype = 'satellite') -> satellite
  
  overlay <-  fortify(x$spatial, region = "MUSYM")
  
  ggmap(satellite) +
    geom_point(x = sites$east[i] + w, y = sites$south[i] + w,
               shape = 4, size = 10)  +
    geom_polygon(data = overlay,
                 aes(x=long, y=lat, group = group, fill =  id),
                 alpha = 0.5, color = "black")  +
    scale_fill_hue(l = 40) +
    coord_equal() +
    labs(title = paste0(sites$abv[i]," (",sites$name[i],")"),
         fill= "MUSYM") +
    ggthemes::theme_foundation() +
    theme(plot.background = element_blank()) -> fig
  
  ggsave(filename =  paste0("../_figs/",sites$abv[i],".jpeg"),
         plot = fig, width = 8, height = 6, dpi = 300)
  
  
  ### Extract tabular data ##################################################
  component <- x$tabular$component
  chorizon <- x$tabular$chorizon
  mapunit <-x$tabular$mapunit
  muaggatt <- x$tabular$muaggatt
  
  ### Data trasformations ###################################################
  
  mapunit %>% 
    mutate(abv = sites$abv[i]) %>%
    dplyr::select(abv,musym,mukey) %>%
    left_join(component %>% 
                dplyr::select(compname,taxclname,drainagecl,
                              mukey,cokey,slope.r,hydgrp,majcompflag)) %>%
    filter(majcompflag == "Yes") %>%
    mutate(musym = as.character(musym))->  majcomp
  
  majcomp %>%
    left_join(chorizon) %>%
    left_join(area) %>%
    arrange(musym,compname,mukey,hzdept.r) %>%
    group_by(abv,musym,compname,mukey,cokey,area) %>% 
    dplyr::select(c("slope.r", "hydgrp","hzname",
                    "hzdept.r","hzdepb.r","hzthk.r",
                    "sandtotal.r","claytotal.r",
                    "dbthirdbar.r","dbovendry.r", 
                    "om.r","ksat.r",
                    "wfifteenbar.r","wthirdbar.r","ph1to1h2o.r")) %>%
    `names<-`(c("abv","musym","compname","mukey","cokey","area",
                "slope","hydrogroup","horizonName","top","bottom",
                "thick","sand","clay","wetbd",
                "drybd", "om","ksat","ll","dul","ph")) %>%
    mutate(slope_code = .bincode(slope, breaks=c(0,2,5,10,100))) %>% # Bin slope groups
    left_join(readRDS("../_dataIn/hydrogroup.rds")) %>%
    mutate(hzC = min(ifelse(grepl("C",horizonName),top,NA),na.rm = T),
           hzB = min(ifelse(grepl("B",horizonName),top,NA),na.rm = T),
           hzWT = min(ifelse(grepl("g",horizonName),top,NA),na.rm = T),
           thick = ifelse(is.na(thick),bottom - top, thick),
           center = trunc(top + thick/2)) %>%
    left_join(muaggatt %>%
                mutate(watertable =  as.numeric(wtdepaprjunmin)) %>%
                dplyr::select(mukey,watertable)) %>%
    bind_rows(out) -> out
}

# Write data files #############################################################
saveRDS(out,"../_dataOut/soilData_raw.rds")
write.csv(out,"../_dataOut/soilData_raw.csv", row.names = F)
