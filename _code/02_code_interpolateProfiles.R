#############################
##
## Date: March 28 2018
## Purpose: Interpolate soil info to deeper depths for APSIM
##
## Author: Rafa, modified by Gina
##
## Inputs: Output from downloadData script, in _dataOut
##           soilData_raw.rds
##
## Outputs: _figs: interpolated_profiles.jpeg
##          _dataOut: soilData_interpolated_profiles.rds
##          _dataOut: soilData_interpolated_profiles.csv
############################################################

## NOTE: Restart R before running. Some functions get confused if you don't. 

# Set working directory to wherever active document is kept ######################################################
#
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

rm(list=ls()) #Clears environment, not always necessary

library(tidyverse)


### Read raw data #########################################################
readRDS("../_dataOut/soilData_raw.rds")  %>%
  filter(!(compname %in% c("Water","Urban land"))) %>%
  filter(!is.na(sand))-> h

length(unique(h$abv))

soilLayer_breaks <- c(1,5,10,15,30,45,60,80,100,120,150,200,250,300,
                      350,400,450,500,550,600) # in cm
saveRDS(soilLayer_breaks,"../_dataOut/soilLayer_breaks.rds")

### Helper functions ######################################################
source("XX_code_helper.R")

### Add extra layer to the bottom if not depth enough ##################### 

h %>%
  rbind(h %>% 
          group_by(mukey) %>%
          mutate(maxdepth = max(bottom)) %>%
          filter(bottom == maxdepth) %>%
          filter(maxdepth < max(soilLayer_breaks)) %>%
          mutate(top = bottom + 1,
                 bottom = max(soilLayer_breaks),
                 thick = bottom - top)) %>%
  select(-maxdepth) %>%
  group_by() -> h

### Expand by each cm #####################################################

x <- data.frame(expand.grid(compname = unique(h$compname),
                            mukey = unique(h$mukey),
                            cokey = unique(h$cokey),
                            center = 1:max(soilLayer_breaks)))

h %>% 
  select(abv:hydrogroup,watertable)%>%
  unique() %>%
  full_join(x %>%
              left_join(h %>% 
                          select(c("abv","compname","musym","cokey","slope","hydrogroup","top",
                                   "bottom","thick","slope_code","CN2","hzB",
                                   "hzC","hzWT"))) %>%
              left_join(h) %>%
              filter(center >= top & center <= bottom) %>%
              arrange(abv,mukey,cokey,center) %>%
              group_by(abv,mukey,cokey,center) %>%
              summarise_at(c("CN2","hzB","hzC","hzWT","sand","clay","wetbd",
                             "drybd","om","ksat","ll","dul","ph"),
                           mean) %>%
              arrange(mukey,cokey,center) %>%
              mutate(sand = Int(sand,center),
                     clay = Int(clay,center),
                     wetbd = Int(wetbd,center),
                     drybd = Int(drybd,center),
                     om = Int(om,center),
                     ksat = Int(ksat,center),
                     ll = Int(ll,center),
                     dul = Int(dul,center),
                     ph = Int(ph,center))) %>%
  filter(!is.na(musym)) -> horizon

horizon %>%
  gather(variable,value,CN2:ph,watertable) %>%
  ggplot(aes(x=center,y=value,colour = as.factor(paste(abv,musym)))) +
  geom_point() + 
  #geom_line() + 
  scale_x_reverse() +
  coord_flip() + 
  facet_wrap(~variable, scales = "free_x")+ 
  theme(legend.position = "none") -> fig

ggsave("../_figs/interpolated_profiles.jpeg", plot = fig, width = 10, height = 8)

### Save expanded data #########################################################

saveRDS(horizon, "../_dataOut/soilData_interpolated_profiles.rds")
write.csv(horizon, "../_dataOut/soilData_interpolated_profiles.csv")
