require(tidyverse)
require(XML)

### Helper functions ######################################################
source("../Rscripts/helper.R")

# Build Soils ##############################################################################################################

#name of APSIM xml

xmlname = "KeithSites.xml"

soilLayer_breaks <- readRDS("../data/soilLayer_breaks.rds")

drainage_parms = list(boundary_layer = TRUE,
                      DrainDepth  = 1150, #Depth of subsurface drain (mm)
                      DrainSpacing = 13000, #Distance between subsurface drains (mm)
                      DrainRadius  = 100, #Radius of each subsurface drain (mm)
                      Klat  = 1000 #Lateral saturated soil water conductivity (mm/d)
)

readRDS("../data/soilData_interpolated_profiles.rds") %>%
   mutate(
     ### HERE IS WHERE YOU CAN MAKE CHANGES ###################################
     bd = ifelse( wetbd < 0.9, 0.9, ifelse(wetbd > 1.8, 1.8,wetbd)),
     ksat =pmin(ksat*100/1.157,500),
     #SaxtonRawls(pSand=sand ,pClay=clay, pOM=om)$KSAT*24), # mm/day
     sat =SaxtonRawls(pSand=sand ,pClay=clay, pOM=om)$SAT/100,
     PO =1-bd/2.65,
     Salb = 0.15, # Bare soil albedo
     MWCON = 1, #(0-1)
     dul =dul/100,
     ll =ll/100,
     SWCON =(PO-dul)/PO,
     AirDry =ifelse(center<=15,0.9,ifelse(center<=30,0.95,1))*ll,
     # Sat and ksat calculations
     #sat = ifelse(center > hzC, (sat - dul)*exp((center- hzC)*-0.01) + dul,sat),
     #ksat = ifelse(center > hzC, (ksat)*exp((center- hzC)*-0.03),ksat),
     #ksat = ifelse(ksat < 0.01, 0.01, ksat),
     #ksat = ifelse(ksat > 300, 300, ksat),
     #
     U =ifelse(clay<=20,5+0.175*clay,
               ifelse(clay<=40,7.5+0.05*clay,
                      ifelse(clay<=50,11.5-0.05*clay,
                             ifelse(clay<=70,12.75-0.075*clay,
                                    ifelse(clay<=80,11-0.05*clay,0))))),# mm
     cona =ifelse(clay<=30,0.025*clay+3.25,
                  ifelse(clay<=50,4,
                         ifelse(clay<=70,-0.025*clay+5.25,
                                ifelse(clay<=80,3.5,0)))), # mm/d^5
     DiffusConst =40,
     DiffusSlope =16,
     CN2 =ifelse(is.na(CN2),80,CN2) - 5,
     CNRed = 20, #residue runoff,
     CNCov	 = 0.8,
     EnrAcoeff = 7.4,
     EnrBcoeff = 0.2,
     XF_maize = ifelse(center<200,1,0),
     KL_maize = ifelse(center<=20,0.08,0.09*exp(-0.007*center)),
     e	= 0.5,  #ifelse(F4=$BC$3,0.07,IF(F4=$BC$4,0.03,0.05))
     # Soil chemical properties
     ph =0.52+1.06*ph, #pH 1:5
     OC =om/1.72, # %
     # OC adjustments
     OCdiff = c(NA, round(diff(OC),1)),
     OCthr = ifelse(is.finite(hzC),hzC, 120),
     OC = ifelse(center >= OCthr & OCdiff == 0,
                 OC*exp((center- OCthr)*-0.035),
                 OC),
     #
     FInert =ifelse(center<=1,0.4,
                    ifelse(center<=10,0.4,
                           ifelse(center<60,0.008*center+0.32,
                                  ifelse(center<=120,0.8,
                                         ifelse(center<180,0.0032*center+0.42,
                                                ifelse(center<=300,0.99,0)))))), #(0-1)
     FInert = ifelse(center > hzC, 0.99, FInert),
     FBiom =ifelse(center<=10,0.04,
                   ifelse(center<=20,0.055-0.0015*center,
                          ifelse(center<=30,0.03-0.0005*center,
                                 ifelse(center<60,0.0216-0.0002*center,
                                        ifelse(center<=300,0.01,0)))))*2, #(0-1)
     FBiom = ifelse(center > hzC, 0, FBiom),
     RootCN =45,
     SoilCN =13,
     RootWt =1000,
     # SWIM Specific
     KDul =0.1,
     PSIDul =-100,  #changing to -100 from -300, PME
     VC=TRUE,
     DTmin =0,
     DTmax =1440,
     MaxWaterIncrement =10,
     SpaceWeightingFactor =0,
     SoluteSpaceWeightingFactor =0,
     Diagnostics =FALSE,
     Dis  =15,
     Disp  =1,
     A  =1,
     DTHC  =0,
     DTHP  =1,
     WaterTableCl  =0,
     WaterTableNO3  =0,
     WaterTableNH4  =0,
     WaterTableUrea  =0,
     WaterTableTracer  =0,
     WaterTableMineralisationInhibitor  =0,
     WaterTableUreaseInhibitor  =0,
     WaterTableNitrificationInhibitor  =0,
     WaterTableDenitrificationInhibitor  =0,
     NO3Exco =0,
     NO3FIP = 1,
     NH4Exco =100,
     NH4FIP  =1,
     UreaExco  =0,
     UreaFIP =1,
     ClExco  =0,
     ClFIP =1,
     DrainDepth  =drainage_parms$DrainDepth,
     DrainSpacing  =drainage_parms$DrainSpacing,
     DrainRadius  =drainage_parms$DrainRadius,
     Klat  =drainage_parms$Klat,
     ImpermDepth  = 2500, 
     WaterTableDepth  =  ifelse(!is.finite(watertable),(hzWT)*10,(watertable)*10),
     WaterTableDepth  =  ifelse(!is.finite(WaterTableDepth),2500,WaterTableDepth),
     # Inital conditons
     initNO3 = OC,
     initNH4 = initNO3*0.5,
     initSW = dul, #ifelse(center < DrainDepth/10,dul,sat),
     # Layer codes
     layer =.bincode(center, breaks=c(-1,soilLayer_breaks[soilLayer_breaks>0]))
     ############################################################################
    ) %>%
  group_by(abv,musym,compname,mukey,area,slope,hydrogroup,layer) %>%
  mutate(thick = 10 + (max(center)-min(center))*10) %>%
  summarise_all("mean") %>%
  group_by() %>%
  filter(!is.na(sand)) %>%
  mutate(compname = paste(abv,compname,musym,sep = "_")) -> horizon

# Make soil files #######################################################

SSURGO2APSIM(horizon %>%
               as.data.frame(),
             crops = c("maize","soybean","wheat","bambatsi"),
             include_Swim = TRUE, include_tile = FALSE,
             filename = xmlname)

# View data #############################################################

horizon %>%
  ggplot(aes(x=center)) +
  geom_vline(aes(xintercept = WaterTableDepth/10), linetype = "dotted") + 
  geom_text(aes(y = 0.50, x= ifelse(is.finite(WaterTableDepth),WaterTableDepth/10 + 10,NA), label =  "WT")) +
  geom_ribbon(aes(ymin = ll, ymax=dul),fill = "blue", alpha = 0.5) + 
  geom_line(aes(y = sat), linetype = 2) + 
  #geom_text(aes(x=-10, y = mean(sat)),label = "Sat") + 
  annotate("text",x=240, y = 0.45,label = "Sat") + 
  annotate("text",x=240, y = 0.30,label = "DUL") +
  annotate("text",x=240, y = 0.18,label = "LL") +
  scale_x_reverse() +
  coord_flip() + 
  facet_wrap(~compname, ncol = 6, dir = "h")+ 
  theme(legend.position = "none") + 
  labs(y = "% moisture",
       x = "Depth (cm)") +
  ggthemes::theme_base() -> fig

ggsave("../data/calculated_profiles_hydr.jpeg", plot = fig, width = 16, height = 32)

horizon %>%
  mutate(FBiom = OC*FBiom,
         FInert = OC*FInert,
         FHum = OC - FBiom - FInert) %>%
  gather(pool, value, FBiom, FInert, FHum) %>%
  ggplot(aes(x=center)) +
  geom_area(aes(y = value, fill= pool), colour= "black") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_vline(aes(xintercept = hzB), linetype = "dotted") + 
  geom_vline(aes(xintercept = hzC), linetype = "dotted") + 
  geom_text(aes(y = 4.5, x= ifelse(is.finite(hzC),10,NA), label =  "A")) +
  geom_text(aes(y = 4.5, x= ifelse(is.finite(hzB),hzB,NA)+ 10, label =  "B")) +
  geom_text(aes(y = 4.5, x= ifelse(is.finite(hzC),hzC,NA) + 10, label =  "C")) +
  scale_x_reverse() +
  coord_flip() + 
  facet_wrap(~compname, ncol = 6, dir = "h", scales = "free")+ 
  theme(legend.position = "none") + 
  labs(y = "% OC",
       x = "Depth (cm)",
       fill = "Pool:") +
  ggthemes::theme_base() + theme(legend.position = "top") + 
  scale_fill_manual(values = c("red3","orange1","wheat4")) -> fig

ggsave("../data/calculated_profiles_SOM.jpeg", plot = fig, width = 16, height = 32)
