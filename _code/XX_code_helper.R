
# Interpolate function ################################################

Int <- function (Y,X,method = "linear" ){
  
  if(all(is.na(X)) | all(is.infinite(X)) | all(is.na(Y)) | all(is.infinite(Y))){
    
    rep(NA, length(X))
  
    }
  
  else { 
    
    if(sum(!is.na(Y)) == 1) {
      
      rep(unique(Y[!is.na(Y)]), length(X))
      
    } else {
      
      approx(X, Y, 
             rule=2,
             xout=1:length(X)-1,
             method=method)$y
    }
    
    }
}



# Saxton and Rawls function ############################################

SaxtonRawls <-function(pSand, pClay, pOM){
  
  pSand <- pSand/100
  pClay <- pClay/100
  pOM <- pOM/100
  
  # calc LL15 (theta_1500)
  theta_1500t = -0.024*pSand + 0.487*pClay + 0.006*pOM +
    0.005*pSand*pOM - 0.013*pClay*pOM + 0.068*pSand*pClay + 0.031
  LL15 = theta_1500t + (0.14*theta_1500t - 0.02)
  LL15 = round(pmax(0.01, pmin(0.99, LL15)),3)
  
  # calc DUL (theta_33)
  theta_33t = -0.251*pSand + 0.195*pClay + 0.11*pOM +
    0.006*pSand*pOM - 0.027*pClay*pOM + 0.452*pSand*pClay + 0.299
  #DUL = theta_33t + (1.283*theta_33t^2 - 0.374*theta_33t - 0.015)
  DUL = theta_33t + (1.283*theta_33t^2 - 0.374*theta_33t - 0.015)
  DUL = round(pmax(0.01, pmin(0.99, DUL)),3)
  
  # calc SAT-33 KPa moisture
  theta_sat33t = 0.278*pSand + 0.034*pClay + 0.022*pOM -
    0.018*pSand*pOM - 0.027*pClay*pOM - 0.584*pSand*pClay + 0.078
  theta_sat33 = theta_sat33t + (0.636*theta_sat33t - 0.107)
  
  # calc SAT
  SAT = DUL + theta_sat33 - 0.097*pSand + 0.043
  SAT = round(pmax(0.01, pmin(0.99, SAT)),3)
  
  # calc BD
  BD = (1 - SAT)*2.65
  BD = round(pmax(1.0, pmin(2.1, BD)),3)
  
  # calc ksat (saturated water conductivity)
  lambda = (log(DUL)-log(LL15)) / (log(1500)-log(33))
  ksat = 1930*((SAT-DUL)^(3-lambda))
  SWCON = round(0.15 + pmin(ksat,75)/100, 3)
  
  res <- list(BD, LL15*100, DUL*100, SAT*100, SWCON*100,ksat)
  names(res) <- c("BD", "LL15", "DUL", "SAT", "SWCON","KSAT")
  
  return(res)
}


# xlmcomplie Function #########################################

xmlCompile <- function(x, crops = c("maize","soybean","wheat","bambatsi"), 
                       include_Swim = TRUE, 
                       include_tile = TRUE,
                       inifile = NULL) {
  
  require(XML)
  require(Hmisc)
  require(lubridate)
  
  #<folder version="35" creator="Apsim 7.8-r3899" name="Soils">
  folder <- newXMLNode("folder", attrs = list(version="37",
                                              creator="Apsim 7.9-r4065",
                                              name="IL Water Table Sites"))

  x <- split(x,x$compname)
  
  for(i in 1:length(x)){
    
    print(names(x)[i])
    data <- x[[i]]
    
    ## <Soil name="Default">
    Soil <- newXMLNode("Soil", 
                       attrs = list(name=unique(data$compname)),
                       parent = folder)
    
    ### <RecordNumber>0</RecordNumber>
    RecordNumber <- newXMLNode("RecordNumber",parent = Soil)
    xmlValue(RecordNumber) <-  1
    
    ### <SoilType>Nicollet</SoilType>
    SoilType <- newXMLNode("SoilType",parent = Soil)
    xmlValue(SoilType) <- unique(data$compname)
    
    ### <Region>Story</Region>
    county <- unique(data$compname) #capitalize(latlong2county(lat=data$coordinates[1], long=data$coordinates[2]))
    Region <- newXMLNode("Region",parent = Soil)
    xmlValue(Region) <- county
    
    ### <State>Iowa</State>
    State <- newXMLNode("State",parent = Soil)
    xmlValue(State) <- "IL" #county[1]
    
    ### <Country>US</Country>
    Country <- newXMLNode("Country",parent = Soil)
    xmlValue(Country) <- "USA"
    
    ### <ApsoilNumber>1</ApsoilNumber>
    ApsoilNumber <- newXMLNode("ApsoilNumber",parent = Soil)
    xmlValue(ApsoilNumber) <- 1
    
    ### <Latitude>0</Latitude>
    Latitude <- newXMLNode("Latitude",parent = Soil)
    xmlValue(Latitude) <- 0 #round(data$coordinates[1],2)
    
    ### <Longitude>-0</Longitude>
    Longitude <- newXMLNode("Longitude",parent = Soil)
    xmlValue(Longitude) <- 0 #round(data$coordinates[2],2)
    
    ### <YearOfSampling>0</YearOfSampling>
    YearOfSampling <- newXMLNode("YearOfSampling",parent = Soil)
    xmlValue(YearOfSampling) <- as.character(year(Sys.Date()))
    
    ### <DataSource>ssurgo2apsim</DataSource>
    DataSource <- newXMLNode("DataSource",parent = Soil)
    xmlValue(DataSource) <- "ssurgo2apsim"
    
    ### <Comments></Comments>
    Comments <- newXMLNode("Comments",parent = Soil)
    xmlValue(Comments) <- "None"
    
    
    if(!is.null(inifile)){
    ### <ini>
      ini <- newXMLNode("ini",parent = Soil)
      filename <- newXMLNode("filename",parent = ini,
                             attrs = list(input="yes"))
      xmlValue(filename) <- inifile
      
    }
    
    
    ### <Water>
    Water <- newXMLNode("Water",parent = Soil)
    
    #### <Thickness>
    Thickness <- newXMLNode("Thickness",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = Thickness)
      xmlValue(double) <- data$thick[j]
    }
    
    #### <BD>
    BD <- newXMLNode("BD",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = BD)
      xmlValue(double) <- round(data$bd[j],2)
    }
    
    #### <AirDry>
    AirDry <- newXMLNode("AirDry",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = AirDry)
      xmlValue(double) <- round(data$AirDry[j],2)
    }
    
    #### <LL15>
    LL15 <- newXMLNode("LL15",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = LL15)
      xmlValue(double) <- round(data$ll[j],2)
    }
    
    #### <DUL>
    DUL <- newXMLNode("DUL",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = DUL)
      xmlValue(double) <- round(data$dul[j],2)
    }
    
    #### <SAT>
    SAT <- newXMLNode("SAT",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = SAT)
      xmlValue(double) <- round(data$sat[j],2)
    }
    
    #### <KS>
    KS <- newXMLNode("KS",parent = Water)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = KS)
      xmlValue(double) <- round(data$ksat[j],1)
    }
    
    #### <SoilCrop>
    for(k in 1:length(crops)){
      SoilCrop <- newXMLNode("SoilCrop", attrs = list(name = crops[k]),parent = Water)
      
      ##### <Thickness>
      Thickness <- newXMLNode("Thickness",parent = SoilCrop)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = Thickness)
        xmlValue(double) <- data$thick[j]
      }
      
      ##### <LL>
      LL <- newXMLNode("LL",parent = SoilCrop)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = LL)
        xmlValue(double) <- round(data$ll[j],2)
      }
      
      ##### <KL>
      KL <- newXMLNode("KL",parent = SoilCrop)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = KL)
        xmlValue(double) <- round(data$KL_maize[j],3)
      }
      
      ##### <XF>
      XF <- newXMLNode("XF",parent = SoilCrop)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = XF)
        xmlValue(double) <- round(data$XF_maize[j],3)
      }
    }
    
    if(include_Swim) {
      
      ### <Swim>
      Swim <- newXMLNode("Swim",parent = Soil)
      
      #### <Salb>
      Salb <- newXMLNode("Salb",parent = Swim)
      xmlValue(Salb) <- unique(data$Salb)
      
      #### <CN2Bare>
      CN2Bare <- newXMLNode("CN2Bare",parent = Swim)
      xmlValue(CN2Bare) <- (round(unique(data$CN2),0))
      
      #### <CNRed>
      CNRed <- newXMLNode("CNRed",parent = Swim)
      xmlValue(CNRed) <- (unique(data$CNRed))
      
      #### <CNCov>
      CNCov <- newXMLNode("CNCov",parent = Swim)
      xmlValue(CNCov) <- (unique(data$CNCov))
      
      #### <KDul>
      KDul <- newXMLNode("KDul",parent = Swim)
      xmlValue(KDul) <- (unique(data$KDul))
      
      #### <PSIDul>
      PSIDul <- newXMLNode("PSIDul",parent = Swim)
      xmlValue(PSIDul) <- (unique(data$PSIDul))
      
      #### <VC>
      VC <- newXMLNode("VC",parent = Swim)
      xmlValue(VC) <- tolower(as.character(as.logical(unique(data$VC))))
      
      #### <DTmin>
      DTmin <- newXMLNode("DTmin",parent = Swim)
      xmlValue(DTmin) <- (unique(data$DTmin))
      
      #### <DTmax>
      DTmax <- newXMLNode("DTmax",parent = Swim)
      xmlValue(DTmax) <- (unique(data$DTmax))
      
      #### <MaxWaterIncrement>
      MaxWaterIncrement <- newXMLNode("MaxWaterIncrement",parent = Swim)
      xmlValue(MaxWaterIncrement) <- (unique(data$MaxWaterIncrement))
      
      #### <SpaceWeightingFactor>
      SpaceWeightingFactor <- newXMLNode("SpaceWeightingFactor",parent = Swim)
      xmlValue(SpaceWeightingFactor) <- (unique(data$SpaceWeightingFactor))
      
      #### <SoluteSpaceWeightingFactor>
      SoluteSpaceWeightingFactor <- newXMLNode("SoluteSpaceWeightingFactor",parent = Swim)
      xmlValue(SoluteSpaceWeightingFactor) <- (unique(data$SoluteSpaceWeightingFactor))
      
      #### <Diagnostics>
      Diagnostics <- newXMLNode("Diagnostics",parent = Swim)
      xmlValue(Diagnostics) <- tolower(as.character(as.logical(unique(data$Diagnostics))))
      
      #### <SwimSoluteParameters>
      SwimSoluteParameters <- newXMLNode("SwimSoluteParameters",parent = Swim)
      
      #####<Dis>
      Dis <- newXMLNode("Dis",parent = SwimSoluteParameters)
      xmlValue(Dis) <- (unique(data$Dis))
      
      #####<Disp>1</Disp>
      Disp <- newXMLNode("Disp",parent = SwimSoluteParameters)
      xmlValue(Disp) <- (unique(data$Disp))
      
      #####<A>1</A>
      A <- newXMLNode("A",parent = SwimSoluteParameters)
      xmlValue(A) <- (unique(data$A))
      
      #####<DTHC>0</DTHC>
      DTHC <- newXMLNode("DTHC",parent = SwimSoluteParameters)
      xmlValue(DTHC) <- (unique(data$DTHC))
      
      #####<DTHP>1</DTHP>
      DTHP <- newXMLNode("DTHP",parent = SwimSoluteParameters)
      xmlValue(DTHP) <- (unique(data$DTHP))
      
      #####<WaterTableCl>
      WaterTableCl <- newXMLNode("WaterTableCl",parent = SwimSoluteParameters)
      xmlValue(WaterTableCl) <- (unique(data$WaterTableCl))
      
      #####<WaterTableNO3>
      WaterTableNO3 <- newXMLNode("WaterTableNO3",parent = SwimSoluteParameters)
      xmlValue(WaterTableNO3) <- (unique(data$WaterTableNO3))
      
      #####<WaterTableNH4>
      WaterTableNH4 <- newXMLNode("WaterTableNH4",parent = SwimSoluteParameters)
      xmlValue(WaterTableNH4) <- (unique(data$WaterTableNH4))
      
      #####<WaterTableUrea>
      WaterTableUrea <- newXMLNode("WaterTableUrea",parent = SwimSoluteParameters)
      xmlValue(WaterTableUrea) <- (unique(data$WaterTableUrea))
      
      #####<WaterTableTracer>
      WaterTableTracer <- newXMLNode("WaterTableTracer",parent = SwimSoluteParameters)
      xmlValue(WaterTableTracer) <- (unique(data$WaterTableTracer))
      
      #####<WaterTableMineralisationInhibitor>
      WaterTableMineralisationInhibitor <- newXMLNode("WaterTableMineralisationInhibitor",parent = SwimSoluteParameters)
      xmlValue(WaterTableMineralisationInhibitor) <- (unique(data$WaterTableMineralisationInhibitor))
      
      #####<WaterTableUreaseInhibitor>
      WaterTableUreaseInhibitor <- newXMLNode("WaterTableUreaseInhibitor",parent = SwimSoluteParameters)
      xmlValue(WaterTableUreaseInhibitor) <- (unique(data$WaterTableUreaseInhibitor))
      
      #####<WaterTableNitrificationInhibitor>
      WaterTableNitrificationInhibitor <- newXMLNode("WaterTableNitrificationInhibitor",parent = SwimSoluteParameters)
      xmlValue(WaterTableNitrificationInhibitor) <- (unique(data$WaterTableNitrificationInhibitor))
      
      #####<WaterTableDenitrificationInhibitor>
      WaterTableDenitrificationInhibitor <- newXMLNode("WaterTableDenitrificationInhibitor",parent = SwimSoluteParameters)
      xmlValue(WaterTableDenitrificationInhibitor) <- (unique(data$WaterTableDenitrificationInhibitor))
      
      #####<Thickness>
      
      Thickness <- newXMLNode("Thickness",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = Thickness)
        xmlValue(double) <- data$thick[j]
      }
      
      #### <NO3Exco>
      NO3Exco <- newXMLNode("NO3Exco",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = NO3Exco)
        xmlValue(double) <- data$NO3Exco[j]
      }
      
      #### <NO3FIP>
      NO3FIP <- newXMLNode("NO3FIP",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = NO3FIP)
        xmlValue(double) <- data$NO3FIP[j]
      }
      
      #### <NH4Exco>
      NH4Exco <- newXMLNode("NH4Exco",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = NH4Exco)
        xmlValue(double) <- data$NH4Exco[j]
      }
      
      #### <NH4FIP>
      NH4FIP <- newXMLNode("NH4FIP",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = NH4FIP)
        xmlValue(double) <- data$NH4FIP[j]
      }
      
      #### <UreaExco>
      UreaExco <- newXMLNode("UreaExco",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = UreaExco)
        xmlValue(double) <- data$UreaExco[j]
      }
      
      #### <UreaFIP>
      UreaFIP <- newXMLNode("UreaFIP",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = UreaFIP)
        xmlValue(double) <- data$UreaFIP[j]
      }
      
      #### <ClExco>
      ClExco <- newXMLNode("ClExco",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = ClExco)
        xmlValue(double) <- data$ClExco[j]
      }
      
      #### <ClFIP>
      ClFIP <- newXMLNode("ClFIP",parent = SwimSoluteParameters)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = ClFIP)
        xmlValue(double) <- data$ClFIP[j]
      }
      
      if (include_tile){
        
        #### <SwimSubsurfaceDrain>
        SwimSubsurfaceDrain <- newXMLNode("SwimSubsurfaceDrain",parent = Swim)
        
        #####<DrainDepth>
        DrainDepth <- newXMLNode("DrainDepth",parent = SwimSubsurfaceDrain)
        xmlValue(DrainDepth) <- (unique(data$DrainDepth))
        
        #####<DrainSpacing>
        DrainSpacing <- newXMLNode("DrainSpacing",parent = SwimSubsurfaceDrain)
        xmlValue(DrainSpacing) <- (unique(data$DrainSpacing))
        
        #####<DrainRadius>
        DrainRadius <- newXMLNode("DrainRadius",parent = SwimSubsurfaceDrain)
        xmlValue(DrainRadius) <- (unique(data$DrainRadius))
        
        #####<Klat>
        Klat <- newXMLNode("Klat",parent = SwimSubsurfaceDrain)
        xmlValue(Klat) <- (unique(data$Klat))
        
        #####<ImpermDepth>
        ImpermDepth <- newXMLNode("ImpermDepth",parent = SwimSubsurfaceDrain)
        xmlValue(ImpermDepth) <- (unique(data$ImpermDepth))
        
      }
      
      if(!is.na(unique(data$WaterTableDepth))) {
        
        #### <SwimWaterTable>
        SwimWaterTable <- newXMLNode("SwimWaterTable",parent = Swim)
        
        #####<WaterTableDepth>
        WaterTableDepth <- newXMLNode("WaterTableDepth",parent = SwimWaterTable)
        
        if(include_tile) {
          xmlValue(WaterTableDepth) <- (unique(data$DrainDepth)) - 50
          
        } else {
          
          xmlValue(WaterTableDepth) <- round(unique(data$WaterTableDepth))
        
        }
        
        
      }
      
    } else {
      
      ### <SoilWater>
      SoilWater <- newXMLNode("SoilWater",parent = Soil)
      
      #### <SummerCona>
      SummerCona <- newXMLNode("SummerCona",parent = SoilWater)
      xmlValue(SummerCona) <-(round(data$cona,2)[1])
      
      #### <SummerU>
      SummerU <- newXMLNode("SummerU",parent = SoilWater)
      xmlValue(SummerU) <- (round(data$U,2)[1])
      
      #### <SummerDate>
      SummerDate <- newXMLNode("SummerDate",parent = SoilWater)
      xmlValue(SummerDate) <- "1-jun"
      
      #### <WinterCona>
      WinterCona <- newXMLNode("WinterCona",parent = SoilWater)
      xmlValue(WinterCona) <-as.character(round(data$cona,2)[1])
      
      #### <WinterU>
      WinterU <- newXMLNode("WinterU",parent = SoilWater)
      xmlValue(WinterU) <- as.character(round(data$U,2)[1])
      
      #### <WinterDate>
      WinterDate <- newXMLNode("WinterDate",parent = SoilWater)
      xmlValue(WinterDate) <- "1-nov"
      
      #### <DiffusConst>
      DiffusConst <- newXMLNode("DiffusConst",parent = SoilWater)
      xmlValue(DiffusConst) <-  unique(data$DiffusConst)
      
      #### <DiffusSlope>
      DiffusSlope <- newXMLNode("DiffusSlope",parent = SoilWater)
      xmlValue(DiffusSlope) <-  unique(data$DiffusSlope)
      
      #### <Salb>
      Salb <- newXMLNode("Salb",parent = SoilWater)
      xmlValue(Salb) <- unique(data$Salb)
      
      #### <CN2Bare>
      CN2Bare <- newXMLNode("CN2Bare",parent = SoilWater)
      xmlValue(CN2Bare) <- (round(unique(data$CN2),0))
      
      #### <CNRed>
      CNRed <- newXMLNode("CNRed",parent = SoilWater)
      xmlValue(CNRed) <- (unique(data$CNRed))
      
      #### <CNCov>
      CNCov <- newXMLNode("CNCov",parent = SoilWater)
      xmlValue(CNCov) <- (unique(data$CNCov))
      
      #### <Slope>NaN</Slope>
      #### <DischargeWidth>NaN</DischargeWidth>
      #### <CatchmentArea>NaN</CatchmentArea>
      #### <MaxPond>NaN</MaxPond>
      
      #### <Thickness>
      Thickness <- newXMLNode("Thickness",parent = SoilWater)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = Thickness)
        xmlValue(double) <- data$thick[j]
      }
      
      #### <SWCON>
      SWCON <- newXMLNode("SWCON",parent = SoilWater)
      for(j in 1:length(data$thick)){ 
        double <- newXMLNode("double",parent = SWCON)
        xmlValue(double) <- round(data$SWCON[j],3)
      }
      
    }
    
    
    ### SoilOrganicMatter
    SoilOrganicMatter <- newXMLNode("SoilOrganicMatter",parent = Soil)
    
    #### <RootCN>
    RootCN <- newXMLNode("RootCN",parent = SoilOrganicMatter)
    xmlValue(RootCN) <- unique(round(data$RootCN))
    
    #### <RootWt>
    RootWt <- newXMLNode("RootWt",parent = SoilOrganicMatter)
    xmlValue(RootWt) <- unique(round(data$RootWt))
    
    #### <SoilCN>13</SoilCN>
    SoilCN <- newXMLNode("SoilCN",parent = SoilOrganicMatter)
    xmlValue(SoilCN) <- unique(round(data$SoilCN))
    
    #### <EnrACoeff>
    EnrACoeff <- newXMLNode("EnrACoeff",parent = SoilOrganicMatter)
    xmlValue(EnrACoeff) <- unique(round(data$EnrAcoeff,2))
    
    #### <EnrBCoeff>
    EnrBCoeff <- newXMLNode("EnrBCoeff",parent = SoilOrganicMatter)
    xmlValue(EnrBCoeff) <- unique(round(data$EnrBcoeff,2))
    
    #### <Thickness>
    Thickness <- newXMLNode("Thickness",parent = SoilOrganicMatter)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = Thickness)
      xmlValue(double) <- data$thick[j]
    }
    
    #### <OC>
    OC <- newXMLNode("OC",parent = SoilOrganicMatter)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = OC)
      xmlValue(double) <- round(data$OC[j],2)
    }
    
    #### <FBiom>
    FBiom <- newXMLNode("FBiom",parent = SoilOrganicMatter)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = FBiom)
      xmlValue(double) <- round(data$FBiom[j],4)
    }
    
    #### <FInert>
    FInert <- newXMLNode("FInert",parent = SoilOrganicMatter)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = FInert)
      xmlValue(double) <- round(data$FInert[j],4)
    }
    
    #### <OCUnits>
    OCUnits <- newXMLNode("OCUnits",parent = SoilOrganicMatter)
    xmlValue(OCUnits) <- "Total"
    
    ### <Analysis>
    Analysis <- newXMLNode("Analysis",parent = Soil)
    
    #### <Thickness>
    Thickness <- newXMLNode("Thickness",parent = Analysis)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = Thickness)
      xmlValue(double) <- data$thick[j]
    }
    
    #### <PH>
    PH <- newXMLNode("PH",parent = Analysis)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = PH)
      xmlValue(double) <- round(data$ph[j],2)
    }
    
    #### <ParticleSizeSand>
    ParticleSizeSand <- newXMLNode("ParticleSizeSand",parent = Analysis)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = ParticleSizeSand)
      xmlValue(double) <- round(data$sand[j],1)
    }
    
    #### <ParticleSizeClay>
    ParticleSizeClay <- newXMLNode("ParticleSizeClay",parent = Analysis)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = ParticleSizeClay)
      xmlValue(double) <- round(data$clay[j],1)
    }
    
    #### <ParticleSizeSilt>
    ParticleSizeSilt <- newXMLNode("ParticleSizeSilt",parent = Analysis)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = ParticleSizeSilt)
      xmlValue(double) <- round(100 - data$sand[j] - data$clay[j],1)
    }
    
    ### Sample
    Sample <- newXMLNode("Sample",parent = Soil, attrs = list(name="Intial conditions"))
    
    #### <Thickness>
    Thickness <- newXMLNode("Thickness",parent = Sample)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = Thickness)
      xmlValue(double) <- data$thick[j]
    }
    
    #### <NO3>
    NO3 <- newXMLNode("NO3",parent = Sample)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = NO3)
      xmlValue(double) <- round(data$initNO3[j],2) # same as OC but in ppm
    }
    
    #### <NH4>
    NH4 <- newXMLNode("NH4",parent = Sample)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = NH4)
      xmlValue(double) <- round(data$initNH4[j],2) # same as 1/2 oc but ppm
    }
    
    #### <SW>
    SW <- newXMLNode("SW",parent = Sample)
    for(j in 1:length(data$thick)){ 
      double <- newXMLNode("double",parent = SW)
      xmlValue(double) <- round(data$initSW[j],2)
    }
    
    #### <NO3Units>
    NO3Units <- newXMLNode("NO3Units",parent = Sample)
    xmlValue(NO3Units) <- "ppm"
    
    #### <NH4Units>
    NH4Units <- newXMLNode("NH4Units",parent = Sample)
    xmlValue(NH4Units) <- "ppm"
    
    #### <SWUnits>
    SWUnits <- newXMLNode("SWUnits",parent = Sample)
    xmlValue(SWUnits) <- "Volumetric"
    
  }
  
  return(folder)
  
} 

# SSURGO2APSIM Function #################################################

SSURGO2APSIM <- function(data,
                         filename,
                         crops = c("maize","soybean","wheat","bambatsi"),
                         include_Swim = TRUE, include_tile = TRUE
                         #inifile = "C:\\Users\\rmartine\\Dropbox\\PhD_Martinez-Feria\\_modeling\\paper3\\scenarios\\xlmfiles\\Soil.xml"
                         ) {
  
  out <- xmlCompile(data, crops, include_Swim, include_tile)
  
  writeLines(saveXML(out),filename)
}


