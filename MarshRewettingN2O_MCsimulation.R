
# Supporting files for the code below to work can be downloaded from:
# https://github.com/HollieEmery/RewettingPaper 

library(tidyverse)
library(lubridate)

# Step 1: Data Import ####

# Tide data ####

# Downloaded on 10 June 2018 from University of Hawaii Sea Level Center: 
# https://uhslc.soest.hawaii.edu/datainfo/

## Note: unzip the archive "h741a.zip" available on GitHub to a directory named "UHSLC" for the code below to work.

span <- 1949:2013 # years in this analysis

UHSLC = data.frame()
for(y in span){
  lcode = ifelse(y >= 2000, "i", "h")
  filename = paste0("UHSLC/", lcode, "741a", substr(y,3,4), ".dat")
  temp = read_table(file = filename,
                    na = "9999",
                    skip=1,
                    col_names = c("station", "City", "DateCode", paste(0:11))) %>%
    separate(DateCode, into = c("date","AMPM"), sep = 8) %>%
    type_convert(col_types = cols(date = col_date(format = "%Y%m%d"), AMPM = col_character())) %>%
    gather(paste(0:11), key = "hour", value = "ht_mm",convert=TRUE) %>%
    mutate(hour = ifelse(AMPM==1,hour,hour+12),
           datetime = with_tz(as.POSIXct(paste(date,hour), format=c("%F %H"),tz="UTC"),tzone = "America/New_York"),
           month = month(datetime),
           year = year(datetime)) %>%
    arrange(datetime)
  UHSLC = rbind(UHSLC, temp)
}

# detrend the tide gauge data ####
fit = lm(UHSLC$ht_mm~UHSLC$datetime)

UHSLC <- UHSLC %>%
  mutate(ht_mm_orig = ht_mm,
         ht_mm = as.numeric(ht_mm - fit$coefficients[2] * (datetime - tail(datetime,1))))


# Precipitation data ####

## National Weather Service Hourly Precipitation data for Boston 
### obtained 10 June 2018 by request from www.ncdc.noaa.gov

## Data were manually gap-filled:
### for one 5-hour, three 3-hour, and one 2-hour intervals in 1949-1972, only total accumulation was known, 
### and this total was displayed for the final hourly time point with a note indicating the true interval.
### We divided the total accumulation evenly over all hours in the interval for the "_clean" file.  

nws <- read_csv("NWS_BOS_HPRECIP_clean.csv") %>% # file is available on GitHub
  mutate(year = year(DATE),
         month = month(DATE),
         datetime = force_tz(DATE,tz="America/New_York"),
         # fill in NAs in flag variable
         `Measurement Flag` = if_else(is.na(`Measurement Flag`),"0",`Measurement Flag`),
         # set observations of trace (flag = "T") precipitation to .05 mm
         precip = HPCP 
  ) %>%   
  # remove missing values flagged with 99999, and non-measurements of 0 entered for bookkeeping/trace amnts
  filter(precip < 100 & precip > 0)             


#Step 2: Parameterize with flux data ####

load("MC_params.RDATA")

# Nitrous oxide flux data collected in the field (available at 10.6084/m9.figshare.7037531) 
# and pre-treatment lab flux data (available at 10.6084/m9.figshare.9827294) 
# were used to generate means and standard deviations for nitrous oxide fluxes for each site.
# These are variables with names beginning "baseline_".

# Mean and standard deviations of nitrous oxide fluxes (available at 10.6084/m9.figshare.9827294)
# measured after tide or storm rewetting treatments for each site were calculated for the 
# first hour after rewetting (variables beginning "tide_" and "storm_"), at the five hour mark
# (variables begeinning "tide2_" and "storm2_"), and at the 24-hour mark ("storm3_").



# step 3: Infer rewetting events ####
#Assumptions:
# N2O pulses occur with high tides if it has been at least 5 days since inundation
# storms can induce pulses if it has been 5 days since a storm and 5 days since a high tide

plat <- 4100  # platform height in the marshes corresponds to 4100 mm at the tide gauge

#Initialize:
dats = data.frame(row.names = span, tides = rep(NA,length(span)), 
                  storm1 = NA, storm2 = NA, #storm3 = NA, 
                  nmiss = NA,
                  inund = NA)

#loop through years:
for(y in span){
  # Create flag for triggers caused by tide
  # remove NAs, filter to tide heights that represent inundation in the year y 
  tide = na.omit(UHSLC[UHSLC$year == y & UHSLC$ht_mm > plat, ])
  
  # flag becomes TRUE if it has been 5 days (432000 s) since previous inundation
  tide$pulsetrig = FALSE        # initialize trigger flag...
  nt = length(tide$ht_mm)       # initialize loop length...
  for(i in 2:nt){                                                         
    tide$pulsetrig[i] = tide$datetime[i] > (tide$datetime[i-1] + 432000)
  }
  
  # trim to April-October 
  tide = tide[tide$month > 3 & tide$month < 11,] # inundated hours
  tideyr = UHSLC[UHSLC$year == y & UHSLC$month > 3 & UHSLC$month < 11, ]  # total hours
  
  # calculate tidal inundation fraction
  inund = length(tide$ht_mm) / length(na.omit(tideyr$ht_mm)) # % of time inundated
  nmiss = sum(is.na(tideyr$ht_mm))
  
  # create flag for triggers caused by rain
  # filter to year y and rain events of at least .1 mm
  rain = nws[nws$year == y & nws$HPCP > .1,] 
  
  # flag becomes TRUE if it has been 5 days (432000 s) since previous rainfall
  rain$pulsetrig = FALSE          # initialize trigger flag...
  nt = length(rain$pulsetrig)     # initialize loop length...
  for(i in 2:nt){
    rain$pulsetrig[i] = rain$datetime[i] > (rain$datetime[i-1] + 432000)
  }
  
  rain$pulsetrig2 = FALSE
  
  gg = which(rain$pulsetrig==TRUE & rain$month > 3 & rain$month < 11)
  for(i in gg){
    rain$pulsetrig2[i] = max(
      tideyr[tideyr$datetime > rain$datetime[i]-86400*5 & tideyr$datetime < rain$datetime[i],"ht_mm"]
      ) < plat
  }
  rain = rain[rain$month > 3 & rain$month < 11,]
  
  dats[paste(y),] = c(sum(tide$pulsetrig,na.rm = TRUE),
                      sum(rain$pulsetrig,na.rm = TRUE), 
                      sum(rain$pulsetrig2,na.rm = TRUE),
                      nmiss, 
                      inund)
}


# step 4: monte carlo simulation for years 1949-2013 ####

#Assumptions
# N2O flux is zero during inundation
# N2O flux is low (measured field levels) during low tide 

nmc = 20000              # select number of Monte Carlo iterations
peaty = data.frame(matrix(NA, length(span), 3*6), row.names = span) # initialize
names(peaty) = paste0(rep(c("baseline","tidepulse","stormpulse","total","storm_pct","tide_pct"), 3),
                      rep(c(2.5,50,97.5), each=6))
sandy = peaty
for(y in span){
  print(y) #counter for sanity
  storage = matrix(NA, nmc, 3)
  for(i in 1:nmc){
    #how many storms?
    n_storm = dats[paste(y), "storm2"]
    if(n_storm > 0){
      storms = rep(NA, n_storm)
      for(j in 1:n_storm){
        storm1 = rnorm(2, storm_peat, storm_peat_sd) #pulse from storms - first 2 hours
        storm2 = rnorm(6, storm2_peat, storm2_peat_sd) #pulse from storms - hours 3-8
        storms[j] = sum(storm1, storm2)
      }
      stormsum = sum(storms) #total flux due to storm rewetting
    }else{
      stormsum = 0
    }
    #how many tides?
    n_tide = dats[paste(y), "tides"]
    if(n_tide > 0){
      tides = rep(NA, n_tide)
      for(j in 1:n_tide){
        tidepulse1 = rnorm(2, tide_peat, tide_peat_sd) #pulse from high tide - first 2 hours
        tidepulse2 = rnorm(6, tide2_peat, tide2_peat_sd) #pulse from high tide - hours 3-8
        tides[j] = sum(tidepulse1, tidepulse2)
      }
      tidesum = sum(tides) #total flux due to tidal rewetting
    }else{
      tidesum = 0
    }
    
    hrs = 213*24 - 8*(n_tide + n_storm) #213 days between April 1 and Oct 31. Each pulse takes 8 hours.
    basehrs = hrs * (1-dats[paste(y), "inund"]) # hours experiencing baseline flux
    baseline = sum(rnorm(basehrs, baseline_peat, baseline_peat_sd)) #total baseline flux
    
    storage[i,] = c(baseline, tidesum, stormsum)
  }
  totals = as.data.frame(storage)
  names(totals) = c("baseline","tidepulse","stormpulse")
  totals$total = totals$baseline + totals$tidepulse + totals$stormpulse
  totals$storm_pct = 100 * totals$stormpulse / totals$total
  totals$tide_pct = 100 * totals$tidepulse / totals$total
  
  peaty[paste(y),] = c(apply(totals,2,quantile,c(.025,.5,.975))[1,],apply(totals,2,quantile,c(.025,.5,.975))[2,],apply(totals,2,quantile,c(.025,.5,.975))[3,])
}


for(y in span){
  print(y)
  storage = matrix(NA,nmc,3)
  for(i in 1:nmc){
    #how many storms?
    n_storm = dats[paste(y),"storm2"]
    if(n_storm > 0){
      storms = rep(NA, n_storm)
      for(j in 1:n_storm){
        storm1 = rnorm(2, storm_sand, storm_sand_sd) #pulse from storms - first 2 hours
        storm2 = rnorm(6, storm2_sand, storm2_sand_sd) #pulse from storms - hours 3-8
        storm3 = rnorm(16, storm3_sand, storm3_sand_sd) #pulse from storms - hours 8-24
        storms[j] = sum(storm1, storm2, storm3)
      }
      stormsum = sum(storms) #total flux due to storm rewetting
    }else{
      stormsum = 0
    }
    #how many tides?
    n_tide = dats[paste(y), "tides"]
    if(n_tide > 0){
      tides = rep(NA, n_tide)
      for(j in 1:n_tide){
        tidepulse1 = rnorm(2, tide_sand, tide_sand_sd) #pulse from high tide - first 2 hours
        tidepulse2 = rnorm(6, tide2_sand, tide2_sand_sd) #pulse from high tide - hours 3-8
        tides[j] = sum(tidepulse1, tidepulse2)
      }
      tidesum = sum(tides) #total flux due to tide rewetting
    }else{
      tidesum = 0
    }
    
    hrs = 213*24 - 8*n_tide - 24*n_storm # 213 days between April 1 and Oct 31
    basehrs = hrs*(1-dats[paste(y), "inund"])
    baseline = sum(rnorm(basehrs, baseline_sand, baseline_sand_sd))
    
    storage[i,] = c(baseline,tidesum,stormsum)
  }
  totals = as.data.frame(storage)
  names(totals) = c("baseline", "tidepulse", "stormpulse")
  totals$total = totals$baseline + totals$tidepulse + totals$stormpulse
  totals$storm_pct = 100 * totals$stormpulse / totals$total
  totals$tide_pct = 100 * totals$tidepulse / totals$total
  
  sandy[paste(y),] = c(apply(totals,2,quantile,c(.025,.5,.975))[1,],apply(totals,2,quantile,c(.025,.5,.975))[2,],apply(totals,2,quantile,c(.025,.5,.975))[3,])
}
dats = dats[paste(span),]

#write.csv(cbind(dats, peaty),"peaty_mc_ts.csv")
#write.csv(cbind(dats, sandy),"sandy_mc_ts.csv")

#step 5: monte carlo with interannual variability ####
nmc = 20000 # number of monte carlo iterations

#parameterize inundation fraction and storm/tide rewetting incidence from time series data
inund_mean = mean(dats$inund)
inund_sd = sd(dats$inund)
storm_lambda = mean(dats$storm2)
tide_lambda = mean(dats$tides)

storage = matrix(NA, nmc, 3)
for(i in 1:nmc){
  #draw numbers
  n_storm = rpois(1,storm_lambda) # number of storm rewettings
  n_tide = rpois(1,tide_lambda) # number of tidal rewettings
  in_prob = rnorm(1,inund_mean,inund_sd) # inundation fraction 
  
  #baseline flux
  hrs = 213*24 - 8*(n_tide + n_storm)  # 213 days between April 1 and Oct 31
  basehrs = hrs*(1-in_prob)
  baseline = sum(rnorm(basehrs, baseline_peat, baseline_peat_sd))
  
  #storm pulse totals
  if(n_storm > 0){
    storms = rep(NA,n_storm)
    for(j in 1:n_storm){
      storm1 = rnorm(2,storm_peat,storm_peat_sd) #pulse from storms - first 2 hours
      storm2 = rnorm(6,storm2_peat,storm2_peat_sd) #pulse from storms - hours 3-8
      storms[j] = sum(storm1,storm2)
    }
    stormsum = sum(storms)
  }else{
    stormsum = 0
  }
  
  #tide pulse totals
  if(n_tide > 0){
    tides = rep(NA,n_tide)
    for(j in 1:n_tide){
      tidepulse1 = rnorm(2,tide_peat,tide_peat_sd) #pulse from tides - first 2 hours
      tidepulse2 = rnorm(6,tide2_peat,tide2_peat_sd) #pulse from tides - hours 3-8
      tides[j] = sum(tidepulse1,tidepulse2)
    }
    tidesum = sum(tides)
  }else{
    tidesum = 0
  }
  
  #save result
  storage[i,] = c(baseline,tidesum,stormsum)
}
totals = as.data.frame(storage)
names(totals) = c("baseline","tidepulse","stormpulse")
totals$total = totals$baseline + totals$tidepulse + totals$stormpulse
totals$storm_pct = 100 * totals$stormpulse / totals$total
totals$tide_pct = 100 * totals$tidepulse / totals$total

peaty_tot=totals

storage = matrix(NA,nmc,3)
for(i in 1:nmc){
  #draw numbers
  n_storm = rpois(1,storm_lambda)
  n_tide = rpois(1,tide_lambda)
  in_prob = rnorm(1,inund_mean,inund_sd)
  
  #baseline
  hrs = 213*24 - 8*n_tide + 24*n_storm #213 days between April 1 and Oct 31
  basehrs = hrs*(1-in_prob)
  baseline = sum(rnorm(basehrs,baseline_sand,baseline_sand_sd))
  
  #storms
  if(n_storm > 0){
    storms = rep(NA,n_storm)
    for(j in 1:n_storm){
      storm1 = rnorm(2,storm_sand,storm_sand_sd) #pulse from storms - first 2 hours
      storm2 = rnorm(6,storm2_sand,storm2_sand_sd) #pulse from storms - hours 3-8
      storm3 = rnorm(16,storm3_sand,storm3_sand_sd) #pulse from storms - hours 9-24
      storms[j] = sum(storm1,storm2,storm3)
    }
    stormsum = sum(storms)
  }else{
    stormsum = 0
  }
  
  #tides
  if(n_tide > 0){
    tides = rep(NA,n_tide)
    for(j in 1:n_tide){
      tidepulse1 = rnorm(2,tide_sand,tide_sand_sd) #pulse from tides - first 2 hours
      tidepulse2 = rnorm(6,tide2_sand,tide2_sand_sd)#pulse from tides - hours 3-8
      tides[j] = sum(tidepulse1,tidepulse2)
    }
    tidesum = sum(tides)
  }else{
    tidesum = 0
  }
  
  storage[i,] = c(baseline,tidesum,stormsum)
}
totals = as.data.frame(storage)
names(totals) = c("baseline","tidepulse","stormpulse")
totals$total = totals$baseline + totals$tidepulse + totals$stormpulse
totals$storm_pct = 100 * totals$stormpulse / totals$total
totals$tide_pct = 100 * totals$tidepulse / totals$total

sandy_tot=totals

#write.csv(peaty_tot,"peaty_mc_output.csv",row.names = FALSE)
#write.csv(sandy_tot,"sandy_mc_output.csv",row.names = FALSE)


