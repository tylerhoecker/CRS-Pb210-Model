# R version of Philip Higuera's Matlab function 'CRSModel'
# Implements the constant-rate-of-supply (CRS) model described in 
# Binford et al. 1990. 
#-------------------------------------------------------------------------------
library(tidyverse)
library(reshape2)
rm(list = ls())
############## USER INPUTS FOR 'lake' AND 'supported_method' ###################
# Enter lake name here:
lake = "TL15"
# Define the number of samples to use for estimate background?
# Enter number of samples from bottom to use, otherwise make FALSE
supported_method = FALSE
# Use Radium-226 data to estimate background? 
# Make TRUE if so and include in input data, otherwise make FALSE
manual_background = FALSE 
# Save plots? 
SAVE = TRUE; library(ggplot2) ; library(cowplot)
# Create BACON input csv?
BACON = TRUE
#-------------------------------------------------------------------------------
# Select all and run <><><><><><><><><><><><><><><><><>

#-------------------------------------------------------------------------------
# Set working directory
setwd(paste("L:/1_projectsData/AK_PalEON/Analysis/Lakes/",lake,"/210Pb",sep=""))
# Read in Pb210 data
data = read.csv(paste(lake,"_Pb210_data.csv",sep=''))
Pb210_data = data[data$type=='Pb210',]
Ra226_data = data[data$type=='Ra226',]
#-------------------------------------------------------------------------------
##  Redefine data into Binford et al. 1990 terms
# [cm] top of analyzed interval
INTTOP = Pb210_data[,"topDepth"] 
# [cm] bottom of analyzed interval
INTBOT = Pb210_data[,"botDepth"]    
# [DPM g^1] total 210Pb activity
TOTACT = Pb210_data[,"DPMg"]    
# [DPM g^1] std on total 210Pb activity, i.e.  210Pb counting error
SDACT  = Pb210_data[,"error"]                 
# [g/cm3] bulk density                                       
RHO   = Pb210_data[,"bulkDens"]
# number of samples
n = length(INTTOP)                  

#-------------------------------------------------------------------------------
## Set up parameters
# decay constant for 210Pb
k = 0.03114   
# parameter 'p' for smoothing spline
p_210 = 1/(1+(mean(diff(INTTOP)))^3/60)   

#-------------------------------------------------------------------------------
# or define background
## Estimate background 210Pb activity
# total activity i^th sample to the end

if (supported_method == 0 & manual_background == FALSE){
  for (i in (n-2):2){
    bkgACT = TOTACT[i:length(TOTACT)]
    bkgMeanACT= mean(bkgACT)
    bkgStdACT = sd(bkgACT)
    if ((bkgMeanACT + bkgStdACT) < TOTACT[i-1]) {
      # the ith sample and below are considered "background"
      BACKGROUND_in = i
      break
      } else if (i == 2) {
      print("Background not distinct using Binford 1990 criteria (p. 257)...")
    } 
  }
}

if (supported_method > 0) BACKGROUND_in = n - supported_method + 1


# [DPM g^1] supported 210Pb activity
if (manual_background == FALSE){
  BACKGROUND = mean(TOTACT[BACKGROUND_in:length(TOTACT)])
  SDSUP = sd(TOTACT[BACKGROUND_in:length(TOTACT)])
  if (is.na(SDSUP)) SDSUP = Pb210_data[BACKGROUND_in,"error"]
}


if (manual_background == TRUE){
  BACKGROUND = mean(Ra226_data[,'DPMg'])
  SDSUP = sd(Ra226_data[,'DPMg'])
  BACKGROUND_in = which(TOTACT < BACKGROUND)[1]
  if (is.na(BACKGROUND_in)) BACKGROUND_in = n
}

# UNSUPACT = Unsupport 210Pb activity [DPM / g]
UNSUPACT = TOTACT - BACKGROUND
UNSUPACT[BACKGROUND_in:length(UNSUPACT)] = 0 # Background samples set to 0 
UNSUPACT[UNSUPACT<0] = 0 # negative values reset to 0

# RHOXCEE = 210Pb in 1 cm3 of wet [DPM / cm^3]  
RHOXCEE = RHO * UNSUPACT  

# INTEGSEG = 210Pb activity integreated over the thickness of the analyzed
# intervals. [DPM / cm^2] 
INTEGSEG = RHOXCEE * (INTBOT - INTTOP) 

# INTEGINT = Interpolated 210Pb activity integreated over the thickness 
# of the non-analyzed intervals. [DPM / cm^2] 
#
# The interpolation is done by the log rule of integration. This value 
# will be zero in cores with analyses done on contiguous intervals, which 
# is the recommended procedure for use with the CRS model. 

# Note that INTEGINT(i) corresponds to the intervel BELOW the ith sample. 
# I.e. INTEGINT(BACKGROUND_in:end) should be zero, because
# UNSUPACT(BACKGROUND_in:end) is zero by definition   
INTEGINT = ( (RHOXCEE[2:length(RHOXCEE)] - RHOXCEE[1:(length(RHOXCEE)-1)]) / 
             log( RHOXCEE[2:length(RHOXCEE)] / RHOXCEE[1:(length(RHOXCEE)-1)] ))*
                  (INTTOP[2:length(RHOXCEE)] - INTBOT[1:(length(RHOXCEE)-1)]) 

# Behavior of this expression is erratic near BACKGROUND_in due to
# undefined terms (log(0), x/0). Below, I explicitly define INGEGINT=0 
# between background samples, and use a simple linear average to derive 
# INTEGINT for the unsampled sediment just above background. 
INTEGINT[BACKGROUND_in:n] = 0
INTEGINT[BACKGROUND_in-1] = (INTTOP[BACKGROUND_in] - INTBOT[BACKGROUND_in-1])* 
  RHOXCEE[BACKGROUND_in-1]/2 

# note this last term is an average of RHOXCEE(BACKGROUND_in-1) and RHOXCEE(BACKGROUND_in)=0

# SUMBOT, SUMTOP = Total 210Pb activity integrated from the bottom of the 
# core to the bottom or top (respectively) of the interval. [DPM / cm^2] 

# Requires iterative calculation
SUMBOT = c(rep(0,n))
SUMTOP = c(rep(0,n))
for (i in (BACKGROUND_in-1):1){
  SUMBOT[i] = SUMTOP[i+1] + INTEGINT[i] 
  SUMTOP[i]= SUMBOT[i] + INTEGSEG[i] 
  
}

# TTOP, TBOT = Time since core was taken [yr] at the top or bottom of the
# interval (respectively)
TTOP = log(SUMTOP[1] / SUMTOP) / k
TTOP[BACKGROUND_in:length(TTOP)] = NA # Age not calculated in background
TBOT = log(SUMBOT[1] / SUMBOT) / k
TBOT[BACKGROUND_in:length(TBOT)] =  NA # Age not calculated in background

# YBPTOP = ages in calander years before CE 1950 [cal ybp] 
YBPTOP = TTOP - (Pb210_data$coreDate[1]-1950)   
YBPBOT = TBOT - (Pb210_data$coreDate[1]-1950)

# SEDRATE = Accumulation rate of bulk sediment [g / cm^2 / yr]
SEDRATE = k * ((SUMTOP + SUMBOT)/2) / UNSUPACT
SEDRATE[BACKGROUND_in:length(SEDRATE)] = NA # Age not calculated in background


## Error propagation
# In this section, error propagation calculations are applied in order to
# estimate the age uncertainty from uncertainty inherent in the 210Pb
# measurements. 

# SDRHO = Standard deviation of RHO. Assumed = 0 as in Binford 1990.
SDRHO = c(rep(0,n))

# SDUNSUP = Standard deviation of unsupported 210Pb activity
SDUNSUP = (SDACT^2 + SDSUP^2)^0.5 

# SDRHOCEE = Standard deviation of variable RHOCEE.
SDRHOCEE = (SDRHO^2 * UNSUPACT^2  +  SDUNSUP^2 * RHO^2)^0.5

# SDINTSEG = Standard deviation of the vairiable INTSEG.
# Note that Binford has unneeded squares and square roots that cancel
SDINTSEG = SDRHOCEE * (INTBOT-INTTOP) 

# SDINTINT = Standard deviation of the variable INTEGINT.
# I modified the error calcs slightly. Based on my own error propagation
# calculations, the correct coefficient on this expression is (2^-0.5).
# Binford uses 0.5, but I believe his equation is wrong because it
# doesn't match his own results. Anyway, the expression I use below is
# more conservative than either (coefficient simply set to 1). 
SDINTINT2 = SDINTSEG[1:(length(SDINTSEG)-1)] * 
           ( INTTOP[2:length(INTTOP)] - INTBOT[1:(length(INTBOT)-1)] )

# SDSUMBOT, SDSUMTOP = Standard deviation of the variables SUMBOT, SUMTOP
# Sets values at and below BACKGROUND_in at 0, which is sensible since all
# variables entering the variance calculations for those samples are
# defined theoretically (i.e., background has 0 unsupported activity), so 
# have variance 0 themselves.
SDSUMBOT = c(rep(0,n))
SDSUMTOP = c(rep(0,n)) 
for (i in (BACKGROUND_in-1):1){
  SDSUMBOT[i] = (SDINTINT2[i]^2 + SDSUMTOP[i+1]^2)^0.5 
  SDSUMTOP[i] = (SDINTSEG[i]^2 + SDSUMBOT[i]^2)^0.5 
  
}

# SDTTOP, SDTBOT, SDSEDRT = Standard deviation of corresponding variables
SDTTOP = ((k^-2 * SUMTOP[1]^-2 * SDSUMTOP[1]^2) + 
            (k^-2 * SUMTOP^-2 * SDSUMTOP^2))^0.5
SDTBOT = ((k^-2 * SUMTOP[1]^-2 * SDSUMTOP[1]^2) + 
            (k^-2 * SUMBOT^-2 * SDSUMBOT^2))^0.5 
SDSEDRT = (k^2 * (SUMTOP/UNSUPACT) * 
             ((SDSUMTOP/SUMTOP)^2 + (SDUNSUP/UNSUPACT)^2))^0.5

# Since ages are not calculated for samples with 0 unsupported activity,
# standard deviations for those values are undefined and therefore
# removed entirely.
SDTTOP[BACKGROUND_in:length(SDTTOP)] = NA # Age not calculated in background
SDTBOT[BACKGROUND_in:length(SDTTOP)] = NA # Age not calculated in background
SDSEDRT[BACKGROUND_in:length(SDTTOP)] = NA # Age not calculated in background          

# Combine relevant results and output to csv file
output = data.frame('topDepth'=INTTOP, 'botDepth'=INTBOT,
                    'ageTop'=YBPTOP, 'StdTop'=SDTTOP, 
                    'ageBot'= YBPBOT, 'StdBot' = SDTBOT,
                    'TotalPb210'=TOTACT, 'SdPb210'=SDACT,
                    'bkgrnd' = c(rep(NA,BACKGROUND_in),rep(1,(n-BACKGROUND_in))))
output

# Plot
xlimit = c(max(output$botDepth[!is.na(output$ageTop)])+.5,min(output$topDepth)-.25)
agePlot = 
  ggplot(output, aes(x= topDepth, y= ageTop)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax = ageTop + 1.96*StdTop, ymin = ageTop - 1.96*StdTop), width=0.1) +
  xlim(xlimit) + scale_y_reverse() +
  xlab("Sample Bottom Depth (cm)") + ylab("Age (years BP)") +
  ggtitle(expression(paste(""^210,Pb,'-CRS-based age estimates'))) +
  theme_bw(base_size = 14)


bkgLabel = paste('est. supported activity:',round(BACKGROUND,2), '(+/-',round(SDSUP,2),')')
if (manual_background == TRUE){
  pb210Plot=
    ggplot(output, aes(x= topDepth)) +
    geom_hline(aes(yintercept = BACKGROUND),color="red2") +
    geom_hline(aes(yintercept = BACKGROUND + SDSUP),color="red2",linetype= 2) +
    geom_hline(aes(yintercept = BACKGROUND - SDSUP),color="red2",linetype= 2) +
    geom_point(aes(y= TotalPb210), size= 2) +
    geom_errorbar(aes(ymax = TotalPb210 + SdPb210, ymin = TotalPb210 - SdPb210), width=0.25) +
    geom_point(data= Ra226_data, 
               aes(y=DPMg), size = 2, color ='turquoise4') +
    geom_errorbar(data= Ra226_data, color ='turquoise4', width=0.25,
                  aes(ymax= DPMg + error, ymin= DPMg - error)) +
    geom_text(aes(mean(output$topDepth),BACKGROUND, label = bkgLabel, vjust = -1), color = 'red2') + 
    coord_trans(y="log10") +
    xlab("Sample Bottom Depth (cm)") + ylab(expression(paste(""^210,Pb,' activity (DPM/g)'))) +
    ggtitle(expression(paste("Total "^210,Pb,' activity'))) +
    theme_bw(base_size = 14)
} else {
  bkgPoints = output[is.na(output[,'ageTop']),] 
  pb210Plot=
    ggplot(output, aes(x= topDepth)) +
    geom_hline(aes(yintercept = BACKGROUND),color="red2") +
    geom_hline(aes(yintercept = BACKGROUND + SDSUP),color="red2",linetype= 2) +
    geom_hline(aes(yintercept = BACKGROUND - SDSUP),color="red2",linetype= 2) +
    geom_point(aes(y= TotalPb210), size= 2) +
    geom_errorbar(aes(ymax = TotalPb210 + SdPb210, ymin = TotalPb210 - SdPb210), width=0.25) +
    
    geom_point(data = bkgPoints, aes(y= TotalPb210), size= 2, color = 'red2') +
    geom_errorbar(data = bkgPoints, aes(ymax = TotalPb210 + SdPb210, ymin = TotalPb210 - SdPb210), 
                  color = 'red2', width=0.25) +
    
    geom_text(aes(mean(output$topDepth),BACKGROUND, label = bkgLabel, vjust = -1), color = 'red2') + 
    coord_trans(y="log10") +
    xlab("Sample Bottom Depth (cm)") + ylab(expression(paste(""^210,Pb,' activity (DPM/g)'))) +
    ggtitle(expression(paste("Total "^210,Pb,' activity'))) +
    theme_bw(base_size = 14)
  }

bothPlots = plot_grid(pb210Plot, agePlot, align = "v", nrow = 1,labels = paste(lake))


if (SAVE == TRUE){
  # Save to .csv file in lake's Pb210 folder
  write.csv(output, paste(lake,"_Pb210_crsR_Results.csv",sep=''), row.names = F)
  
  # Save plot
  pdf(file = paste(lake,"_Pb210_crsR_Results.pdf",sep=''), width = 10, height = 5)
  print(bothPlots)
  dev.off()
}

if (BACON == TRUE){
  len = length(YBPTOP[!is.na(YBPTOP)])
  bcn = data.frame('labID'= c('surface',rep('Flett',len)),
                   'age'= c(-(Pb210_data$coreDate[1]-1950),round(YBPBOT[1:len])), 
                   'error'= c(0,round(SDTBOT[1:len])),
                   'depth'= c(0,INTBOT[1:len]),
                   'cc'= rep(0,len+1))
  write.csv(bcn, 
            paste('L:/1_projectsData/AK_PalEON/Analysis/Lakes/z_Combined_Analyses/Chronologies/BACON/Cores/',lake,'/',
            lake,'.csv',sep=''), row.names = F)
}

dev.off()
bothPlots



