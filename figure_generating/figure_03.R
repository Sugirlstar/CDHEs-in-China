# Note:
# 1. This script was written in Rstudio (Version 1.4.1717), 
#    but was recommended to be ran in the original R for saving memory space.
# 
# 2. We use foldable code sections in Rstudio 
#    to break the large source file into a set of discrete regions 
#    for easy navigation between them. 
#    Each section begins with a group of pound signs (#) 
#    following the title and several trailing dashes (-), 
#    and ends with pound signs (be like: ##End##)
#   
# 3. Comments explaining what the code is doing are above the lines or instantly follow the code.
# 
# 4. To reproduce, just change the "pathname" (line 21) to where you store the file we offered.
#
# * The processes of figure3 and table1 are the same among drought, heatwave and CDHE,
#   here we take drought as the example

# 00 Background setting -------------

pathname = "D:/CDHEs_coding"

library(tcltk)
library(raster)
library(sp)
library(maptools)
library(rgdal)
library(RColorBrewer) 
library(trend)
library(lattice)

# load necessary workspace and set location
load(paste0(pathname,"/SPI.RData"))
load(paste0(pathname,"/STI.RData"))
load(paste0(pathname,"/P_China_1dy_0.5deg_1961_2020.RData")) # precipitation
load(paste0(pathname,"/Tx_China_1dy_0.5deg_1961_2020.RData")) # temperature
load(paste0(pathname,"/Metrics_drought.RData"))
load(paste0(pathname,"/supplement_materials/crs.RData")) # the coordinate reference system
setwd(pathname)

# fundamental parameters
{
  xx=72
  yy=128
  days=length(SPIarray[1,1,])
  st=90 
  da = 10 
  da2 = 5
  thr1 = -1
  thr2 = 0 
  s=1961
  e=2020
  LL=e-s+1 
}
##00End##

# 01 Functions praparing -------------
figure2out = function(objx,ylabname,aty) 
{
  atylabel=as.character(aty)
  ylimget = c(min(aty,na.rm=T),max(aty,na.rm=T))
  
  objx=-DRSTGh
  rst=bp1(s,e,objx)
  #02-01-01-a
  if(rst[1]!=0) 
  {
    FREt=c(s,rst[1,(1:(ncol(rst)-1))],e)
    FREy=FREt
    FREy[1]=FREt[1]*rst[2,1]+rst[3,1]
    for(i in 2: length(FREt) )
      FREy[i]=FREt[i]*rst[2,i-1]+rst[3,i-1]
    FREp1=rst[4,1] ###
    FREp2=rst[4,2] ###
  }
  #02-01-01-b
  if(rst[1]==0) 
  {
    tend=tendency(objx)
    tendx=c(s:e)
    tendy=tendx*tend[[1]]+tend[[3]] 

    mktau=mk.test(objx)$estimates[3]
    mkpvalue=mk.test(objx)$p.value 
  }
  
  
  atx = seq(1960,2020,10)
  atxlabel = as.character(atx)
  plot(c(s:e),objx,type="o",
       xlim = c(1960,2020),
       ylim = ylimget,
       ylab = ylabname,
       pch=20, xlab="Year",
       xaxs = "i",yaxs = "i", xaxt = "n", yaxt = "n",
       lwd=0.8,
       mgp=c(2.4,0.4,0),cex.lab=1.2) 
  
  axis(1,at=atx,atxlabel)
  
  axis(2,at=aty,atylabel,las = 2)
  
  
  #Linear trends
  if(rst[1]!=0)
  {
    lines(FREt, FREy, lwd=2,col="red")
    points(FREt[2],FREy[2],pch=17,cex=3)
  }
  #M-K
  if(rst[1]==0)
    lines(tendx,tendy,lwd=2,col="blue")
  
}
##01End##

# 02 Find a suitable ylim (run before generating figures) -------------
# objx=DRper*100 #objx=DRFREh
# ydur = max(objx,na.rm=T) - min(objx,na.rm=T)
# getylim1 = c(min(objx,na.rm=T), max(objx,na.rm=T)+ydur*0.2)
# getylim1;ydur/4

##02End##

# 03 putout ------------- 
dev.new(width=1300,height=3000,res=72)   
par( mfrow = c(4,1),mar=c(4,4.8,0.5,1.2) )       
figure2out(DRFREh,"Frequency (times)",seq(0.6,2.2,0.4))        
figure2out(DRDURh,"Duration (days)",seq(20,100,20))        
figure2out(-DRSTGh,"Severity",seq(0.6,1.6,0.2))        
figure2out(DRper*100,"Coverage (%)",seq(40,100,10)) 
#manually save as .eps

  # write in csv
  rst1=bp1(s,e,DRFREh)
  rst2=bp1(s,e,DRDURh)
  rst3=bp1(s,e,DRSTGh)
  rst4=bp1(s,e,DRper)
  csvout = cbind(rst1,rst2,rst3,rst4)
  write.table (csvout,sep=",",file ="dr_char.csv") 

##03End##    
          
             
                
                
