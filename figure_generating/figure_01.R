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
# 4. To reproduce, just change the "pathname" (line 19) to where you store the file we offered.


# 00 Background setting -------------

pathname = "D:/CDHEs_coding"

library("sp")
library("maptools")
library("rgdal")
library("raster")
library("tcltk")
library("RColorBrewer") 

# load necessary workspace and set location
load(paste0(pathname,"/SPI.RData"))
load(paste0(pathname,"/STI.RData"))
load(paste0(pathname,"/P_China_1dy_0.5deg_1961_2020.RData")) # precipitation
load(paste0(pathname,"/Tx_China_1dy_0.5deg_1961_2020.RData")) # temperature
load(paste0(pathname,"/supplement_materials/crs.RData")) # the coordinate reference system
setwd(pathname)
# read shp files
Chinasp = readShapePoly(paste0(pathname,"shpfile/china/China")) # not given in the repository
provincesp = readShapePoly(paste0(pathname,"shpfile/china/PROVINCE_region")) # not given in the repository

# fundamental parameters
{
  xx=72
  yy=128
  Longtirange = seq(72,136,0.5)
  Lattirange = seq(18,54,0.5)
  
  # coordinates showing on figures
  xL = seq(72, 136, 5)
  xlabel = paste(xL, "°", sep = "") 
  yL = seq(18, 54, 5)
  ylabel = paste(yL, "°", sep = "")
}
##00End##


# 01 station information ----------
stationInfo = read.table(paste0(pathname,"/supplement_materials/SURF_CHN_MUL_HOR_STATION.csv"), header=F, sep=",") #not given in the repository
stationID = stationInfo[,2]
stationLatti = stationInfo[,4]
stationLongti = stationInfo[,5]
snum = length(stationLatti)

# convert into real value
conv = function(x)
{
  xc = as.character(x)
  clen = nchar(xc)
  if (clen == 5)
  {
    x_min = as.numeric(substr( xc, 1, 3))
    x_sec = as.numeric(substr( xc, 4, 5))
  }else
  {
    x_min = as.numeric(substr( xc, 1, 2))
    x_sec = as.numeric(substr( xc, 3, 4))
  }
  x_out = x_min+x_sec/60  
  return(x_out)  
}

LattiValue = rep(NA,snum)
LongtiValue = rep(NA,snum)
for ( i in 1:snum )
{
  LattiValue[i] = conv(stationLatti[i])
  LongtiValue[i] = conv(stationLongti[i])
}
  
# judge index 
foundindex = function(vv)
{
  if( 0 %in% vv )
    ind = which(vv==0)
  else
  {
    groups = embed(vv,2)
    mutres = groups[,1]*groups[,2]
    ind = which(mutres < 0)
  }
  return(ind)
}

staLatti = rep(NA,snum)
staLongti = rep(NA,snum)
for ( i in 1:snum )
{
  Longtigap = Longtirange - LongtiValue[i]
  staLongti[i] = foundindex(Longtigap)
  Lattigap = Lattirange - LattiValue[i]
  staLatti[i] = foundindex(Lattigap)
}

##00End##

# 02 plotting ----------

stationmap = array(0,dim=c(xx,yy)) # set blank map

for ( i in 1:snum )
{
  x = 73 - staLatti[i]
  y = staLongti[i]
  stationmap[x,y] = stationmap[x,y]+1
}

# classify
classifyType = function(v,a1,a2,a3,a4)
{
    c=NA
      if(v==1)
        c=a1 else 
          if(v==2)
            c = a2 else
              if(v==3)
                c = a3 else
                  if(v>=4)
                    c = a4
              
        return(c)
}

a1 = "#006d2c"
a2 = "#08519c"
a3 = "#810f7c"
a4 = "#bd0026"


sigmap = array(NA,dim=c(xx,yy))
for(x in 1:xx)
  for(y in 1:yy)
      sigmap[x,y] = classifyType(stationmap[x,y],a1,a2,a3,a4)
  

yind=xind=sigindex=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) 
    xind=c(xind,54.25-0.5*x) 
    sigindex=c(sigindex,sigmap[x,y])
  }


par(mar=c(1,1,1,1))
mapr = raster(SPIarray[,,1]) 
extent(mapr) = c(72, 136, 18, 54)
plot(mapr,
     ylim=c(17,55),xlim=c(72,136),
     xaxs="i",yaxs="i",
     xaxt="n",yaxt="n",
     xlab="",ylab="",
     mgp=c(0,0,0))

lines(provincesp)
points(yind,xind,pch=22, col=sigindex, bg=sigindex, cex=0.7,) 
