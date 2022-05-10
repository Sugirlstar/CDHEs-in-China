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
# * The processes of figure4 are the same among drought, heatwave and CDHE,
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

# read shp files
Chinasp = readShapePoly(paste0(pathname,"shpfile/china/China")) # not given in the repository
provincesp = readShapePoly(paste0(pathname,"shpfile/china/PROVINCE_region")) # not given in the repository

# fundamental parameters
{
  xx = 72
  yy = 128
  days = length(P[1,1,])
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

# 01 cuculate linear trend -------------
objx = DRfre
ck = DRFREh
rst=bp1(s,e,ck)
bp=rst[1]-s+1 
tend2=array(dim=c(xx,yy))
pfre2=array(dim=c(xx,yy))

for(x in 1:xx)
  for(y in 1:yy)
    if (TRUE %in% (objx[x,y,]>=0))
    {
      tend2[x,y]=tendency(objx[x,y,(bp+1):(e-s+1)])[[1]]
      pfre2[x,y]=classtype(tendency(objx[x,y,(bp+1):(e-s+1)])[[2]])
    }

tend2r = raster(tend2) 
extent(tend2r) = c(72, 136, 18, 54) 
crs(tend2r) = crs(r) 

yind=xind=pfre1ind=pfre2ind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) 
    xind=c(xind,54.25-0.5*x) 
    pfre2ind=c(pfre2ind,pfre2[x,y])
  }
##01End##

# 02 Putout_linrar trend -------------
max(tend2,na.rm=T)
min(tend2,na.rm=T)

colr1=rev(colorRampPalette(brewer.pal(9,"GnBu")[1:9])(150)) 
colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(150)
#tend2col = c(colr1,colr2)
tend2col = cf(c(colr1,colr2),whitesite="#FFF7EC",
              c(-0.060, 0.08,0),z0=0) 

png(
  file = "fig4_tend_drought.png",
  width = 2000,
  height = 1300,
  res = 72 * 3
)
par(mar=c(1,1,1,2))
plot(
  tend2r,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=tend2col,
  legend.width=1,
  zlim=c(-0.060, 0.08),
  axis.args=list(at=seq(-0.06, 0.08, 0.02),
                 labels=seq(-0.06, 0.08, 0.02),
                 cex.axis=1)
)

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) 
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)

points(yind,xind,pch=pfre2ind,cex=0.5) 
points(124,19,pch=4,cex=1.5) 

dev.off()

##02End##

# 03 caculate M-K trend -------------
zz = MK.raster(objx,"year") 
z <- raster(zz[,,1]) 
extent(z)<- c(72, 136, 18, 54)
res(z)<-res(r) 
projection(z)=projection(r)
sigz=array(dim=c(xx,yy))
for(x in 1:xx)
  for(y in 1:yy)
    if (TRUE %in% (zz[x,y,]>=0))
      sigz[x,y]=classtype(zz[x,y,2])

yind=xind=sigind=NULL
for(x in 1:xx)
  for(y in 1:yy)
  {
    yind=c(yind,71.75+0.5*y) 
    xind=c(xind,54.25-0.5*x)
    sigind=c(sigind,sigz[x,y])
  }

##03End##

# 04 Putout_M-K trend -------------
colr=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
MKcol =  cf(colr,whitesite=colr[50],c(-0.5,0.5,0),z0=0) 

png(
  file = "fig4_MK_drought.png",
  width = 2000,
  height = 1300,
  res = 72 * 3
)
par(mar=c(1,1,1,2))
plot(
  z,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=MKcol,
  legend.width=1,
  zlim=c(-0.5,0.5),
  axis.args=list(at=seq(-0.5, 0.5, 0.2),
                 labels=seq(-0.5, 0.5, 0.2),
                 cex.axis=1)
)

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) 
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)

points(yind,xind,pch=sigind,cex=0.5) 
points(124,19,pch=4,cex=1.5) 
#text(130,19.1,"α=0.05",cex=1.2) 

dev.off()

##04End##
          
          
          
          
