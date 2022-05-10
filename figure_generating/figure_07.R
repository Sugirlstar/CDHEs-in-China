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

library(tcltk)
library(raster)
library(sp)
library(maptools)
library(rgdal)
library(RColorBrewer)
library(trend)
library(lattice)
library(VineCopula)
library(copula)

# load necessary workspace and set location
load(paste0(pathname,"/SPI.RData"))
load(paste0(pathname,"/STI.RData"))
load(paste0(pathname,"/P_China_1dy_0.5deg_1961_2020.RData")) # precipitation
load(paste0(pathname,"/Tx_China_1dy_0.5deg_1961_2020.RData")) # temperature
load(paste0(pathname,"/supplement_materials/crs.RData")) # the coordinate reference system
load(paste0(pathname,"/Metrics_CDHE.RData"))
load(paste0(pathname,"/phase results/droughtFlag.RData"))
load(paste0(pathname,"/phase results/heatwaveFlag.RData"))
setwd(pathname)

# read shp files
Chinasp = readShapePoly(paste0(pathname,"shpfile/china/China")) # not given in the repository
provincesp = readShapePoly(paste0(pathname,"shpfile/china/PROVINCE_region")) # not given in the repository

# fundamental parameters
{
  xx = 72
  yy = 128
  days = length(P[1,1,])
  
}

##00End##

# 01 Functions praparing -------------
findNestOne = function(v,dot1) 
{
  v[which(is.na(v)==TRUE)] = 0
  i = dot1 - 1
  outputIndex = i
  
  if(dot1 == 1)
  {
    print("dot1 > 1")
  }else
  {
    while( v[i] != 1 )
    {
      i=i-1
      outputIndex = i
      if( i == 0 )
      {
        outputIndex = NA
        break
      }
    }
    return(outputIndex)
  }
}

##01End##

# 02 Calculate events with drought leading -------------

leadTime = array( dim=dim(Flag3) ) 
multiEvent = array( dim=dim(Flag3) )

pb = tkProgressBar(title="Progress",label="Completed %", 
                   min=0, max=100, initial = 0, width = 300)
for(x in 1:xx)
{
  for(y in 1:yy)
    if(TRUE %in% (Flag[x,y,] >= -999))
    {
      CDHindex = which( Flag3[x,y,] == 1 ) 
      
      if( length(CDHindex) == 0 )
        next
      
      for( i in 1:length(CDHindex) )
      {
        dot = CDHindex[i]
        
        if( is.na(FD3[ x,y,dot ]) == TRUE )
          FDdotindex = findNestOne( FD3[x,y,], dot ) else
            if( FD3[ x,y,dot ] ==  0)
              FDdotindex = findNestOne( FD3[x,y,], dot ) else 
                FDdotindex = dot
              
              if( is.na(FH3[ x,y,dot ]) == TRUE )
                FHdotindex = findNestOne( FH3[x,y,], dot ) else
                  if( FH3[ x,y,dot ] ==  0)
                    FHdotindex = findNestOne( FH3[x,y,], dot ) else
                      FHdotindex = dot
                    
                    leadTime[x,y,dot] = FHdotindex - FDdotindex  
                    
      }
      
      if( length(CDHindex)>1 )
      {
        test = CDHindex
        
        for(i in 1:(length(test)-1) )
        {
          intervalDays = (test[i]+1):(test[i+1]-1)
          
          t1 = FD[ x,y,intervalDays]
          t2 = FH[ x,y,intervalDays]
          
          if( length(which(t1>0)) == length(t1) ) 
          {
            leadTime[ x,y, test[i+1] ] = length(intervalDays)+1
            multiEvent[x,y,test[i+1] ] = "D"
          }
          if( length(which(t2>0)) == length(t2) ) 
          {
            leadTime[ x,y, test[i+1] ] = -(length(intervalDays)+1)
            multiEvent[x,y,test[i+1]] = "H"
          }
          
        }
      }
      
    }
  info = sprintf("Completed %d%%", round(x*100/xx))  
  setTkProgressBar(pb, value = x*100/xx, 
                   title = sprintf("Progress (%s)",info),label = info) 
  
}
close(pb)  

DorHarray = array( dim = c(xx,yy) )
Dmean = array( dim = c(xx,yy) )
multis = array( dim = c(xx,yy) )
multisType = array( dim = c(xx,yy) )
multisD = array( dim = c(xx,yy) )
multisH = array( dim = c(xx,yy) )

pb = tkProgressBar(title="Progress",label="Completed %", 
                   min=0, max=100, initial = 0, width = 300)
for(x in 1:xx)
{
  for(y in 1:yy)
    if(TRUE %in% (Flag[x,y,] >= -999))
    {
      
      CDHindex = which( Flag3[x,y,] == 1 ) 
      if( length(CDHindex) == 0 )
      {
        DorHarray[x,y] = NA
        Dmean[x,y] = NA
        next
      }
      
      DorHarray[x,y] = length(which(leadTime[x,y,CDHindex]>0)) / length(CDHindex) 
      Dmean[x,y] = mean(leadTime[x,y,CDHindex],na.rm=TRUE) 
      multis[x,y] = length(which( multiEvent[x,y,CDHindex] == "D" ) ) + 
        length(which( multiEvent[x,y,CDHindex] == "H" ) )
      multisD[x,y] = length(which( multiEvent[x,y,CDHindex] == "D" ) )
      multisH[x,y] = length(which( multiEvent[x,y,CDHindex] == "H" ) )
      
      if( multis[x,y] > 0 )
        multisType[x,y] = length(which( multiEvent[x,y,CDHindex] == "D" ) )/multis[x,y]
      
    }
  
  info = sprintf("Completed %d%%", round(x*100/xx))  
  setTkProgressBar(pb, value = x*100/xx, 
                   title = sprintf("Progress (%s)",info),label = info) 
  
}
close(pb)  

DorHarray = DorHarray *100 
DorHarrayr = raster(DorHarray) 
extent(DorHarrayr) = c(72, 136, 18, 54)

##02End##
 
# 03 Distribution plotting -------------
colr = colorRampPalette(brewer.pal(11, "Set3"))(11)[c(4,6,8,11,7)] 

# Manually save as .eps
dev.new(width = 2000,height = 1300)
par(mar=c(1,1,1,2))
plot(
  DorHarrayr,
  ylim=c(17,55),xlim=c(72,136),
  xaxs="i",yaxs="i",
  xaxt="n",yaxt="n",
  xlab="",ylab="",
  mgp=c(0,0,0),
  col=colr,
  legend.width=1,
  zlim=c(0,100),breaks = seq(0, 100, 20),
  axis.args=list(at=seq(0, 100, 20),
                 labels=seq(0, 100, 20),
                 cex.axis=1,las=0)
) 

xL = seq(72, 136, 5)
xlabel = rep("",length(xL)) 
yL = seq(18, 54, 5)
ylabel = rep("",length(yL))

axis(1, xL, xlabel)
axis(2, yL, ylabel,las=2)
lines(Chinasp)
lines(provincesp)
##03End##  
   
# 04 Process plotting in Yunnan -------------
x=59
y=57
testr=r
testr[x,y]=10000
plot(testr)

indx = c(19000:19080)  #20130107-20130328
dimnames(TX)[[3]][indx]

hstime=indx
hstime[which(FH[x,y,indx]!=1)]=NA
dstime=indx
dstime[which(FD[x,y,indx]!=1)]=NA

dev.new(width=1300,height=800,res=72) 
par(mar=c(3,4.5,0.5,4.5))

plot(indx,SPIarray[x,y,indx], type = "l", xaxt = "n", yaxt = "n",pch=20,
     ylab = "", xlab = "Date",col="deepskyblue", ylim=c(-4,4),
     cex=1.2,lwd=2,mgp=c(1.5,0.8,0),
     xaxs = "i",yaxs = "i")

lines(indx,TSI[x,y,indx], type = "l", pch=20, col="brown1",lwd=2,cex=1.2)
abline(h=-1,lty=2)
abline(h=1,lty=2)

points(dstime,SPIarray[x,y,dstime],pch=20,col="deepskyblue",cex=1.5)
points(hstime,TSI[x,y,hstime],pch=18,col="brown1",cex=1.5)


axis(at=seq(-4,4,1)
     ,label=as.character(seq(-4,4,1)),
     side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
     las=2)
axis(at= seq(19000,19080,20), 
     label=c("7  Jan 2013","27 Jan 2013",
             "16 Feb 2013","8  Mar 2013","28 Mar 2013"),
     side=1,mgp = c(2.8, 0.4, 0))

mtext("SPI/-STI",
      side = 2, line = 2.2, cex.lab = 1, 
      las = 0,cex=0.8)

par(new = TRUE)
plot(indx,CDHI[x,y,indx], type = "o", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", pch=17, col="purple2",lwd=2,xaxs = "i",yaxs = "i",
     ylim=c(0.94,1.04))

axis(at=seq(0.980,1,0.01),
     label=as.character(seq(0.980,1,0.01)),
     side = 4, mgp=c(1,0.8,0),
     col.axis="black",cex.axis=0.8,las=2)

mtext("CDHId" ,
      side = 4, line = 2.2, cex.lab = 1,
      las = 0, cex=0.8)

legend(x="top",
       c("SPI","STI","CDHId"),
       lty=c(1,1,1),pch=c(NA,NA,17), lwd=c(2,2,2),
       col=c("deepskyblue","brown1","purple2"),horiz=T,
       bty="n",cex=1.5, x.intersp=c(0.2,0.2,0.2), inset=c(0.01,0) )

##04End##

# 05 Process plotting in NE -------------
x=13
y=103
testr=r
testr[x,y]=10000
plot(testr)

indx = c(21220:21300)  #20130107-20130328
dimnames(TX)[[3]][indx]

hstime=indx
hstime[which(FH[x,y,indx]!=1)]=NA
dstime=indx
dstime[which(FD[x,y,indx]!=1)]=NA

dev.new(width=1300,height=800,res=72)
par(mar=c(3,4.5,0.5,4.5))

plot(indx,SPIarray[x,y,indx], type = "l", xaxt = "n", yaxt = "n",pch=20,
     ylab = "", xlab = "Date",col="deepskyblue", ylim=c(-3,4),
     cex=1.2,lwd=2,mgp=c(1.5,0.8,0),
     xaxs = "i",yaxs = "i")

lines(indx,TSI[x,y,indx], type = "l", pch=20, col="brown1",lwd=2,cex=1.2)
abline(h=-1,lty=2)
abline(h=1,lty=2)

points(dstime,SPIarray[x,y,dstime],pch=20,col="deepskyblue",cex=1.5)
points(hstime,TSI[x,y,hstime],pch=18,col="brown1",cex=1.5)


axis(at=seq(-4,4,1)
     ,label=as.character(seq(-4,4,1)),
     side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
     las=2)
axis(at= seq(21220,21300,20), 
     label=c("5  Feb 2019","25 Feb 2019",
             "17 Mar 2019","6  Apr 2019","26 Apr 2013"),
     side=1,mgp = c(2.8, 0.4, 0))

mtext("SPI/-STI",
      side = 2, line = 2.2, cex.lab = 1, 
      las = 0,cex=0.8)

par(new = TRUE)
plot(indx,CDHI[x,y,indx], type = "o", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", pch=17, col="purple2",lwd=2,xaxs = "i",yaxs = "i",
     ylim=c(0.2,1.5))

axis(at=seq(0.6,1,0.1),
     label=as.character(seq(0.6,1,0.1)),
     side = 4, mgp=c(1,0.8,0),
     col.axis="black",cex.axis=0.8,las=2)

mtext("CDHId" ,
      side = 4, line = 2.2, cex.lab = 1,
      las = 0, cex=0.8)

legend(x="top",
       c("SPI","STI","CDHId"),
       lty=c(1,1,1),pch=c(NA,NA,17), lwd=c(2,2,2),
       col=c("deepskyblue","brown1","purple2"),horiz=T,
       bty="n",cex=1.5, x.intersp=c(0.2,0.2,0.2), inset=c(0.01,0) )

##05End##
  
# 06 Process plotting in Yunnan(2) -------------       
x=63
y=59
testr=r
testr[x,y]=10000
plot(testr)

indx = c(21250:21350)  #20130107-20130328
dimnames(TX)[[3]][indx]

hstime=indx
hstime[which(FH[x,y,indx]!=1)]=NA
dstime=indx
dstime[which(FD[x,y,indx]!=1)]=NA

dev.new(width=1300,height=800,res=72)
par(mar=c(3,4.5,0.5,4.5))

plot(indx,SPIarray[x,y,indx], type = "l", xaxt = "n", yaxt = "n",pch=20,
     ylab = "", xlab = "Date",col="deepskyblue", ylim=c(-4,4),
     cex=1.2,lwd=2,mgp=c(1.5,0.8,0),
     xaxs = "i",yaxs = "i")

lines(indx,TSI[x,y,indx], type = "l", pch=20, col="brown1",lwd=2,cex=1.2)
abline(h=-1,lty=2)
abline(h=1,lty=2)

points(dstime,SPIarray[x,y,dstime],pch=20,col="deepskyblue",cex=1.5)
points(hstime,TSI[x,y,hstime],pch=18,col="brown1",cex=1.5)


axis(at=seq(-4,4,1)
     ,label=as.character(seq(-4,4,1)),
     side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
     las=2)
axis(at= seq(21250,21350,25), 
     label=c("7  Mar 2019","1  Apr 2019",
             "26 Apr 2019","21 MAy 2019","15 June 2013"),
     side=1,mgp = c(2.8, 0.4, 0))

mtext("SPI/-STI",
      side = 2, line = 2.2, cex.lab = 1, 
      las = 0,cex=0.8)


par(new = TRUE)
plot(indx,CDHI[x,y,indx], type = "o", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", pch=17, col="purple2",lwd=2,xaxs = "i",yaxs = "i",
     ylim=c(0,1.6))

axis(at=seq(0.6,1,0.1),
     label=as.character(seq(0.6,1,0.1)),
     side = 4, mgp=c(1,0.8,0),
     col.axis="black",cex.axis=0.8,las=2)

mtext("CDHId" ,
      side = 4, line = 2.2, cex.lab = 1,
      las = 0, cex=0.8)

legend(x="top",
       c("SPI","STI","CDHId"),
       lty=c(1,1,1),pch=c(NA,NA,17), lwd=c(2,2,2),
       col=c("deepskyblue","brown1","purple2"),horiz=T,
       bty="n",cex=1.5, x.intersp=c(0.2,0.2,0.2), inset=c(0.01,0) )        

##06End##
    
  