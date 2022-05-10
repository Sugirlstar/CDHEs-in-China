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

library(sp)
library(maptools)
library(rgdal)
library(raster)
library(tcltk)
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
setwd(pathname)

# read shp files
Chinasp = readShapePoly(paste0(pathname,"shpfile/china/China")) # not given in the repository
provincesp = readShapePoly(paste0(pathname,"shpfile/china/PROVINCE_region")) # not given in the repository

# fundamental parameters
{
  xx = 72
  yy = 128
  days = length(P[1,1,])
  spithr1 = -1 
  hwthr = -1 
  da=10 #最低连续时间
  st=90 
  s=1961
  e=2020
  LL=e-s+1
  # coordinates showing on figures
  xL = seq(72, 136, 5)
  xlabel = paste(xL, "°", sep = "") #设置坐标轴
  yL = seq(18, 54, 5)
  ylabel = paste(yL, "°", sep = "")
}
##00End##

# 01 Functions praparing -------------

# 01-01 
fenqi=function(v,a) 
{
  v[is.na(v)==TRUE]=-999
  z=list()
  for(j in 1:length(a) )
  {
    i=j
    inx=a[i]
    b=NULL
    while( v[inx]==1 )
    {
      if(inx == length(v))
      {
        b=c(b,inx)
        break
      }
      b=c(b,inx)
      inx=inx+1
    }
    z[[j]]=b
  }
  return(z)
}

# 01-02 Counting the length of consecutive 1
sumfun = function(v,i) #给参:v为某格点逐日spi,起算位置i
{
  a=0
  while(v[i]==1)
  {
    if( i == length(v) ) #放在i=i+1前面
    {a=a+1; break}
    a=a+1
    i=i+1
  }
  return(a)
}

# 01-03 Counting the length of consecutive 0
sumfun2 = function(v,i) #给参:v为某格点逐日spi,起算位置i
{
  a=0
  while(v[i]==0)
  {
    if( i == length(v) )
    {a=a+1; break}
    a=a+1
    i=i+1
  }
  return(a)
}

# 01-04 Linear regression
tendency = function(v) #默认:由1961开始
{
  long=length(v)
  fit=lm(v~c(1961:(1961+long-1)))
  ten=fit$coefficients[2]
  intercept=fit$coefficients[1]
  pvalue=summary(fit)$coefficients[,4][2]
  tp=list(ten,pvalue,intercept)
  return(tp)
}

# 01-05 FLMIP function (only one breakpoint)
bp1 = function(s,e,STGh)
{
  # check whether adjacent slopes have opposite symbols
  muti=function(v)  
  {
    x=length(v)-1
    z=0
    for( i in 1:x )
    { p=v[i]*v[i+1]
    if(p>0)
      break
    z=z+1
    }
    return(x==z)
  }
  
  m=e-s+1 
  Y=matrix(STGh,m,1)  
  T=s:e  
  
  SSR=10000000
  RST=0
  
  bp=(s+9):(e-9)   
  b=bp-s+1
  for( j in 1:length(bp) )  
  {
    A0=matrix(0,m,2)  
    A0[1: b[j],1]=s:bp[j]
    A0[ (b[j]+1):m ,1 ]=bp[j]
    A0[ (b[j]+1):m, 2 ]=1:(e-bp[j])
    c=matrix(rep(1,m),m,1)  
    A=cbind(A0,c)
    S=solve(t(A)%*%A)%*%t(A)%*%Y
    
    if(muti(S[1:2]))   
    {
      rst=matrix(0,4,2)  
      rownames(rst)=c("year","a","c","p of a") # y=a*x+c
      rst[1,1]=bp[j]
      rst[2,]=S[1: 2]
      rst[3,1]=S[3]
      rst[3,2]=rst[3,1]+(S[1]-S[2])*bp[j]
      
      # Significance evaluate
      
      rdf1=rst[1,1]-s+1-2 
      rdf2=e-rst[1,1]+1-2
      
      ssr1=0 
      for(u in 1:b[j])
        ssr1=ssr1+(Y[u]-(rst[2,1]*T[u]+rst[3,1]))^2
      
      ssr2=0 
      for( u in (b[j]+1):m )
        ssr2=ssr2+(Y[u]-(rst[2,2]*T[u]+rst[3,2]))^2  
      
      vv1=sqrt(ssr1/rdf1) #Residual standard error
      vv2=sqrt(ssr2/rdf2) #Residual standard error
      
      stderr1=vv1/sqrt(sum( (c(s:rst[1,1])-mean(c(s:rst[1,1])))^2))
      stderr2=vv2/sqrt(sum( ( c((rst[1,1]+1):e)-mean(c((rst[1,1]+1):e)) )^2))
      
      tval1=rst[2,1]/stderr1
      tval2=rst[2,2]/stderr2
      
      pr1=2*pt(abs(tval1),rdf1,lower.tail=FALSE)
      pr2=2*pt(abs(tval2),rdf2,lower.tail=FALSE)
      
      rst[4,1]=pr1
      rst[4,2]=pr2
      
      ssr=ssr1+ssr2
      
      if( pr1>0.05 & pr2>0.05 )
        ssr=1000000000000  
      
      if(ssr<SSR)
      { SSR=ssr
      RST=rst }
    }
  }
  return(RST)
}

# 01-06 M-K test in each grid
MK.raster = function(xraster, type="year", month0=1)
{
  library(Kendall)
  library(raster)
  x = as.array(xraster) 
  year0=1961
  D = dim(x)
  MK.xraster = array(data=NA,dim=c(D[1],D[2],3)) 
  if (type == "year"){
    for (i in 1:D[1])
      for (j in 1:D[2])
        if (TRUE %in% (x[i, j, ] >= -9999))
        {
          if( length(which(x[i,j,]>-9999))>2 ) # require at least 3 values
          {
            xts = ts(x[i,j,],start=year0,frequency=1)
            z = MannKendall(xts)
            MK.xraster[i,j,1:2] = c(z$tau,z$sl)
          }else
            MK.xraster[i,j,1:2]=NA 
        }
  }
  return(MK.xraster)
}

# 01-07 Significance identifier classification (for figure)
classtype = function(v)
{
  if(is.na(v)==TRUE)
    ct=NA else
      if(v<=0.05)  
        ct=4 else 
          ct=NA
        return(ct)
}

# 01-08 Relocating the color bar
cf = function(colr,whitesite="white",z,z0=0) 
{
  # colr: the values of color bar
  # whitesite: the color number (or name) you wish to settle,
  #            must contain in the "colr"
  # z: the raster or array you wish to plot
  # z0: the value you wish to aligned with "whitesite"
  # return to the color bar where "whitesite" is aligned with z0
  
  zz=as.matrix(z)
  z1=min(zz,na.rm=TRUE) 
  z2=max(zz,na.rm=TRUE) 
  zL1=(z0-z1)/(z2-z1) 
  zL2=(z2-z0)/(z2-z1) 
  cL=length(colr)
  if(whitesite=="mid")
    c0=round(cL/2)  else
      c0=which(colr==whitesite) 
  cL1=(c0-1)/cL 
  cL2=(cL-c0+1)/cL 
  
  if(z0<z1) 
  {
    x=round((z1-z0)/(z2-z0)*(cL-c0)+c0)
    colr_result=colr[x:cL]
  }else
    if(z0>z2)
    {
      x=round((z2-z1)/(z0-z1)*c0)
      colr_result=colr[1:x]
    } else
      if(zL1==0)
        colr_result=colr[c0:cL] else
          if(zL2==0)
            colr_result=colr[c0:cL]=colr[1:c0] else
              if(zL1>cL1)
              {
                x=round(c0/zL1) 
                colr_result=colr[1:x]
              }else
                if(zL1<cL1)
                {
                  x=round((zL2*cL+c0-cL)/zL2) 
                  colr_result=colr[x:cL]
                }else
                  colr_result=colr
  
  return(colr_result)        
}

##01End##

# 02 Calculate the joint probability -------------
vindex = dimnames(TX)[[3]]
mds=levels(as.factor(substring(vindex,5,8)))
flowwindows=rep(-15:15,each=LL)
iflagarray = array(dim=dim(SPIarray))
for(i in 1:days)
  iflagarray[,,i] = (SPIarray[,,i] <= (spithr1) ) & (-TSI[,,i] <= (hwthr) )

pCDHdi = array(dim=dim(SPIarray))

pb = tkProgressBar(title="Progress",label="Completed %", 
                   min=0, max=100, initial = 0, width = 300)
for (x in 1:xx)
{
  for (y in 1:yy)
    if (TRUE %in% (SPIarray[x, y, ] > -9999))
    {
      
      for(i in 1:length(mds))
      {
        nx=which(substring(dimnames(TX)[[3]],5,8)==mds[i]) 
        valueindex=nx
        nxs=nx+flowwindows 
        nxs=subset(nxs,nxs<=days & nxs>0) 
      
        pSPI=pnorm( c( spithr1, SPIarray[x,y,nxs] ) )  
        pTSI=pnorm( c( hwthr, -TSI[x,y,nxs] ) )
        
        iflags = as.vector( c(0, iflagarray[x,y,nxs]) ) 
        iflag = which( iflags == 1  ) 
        if(length(iflag) < 5)
          next
        
        valueind = which(iflagarray[x,y,nx] == 1) 
        
        if( length( valueind ) > 0 ) 
        {
  
          finalIndex = match( pnorm(SPIarray[x,y,nx[valueind]]), pSPI )
          a1 = pSPI
          a2 = pTSI
          arx = cbind( a1 , a2 ) 
          
          frankCopula = BiCopSelect(a1,a2,familyset=5,rotations=F) 
          parm = frankCopula$par 
          mycopula = frankCopula(parm, dim = 2) 
          
          ybpcdhi = pCopula(arx, mycopula) 
          
          pCDHdi[x,y,nx[valueind]] = ybpcdhi[finalIndex]/ybpcdhi[1]
          
          
        }
      }
      
    }
  info = sprintf("Completed %d%%", round(x*100/xx))  
  setTkProgressBar(pb, value = x*100/xx, 
                   title = sprintf("Progress (%s)",info),label = info) 
}
close(pb)

CDHdi = pCDHdi # the joint probability

##02End##


# 03a Clssification (jointP) -------------
x=30
y=50
longtt=71.75+0.5*y #96.75
latt=54.25-0.5*x #39.25
cbind(longtt,latt)

ind0= which( CDHdi[x,y,] < 0.3 & CDHdi[x,y,] >= 0.2 ) # -0.5
ind1= which( CDHdi[x,y,] < 0.2 & CDHdi[x,y,] >= 0.1 ) # -0.8
ind2= which( CDHdi[x,y,] < 0.1 & CDHdi[x,y,] >= 0.05 ) # -1.3
ind3= which( CDHdi[x,y,] < 0.05 & CDHdi[x,y,] >= 0.02 ) # -1.6
ind4= which( CDHdi[x,y,] < 0.02 ) # -2

indx = c(20331:20369)  #19791113-19800103 CDHE
indx5 = c(16828:16837) #20131129-20131211 CDHE

##03aEnd##

# 04a distribution plotting (jointP) -------------
colr=rev(colorRampPalette(brewer.pal(9,"Blues")[1:9])(12))
#03 Putout
dev.new(width=1700,height=1700,res=72) 

xxx=SPIarray[x,y,ind0]
yyy=-TSI[x,y,ind0]
zzz=CDHdi[x,y,ind0]
plot(xxx,yyy,type="p",col=colr[9],
     xlab="SPI",ylab="-STI",pch=16,mgp = c(1.7, 0.7, 0),
     xlim=c(-4,4),ylim=c(-4,4), cex=1.3,
     cex.axis=1.3, cex.lab=1.3,las=1)#

xxx=SPIarray[x,y,ind1]
yyy=-TSI[x,y,ind1]
zzz=CDHdi[x,y,ind1]
points(xxx,yyy,type="p",col=colr[7],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,ind2]
yyy=-TSI[x,y,ind2]
zzz=CDHdi[x,y,ind2]
points(xxx,yyy,type="p",col=colr[5],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,ind3]
yyy=-TSI[x,y,ind3]
zzz=CDHdi[x,y,ind3]
points(xxx,yyy,type="p",col=colr[3],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,ind4]
yyy=-TSI[x,y,ind4]
zzz=CDHdi[x,y,ind4]
points(xxx,yyy,type="p",col=colr[1],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,indx5]
yyy=-TSI[x,y,indx5]
zzz=CDHdi[x,y,indx5]
points(xxx,yyy,type="p",col="darkorange",xlab="SPI",ylab="-STI",pch=3,cex=1.5)

xxx=SPIarray[x,y,indx]
yyy=-TSI[x,y,indx]
zzz=CDHdi[x,y,indx]
points(xxx,yyy,type="p",col="darkorange",xlab="SPI",ylab="-STI",pch=4, cex=1.5)

text(2.1,3.7,"  previous index \n ")

points(2.1,3.2, pch=16,cex=1.8,col=colr[9])
text(2.7,3.2," <0.3")

points(2.1,2.9, pch=16,cex=1.8,col=colr[7])
text(2.7,2.9, "<0.2" )

points(2.1,2.6, pch=16,cex=1.8,col=colr[5])
text(2.7,2.6, "<0.1" )

points(2.1,2.3, pch=16,cex=1.8,col=colr[3])
text(2.7,2.3, "<0.05" )

points(2.1,2.0, pch=16,cex=1.8,col=colr[1])
text(2.7,2.0, "<0.02" )

abline(h=-1,lty=2)
abline(v=-1,lty=2)


points(1.8,-3.6, pch=4,cex=1.8,col="darkorange")
text(3.1,-3.6,"20160830-20161007",cex=1.1)

points(1.8,-4, pch=3,cex=1.8,col="darkorange")
text(3.1,-4,"20070127-20070205",cex=1.1)

##04aEnd##

# 05a process plotting (jointP) -------------
indx = c(20331:20369)  #20160830-20161007 CDHE
indx5 = c(16828:16837) #20070127-20070205 CDHE
{
  ## 05a-01 process1 ----
  dev.new(width=1300,height=800,res=72) 
  par(mar=c(3,4.8,0.5,3))
  plot(indx,SPIarray[x,y,indx], type = "o", xaxt = "n", yaxt = "n",pch=20,
       ylab = "", xlab = "Date",col="darkblue", ylim=c(-3.5,3.5),
       cex=0.8,lwd=2,mgp=c(1.5,0.8,0),
       xaxs = "i",yaxs = "i")
  
  lines(indx,-TSI[x,y,indx], type = "o", pch=20, col="brown",lwd=2)
  abline(h=-1,lty=2)
  
  axis(at=seq(-3.5,3.5,1)
       ,label=as.character(seq(-3.5,3.5,1)),
       side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
       las=2)
  axis(at=c(seq(20331,20369,10),20369), 
       label=c("27 Aug 2016","9 Sep 2016",
               "19 Sep 2016","29 Sep 2016","7 Oct 2016"),
       side=1,mgp = c(2.8, 0.4, 0))
  
  mtext("SPI/-STI",
        side = 2, line = 2.2, cex.lab = 1, 
        las = 0, col="darkblue",cex=0.8)
  
  legend(x="bottomright",
         c("SPI","-STI"),
         lty=c(1,1),pch=c(20,20), lwd=c(2,2),col=c("darkblue","brown"),
         text.col=c("darkblue","brown"),
         bty="n",cex=1.5, x.intersp=c(0.2,0.2), inset=c(0.01,0) )
  
  ##05-01End##
  
  ## 05a-02 process2 ----
  dev.new(width=1300,height=800,res=72) 
  par(mar=c(3,4.8,0.5,3))

  plot(indx5,SPIarray[x,y,indx5], type = "o", xaxt = "n", yaxt = "n",pch=20,
       ylab = "", xlab = "Date",col="darkblue", ylim=c(-3.5,3.5),
       cex=0.8,lwd=2,mgp=c(1.5,0.8,0),
       xaxs = "i",yaxs = "i")
  
  lines(indx5,-TSI[x,y,indx5], type = "o", pch=20, col="brown",lwd=2)
  abline(h=-1,lty=2)
  
  axis(at=seq(-3.5,3.5,1)
       ,label=as.character(seq(-3.5,3.5,1)),
       side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
       las=2)
  
  axis(at=seq(16828,16837,3), 
       label=c("27 Jan 2007","30 Jan 2007",
               "2 Feb 2007","5 Feb 2007"),
       side=1,mgp = c(2.8, 0.4, 0))
  
  mtext("SPI/-STI",
        side = 2, line = 2.2, cex.lab = 1, 
        las = 0, col="darkblue",cex=0.8)
  
  legend(x="bottomright",
         c("SPI","-STI"),
         lty=c(1,1),pch=c(20,20), lwd=c(2,2),col=c("darkblue","brown"),
         text.col=c("darkblue","brown"),
         bty="n",cex=1.5, x.intersp=c(0.2,0.2), inset=c(0.01,0) )
  
  ##05-02End##
  
}
##05aEnd##

# 03b Classification (CDHId) -------------
load(paste0(pathname,"/Metrics_CDHE.RData")) # note that the variable "CDHdi" will be covered

x = 30
y = 50
longtt=71.75+0.5*y 
latt=54.25-0.5*x 

#0.3 0.2 0.1 0.05 0.02
ind0=which( CDHdi[x,y,] > 0 & CDHdi[x,y,]<=0.3 )
ind1=which(CDHdi[x,y,]> 0.3 & CDHdi[x,y,]<=0.6 )
ind2=which(CDHdi[x,y,]> 0.6 & CDHdi[x,y,]<=0.9 )
ind3=which(CDHdi[x,y,]> 0.9 & CDHdi[x,y,]<=0.98 )
ind4=which(CDHdi[x,y,]> 0.98)
indxx = which( CDHdi[x,y,] > 0)
indx=indxx[which(Flag2[x,y,indxx]==1)] 
##03bEnd##


# 04b distribution plotting (CDHId) -------------
colr=rev(colorRampPalette(brewer.pal(9,"Blues")[1:9])(12))

dev.new(width=1700,height=1700,res=72)  
xxx=SPIarray[x,y,ind0]
yyy=-TSI[x,y,ind0]
zzz=CDHI[x,y,ind0]
plot(xxx,yyy,type="p",col=colr[9],
     xlab="SPI",ylab="-STI",pch=16,mgp = c(1.7, 0.7, 0),
     xlim=c(-4,4),ylim=c(-4,4), cex=1.3,
     cex.axis=1.3, cex.lab=1.3,las=1)#

xxx=SPIarray[x,y,ind1]
yyy=-TSI[x,y,ind1]
zzz=CDHI[x,y,ind1]
points(xxx,yyy,type="p",col=colr[7],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,ind2]
yyy=-TSI[x,y,ind2]
zzz=CDHI[x,y,ind2]
points(xxx,yyy,type="p",col=colr[5],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,ind3]
yyy=-TSI[x,y,ind3]
zzz=CDHI[x,y,ind3]
points(xxx,yyy,type="p",col=colr[3],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,ind4]
yyy=-TSI[x,y,ind4]
zzz=CDHI[x,y,ind4]
points(xxx,yyy,type="p",col=colr[1],xlab="SPI",ylab="-STI",pch=16)

xxx=SPIarray[x,y,indx]
yyy=-TSI[x,y,indx]
zzz=CDHI[x,y,indx]
points(xxx,yyy,type="p",col="orangered",xlab="SPI",ylab="-STI",pch=4,cex=1.5)


text(2.1,3.7,"  CDHId \n ")
points(2.1,3.2, pch=16,cex=1.8,col=colr[9])
text(2.7,3.2," >0")
points(2.1,2.9, pch=16,cex=1.8,col=colr[7])
text(2.7,2.9, ">0.3" )
points(2.1,2.6, pch=16,cex=1.8,col=colr[5])
text(2.7,2.6, ">0.6" )
points(2.1,2.3, pch=16,cex=1.8,col=colr[3])
text(2.7,2.3, ">0.9" )
points(2.1,2.0, pch=16,cex=1.8,col=colr[1])
text(2.7,2.0, ">0.98" )

abline(h=-1,lty=2)
abline(v=-1,lty=2)

points(1.8,-3.6, pch=4,cex=1.8,col="orangered")
text(3.1,-3.6,"19950526-19950604",cex=1.1)
      
##04bEnd##

# 05b process plotting (CDHId) -------------
dev.new(width=1300,height=800,res=72) 
par(mar=c(3,4.8,0.5,3))

plot(indx,SPIarray[x,y,indx], type = "o", xaxt = "n", yaxt = "n",pch=20,
     ylab = "", xlab = "Date",col="darkblue", ylim=c(-3.5,-0.5),
     cex=0.8,lwd=2,mgp=c(1.5,0.8,0),
     xaxs = "i",yaxs = "i")

lines(indx,-TSI[x,y,indx], type = "o", pch=20, col="brown",lwd=2)
abline(h=-1,lty=2)

axis(at=seq(-3.5,-0.5,1)
     ,label=as.character(seq(-3.5,-0.5,1)),
     side = 2,col.axis="black",cex.axis=0.8,mgp=c(1,0.8,0),
     las=2)
axis(at=seq(12564,12573,3), 
     label=c("26 May 1995","29 May 1995",
             "1 June 1995","4 June 1995"),
     side=1,mgp = c(2.8, 0.4, 0))

mtext("SPI/-STI",
      side = 2, line = 2.2, cex.lab = 1, 
      las = 0, col="black",cex=0.8)

legend(x="bottomright",
       c("SPI","-STI"),
       lty=c(1,1),pch=c(20,20), lwd=c(2,2),col=c("darkblue","brown"),
       text.col=c("darkblue","brown"),
       bty="n",cex=1.5, x.intersp=c(0.2,0.2), inset=c(0.01,0) )

##05bEnd##    
    
  
  

