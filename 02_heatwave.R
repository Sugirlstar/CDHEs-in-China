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
# 4. To reproduce, just change the "pathname" (line 23) to where you store the file we offered.
#
# 5. Step05 and Step06 are not the main code, 
#    but just for preview the result.
#    The formal figure generating code please refer to the folder "EF_figure_generating"


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

# load necessary workspace and set file location
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
  days = length(TX[1,1,])
  st = 2 # the start day of daily STI (same as scale)
  da = 5 # the minimum duration
  da2 = 3 # the maximum interval
  thr1 = -1 # first threshold
  thr2 = 0 # second threshold for merging
  s=1961
  e=2020
  LL=e-s+1
  # coordinates showing on figures
  xL = seq(72, 136, 5)
  xlabel = paste(xL, "°", sep = "") #设置坐标轴
  yL = seq(18, 54, 5)
  ylabel = paste(yL, "°", sep = "")
  # name with the parameters
  # format: extreme_index_thr1_da_thr2_da2
  ctype = "heatwave--1-5days-0-3days"
}

##00End##


# 01 Functions praparing -------------

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

# 02 STI calculating -------------
TSI=array(dim=c(dim(TX)))
dates=substring(dimnames(TX)[[3]],5,8)
fdates=levels(factor(dates)) 
nday = 7 # a 15-days window

pb = tkProgressBar(title="Progress",label="Completed %", 
                   min=0, max=100, initial = 0, width = 300)
for(x in 1:xx)
  for(y in 1:yy)
    if (TRUE %in% (TX[x, y, ] >= -999))
      for(i in 1:366)
      {
        today = fdates[i] 
        f <- which(dates == today) 
        ys=TX[x,y,f]
        tdays = f + rep(-nday:nday,each=LL)
        ff=subset(tdays,tdays<=days & tdays>0) 
        ymean=mean(TX[x,y,ff],na.rm=TRUE)
        ysd=sd(TX[x,y,ff],na.rm=TRUE)
        TSI[x,y,f]=(ys-ymean)/ysd
        
        info = sprintf("Completed %d%%", round(x*100/xx))  
        setTkProgressBar(pb, value = x*100/xx, 
                         title = sprintf("Progress (%s)",info),label = info) 
      }
close(pb)

save(TSI, file = "STI.RData")

##02End##

# 03 Heatwave identify (daily) -------------
{
  ## 03-01 judge using thr1 ----
  Flag <- array(NA,dim = c(xx,yy,days)) 
  for(i in st:days)
    Flag[,,i] = ( -TSI[,,i] <= thr1 ) 
  ##03-01End##
  
  ## 03-02 identify consecutive 1 for at least da days ----
  Flag[which(is.na(Flag)==TRUE)] = -9999 #空值置为-9999，最后调回空值
  pb = tkProgressBar(title="Progress",label="Completed %", 
                     min=0, max=100, initial = 0, width = 300)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if (TRUE %in% (Flag[x,y,] >= -999))
      {
        for( i in st:days ) #st=2
          if(Flag[x,y,i] == 1 & Flag[x,y,i-1] != 1)
          {
            duration=sumfun(Flag[x,y,],i)
            if(duration < da)
              Flag[x,y,i:(i+duration-1)] <- rep(0,duration) #注意(i+duration-1)括号别漏
          }
      }
    info = sprintf("Completed %d%%", round(x*100/xx))  
    setTkProgressBar(pb, value = x*100/xx, 
                     title = sprintf("Progress (%s)",info),label = info) 
  }
  close(pb)
  ##03-02End##
  
  ## 03-03 Merging ----
  Flag2=Flag
  pb = tkProgressBar(title="Progress",label="Completed %", 
                     min=0, max=100, initial = 0, width = 300)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if(TRUE %in% (Flag2[x,y,] >= -999))
      {
        Flag2[x,y,st-1]=0
        for( i in st:length(Flag2[x,y,]) )
          if( Flag2[x,y,i] != 1 & Flag2[x,y,i-1] == 1)
          {
            duration=sumfun2(Flag2[x,y,],i)
            if( duration <= da2 & all( -TSI[x,y,i:(i+duration-1)] <= thr2) )
              Flag2[x,y,i:(i+duration-1)] = rep(1,duration)
          }
      }
    
    info = sprintf("Completed %d%%", round(x*100/xx))  
    setTkProgressBar(pb, value = x*100/xx, 
                     title = sprintf("Progress (%s)",info),label = info) 
  }
  close(pb)
  ##03-03End##
  
  ## 03-04 Counting----
  # set Flag3 as the times, i.e., each heatwave event only the first day is 1, others 0
  Flag3= array(dim = c(xx, yy, dim(Flag2)[3]))
  for (x in 1:xx)
    for (y in 1:yy)
      if (TRUE %in% (Flag2[x, y,] >= -999))
        for(i in days:st)
          if( Flag2[x,y,i] == 1 & Flag2[x,y,i-1] != 1 )
            Flag3[x,y,i] = 1
  
  Flag2[which(Flag2 == -9999)] = NA
  Flag[which(Flag == -9999)] = NA
  ##03-04End##
}
##03End##

# 04 Metrics calculating -------------
{
  ## 04-01 initialization ----
  SPId=TSI*Flag2
  SPId[which(SPId==0)]=NA 
  Hy = substring(dimnames(TX)[[3]], 1, 4)
  fhy = factor(Hy) 
  DRfre = array(dim = c(xx, yy, length(levels(fhy)))) 
  DRdur = array(dim = c(xx, yy, length(levels(fhy)))) 
  DRstg = array(dim = c(xx, yy, length(levels(fhy))))
  dimnames(DRfre)[[3]] = levels(fhy)
  dimnames(DRdur)[[3]] = levels(fhy) 
  dimnames(DRstg)[[3]] = levels(fhy) 
  ##04-01End##
  
  ## 04-02 annual value at each grid calculating ----
  for (x in 1:xx)
    for (y in 1:yy)
      if (TRUE %in% (Flag[x, y,] >= -999))
      {
        DRfre[x, y,] = tapply(Flag3[x, y,], fhy, sum, na.rm = TRUE) 
        DRdur[x, y,] = tapply(Flag2[x, y,], fhy, sum, na.rm = TRUE)
        DRstg[x, y,] = tapply(SPId[x, y,], fhy, mean, na.rm = TRUE)
      }
  ##04-02End##
  ## 04-03 multi-year mean value at each grid calculating ----
  DRFRE = array(dim = c(xx, yy))
  DRDUR = array(dim = c(xx, yy))
  DRSTG = array(dim = c(xx, yy))
  for (x in 1:xx)
    for (y in 1:yy)
    {
      DRFRE[x, y] = mean(DRfre[x, y,], na.rm = TRUE)
      DRDUR[x, y] = mean(DRdur[x, y,], na.rm = TRUE)
      DRSTG[x, y] = mean(DRstg[x, y,], na.rm = TRUE)
    }
  # convert to raster to calculate area weight
  DRFREr = raster(DRFRE)
  DRDURr = raster(DRDUR)
  DRSTGr = raster(DRSTG)
  extent(DRFREr) = c(72, 136, 18, 54) 
  extent(DRDURr) = c(72, 136, 18, 54) 
  extent(DRSTGr) = c(72, 136, 18, 54)
  crs(DRFREr) = crs(r) 
  crs(DRDURr) = crs(r) 
  crs(DRSTGr) = crs(r) 
  ##04-03End##
  
  ## 04-04 calculate national multi year average (consider area weight)----
  r1 = raster(DRfre[, , 1]) 
  extent(r1) = c(72, 136, 18, 54)
  crs(r1) = crs(r)
  w = area(r1, weights = TRUE, na.rm = TRUE) 
  w = as.matrix(w) 
  DRFREh = array(dim = c(LL))
  DRDURh = array(dim = c(LL))
  DRSTGh = array(dim = c(LL))
  for (i in 1:LL)
  {
    y1 = DRfre[, , i] * w  
    DRFREh[i] = sum(y1, na.rm = T)  
    y2 = DRdur[,,i] * w
    DRDURh[i] = sum(y2, na.rm = T)
    y3 = DRstg[,,i] * w
    DRSTGh[i] = sum(y3, na.rm = T)
  }
  # calculate the coverage
  DR1=DRfre[,,1]*0+1 #每个格点定为1
  DRper=NULL
  for( i in 1:length(DRfre[1,1,]) )
  {
    pt=which(DRfre[,,i]>0) #有事件发生的点位置
    DRper=c(DRper,sum(DR1[pt]*w[pt],na.rm=TRUE)) #发生点位*面积权重
  }
  names(DRper)=levels(fhy)
  ##04-04End##
}
save.image(paste0(pathname,"/Metrics_heatwave.RData"))
##04End##

# the rest... ----------------
# the rest is the same as in 01_drought

