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

# load necessary workspace and set location
load(paste0(pathname,"/P_China_1dy_0.5deg_1961_2020.RData")) # precipitation
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
  scale = 90 
  st = 90 # the start day of daily SPI (same as scale)
  da = 10 # the minimum duration
  da2 = 5 # the maximum interval
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
  ctype = "drought_spi90_-1_10days_0_5days" 
}

##00End##

  
# 01 Functions praparing -------------

# 01-01 Calculating daily SPI using accumulated precipitation
SPIh = function(p)  
{
  p.NA<-which(is.na(p))
  y<-p[which(! is.na(p))] 
  N<-length(y)

  x<-y[which(y>0)] 
  A<-log(mean(x))-sum(log(x))/length(x)
  a<-(1+sqrt(1+4*A/3))/(4*A) #estimation of alpha
  b<-mean(x)/a #estimation of beta
  
  g_x<-function(x) {(x^(a-1))*exp(-x/b)/(b^a * gamma(a))}
  G_x<-numeric(N)
  for(i in 1:N) 
    ifelse(y[i]>0, G_x[i]<-integrate(g_x,0,y[i])$value, G_x[i]<-0)
  q<-length(y[which(y==0)])/N
  H_x<-q+(1-q)*G_x
  
  c0<-2.515517
  c1<-0.802853
  c2<-0.010328
  d1<-1.432788
  d2<-0.189269
  d3<-0.001308
  
  Z=rep(NA,length(p))    
  zx=NULL  
  for(i in 1:length(H_x))  
    if( (H_x[i]>0)&(H_x[i]<=0.5) )
    {
      T<-sqrt(log(1/(H_x[i])^2))
      zx = c(zx, -T+(c0+c1*T+c2*T^2)/(1+d1*T+d2*T^2+d3*T^3) ) 
    } else if( (H_x[i]>0.5)&(H_x[i]<1) )
    {
      T<-sqrt(log(1/(1-H_x[i])^2))
      zx = c(zx, T-(c0+c1*T+c2*T^2)/(1+d1*T+d2*T^2+d3*T^3) )
    }
  Z[ which(is.na(p)==FALSE) ] = zx
  
  return(Z)
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



# 02 SPI calculating -------------
{
  ## 02-01 Calculate accumulated P (90days) ----
  acump=array(dim=c(xx,yy,days))
  pb = tkProgressBar(title="Progress",label="Completed %", 
                     min=0, max=100, initial = 0, width = 300)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if (TRUE %in% (P[x, y, ] >= 0))
        for(i in scale:days)
          acump[x,y,i]=sum(P[x,y,(i-scale+1):i],na.rm=TRUE)
    
    info = sprintf("Completed %d%%", round(x*100/xx))  
    setTkProgressBar(pb, value = x*100/xx, 
                     title = sprintf("Progress (%s)",info),label = info) 
  }
  close(pb)
  ##02-01End##
  
  ## 02-02 Calculate daily SPI ----
  SPIarray=array(dim=c(xx,yy,days))
  dates=substring(dimnames(P)[[3]],5,8) 
  fdates=levels(factor(dates)) 
  t0228 = fdates[59]
  f0228 = which(dates == t0228) 
  
  pb = tkProgressBar(title="Progress",label="Completed %", 
                     min=0, max=100, initial = 0, width = 300)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if (TRUE %in% (P[x, y, ] >= 0))
        for(i in 1:366)
        {
          today = fdates[i] 
          index = which(dates == today)
          
          if(i==60) # day0229
          {
            runnian=which( (f0228+1) %in% index ) 
            f0228[runnian]=index 
            index=f0228
          }
          
          SPIarray[x,y,index] = SPIh(acump[x,y,index]) 
          SPIarray[x,y,which(is.nan(SPIarray[x,y,]) == TRUE)] = NA  
          
        }
    info = sprintf("Completed %d%%", round(x*100/xx))  
    setTkProgressBar(pb, value = x*100/xx, 
                     title = sprintf("Progress (%s)",info),label = info) 
  }
  close(pb)
  ##02-02End##  
}

save(TSI, file = "SPI.RData")

##02End##

# 03 Drought identify (daily) -------------
{
  ## 03-01 judge using thr1 ----
  Flag <- array(NA,dim = c(xx,yy,days)) # set judge flag
  for(i in st:days)
    Flag[,,i] = ( SPIarray[,,i] <= thr1 ) 
  ##03-01End##
  
  ## 03-02 identify consecutive 1 for at least da days ----
  Flag[which(is.na(Flag)==TRUE)] = -9999 
  pb = tkProgressBar(title="Progress",label="Completed %", 
                     min=0, max=100, initial = 0, width = 300)
  for(x in 1:xx)
  {
    for(y in 1:yy)
      if (TRUE %in% (Flag[x,y,] >= -999))
      {
        Flag[ x,y, st-1 ]=0 
        for( i in st:days )
          if(Flag[x,y,i] == 1 & Flag[x,y,i-1] != 1)
          {
            duration=sumfun(Flag[x,y,],i)
            if(duration < da) 
              Flag[x,y,i:(i+duration-1)] <- rep(0,duration) 
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
        for( i in st:length(Flag2[x,y,]) )
          if( Flag2[x,y,i] != 1 & Flag2[x,y,i-1] == 1)
          {
            duration=sumfun2(Flag2[x,y,],i)
            if( duration <= da2 & all(SPIarray[x,y,i:(i+duration-1)] < thr2) )
              Flag2[x,y,i:(i+duration-1)] = rep(1,duration)
          }
    info = sprintf("Completed %d%%", round(x*100/xx))  
    setTkProgressBar(pb, value = x*100/xx, 
                     title = sprintf("Progress (%s)",info),label = info) 
  }
  close(pb)
  ##03-03End##
  
  ## 03-04 Counting----
  # set Flag3 as the times, i.e., each drought event only the first day is 1, others 0
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
  SPId=SPIarray*Flag2
  SPId[which(SPId==0)]=NA 
  Hy = substring(dimnames(TX)[[3]], 1, 4)
  fhy = factor(Hy) # convert to factors of year
  DRfre = array(dim = c(xx, yy, length(levels(fhy)))) # yearly frequency
  DRdur = array(dim = c(xx, yy, length(levels(fhy)))) # yearly duration
  DRstg = array(dim = c(xx, yy, length(levels(fhy)))) # yearly stength
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
save.image(paste0(pathname,"/Metrics_drought.RData"))
##04End##

# 05 Defining Class and Methods for plotting (just for preview) -------------
pbt = function(x,...) UseMethod("pbt")

#05-01 Spatial distribution (named with "SD")
pbt.SD = function(objx,ctype)
{
  colr = rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(100)) 
  vname = as.character(substitute(objx))
  
  png(
    file = paste("01SD_",ctype, "_", vname, ".png", sep = ""),
    width = 1800,
    height = 1400,
    res = 72 * 3
  )
  plot(
    objx,
    xaxt = "n",
    yaxt = "n",
    main = paste("SD_", ctype, "_", vname, sep = ""),
    col = colr
  )
  axis(1, xL, xlabel)
  axis(2, yL, ylabel)
  lines(Chinasp)
  lines(provincesp)
  dev.off()
}

#05-02 Temporal series (named with "TS")
pbt.TS = function(objx,ctype) 
{
  vname = as.character(substitute(objx))
  #02-01 
  rst=bp1(s,e,objx)
  #02-01-01-a
  if(rst[1]!=0) 
  {
    FREt=c(s,rst[1,(1:(ncol(rst)-1))],e)
    FREy=FREt
    FREy[1]=FREt[1]*rst[2,1]+rst[3,1]
    for(i in 2: length(FREt) )
      FREy[i]=FREt[i]*rst[2,i-1]+rst[3,i-1]
    FREp1=rst[4,1] 
    FREp2=rst[4,2]
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
  
  #02-01-02
  png(
    file = paste("02TS_",ctype, "_",vname,  ".png", sep = ""),
    width = 2500,
    height = 1800,
    res = 72 * 3
  )
  
  plot(c(s:e),objx,type="o",pch=20,xlab="",main= paste("TS_",ctype,"_",vname,sep=""))
    
      legend("bottomleft", 
             paste0("avg,max,min:", 
                    round(mean(objx,na.rm = TRUE),4), ",",
                    round(max(objx,na.rm = TRUE),4), "(", which.max(objx)+s-1,"),",
                    round(min(objx,na.rm = TRUE),4), "(", which.min(objx)+s-1,")" ) 
             ,bty="n" )
  
  #02-01-02-a
  if(rst[1]!=0)
  {
    lines(FREt, FREy, lwd=3,col="red")
    points(FREt[2],FREy[2],pch=17,cex=2)
      legend("bottom",as.character(rst[1,1]),bty="n")
      legend("topright", 
             paste0(round(rst[2,1],4),",",round(rst[2,2],4)),bty="n")
  }
  #02-01-02-b
  if(rst[1]==0)
  {
    lines(tendx,tendy,lwd=3,col="brown")
    text(1990,max(objx),
         paste("linear Fitting: y=",round(tend[[1]],5),
               "x+",round(tend[[3]],5), ",p-value=",round(tend[[2]],5),sep=""))
    text(1990,quantile(objx,0.97),
         paste("M-K Trend Test: tau=",round(mktau,5),
               ",p-value=",round(mkpvalue,5),sep=""))
  }
  dev.off()
}

#05-03 Spatial temporal distribution (named with "TSTD")
pbt.TSTD = function(objx,ck,ctype) 
{
  vname = as.character(substitute(objx)) 
  rst=bp1(s,e,ck)
  if(rst[1]!=0)
  {
    bp=rst[1]-s+1 
    tend1=tend2=array(dim=c(xx,yy))
    pfre1=pfre2=array(dim=c(xx,yy))
    for(x in 1:xx)
      for(y in 1:yy)
        if (TRUE %in% (objx[x,y,]>=0))
        {
          tend1[x,y]=tendency(objx[x,y,1:bp])[[1]]
          pfre1[x,y]=classtype(tendency(objx[x,y,1:bp])[[2]])
        }
    for(x in 1:xx)
      for(y in 1:yy)
        if (TRUE %in% (objx[x,y,]>=0))
        {
          tend2[x,y]=tendency(objx[x,y,(bp+1):(e-s+1)])[[1]]
          pfre2[x,y]=classtype(tendency(objx[x,y,(bp+1):(e-s+1)])[[2]])
        }

    xposfre1 = (tend1 > 0) &  (pfre1 > 0)
    xnegfre1 = (tend1 < 0) &  (pfre1 > 0) 
    xposfre2 = (tend2 > 0) &  (pfre2 > 0)
    xnegfre2 = (tend2 < 0) &  (pfre2 > 0) 

    xposp1 = round( sum( xposfre1*w , na.rm=T ) , 4 )
    xnegp1 = round( sum( xnegfre1*w , na.rm=T ) , 4 )
    xposp2 = round( sum( xposfre2*w , na.rm=T ) , 4 )
    xnegp2 = round( sum( xnegfre2*w , na.rm=T ) , 4 )

    posfre1 = (tend1 > 0) 
    negfre1 = (tend1 < 0) 
    posfre2 = (tend2 > 0)
    negfre2 = (tend2 < 0) 
    
    posp1 = round( sum( posfre1*w , na.rm=T ) , 4 )
    negp1 = round( sum( negfre1*w , na.rm=T ) , 4 )
    posp2 = round( sum( posfre2*w , na.rm=T ) , 4 )
    negp2 = round( sum( negfre2*w , na.rm=T ) , 4 )
    
    tend1r = raster(tend1) 
    extent(tend1r) = c(72, 136, 18, 54) 
    crs(tend1r) = crs(r) 
    
    tend2r = raster(tend2) 
    extent(tend2r) = c(72, 136, 18, 54) 
    crs(tend2r) = crs(r) 
    
    yind=xind=pfre1ind=pfre2ind=NULL
    for(x in 1:xx)
      for(y in 1:yy)
      {
        yind=c(yind,71.75+0.5*y) 
        xind=c(xind,54.25-0.5*x) 
        pfre1ind=c(pfre1ind,pfre1[x,y])
        pfre2ind=c(pfre2ind,pfre2[x,y])
      }
    
    colr1=rev(colorRampPalette(brewer.pal(7,"GnBu")[1:7])(100)) 
    colr2=colorRampPalette(brewer.pal(9,"OrRd")[1:9])(100)
    
    tend1col = cf(c(colr1,"white",colr2),whitesite="white",tend1r,z0=0) 
    tend2col = cf(c(colr1,"white",colr2),whitesite="white",tend2r,z0=0) 
    
    #tend1
    png(
      file = paste("03TSTD_tend1_",
                   rst[1],"_",ctype,"_",vname,
                   ",(+)",posp1,"(-)",negp1,
                   "(++)",xposp1,"(--)",
                   xnegp1,".png", sep = ""),
      width = 1800,
      height = 1400,
      res = 72 * 3
    )
    plot(
      tend1r,
      xaxt = "n",
      yaxt = "n",
      main = paste("TSTD_tend1_",
                   rst[1],"_",ctype,"_", vname, sep = ""),
      col=tend1col
    )
    axis(1, xL, xlabel)
    axis(2, yL, ylabel)
    lines(Chinasp)
    lines(provincesp)
    #
    points(yind,xind,pch=pfre1ind,cex=0.4) 
    points(126,19,pch=4,cex=1.3) 
    text(130,19.1,"α=0.05",cex=1.2) 
    dev.off()
    
    #tend2
    png(
      file = paste("03TSTD_tend2_",
                   rst[1],"_",ctype,"_", vname,
                   ",(+)",posp2,"(-)",negp2,
                   "(++)",xposp2,"(--)",
                   xnegp2, ".png", sep = ""),
      width = 1800,
      height = 1400,
      res = 72 * 3
    )
    plot(
      tend2r,
      xaxt = "n",
      yaxt = "n",
      main = paste("TSTD_tend2_",rst[1],
                   "_",ctype,"_", vname, sep = ""),
      col=tend2col
    )
    axis(1, xL, xlabel)
    axis(2, yL, ylabel)
    lines(Chinasp)
    lines(provincesp)
    #
    points(yind,xind,pch=pfre2ind,cex=0.4) 
    points(126,19,pch=4,cex=1.3) 
    text(130,19.1,"α=0.05",cex=1.2)
    dev.off()
  } else
    print("No break point")
}

#05-04 M-K spatial distriburion (named with "MKSD")
pbt.MKSD = function(objx,ctype)
{
  vname = as.character(substitute(objx)) 
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
 
  xposfre = (zz[,,1] > 0) &  (sigz > 0)
  xnegfre = (zz[,,1] < 0) &  (sigz > 0) 
  xposp = round( sum( xposfre*w , na.rm=T ) , 4 )
  xnegp = round( sum( xnegfre*w , na.rm=T ) , 4 )
  posfre = (zz[,,1] > 0) 
  negfre = (zz[,,1] < 0) 
  posp = round( sum( posfre*w , na.rm=T ) , 4 )
  negp = round( sum( negfre*w , na.rm=T ) , 4 )
  
  yind=xind=sigind=NULL
  for(x in 1:xx)
    for(y in 1:yy)
    {
      yind=c(yind,71.75+0.5*y)
      xind=c(xind,54.25-0.5*x)
      sigind=c(sigind,sigz[x,y])
    }

  colr=rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100)) 
  colr[50]="white"
  MKcol =  cf(colr,whitesite="white",z,z0=0) 
  
  png(
    file = paste("08MKSD_",ctype, "_",
                 vname,",(+)",posp,"(-)",negp,
                 "(++)",xposp,"(--)",
                 xnegp, ".png", sep = ""),
    width = 1800,
    height = 1400,
    res = 72 * 3
  )
  plot(
    z,
    xaxt = "n",
    yaxt = "n",
    main = paste("MKSD_", ctype, "_",vname ,sep = ""),
    col=MKcol 
  )
  axis(1, xL, xlabel)
  axis(2, yL, ylabel)
  lines(Chinasp)
  lines(provincesp)
  
  points(yind,xind,pch=sigind,cex=0.4)
  points(126,19,pch=4,cex=1.3)
  text(130,19.1,"α=0.05",cex=1.2) 

  dev.off()
}  
  
##05End##

# 06 Figure preview -------------
{
  # SD
  pbt.SD(DRFREr,ctype)
  pbt.SD(DRDURr,ctype)
  pbt.SD(DRSTGr,ctype)
  # TS
  pbt.TS(DRFREh,ctype)
  pbt.TS(DRDURh,ctype)
  pbt.TS(DRSTGh,ctype)
  pbt.TS(DRper,ctype)
  # TSTD
  pbt.TSTD(DRfre,DRFREh,ctype)
  pbt.TSTD(DRdur,DRDURh,ctype)
  pbt.TSTD(DRstg,DRSTGh,ctype)
  # MKSD
  pbt.MKSD(DRfre,ctype)
  pbt.MKSD(DRdur,ctype)
  pbt.MKSD(DRstg,ctype)
  
}

save.image(paste0(pathname,"/Drought_final.RData"))

##06End##


