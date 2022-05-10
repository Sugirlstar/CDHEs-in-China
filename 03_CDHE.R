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
#    The formal figure generating code please refer to the folder "figure_generating"


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
  # name with the parameters
  # format: compound_spithr1_hwthr_da
  ctype = "compound_-1-1-10days"
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

# 02 Correlation calculating -------------
 # * in script "/figure_generating/figure_02.R"
###

##02End##


# 03 Copula fitting -------------
vindex = dimnames(TX)[[3]]
mds=levels(as.factor(substring(vindex,5,8)))
flowwindows=rep(-15:15,each=LL)
iflagarray = array(dim=dim(SPIarray)) 
for(i in 1:days)
  iflagarray[,,i] = (SPIarray[,,i] <= (spithr1) ) & (-TSI[,,i] <= (hwthr) )

coptype = array(dim=c(xx,yy,length(mds)))
pCDHdi = array(dim=dim(SPIarray))
qCDHdi = array(dim=dim(SPIarray))
Q = array(dim=c(xx,yy,length(mds)))

pb = tkProgressBar(title="Progress",label="Completed %", 
                   min=0, max=100, initial = 0, width = 300)
for (x in 1:xx)
{
  for (y in 1:yy)
    if (TRUE %in% (SPIarray[x, y, ] > -9999))
    {
      xs=c(x,ifelse( x+1<xx & x+1>0 , x+1 , NA),ifelse( x-1<xx & x-1>0 , x-1 , NA)   )
      xs=subset(xs,xs!="NA")
      ys=c(y,ifelse( y+1<yy & y+1>0 , y+1 , NA),ifelse( y-1<yy & y-1>0 , y-1 , NA)   )
      ys=subset(ys,ys!="NA")
      
      for(i in 1:length(mds))
      {
        nx=which(substring(dimnames(TX)[[3]],5,8)==mds[i]) 
        nxs=nx+flowwindows 
        nxs=subset(nxs,nxs<=days & nxs>0) 
        valueindex=nxs 
        
        pSPI=pnorm( as.vector(SPIarray[xs,ys,valueindex]) ) 
        pTSI=pnorm( as.vector(-TSI[xs,ys,valueindex]) ) 
        iflags = as.vector( iflagarray[xs,ys,valueindex] )
        
        iflag = which( iflags == 1  ) 
        if(length(iflag) < 5)
          next
        
        q0 = (length(iflag))/length(iflags) 
        Q[x,y,i] = q0 
        
        valueind = which(iflagarray[x,y,nx] == 1)
        
        if( length( valueind ) > 0 ) 
        {
          finalIndex = match( SPIarray[x,y,nx[valueind]] ,
                              (as.vector(SPIarray[xs,ys,valueindex]))[iflag] )
          a1 = pSPI[iflag]/pnorm(spithr1)
          a2 = pTSI[iflag]/pnorm(spithr1)
          arx = cbind( a1 , a2 ) 
          
          selectedCopula = BiCopSelect(a1,a2,familyset=c(3,4,5),rotations=F) #3,4,5 : "Clayton","Gumbel","Frank"
          coptype[x,y,i]=selectedCopula[[1]] 
          
          frankCopula = BiCopSelect(a1,a2,familyset=5,rotations=F) 
          parm = frankCopula$par 
          mycopula = frankCopula(parm, dim = 2) 
          
          ybpcdhi = pCopula(arx, mycopula) 
          pcdhi = ybpcdhi[finalIndex]
          
          qCDHdi[x,y,nx[valueind]] = pcdhi*q0
          pCDHdi[x,y,nx[valueind]] = pcdhi
        }
      }
    }
   
  
  info = sprintf("Completed %d%%", round(x*100/xx))  
  setTkProgressBar(pb, value = x*100/xx, 
                   title = sprintf("Progress (%s)",info),label = info) 
}
close(pb)

#y=1-x
CDHdi = 1 - pCDHdi
##03End##

# 04 CDHE identify (daily) -------------
{
  ## 04-01 first judge ----
  Flag <- array(NA,dim = c(xx,yy,days)) # set judge flag
  for(i in st:days)
    Flag[,,i] = ( CDHdi[,,i] >= 0 )
  ##04-01End##
  
  ## 04-02 identify consecutive 1 for at least da days ----
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
          if( (Flag[x,y,i] == 1 & Flag[x,y,i-1] != 1) )
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
  
  Flag2=Flag
  ##04-02End##
  
  ## 04-03 Counting----
  Flag2[ which( is.na(Flag2) == TRUE ) ] = -9999
  Flag3= array(dim = c(xx, yy, dim(Flag2)[3]))
  for (x in 1:xx)
    for (y in 1:yy)
      if (TRUE %in% (Flag2[x, y,] >= -999))
        for(i in days:st)
          if( Flag2[x,y,i] == 1 & Flag2[x,y,i-1] != 1 )
            Flag3[x,y,i] = 1
  
  Flag[which(Flag == -9999)] = NA
  Flag2[which(Flag2 == -9999)] = NA
  ##04-03End##
  
}
##04End##


# 05 Metrics calculating -------------
{
  ## 05-01 initialization ----
  CDHI=CDHdi*Flag2 
  CDHI[which(CDHI==0)]=NA 
  Hy = substring(dimnames(TX)[[3]], 1, 4)
  fhy = factor(Hy) 
  CDHfre = array(dim = c(xx, yy, length(levels(fhy)))) 
  CDHdur = array(dim = c(xx, yy, length(levels(fhy))))
  CDHstg = array(dim = c(xx, yy, length(levels(fhy)))) 
  CDHmag = array(dim = c(xx, yy, length(levels(fhy)))) 
  dimnames(CDHfre)[[3]] = levels(fhy)
  dimnames(CDHdur)[[3]] = levels(fhy) 
  dimnames(CDHstg)[[3]] = levels(fhy) 
  dimnames(CDHmag)[[3]] = levels(fhy) 
  ##05-01End##
  
  ## 05-02 annual value at each grid calculating ----
  for (x in 1:xx)
    for (y in 1:yy)
      if (TRUE %in% (SPIarray[x, y,] >= -999))
      {
        CDHfre[x, y,] = tapply(Flag3[x, y,], fhy, sum, na.rm = TRUE) 
        CDHdur[x, y,] = tapply(Flag2[x, y,], fhy, sum, na.rm = TRUE) 
        CDHstg[x, y,] = tapply(CDHI[x, y,], fhy, mean, na.rm = TRUE) 
        CDHmag[x, y,] = tapply(CDHI[x, y,], fhy, sum, na.rm = TRUE) 
      }
  ##04-02End##
  
  ## 05-03 multi-year mean value at each grid calculating ----
  CDHFRE = array(dim = c(xx, yy))
  CDHDUR = array(dim = c(xx, yy))
  CDHSTG = array(dim = c(xx, yy))
  CDHMAG = array(dim = c(xx, yy))
  for (x in 1:xx)
    for (y in 1:yy)
    {
      CDHFRE[x, y] = mean(CDHfre[x, y,], na.rm = TRUE)
      CDHDUR[x, y] = mean(CDHdur[x, y,], na.rm = TRUE)
      CDHSTG[x, y] = mean(CDHstg[x, y,], na.rm = TRUE)
      CDHMAG[x, y] = mean(CDHmag[x, y,], na.rm = TRUE)
    }
  
  # convert to raster to calculate area weight
  CDHFREr = raster(CDHFRE)
  CDHDURr = raster(CDHDUR)
  CDHSTGr = raster(CDHSTG)
  CDHMAGr = raster(CDHMAG)
  extent(CDHFREr) = c(72, 136, 18, 54) 
  extent(CDHDURr) = c(72, 136, 18, 54) 
  extent(CDHSTGr) = c(72, 136, 18, 54) 
  extent(CDHMAGr) = c(72, 136, 18, 54) 
  crs(CDHFREr) = crs(r) 
  crs(CDHDURr) = crs(r)
  crs(CDHSTGr) = crs(r) 
  crs(CDHMAGr) = crs(r) 
  ##05-03End##
  
  ## 05-04 calculate national multi year average (consider area weight)----
  r1 = raster(CDHfre[, , 1])  
  extent(r1) = c(72, 136, 18, 54)
  crs(r1) = crs(r)
  w = area(r1, weights = TRUE, na.rm = TRUE) 
  w = as.matrix(w) 
  
  CDHFREh = array(dim = c(LL))
  CDHDURh = array(dim = c(LL))
  CDHSTGh = array(dim = c(LL))
  CDHMAGh = array(dim = c(LL))
  
  for (i in 1:LL)
  {
    y1 = CDHfre[, , i] * w 
    CDHFREh[i] = sum(y1, na.rm = T) 
    y2 = CDHdur[,,i] * w
    CDHDURh[i] = sum(y2, na.rm = T)
    y3 = CDHstg[,,i] * w
    CDHSTGh[i] = sum(y3, na.rm = T)
    y4 = CDHmag[,,i] * w
    CDHMAGh[i] = sum(y4, na.rm = T)
  }
  # calculate the coverage
  CDH1=CDHfre[,,1]*0+1 #每个格点定为1
  CDHper=NULL
  for( i in 1:length(CDHfre[1,1,]) )
  {
    pt=which(CDHfre[,,i]>0) #有事件发生的点位置
    CDHper=c(CDHper,sum(CDH1[pt]*w[pt],na.rm=TRUE)) #发生点位*面积权重
  }
  names(CDHper)=levels(fhy)
  ##05-04End##
  
}
save.image(paste0(pathname,"/Metrics_CDHE.RData"))
##05End##


# the rest... ----------------
# the rest is the same as in 01_drought

