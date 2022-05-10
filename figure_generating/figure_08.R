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

# fundamental parameters
{
  s=1961
  e=2020
}
##00End##

# 01 Functions praparing -------------
tendency=function(v) 
{
  long=length(v)
  fit=lm(v~c(1961:(1961+long-1)))
  ten=fit$coefficients[2]
  intercept=fit$coefficients[1]
  pvalue=summary(fit)$coefficients[,4][2]
  tp=list(ten,pvalue,intercept)
  return(tp)
}

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

figure3out = function(objx,ylabname,aty,spoint) 
{
  atylabel=as.character(aty)
  ylimget = c(min(aty,na.rm=T),max(aty,na.rm=T))
  
  rst=bp1(s,e,objx)

  if(rst[1]!=0) #
  {
    FREt=c(s,rst[1,(1:(ncol(rst)-1))],e)
    FREy=FREt
    FREy[1]=FREt[1]*rst[2,1]+rst[3,1]
    for(i in 2: length(FREt) )
      FREy[i]=FREt[i]*rst[2,i-1]+rst[3,i-1]
    FREp1=rst[4,1] ###
    FREp2=rst[4,2] ###
  }

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
  
  if(  !(is.na(spoint[1]) == TRUE) )
    points(c(1961:2020),spoint,pch=".")
    
}

##01End###


#02 read data -------------
ctype = "GLOTI"
df = read.table(paste0(pathname,"supplement_materials/GLOTI_1961_2020.csv"), header=F, sep=",")
GLOTI6120 = df[,2]
#02
ctype = "PDO"
df = read.table(paste0(pathname,"supplement_materials/PDO61_20.csv"), header=F, sep=",")
PDO6120 = df[,2]

##02End###


#03 Plotting -------------
objx=PDO6120 
ydur = max(objx,na.rm=T) - min(objx,na.rm=T)
getylim1 = c(min(objx,na.rm=T), max(objx,na.rm=T)+ydur*0.2)
getylim1;ydur/4

s.res = pettitt.test(GLOTI6120)
print(s.res)
n = s.res$nobs
i = s.res$estimate
s.1 = mean(GLOTI6120[1:i])
s.2 = mean(GLOTI6120[(i+1):n])
spoint = c(rep(s.1,i), rep(s.2,(n-i)))

# Manually save as .eps
dev.new(width=1300,height=3000,res=72) 
par( mfrow = c(2,1),mar=c(4,4.8,0.5,1.2) )   
figure3out(PDO6120,"PDO",seq(-2.0,1.5,0.5),spoint=NA)
figure3out(GLOTI6120,"Temperature anomaly",seq(-0.3,1.2,0.3),spoint)        

##03End###

