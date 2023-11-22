library(raster)
library(ncdf4)

pathtofile = "H:/LHdaily"
file_list = list.files(pathtofile,pattern="\\.nc$",full.names=T)
file_list = sort(file_list) #order the file by names

#####
#transfer the nc to matrix and give the name
opnames = c("LH00","LH12")
for(i in 1:length(file_list))
{
  input_nc = file_list[i]
  ncfile = nc_open(input_nc)
  rstack = stack(input_nc, varname = names(ncfile$var)[1] )
  formattedD = as.vector(sapply(names(rstack),function(data) gsub("X|\\.","",data)))
  vals = values(rstack)
  arraytest = array(vals,dim=c(dim(rstack)[2],dim(rstack)[1],dim(rstack)[3]))
  arrayfinal = aperm(arraytest, c(2,1,3))
  dimnames(arrayfinal)[[3]] = substring(formattedD,1,8)
  assign(opnames[i],arrayfinal)
}
#读入：LH00, LH12
#####

#####
#read in the raw data, give a mask
pathname = "D:/hyj/CDHE/CDHEs_20230611"
setwd(pathname)
load(paste0(pathname, "/TX_China_1dy_0.5deg_1961_2020.RData")) 
{
  xx = 72
  yy = 128
  days = length(TX[1,1,])
  s=1961
  e=2020
  LL=e-s+1
  daynames = dimnames(TX)[[3]]
  
  r1 = raster(TX[, , 1]) 
  extent(r1) = c(72, 136, 18, 54)
  crs(r1) = CRS("+proj=longlat +ellps=WGS84
                      +datum=WGS84 +no_defs")
  wt = raster::area(r1, weights = TRUE, na.rm = TRUE) 
  w = as.matrix(wt) 
  
  CHsp = (TX[,,1300]+1000)/(TX[,,1300]+1000) 
  length(which(CHsp==1)) #3825
}

slhf_00_Dp1 = LH00[,,match(daynames,dimnames(LH00)[[3]])+1]
slhf_12_Dp1 = LH12[,,match(daynames,dimnames(LH12)[[3]])+1]
slhf_12_D = LH12[,,match(daynames,dimnames(LH12)[[3]])]

for(i in 1:days)
{
  slhf_00_Dp1[,,i] = slhf_00_Dp1[,,i] * CHsp
  slhf_12_Dp1[,,i] = slhf_12_Dp1[,,i] * CHsp
  slhf_12_D[,,i] = slhf_12_D[,,i] * CHsp
}

slhf_6120 = slhf_00_Dp1 + slhf_12_Dp1 - slhf_12_D


# standadized the unit


# LH (downwards are positive
LH = -slhf_6120/86400

a = ls()
keepv = c('LH')
removeones = setdiff(a,keepv)
rm(list=removeones)

save.image("ERA5Land_LandAtms_add_LH.RData")
#####





