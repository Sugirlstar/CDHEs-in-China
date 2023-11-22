library(raster)
library(ncdf4)

pathtofile = "H:/"
file_list = list.files(pathtofile,pattern="\\.nc$",full.names=T)
file_list = sort(file_list) #order the file by names

#####
#transfer the nc to matrix and give the name
for(i in 1:length(file_list))
{
  input_nc = file_list[i]
  ncfile = nc_open(input_nc)
  rstack = stack(input_nc, varname = names(ncfile$var)[1] )
  formattedD = as.vector(sapply(names(rstack),function(data) gsub("X|\\.","",data)))
  vals = values(rstack)
  arraytest = array(vals,dim=c(dim(rstack)[2],dim(rstack)[1],dim(rstack)[3]))
  arrayfinal = aperm(arraytest, c(2,1,3))
  dimnames(arrayfinal)[[3]] = formattedD
  assign(names(ncfile$var)[1],arrayfinal)
}
#读入：d2m, t2m, slhf, sshf, swvl1, swvl2, swvl3
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

d2m_6120 = d2m[,,match(daynames,dimnames(d2m)[[3]])]
t2m_6120 = t2m[,,match(daynames,dimnames(t2m)[[3]])]
slhf_6120 = slhf[,,match(daynames,dimnames(slhf_6120)[[3]])]
swvl1_6120 = swvl1[,,match(daynames,dimnames(swvl1_6120)[[3]])]
swvl2_6120 = swvl2[,,match(daynames,dimnames(swvl2_6120)[[3]])]
swvl3_6120 = swvl3[,,match(daynames,dimnames(swvl3_6120)[[3]])]

for(i in 1:days)
{
  d2m_6120[,,i] = d2m_6120[,,i] * CHsp
  t2m_6120[,,i] = t2m_6120[,,i] * CHsp
  slhf_6120[,,i] = slhf_6120[,,i] * CHsp
  swvl1_6120[,,i] = swvl1_6120[,,i] * CHsp
  swvl2_6120[,,i] = swvl2_6120[,,i] * CHsp
  swvl3_6120[,,i] = swvl3_6120[,,i] * CHsp
}

length(which(d2m_6120[,,1]>-9999)) #3795

#####

#####
# standadized the unit
d2m_6120 = d2m_6120-273.15
t2m_6120 = t2m_6120-273.15
# 
es = 6.112 * exp((17.67*d2m_6120) / (d2m_6120+243.5)) 
er = 6.112 * exp((17.67*t2m_6120) / (t2m_6120+243.5))
VPD = er-es
# 
SM100 = swvl1_6120*7/100 + swvl2_6120*21/100 + swvl3_6120*72/100
SM28 = swvl1_6120*7/28 + swvl2_6120*21/28

a = ls()
keepv = c('VPD','SM100','SM28')
removeones = setdiff(a,keepv)
rm(list=removeones)

save.image("ERA5Land_LandAtms_add.RData")
#####





