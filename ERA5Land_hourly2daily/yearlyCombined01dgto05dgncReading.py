import os
import netCDF4 as nc
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# 设置目录路径：为当前工作目录
current_directory = os.getcwd()
# Get the file path
filedir = os.path.join(current_directory,"ERA5Land_China_dailyData_1960_2022_timeshifted")

filepathName = ["surface_latent_heat_flux_00","surface_latent_heat_flux_12"]

for fn in filepathName:
    filepath = os.path.join(filedir,fn)
    files = [f for f in os.listdir(filepath) if f.endswith('.nc')]
    files = sorted(files)
    datasets = []
    for file in files:
        # 读取每年的数据
        with xr.open_dataset(os.path.join(filepath,file)) as ds:
            varname = list(ds.variables.keys())[3]
            # aggragate
            new_lon = np.arange(72,136,0.5)
            new_lat = np.arange(18,54,0.5)
            dsnew = ds.interp(longitude=new_lon, latitude=new_lat)
            annual_data = dsnew[varname]
            datasets.append(annual_data)

        print(f"{file} ReadingIn")   
    print("------ALL read in------")
    combined = xr.concat(datasets, dim="time")
    combined.to_netcdf(os.path.join(current_directory,f"{fn}_China_1961-2022_daily_ERA5Land_05dg_shifted.nc"))
    print(f"{fn}_combined successfully")

print("ALL_done")


# def plot_spatial_distribution(new_A, t):
#                 # Extracts the spatial distribution at a given point in time
#                 spatial_data = new_A[t, :, :]

#                 # plot
#                 plt.figure(figsize=(10,6))
#                 plt.imshow(spatial_data, cmap='viridis', origin='lower')
#                 plt.colorbar(label="Value")
#                 plt.title(f"Spatial Distribution at Timestep {t}")
#                 plt.xlabel("Longitude")
#                 plt.ylabel("Latitude")
#                 plt.show()

#             test = np.array(annual_data)
#             plot_spatial_distribution(annual_data,1)

#             datasets.append(annual_data)




