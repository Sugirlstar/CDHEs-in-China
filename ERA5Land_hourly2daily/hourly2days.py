import xarray as xr
import pandas as pd
import os

def instHour2Day(input_dir,output_dir):

    # region -----------------------Step 1---------------------------------
    folderName = os.path.basename(input_dir)  # get the name of the folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith('.nc')]
    files.sort()  # order by the file name
    # endregion END

    # region -----------------------Step 2---------------------------------
    for file in files:

        pure_filename = file.rsplit('.', 1)[0]
        year = pure_filename[-4:]
        outputFname = os.path.join(output_dir, f'{folderName}_{year}_daily.nc')

        if not os.path.exists(outputFname):
            # read annual data
            ds = xr.open_dataset(os.path.join(input_dir, file))
            # get the vairable name
            variable_name = [var for var in ds.data_vars][0]
            # daily mean
            # move 12 hours backward to get Beijing Time
            ds['time'] = pd.to_datetime(ds['time'].values)+pd.Timedelta(hours=12)
            daily_mean = ds[variable_name].resample(time='1D').mean(dim='time')
            
            # write to the nc file
            daily_mean.to_netcdf(outputFname)
            print(f"******dailyAggragated: {outputFname}******")

            # delete the var to save space
            del ds, daily_mean
        else:
            print(outputFname," file already exist, skip")
            continue
    # endregion END


def acumHour2Day(input_dir,output_dir):
    # region -----------------------Step 1---------------------------------
    folderName = os.path.basename(input_dir)  # get the name of the folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith('.nc')]
    files.sort()  # order by the file name
    # endregion END

    # region -----------------------Step 2---------------------------------
    for file in files:

        pure_filename = file.rsplit('.', 1)[0]
        year = pure_filename[-4:]
        outputFname1 = os.path.join(output_dir, f'{folderName}_{year}_daily00.nc')
        outputFname2 = os.path.join(output_dir, f'{folderName}_{year}_daily12.nc')
        
        if not os.path.exists(outputFname1):
            #  read annual data
            ds = xr.open_dataset(os.path.join(input_dir, file))
            # get the vairable name
            variable_name = [var for var in ds.data_vars][0]
            data = ds[variable_name]
            # Extract data at 00:00 and 12:00 per day
            daily = data.sel(time=data.time.dt.hour == 0)
            daily12 = data.sel(time=data.time.dt.hour == 12)

            # write to the nc file
            daily.to_netcdf(outputFname1)
            daily12.to_netcdf(outputFname2)
            print(f"******dailyAggragated: {outputFname2}******")

            # delete the var to save space
            del data, daily
        else:
            print(outputFname1," file already exist, skip")
    # endregion END

accumtype = pd.read_csv("ERA5Land_accumulationVar.csv")
# Extract the 'varname' column
varname_values = accumtype['varname'].values
# 设置目录路径：为当前工作目录
current_directory = os.getcwd()
# Get all items in the directory
all_items = os.listdir(
    current_directory+'ERA5Land_China_hourlyData_1960_2022')

for item in all_items:
    filepath = os.path.join(
        current_directory, 'ERA5Land_China_hourlyData_1960_2022', item)
    outputpath = os.path.join(
        current_directory, 'ERA5Land_China_dailyData_1960_2022_timeshifted', item)
    if os.path.isdir(filepath) and item in varname_values:
        acumHour2Day(filepath,outputpath)  
    elif os.path.isdir(filepath):
        instHour2Day(filepath,outputpath) 

print("all done")
