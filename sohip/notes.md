
# SOHIP Observation date/time

```shell
2023-06-12 18:30 - South West Indian ocean
2023-06-13 04:00 - South East Pacific
2023-06-15 03:30 - South Central Pacific
2023-06-19 21:00 - Equatorial Indian Ocean
2023-06-21 14:30 - South Central Indian Ocean
2023-06-21 21:00 - Equatorial Indian Ocean
```

# Curtain Interpolated Data

```shell
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-eq-ind-v1.2023-06-19_21_00.curtain.ncells_2.len_1200.spc_2.nc
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-eq-ind-v1.2023-06-21_21_00.curtain.ncells_2.len_1200.spc_2.nc
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-ptgnia-v1.2023-06-14_02_00.curtain.ncells_2.len_1200.spc_2.nc
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-sc-ind-v1.2023-06-21_15_00.curtain.ncells_2.len_1200.spc_2.nc
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-sc-pac-v1.2023-06-15_04_00.curtain.ncells_2.len_1200.spc_2.nc
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-se-pac-v1.2023-06-13_04_00.curtain.ncells_2.len_1200.spc_2.nc
/global/cfs/cdirs/m4842/whannah/curtain_data/256x2-sw-ind-v1.2023-06-12_19_00.curtain.ncells_2.len_1200.spc_2.nc
```

# ERA5 validation data

```shell
# Dates for RRM grids - init ~12 hour prior to observation
2023-06-14 01:30 - init 2023-06-13 19Z - Patagonia / South Atlantic
2023-06-12 18:30 - init 2023-06-12 06Z - South West Indian ocean
2023-06-13 04:00 - init 2023-06-12 16Z - South East Pacific
2023-06-15 03:30 - init 2023-06-14 15Z - South Central Pacific
2023-06-19 21:00 - init 2023-06-19 09Z - Equatorial Indian Ocean
2023-06-21 14:30 - init 2023-06-21 02Z - South Central Indian Ocean
2023-06-21 21:00 - init 2023-06-21 09Z - Equatorial Indian Ocean

# commands to get initialization data - copied from HICCUP notes
DATE=20230613; HR=19; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230612; HR=06; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230612; HR=16; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230614; HR=15; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230619; HR=09; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230621; HR=02; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &
DATE=20230621; HR=09; nohup python -u ~/HICCUP/get_hindcast_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --output-root=/global/cfs/projectdirs/m4842/whannah/HICCUP > ${HOME}/E3SM_grid_support/2025-SOHIP-RRM/logs_hiccup/log.2025-SOHIP-get-ERA5.${DATE}.${HR} &

# same commands uising get_validation_data
DATE=20230613; HR=19; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230612; HR=06; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230612; HR=16; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230614; HR=15; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230619; HR=09; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230621; HR=02; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
DATE=20230621; HR=09; python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=${DATE} --start-hour=${HR} --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5

# get all dates at once:
python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230611 --final-date=20230622 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5


# get all dates - one at a time:
python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230611 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230612 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230613 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5

python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230614 --final-date=20230616 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230614 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230615 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230616 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5

python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230617 --final-date=20230622 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230617 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230618 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230619 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230620 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230621 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5
# python -u ~/HICCUP/get_validation_data.ERA5.py --start-date=20230622 --data-freq=1h --output-root=/global/cfs/projectdirs/m4842/whannah/ERA5


```

# IMERG grid file

```shell
ncremap -G ttl='Equi-Angular grid 0.1x0.1 deg'#latlon=1800,3600#lat_typ=uni#lon_typ=180_wst -g /pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc
```


# COSMIC2 Data

```shell
#---------------------------------------------------------------------------------------------------
ncdump -h /global/cfs/cdirs/m4842/whannah/COSMIC2/atmPrf/2023/atmPrf_2023_109/atmPrf_C2E2.2023.109.00.01.R21_0001.0001_nc | grep "long_name"
    Lat:long_name = "Latitude of perigee point" ;
    Lon:long_name = "Longitude of perigee point" ;
    MSL_alt:long_name = "Mean sea level altitude of perigee point" ;
    Ref:long_name = "Refractivity" ;
    Azim:long_name = "Azimuth of occultation plane wrt North" ;
    Pres:long_name = "Dry pressure" ;
    Temp:long_name = "Dry temperature" ;
    Bend_ang1:long_name = "L1 bending angle" ;
    Bend_ang2:long_name = "L2 bending angle" ;
    Bend_ang:long_name = "Ionosphere-corrected bending angle" ;
    Impact_height:long_name = "Impact height" ;
    Ref_L1:long_name = "Refractivity " ;
    Bend_ang_stdv:long_name = "Bending angle uncertainty" ;
    Bend_ang_conf:long_name = "Bending angle confidence" ;
#---------------------------------------------------------------------------------------------------
ncdump -h /global/cfs/cdirs/m4842/whannah/COSMIC2/wetPf2/2023/wetPf2_2023_109/wetPf2_C2E2.2023.109.00.01.R21_0001.0001_nc | grep "long_name"
    MSL_alt:long_name = "Geometric height above MSL" ;
    lat:long_name = "Geodetic latitude" ;
    lon:long_name = "Geodetic longitude" ;
    gph:long_name = "Geopotential height above MSL" ;
    Temp:long_name = "Temperature" ;
    Pres:long_name = "Pressure" ;
    Vp:long_name = "Water vapor pressure" ;
    sph:long_name = "Specific humidity" ;
    rh:long_name = "Relative humidity" ;
    ref:long_name = "Refractivity" ;
    temp_dry:long_name = "DRY T derived from variationally optimized BA (before entering 1DVAR)" ;
    pres_dry:long_name = "DRY p derived from variationally optimized BA (before entering 1DVAR)" ;
#---------------------------------------------------------------------------------------------------
```
