
# BEST

```shell
ncremap -G ttl='BEST grid'#latlon==180,360#lat_typ=uni#lon_typ=180_wst -g grid_files/BEST_180x360_scrip.20241001.nc
SRC_GRID_FILE=grid_files/ne30pg2_scrip.nc
DST_GRID_FILE=grid_files/BEST_180x360_scrip.20241001.nc
MAP_FILE=map_files/map_ne30pg2_to_180x360_BEST_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}
```

# HadSST

```shell
SRC_GRID_FILE=grid_files/HadSST_180x360_scrip.20241001.nc
DST_GRID_FILE=grid_files/ne30pg2_scrip.nc
MAP_FILE=map_files/map_180x360_HadSST_to_ne30pg2_traave_20241001.nc
ncremap -G ttl='HadSST Equi-Angular grid 180x360'#latlon=180,360#lat_typ=uni#lat_drc=n2s#lon_typ=180_wst  -g ${SRC_GRID_FILE}
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

SRC_DATA=/global/cfs/cdirs/m3312/whannah/obs_data/HadSST/HadISST_sst.nc
DST_DATA=/global/cfs/cdirs/m3312/whannah/obs_data/HadSST/HadISST_sst.remap_ne30pg2.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
```


# HadCRU data

```shell
ncremap -G ttl='HadCRU grid'#latlon==36,72#lat_typ=uni#lon_typ=grn_ctr -g grid_files/HadCRU_36x72_scrip.20241001.nc
SRC_GRID_FILE=grid_files/ne30pg2_scrip.nc
DST_GRID_FILE=grid_files/HadCRU_36x72_scrip.20241001.nc
MAP_FILE=map_files/map_ne30pg2_to_36x72_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

DATA_ROOT=/global/cfs/cdirs/m3312/whannah/obs_data/HadCRU
SRC_DATA=${DATA_ROOT}/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc
DST_DATA=${DATA_ROOT}/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.remap_ne30pg2.nc
MAP_FILE=map_files/map_36x72_to_ne30pg2_traave_20241001.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
```

# CERES

```shell
SRC_GRID_FILE=grid_files/cmip6_180x360_scrip.20181001.nc
DST_GRID_FILE=grid_files/ne30pg2_scrip.nc
MAP_FILE=map_files/map_180x360_to_ne30pg2_traave_20241001.nc
ncremap -G ttl=Equi-Angular grid 1x1 degree, dimensions 180x360, cell edges on Poles/Equator and Prime Meridian/Date Line#latlon=180,360#lat_typ=uni#lon_typ=grn_wst -g ${SRC_GRID_FILE}
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}
DATA_ROOT=/global/cfs/cdirs/m3312/whannah/obs_data/CERES
SRC_DATA=${DATA_ROOT}/CERES_EBAF-TOA_Ed4.2_Subset_200003-202406.nc
DST_DATA=${DATA_ROOT}/CERES_EBAF-TOA_Ed4.2_Subset_200003-202406.remap_ne30pg2.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
```