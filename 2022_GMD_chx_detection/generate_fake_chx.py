import os, ngl, sys, numba, copy, xarray as xr, numpy as np, re, string, subprocess as sp
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import numba, itertools
import pg_checkerboard_utilities as pg
class tcolor: ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'

mac_file = os.getenv('HOME')+'/Data/Obs/MAC/daily/tlwp1deg_maclwpv1.200501.nc4'

fig_type = 'png'
fig_file = 'generate_fake_chx'
# tmp_file = 'generate_fake_chx'

scripfile_fv_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
scripfile_pg_path = 'scrip_files/ne30pg2_scrip.nc'

map_file = os.getenv('HOME')+'/maps/map_180x360_to_ne30pg2_aave_20210113.nc'

src_file_name = 'fake_chx_fv.nc'
dst_file_name = 'fake_chx_pg.nc'

regenerate_data = False

subset_min_length = 4

# scripfile_fv = xr.open_dataset(scripfile_fv_path)
# print(scripfile_fv)
# exit()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if regenerate_data:
    #---------------------------------------------------------------------------
    # Generate a checkerboard signal on the native MAC grid
    #---------------------------------------------------------------------------
    
    @numba.njit
    def generate_chx(chx,num_lat,num_lon):
        for j in range(num_lat):
            for i in range(num_lon):
                if (i%2): 
                    if (j%2): 
                        chx[j,i] = 0
                    else:
                        chx[j,i] = 1
                else:
                    if (j%2): 
                        chx[j,i] = 1
                    else:
                        chx[j,i] = 0

    mac_ds = xr.open_dataset(mac_file)
    num_lat,num_lon = len(mac_ds['lat']),len(mac_ds['lon'])
    chx = np.zeros([num_lat,num_lon])
    generate_chx(chx,num_lat,num_lon)

    fv_ds = xr.Dataset()
    fv_ds['lat'] = mac_ds['lat']
    fv_ds['lon'] = mac_ds['lon']
    fv_ds['chx'] = (['lat','lon'],chx)

    #---------------------------------------------------------------------------
    # Regrid the fake data to ne30pg2
    #---------------------------------------------------------------------------

    def run_cmd(cmd,indent=' '*4):
        print(indent+tcolor.GREEN + cmd + tcolor.ENDC)
        try:
            sp.check_output(cmd,shell=True,universal_newlines=True)
        except sp.CalledProcessError as error:
            exit(error.output)

    

    fv_ds.to_netcdf(path=src_file_name,mode='w')

    # # deal with fill value to avoid warning - too slow!
    # cmd = run_cmd(f'ncatted -a _FillValue,chx,m,f,1.0e36 {src_file_name} {src_file_name}')

    run_cmd(f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name}')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

fv_ds = xr.open_dataset(src_file_name)
pg_ds = xr.open_dataset(dst_file_name)

data_fv = fv_ds['chx']
data_pg = pg_ds['chx']

#-------------------------------------------------------------------------------
# Pattern detection
#-------------------------------------------------------------------------------

### full set of patterns
rotate_sets = True; sets = pg.all_possible_sets

(num_set,set_len) = sets.shape
set_coord,nn_coord = np.arange(num_set),np.arange(set_len)
sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})

set_labels = pg.get_set_labels(sets)

if True:
    # find neighbors
    scripfile_pg = xr.open_dataset(scripfile_pg_path)
    ( neighbors, bearings ) = pg.get_neighbors_and_bearings(scripfile_pg)

    # print('      calculating anomalies...')
    data_dum = data_pg.expand_dims(dim='time',axis=0)
    # print(); print(data_dum); print()
    nn_anomalies = pg.remove_neighbor_mean_numba( data_dum.values, neighbors, bearings, sort_by_bearing=False)

    # define array of neighbor states
    nn_states = np.sign( nn_anomalies )
    nn_states = xr.where(nn_states<=0,0,nn_states)
    nn_states = nn_states.astype(type(sets.values[0,0]))

    # print('      counting sets...')
    # ncol = np.arange(len(data_pg))
    ncol = len(data_pg)
    cnt_list = pg.count_sets_per_col_numba( nn_states, sets.values, ncol, rotate_sets )

    dims = ('ncol','set')
    cnt_ds = xr.Dataset()
    cnt_ds['cnt'] = ( dims, cnt_list )
    # cnt_ds.coords['ncol'] = ('ncol',data_pg['ncol'].values)
    cnt_ds.coords['set'] = ('set',set_coord)

    # print('      writing to file: '+cnt_tmp_file)
    # cnt_ds.to_netcdf(path=case_tmp_file,mode='w')

### combine sets that contain partial checkerboard

chk_idx,nox_idx = [],[]
for s in range(num_set):
    if pg.is_partial_checkerboard(sets[s,:],subset_length=subset_min_length):
        chk_idx.append(s)
    else:
        nox_idx.append(s)

ds_sum_nox = cnt_ds.isel(set=nox_idx).sum(dim='set')
ds_sum_chk = cnt_ds.isel(set=chk_idx).sum(dim='set')
ds_sum_nox = ds_sum_nox.expand_dims('set').assign_coords(coords={'set':np.array([0])})
ds_sum_chk = ds_sum_chk.expand_dims('set').assign_coords(coords={'set':np.array([1])})
# cnt_ds = xr.concat( [ ds_sum_nox, ds_sum_chk ], 'set' )
cnt_ds = xr.concat( [ ds_sum_chk ], 'set' ) ### only show partial chx

# set_labels = ['no checkerboard','partial checkerboard']
set_labels = ['partial checkerboard'] ### only show partial chx
num_set = len(set_labels)
sort_idx = [c for c in range(num_set)]

#-------------------------------------------------------------------------------
# plot the data on both grids, plus the map of chx occurrence
#-------------------------------------------------------------------------------
wks = ngl.open_wks(fig_type,fig_file)
plot = [None]*3

res = hs.res_contour_fill_map()
res.cnFillMode    = "CellFill"
res.mpOutlineBoundarySets = 'NoBoundaries'
res.mpCenterLonF = 0.
res.mpMinLatF = -60
res.mpMaxLatF =  60
# res.mpMinLonF = -60
# res.mpMaxLonF =  60

res_fv = copy.deepcopy(res)
res_pg = copy.deepcopy(res)


scripfile_fv = xr.open_dataset(scripfile_fv_path)
res_fv.sfXArray      = scripfile_fv['grid_center_lon'].values
res_fv.sfYArray      = scripfile_fv['grid_center_lat'].values
res_fv.sfXCellBounds = scripfile_fv['grid_corner_lon'].values
res_fv.sfYCellBounds = scripfile_fv['grid_corner_lat'].values

scripfile_pg = xr.open_dataset(scripfile_pg_path)
res_pg.sfXArray      = scripfile_pg['grid_center_lon'].values
res_pg.sfYArray      = scripfile_pg['grid_center_lat'].values
res_pg.sfXCellBounds = scripfile_pg['grid_corner_lon'].values
res_pg.sfYCellBounds = scripfile_pg['grid_corner_lat'].values


plot[0] = ngl.contour_map(wks, data_fv.values.ravel(), res_fv) 
plot[1] = ngl.contour_map(wks, data_pg.values, res_pg) 

plot[2] = ngl.contour_map(wks, cnt_ds['cnt'].isel(set=0), res_pg) 

hs.set_subtitles(wks, plot[0], '', '180x360', '', font_height=0.01)
hs.set_subtitles(wks, plot[1], '', 'ne30pg2', '', font_height=0.01)

#-------------------------------------------------------------------------------
# Finalize plot
#-------------------------------------------------------------------------------

pnl_res = hs.setres_panel()
pnl_res.nglPanelYWhiteSpacePercent = 5

# pnl_res.nglPanelFigureStrings            = list(string.ascii_lowercase)
# pnl_res.nglPanelFigureStringsJust        = "TopLeft"
# pnl_res.nglPanelFigureStringsFontHeightF = 0.01

# layout = [1,len(plot)]
layout = [len(plot),1]

ngl.panel(wks,plot,layout,pnl_res)
ngl.end()

hc.trim_png(fig_file)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
