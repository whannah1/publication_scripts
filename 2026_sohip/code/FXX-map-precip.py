from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#-------------------------------------------------------------------------------
case_opts_list = []
case_list = []
def add_case(case_in,**kwargs):
    case_list.append(case_in)
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    if 'g' not in tmp_opts:
        tmp_opts['g'] = get_grid_file(case_in)
    tmp_opts['p'] = case_root
    tmp_opts['s'] = 'run'
    case_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------
imerg_grid = '/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc'
imerg_root = '/global/cfs/cdirs/m4842/whannah/IMERG/hourly_hdf'
# imerg_root = '/global/cfs/cdirs/m4842/whannah/IMERG/hourly_netcdf'

# file_list = sorted(glob.glob(f'{imerg_root}/*'))
# for f in file_list:
#     os.system(f'ncdump -k {f}')
# exit()

# add_case('IMERG',n='IMERG',xtime='2023-06-19 21:00',g=imerg_grid,xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('IMERG',n='IMERG',xtime='2023-06-21 21:00',g=imerg_grid,xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('IMERG',n='IMERG',xtime='2023-06-14 02:00',g=imerg_grid,xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  -51.66,slon=   -28.91)
# add_case('IMERG',n='IMERG',xtime='2023-06-21 15:00',g=imerg_grid,xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('IMERG',n='IMERG',xtime='2023-06-15 04:00',g=imerg_grid,xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('IMERG',n='IMERG',xtime='2023-06-13 04:00',g=imerg_grid,xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
add_case('IMERG',n='IMERG',xtime='2023-06-12 19:00',g=imerg_grid,xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  -51.66,slon=   -28.91)
# add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

alt_subplot_method = True # set to true for viewing multiple model cases with different projections

path_len_km = 1200  # total path distance in meters (dist_km)
path_spc_km = 5      # spacing between interpolated path points
path_ncells = 1      # number of cells to consider nearest to each point (ncll)

# view width/height in degrees
dx,dy = 40,20
# dx,dy = 60,30
# dx,dy = 180,90

htype_2D = 'output.scream.2D.10min.INSTANT.nmins_x10'
htype_3D = 'output.scream.3D.10min.INSTANT.nmins_x10'
# first_file,num_files = -1,1

#-------------------------------------------------------------------------------
var,var_str,plev_list,klev_list = [],[],[],[]
unit_fac_list = []
obs_var_list = []
htype_list = []
def add_var(var_name,vstr=None,plev=None,klev=None,unit_fac=None,obs_var=None,htype=None): 
    var.append(var_name);
    var_str.append(var_name if vstr is None else vstr)
    plev_list.append(plev); klev_list.append(klev)
    unit_fac_list.append(unit_fac)
    obs_var_list.append(obs_var)
    htype_list.append(htype)
#-------------------------------------------------------------------------------

add_var('precip_total_surf_mass_flux', vstr='precip',unit_fac=86400*1e3,htype=htype_2D)

fig_file = 'figs/FXX-map-precip.png'

#-------------------------------------------------------------------------------
print_stats = True
var_x_case = True
plot_diff = False

# num_plot_col = int(np.sqrt(len(case_list)))
# use_common_label_bar = True

#---------------------------------------------------------------------------------------------------
# Set up plot resources

if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case_list)

# set up figure object
(d1,d2) = (num_var,num_case) if var_x_case else (num_case,num_var)

fdx,fdy=15,10;figsize = (fdx*num_case,fdy*num_var) if var_x_case else (fdx*num_var,fdy*num_case)

if alt_subplot_method:
    fig = plt.figure(figsize=figsize,layout="constrained")
else:
    ctr_lon = case_opts_list[0]['xlon']
    proj_plot = ccrs.PlateCarree(central_longitude=ctr_lon)
    fig,axs = plt.subplots( d1,d2,figsize=figsize,layout="constrained",squeeze=False,
                            subplot_kw={'projection':proj_plot})
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
path_lat = [None]*num_case
path_lon = [None]*num_case
for v in range(num_var):
    hapy.print_line()
    print('  var: '+hapy.tclr.MAGENTA+var[v]+hapy.tclr.END)
    data_list,lat_list,lon_list = [],[],[]
    std_list,cnt_list = [],[]
    for c in range(num_case):
        print('    case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
        case_opts = case_opts_list[c]
        case_root = case_opts['p']
        case_sub  = case_opts['s']
        #-------------------------------------------------------------------------
        # read the data
        if case_list[c]=='IMERG':
            file_path = f'{imerg_root}/3B-HHR.MS.MRG.3IMERG.2023*'
            file_list = sorted(glob.glob(file_path))
        else:
            file_path = f'{case_root}/{case_list[c]}/{case_sub}/*{htype_list[v]}*'
            file_list = sorted(glob.glob(file_path))
            if 'first_file' in locals(): file_list = file_list[first_file:]
            if 'num_files'  in locals(): file_list = file_list[:num_files]

        if file_list==[]:
            print(); print('ERROR - file_list is empty!')
            print(); print(f'file_path: {file_path}')
            print()
        #-----------------------------------------------------------------------
        dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
        if case_list[c]=='IMERG':
            target_time = cftime.DatetimeJulian(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
        else:
            target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
        file_list = reduce_file_list_to_target_time(file_list,target_time)
        #-----------------------------------------------------------------------
        for f in file_list: print(f'    {hapy.tclr.YELLOW}{f}{hapy.tclr.END}')
        #-----------------------------------------------------------------------
        kwargs = {}
        if case_list[c]=='IMERG':
            kwargs['group'] = 'Grid'
            kwargs['engine'] = 'netcdf4'
            # kwargs['use_cftime']=True
        ds = ux.open_mfdataset(case_opts['g'], file_list, **kwargs)
        ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        #-----------------------------------------------------------------------
        # if case_list[c]=='IMERG':
        #     actual_ds_time = ds.time_bnds.isel(nv=1).values
        # else:
        #     actual_ds_time = ds.time.values
        # print()
        # print(f'actual_ds_time: {actual_ds_time}')
        # print()
        #-----------------------------------------------------------------------
        tvar = var[v]
        if case_list[c]=='IMERG':
            data = ds['precipitation'].stack(n_face=('lat','lon'))*24. # mm/hr to mm/day
        else:
            data = ds[tvar]#.lon_to_180()
        #-----------------------------------------------------------------------
        if 'lev' in data.dims:
            if klev_list[v] is     None and plev_list[v] is     None: raise ValueError('ERROR - no p/k level specified')
            if klev_list[v] is not None and plev_list[v] is not None: raise ValueError('ERROR - cannot specify both p/k levels')
            if klev_list[v] is not None: data = data.isel(lev=klev_list[v])
            if plev_list[v] is not None:
                data = hapy.vinth2p_simple(data, ds['hyam'], ds['hybm'], plev_list[v], ds['ps'], interp_type='linear')
                if len(data.plev)==1: data = data.isel(plev=0)
        #-----------------------------------------------------------------------
        # adjust units
        if case_list[c]!='IMERG':
            if unit_fac_list[v] is not None: data = data * unit_fac_list[v]
        #-----------------------------------------------------------------------
        # if add_obs:
        #     obs_file_path = f'/global/cfs/projectdirs/m4842/whannah/ERA5/ERA5_validation.atm.*'
        #     obs_file_list = sorted(glob.glob(obs_file_path))
        #     obs_grid_file = f'/global/cfs/projectdirs/m3312/whannah/HICCUP/files_grid/scrip_ERA5_721x1440.nc'
        #     ds = ux.open_mfdataset(obs_grid_file, obs_file_list, data_vars='all')
        #-----------------------------------------------------------------------
        # print stats before time averaging
        if print_stats: hapy.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
        #-----------------------------------------------------------------------
        data_list.append( data )

    #-----------------------------------------------------------------------------------------------
    data_min = np.min([np.nanmin(d) for d in data_list])
    data_max = np.max([np.nanmax(d) for d in data_list])
    #---------------------------------------------------------------------------
    if plot_diff:
        tmp_data = copy.deepcopy(data_list)
        for c in range(num_case): tmp_data[c] = data_list[c] - data_list[diff_base]
        # diff_data_max = np.max([np.nanmax(d) for d in tmp_data])
        # diff_data_min = np.min([np.nanmin(d) for d in tmp_data])
        diff_data_max = np.max([np.nanmax(np.absolute(d)) for d in tmp_data])
        diff_data_min = -1*diff_data_max
        for c in range(num_case):
            if c!=diff_base:
                data_list[c] = data_list[c] - data_list[diff_base]
    #---------------------------------------------------------------------------
    # set color map
    cmap = 'viridis'
    if 'precip' in var[v]                   : cmap = cmocean.cm.rain
    if 'horiz_winds' in var[v]              : cmap = cmocean.cm.balance
    if var[v]=='U'                          : cmap = cmocean.cm.balance
    if var[v]=='T_2m'                       : cmap = cmocean.cm.amp
    if var[v]=='LiqWaterPath'               : cmap = cmocean.cm.rain
    if var[v]=='IceWaterPath'               : cmap = cmocean.cm.rain
    #---------------------------------------------------------------------------
    # set color bar levels
    clev = None
    # if 'precip' in var[v]    : clev = np.arange(1,20+1,1)
    # if 'precip' in var[v]    : clev = np.arange(4,80+4,4)
    if 'precip' in var[v]      : clev = np.logspace( 0, 3, num=30).round(decimals=2)
    # if 'prec' in var[v]      : clev = np.logspace( -2, 2, num=20).round(decimals=2)
    if var[v]=='VapWaterPath'  : clev = np.arange(10,90+1,1)
    if var[v]=='LiqWaterPath'  : clev = np.arange(1,61+3,3)/1e2
    if var[v]=='IceWaterPath'  : clev = np.arange(1,61+3,3)/1e1
    if var[v]=='LiqWaterPath'  : clev = np.logspace( -2, 1.2, num=20).round(decimals=2)
    if var[v]=='IceWaterPath'  : clev = np.logspace( -2, 1.2, num=20).round(decimals=2)
    #---------------------------------------------------------------------------
    # create the path between ISS and tangent positions
    if v==0:
        for c in range(num_case):
            case_opts = case_opts_list[c]
            if  'slat' in case_opts and 'slon' in case_opts \
            and 'tlat' in case_opts and 'tlon' in case_opts:
                slat,slon,tlat,tlon = case_opts['slat'],case_opts['slon'],case_opts['tlat'],case_opts['tlon']
                if  slat is not None and slon is not None\
                and tlat is not None and tlon is not None:
                    slat,slon,tlat,tlon = float(slat),float(slon),float(tlat),float(tlon)
                    # define path outward from a given center location
                    (npts,path_coord,path_lat[c],path_lon[c]) = calculate_path(tlat,tlon,slat,slon,path_len_km,path_spc_km)
    #---------------------------------------------------------------------------
    # Create plot
    for c in range(num_case):
        case_opts = case_opts_list[c]
        ctr_lon   = case_opts['xlon']
        ctr_lat   = case_opts['xlat']

        # proj_data = ccrs.PlateCarree(central_longitude=180)
        proj_plot = ccrs.PlateCarree(central_longitude=ctr_lon)
        # proj_plot = ccrs.Robinson(central_longitude=ctr_lon)

        #-----------------------------------------------------------------------
        img_kwargs = {}
        img_kwargs['origin']    = 'lower'
        img_kwargs['cmap']      = cmap
        # img_kwargs['transform'] = proj_data

        if plot_diff and c!=0:
            img_kwargs['cmap']   = cmocean.cm.balance
            img_kwargs['vmin']   = diff_data_min
            img_kwargs['vmax']   = diff_data_max
            clev = None

        if clev is not None: img_kwargs['norm'] = mcolors.BoundaryNorm(clev, ncolors=256)

        if alt_subplot_method:
            (d1,d2,ip) = (num_var,num_case,v*num_case+c+1) if var_x_case else (num_case,num_var,c*num_var+v+1)
            ax = fig.add_subplot(d1,d2,ip,projection=proj_plot)
        else:
            ax = axs[v,c] if var_x_case else axs[c,v]

        ax.coastlines(linewidth=0.2,edgecolor='white')
        ax.set_title(case_opts['n'],    fontsize=20, loc='left')
        ax.set_title(case_opts['xtime'],fontsize=20, loc='center')
        ax.set_title(var_str[v],        fontsize=20, loc='right')
        
        # ax.set_global()

        ax.set_extent( [ctr_lon-dx, \
                        ctr_lon+dx, \
                        max(-90,ctr_lat-dy), \
                        min( 90,ctr_lat+dy)], \
                        crs=ccrs.PlateCarree() )

        
        raster = data_list[c].to_raster(ax=ax)
        img = ax.imshow(raster, extent=ax.get_xlim()+ax.get_ylim(), **img_kwargs)

        # orientation = 'vertical' if var_x_case else 'horizontal'
        # if c==num_case-1: fig.colorbar(img, ax=ax, fraction=0.02, orientation=orientation)

        cbar = fig.colorbar(img, ax=ax, fraction=0.02, orientation='vertical')
        cbar.ax.tick_params(labelsize=18)

        #-----------------------------------------------------------------------
        # Plot ISS position marker (X)
        if 'slat' in case_opts and 'slon' in case_opts:
            slat,slon = case_opts['slat'],case_opts['slon'] # ISS position
            if slat is not None and slon is not None:
                ax.scatter(slon, slat, 
                                marker='x', 
                                s=100, 
                                c='red', 
                                linewidths=3, 
                                transform=ccrs.PlateCarree(),
                                zorder=5)
        #-----------------------------------------------------------------------
        # Plot tangent point marker (circle)
        if 'tlat' in case_opts and 'tlon' in case_opts:
            tlat,tlon = case_opts['tlat'],case_opts['tlon'] # tangent point of retreival
            if tlat is not None and tlon is not None:
                ax.scatter(tlon, tlat, 
                                marker='o', 
                                s=100, 
                                c='red', 
                                # edgecolors='red', 
                                # facecolors='none', 
                                linewidths=2,
                                transform=ccrs.PlateCarree(),
                                zorder=5)
        #-----------------------------------------------------------------------
        # overlay path markers
        if path_lon[c] is not None and path_lat[c] is not None:
            # ax.scatter(path_lon[c], path_lat[c], 
            #             marker='o', 
            #             s=1, 
            #             c='red', 
            #             linewidths=1,
            #             transform=ccrs.PlateCarree(),
            #             zorder=5)
            ax.plot(path_lon[c], path_lat[c], 
                      color='red', 
                      linewidth=1,
                      transform=ccrs.PlateCarree())

#---------------------------------------------------------------------------------------------------
# Finalize plot
# plt.tight_layout()
fig.savefig(fig_file, dpi=100, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
