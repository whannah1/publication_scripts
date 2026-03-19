from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
#-------------------------------------------------------------------------------
'''
python -i
exec( open( 'code/FXX-curtain-v2.py' ).read() )
'''
#-------------------------------------------------------------------------------
case_opts_list = []
case_list = []
def add_case(case_in,**kwargs):
    case_list.append(case_in)
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    tmp_opts['g'] = get_grid_file(case_in)
    tmp_opts['p'] = case_root
    tmp_opts['s'] = 'run'
    case_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
# add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)

#-------------------------------------------------------------------------------
var_list,var_opts_list = [],[]
def add_var(var_name,**kwargs): 
    var_list.append(var_name);
    tmp_opts = {}
    for k, val in kwargs.items(): tmp_opts[k] = val
    var_opts_list.append(tmp_opts)
#-------------------------------------------------------------------------------
add_var('T_mid',         vstr='Temperature [K]')
# add_var('p_mid',         vstr='Pressure [mb]')
add_var('wind_speed',    vstr='Wind Speed [m/s]')
add_var('omega',         vstr='Omega [Pa/s]')
# add_var('qv',            vstr='Water Vapor')
# add_var('horiz_winds_u', vstr='U wind')
# add_var('horiz_winds_v', vstr='V wind')

#-------------------------------------------------------------------------------

fig_file = 'figs/FXX-curtain-v2.png'

curtain_data_root = '/global/cfs/cdirs/m4842/whannah/curtain_data'

path_len_km = 1200   # total path distance [km]
path_spc_km = 2      # path point spacing [km]
path_ncells = 2      # number of cells to average for each path point

# first_file,num_files = -1,1

#-------------------------------------------------------------------------------

print_stats = False
var_x_case  = False
# plot_diff   = False

# num_plot_col = 2#int(np.sqrt(len(case_list)))
# use_common_label_bar = True

#---------------------------------------------------------------------------------------------------
if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var_list),len(case_list)
#---------------------------------------------------------------------------------------------------
# set up figure object

num_panels = num_case * num_var

if 'num_plot_col' in locals() and num_plot_col is not None:
    nrows = int(np.ceil(num_panels / float(num_plot_col)))
    ncols = num_plot_col
else:
    (nrows,ncols) = (num_var,num_case) if var_x_case else (num_case,num_var)

# fdx,fdy=10,15;figsize = (fdx*num_case,fdy*num_var) if var_x_case else (fdx*num_var,fdy*num_case)
# fig = plt.figure(figsize=figsize)

fig, axs = plt.subplots(nrows, ncols, figsize=(10*ncols, 5*nrows), squeeze=False)
axs_flat = axs.flatten()


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
    var_opts = var_opts_list[v]
    hapy.print_line()
    print('  var: '+hapy.tclr.MAGENTA+var_list[v]+hapy.tclr.END)
    data_list,lat_list,lon_list = [],[],[]
    std_list,cnt_list = [],[]
    for c in range(num_case):
        print('\n'+' '*4+'case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
        case_opts = case_opts_list[c]
        case_root = case_opts['p']
        case_sub  = case_opts['s']
        #-----------------------------------------------------------------------
        # read the data
        f_name = case_opts['n']
        f_time = case_opts['xtime'].replace(' ','_').replace(':','_')
        file_path = f'{curtain_data_root}/{f_name}.{f_time}.curtain.ncells_{path_ncells}.len_{int(path_len_km)}.spc_{int(path_spc_km)}.nc'
        #-----------------------------------------------------------------------
        # print(f'\nfile_path: {file_path}\n')
        #-----------------------------------------------------------------------
        ds = xr.open_dataset(file_path)
        #-----------------------------------------------------------------------
        dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
        target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
        ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        #-----------------------------------------------------------------------
        if var_list[v]=='wind_speed':
            data = np.sqrt( np.square(ds['horiz_winds_u']) + np.square(ds['horiz_winds_u']) )
        else:
            data = ds[var_list[v]]
        #-----------------------------------------------------------------------
        if 'anomaly' in var_opts and var_opts['anomaly']:
            data = data - data.mean(dim='path_coord')
        #-----------------------------------------------------------------------
        # # adjust units
        # if 'unit_fac' in var_opts:
        #     data = data * var_opts['unit_fac']
        #-----------------------------------------------------------------------
        if print_stats: hapy.print_stat(data_interp_hgt,name=var_list[v],stat='naxsh',indent=' '*4,compact=True)
        #-----------------------------------------------------------------------
        data_list.append( data )
    #-----------------------------------------------------------------------------------------------
    data_min = np.min([np.nanmin(d) for d in data_list])
    data_max = np.max([np.nanmax(d) for d in data_list])
    #---------------------------------------------------------------------------
    # set color map
    cmap = 'viridis'
    # cmap = cmocean.cm.rain
    # cmap = cmocean.cm.balance
    # cmocean.cm.amp
    #---------------------------------------------------------------------------
    # set color bar levels
    clev = None
    # if 'precip' in var_list[v]      : clev = np.logspace( 0, 3, num=30).round(decimals=2)
    # if var_list[v]=='VapWaterPath'  : clev = np.arange(10,90+1,1)
    #---------------------------------------------------------------------------
    # Create plot
    for c in range(num_case):
        case_opts = case_opts_list[c]
        #-----------------------------------------------------------------------
        img_kwargs = {}
        # img_kwargs['origin']    = 'lower'
        img_kwargs['cmap']      = cmap
        img_kwargs['levels']    = 60

        # if plot_diff and c!=0:
        #     img_kwargs['cmap']   = cmocean.cm.balance
        #     img_kwargs['vmin']   = diff_data_min
        #     img_kwargs['vmax']   = diff_data_max
        #     clev = None

        if clev is not None:
            img_kwargs['norm'] = mcolors.BoundaryNorm(clev, ncolors=256)

        # (d1,d2,ip) = (num_var,num_case,v*num_case+c+1) if var_x_case else (num_case,num_var,c*num_var+v+1)
        # ax = fig.add_subplot(d1,d2,ip)        

        ip = v*num_case+c if var_x_case else c*num_var+v
        ax = axs_flat[ip]

        fontsize = 15

        ax.set_title(case_opts['n'],    fontsize=fontsize, loc='left')
        ax.set_title(case_opts['xtime'],fontsize=fontsize, loc='center')
        ax.set_title(var_opts['vstr'],  fontsize=fontsize, loc='right')

        data = data_list[c].values
        path = data_list[c]['path_coord'].values
        hgth = data_list[c]['height'].values / 1e3
        contour = ax.contourf(path, hgth, data.T, **img_kwargs)
        ax.set_box_aspect(0.6)

        cbar = fig.colorbar(contour, ax=ax, fraction=0.02, orientation='vertical')
        cbar.ax.tick_params(labelsize=fontsize)

        ax.set_xlabel('Distance from Tanget Point [km]',fontsize=fontsize)
        ax.set_ylabel('Altitude [km]',fontsize=fontsize)

#---------------------------------------------------------------------------------------------------
# Finalize plot
plt.tight_layout(w_pad=1, h_pad=2)
fig.savefig(fig_file, dpi=200, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
