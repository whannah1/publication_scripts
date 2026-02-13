from sohip_methods import *
case_root = '/pscratch/sd/w/whannah/scream_scratch/pm-gpu'
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
var,var_str,klev_list = [],[],[]
unit_fac_list = []
def add_var(var_name,vstr=None,klev=0,unit_fac=None): 
    var.append(var_name);
    var_str.append(var_name if vstr is None else vstr)
    klev_list.append(klev)
    unit_fac_list.append(unit_fac)
#-------------------------------------------------------------------------------
add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-19.09.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-19 21:00',xlat=0,  xlon=  90,tlat= -6.99,tlon=  84.74,slat= 10.24,slon=  94.16)
# add_case('2025-SOHIP-RRM-00.256x2-eq-ind-v1.2023-06-21.02.NN_420',n='256x2-eq-ind-v1',xtime='2023-06-21 21:00',xlat=-5, xlon=  80,tlat= -3.05,tlon=  75.97,slat= 13.56,slon=  85.35)
# add_case('2025-SOHIP-RRM-00.256x2-ptgnia-v1.2023-06-13.19',       n='256x2-ptgnia-v1',xtime='2023-06-14 02:00',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)
# add_case('2025-SOHIP-RRM-00.256x2-sc-ind-v1.2023-06-21.09',       n='256x2-sc-ind-v1',xtime='2023-06-21 15:00',xlat=-50,xlon=  80,tlat=-52.49,tlon=  67.04,slat=-51.03,slon=  98.64)
# add_case('2025-SOHIP-RRM-00.256x2-sc-pac-v1.2023-06-14.15',       n='256x2-sc-pac-v1',xtime='2023-06-15 04:00',xlat=-35,xlon=-135,tlat=-34.73,tlon=-136.73,slat=-43.76,slon=-114.47)
# add_case('2025-SOHIP-RRM-00.256x2-se-pac-v1.2023-06-12.16',       n='256x2-se-pac-v1',xtime='2023-06-13 04:00',xlat=-50,xlon= -95,tlat=-49.60,tlon= -94.45,slat=-51.80,slon= -63.70)
# add_case('2025-SOHIP-RRM-00.256x2-sw-ind-v1.2023-06-12.06',       n='256x2-sw-ind-v1',xtime='2023-06-12 19:00',xlat=-50,xlon=  45,tlat=-49.61,tlon=  45.20,slat=-51.79,slon=  75.97)

# add_case('2025-SOHIP-RRM-00.256x3-ptgnia-v1.2023-06-13.19',       n='256x3-ptgnia-v1',xlat=-50,xlon= -60,tlat=-49.46,tlon= -60.24,slat=  None,slon=   None)

path_len_km = 800e3  # total path distance in meters (dist_km)
path_spc_km = 2      # spacing between interpolated path points
path_ncells = 1      # number of cells to consider nearest to each point (ncll)

# view width/height in degrees
dx,dy = 40,20
# dx,dy = 60,30
# dx,dy = 180,90

htype_2D = 'output.scream.2D.10min.INSTANT.nmins_x10'
htype_3D = 'output.scream.3D.10min.INSTANT.nmins_x10'
# first_file,num_files = -1,1
use_snapshot,ss_t = False,-1

#-------------------------------------------------------------------------------

add_var('T_mid',obs_var='T',vstr='1mb Temperature',plev=1e2,htype=htype_3D)

# add_var('ps')
# add_var('omega')
# add_var('horiz_winds')
# add_var('qv')
# add_var('T_mid')
# add_var('z_mid')

fig_file = 'figs/FXX-curtain-v1.png'

#-------------------------------------------------------------------------------
use_remap,remap_grid = False,'90x180' # 90x180 / 180x360

print_stats = True

var_x_case = False

plot_diff = False

num_plot_col = int(np.sqrt(len(case_list)))

use_common_label_bar = True

if use_snapshot: print(); print(f'{hapy.tclr.RED}WARNING - snapshot mode enabled! (ss_t={ss_t}){hapy.tclr.END}'); print()

#---------------------------------------------------------------------------------------------------
# Set up plot resources

if case_list==[]: raise ValueError('ERROR - case list is empty!')
num_var,num_case = len(var),len(case_list)

# set up figure object
(d1,d2) = (num_var,num_case) if var_x_case else (num_case,num_var)

fdx=10;figsize = (fdx*num_case,fdx*num_var) if var_x_case else (fdx*num_var,fdx*num_case)
# figsize = (30,30)
# figsize = (10,10)

fig = plt.figure(figsize=figsize)

# axs = np.empty([d1,d2])

#---------------------------------------------------------------------------------------------------
def get_comp(case):
    comp = 'scream'
    return comp
#---------------------------------------------------------------------------------------------------
def get_ctr_str(glb_avg=None): return ''
# def get_ctr_str(glb_avg=None):
#    ctr_str = ''
#    if 'lat1' in globals():
#       lat1_str = f'{lat1}N' if lat1>=0 else f'{(lat1*-1)}S'
#       lat2_str = f'{lat2}N' if lat2>=0 else f'{(lat2*-1)}S'
#       ctr_str += f' {lat1_str}:{lat2_str} '
#    if 'lon1' in globals():
#       lon1_str = f'{lon1}E' #if lon1>=0 and lon1<=360 else f'{(lon1*-1)}S'
#       lon2_str = f'{lon2}E' #if lon2>=0 and lon2<=360 else f'{(lon2*-1)}S'
#       ctr_str += f' {lon1_str}:{lon2_str} '
#    # if glb_avg is not None:
#    #    # add logic here t display global average value?
#    return ctr_str
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
path_lat = [None]*num_case
path_lon = [None]*num_case
for v in range(num_var):
    hapy.print_line()
    print('  var: '+hapy.tclr.MAGENTA+var[v]+hapy.tclr.END)
    data_list,lat_list,lon_list = [],[],[]
    glb_avg_list = []
    std_list,cnt_list = [],[]
    for c in range(num_case):
        print('    case: '+hapy.tclr.GREEN+case_list[c]+hapy.tclr.END)
        case_opts = case_opts_list[c]
        case_root = case_opts['p']
        case_sub  = case_opts['s']
        #-----------------------------------------------------------------------
        # read the data
        file_path = f'{case_root}/{case_list[c]}/{case_sub}/*{htype_list[v]}*'
        
        file_list = sorted(glob.glob(file_path))
        if 'first_file' in locals(): file_list = file_list[first_file:]
        if 'num_files'  in locals(): file_list = file_list[:num_files]

        for f in file_list: print(f'    {hapy.tclr.YELLOW}{f}{hapy.tclr.END}')

        if file_list==[]:
            print(); print('ERROR - file_list is empty!')
            print(); print(f'file_path: {file_path}')
            print()
        #-----------------------------------------------------------------------
        dt = datetime.datetime.strptime(case_opts['xtime'], '%Y-%m-%d %H:%M')
        target_time = cftime.DatetimeNoLeap(dt.year, dt.month, dt.day, dt.hour, dt.minute, 0)
        file_list = reduce_file_list_to_target_time(file_list,target_time)
        #-----------------------------------------------------------------------
        for f in file_list: print(f'    {hapy.tclr.YELLOW}{f}{hapy.tclr.END}')
        #-----------------------------------------------------------------------
        ds = ux.open_mfdataset(case_opts['g'], file_list, data_vars='all')
        ds = ds.sel(time=target_time, method='nearest') # select time nearest SOHIP observation
        #-----------------------------------------------------------------------
        data = ds[var[v]]
        #-----------------------------------------------------------------------
        # print(); print(data)
        # exit()
        #-----------------------------------------------------------------------
        # # Get rid of lev dimension
        # if 'lev' in data.dims:
        #    if klev_list[v] is not None:
        #       data = data.isel(lev=klev_list[v])
        #-----------------------------------------------------------------------
        # adjust units
        if unit_fac_list[v] is not None: data = data * unit_fac_list[v]
        #-----------------------------------------------------------------------
        # print stats before time averaging
        if print_stats:
            hapy.print_stat(data,name=var[v],stat='naxsh',indent='    ',compact=True)
        #-----------------------------------------------------------------------
        # # average over time dimension
        # if 'time' in data.dims : 
        #    if use_snapshot:
        #       data = data.isel(time=ss_t)
        #    else:
        #       data = data.mean(dim='time')
        #-----------------------------------------------------------------------
        # Calculate area weighted global mean
        if 'area' in locals() :
            gbl_mean = ( (data*area).sum() / area.sum() ).values 
            print(hapy.tclr.CYAN+f'      Area Weighted Global Mean : {gbl_mean:6.4}'+hapy.tclr.END)
            glb_avg_list.append(gbl_mean)
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
    # cmap = cmocean.cm.rain
    # cmap = cmocean.cm.balance
    # cmocean.cm.amp
    #---------------------------------------------------------------------------
    # set color bar levels
    clev = None
    # if 'precip' in var[v]      : clev = np.logspace( 0, 3, num=30).round(decimals=2)
    # if var[v]=='VapWaterPath'  : clev = np.arange(10,90+1,1)
    #---------------------------------------------------------------------------
    # create the path between ISS and tangent positions
    if v==0:
        for c in range(num_case):
            case_opts = case_opts_list[c]
            slat,slon = float(case_opts['slat']),float(case_opts['slon'])
            tlat,tlon = float(case_opts['tlat']),float(case_opts['tlon'])
            # define path outward from a given center location
            (npts,path_coord,path_lat[c],path_lon[c]) = calculate_path(tlat,tlon,slat,slon,path_len_km,path_spc_km)
    #---------------------------------------------------------------------------
    # Create plot
    for c in range(num_case):
        case_opts = case_opts_list[c]
        ctr_lon   = case_opts['xlon']
        ctr_lat   = case_opts['xlat']

        #-----------------------------------------------------------------------
        img_kwargs = {}
        img_kwargs['origin']    = 'lower'
        img_kwargs['cmap']      = cmap

        if plot_diff and c!=0:
            img_kwargs['cmap']   = cmocean.cm.balance
            img_kwargs['vmin']   = diff_data_min
            img_kwargs['vmax']   = diff_data_max
            clev = None

        if clev is not None: img_kwargs['norm'] = mcolors.BoundaryNorm(clev, ncolors=256)

        (d1,d2,ip) = (num_var,num_case,v*num_case+c+1) if var_x_case else (num_case,num_var,c*num_var+v+1)

        ax = fig.add_subplot(d1,d2,ip)

        ax.set_title(case_opts['n'],  fontsize=15, loc='left')
        ax.set_title(var_str[v],      fontsize=15, loc='right')

        contour = ax.contour(X, Y, Z, levels=10, colors='black', **img_kwargs)

        cbar = fig.colorbar(contour, ax=ax, fraction=0.02, orientation='vertical')
        cbar.ax.tick_params(labelsize=15)

#---------------------------------------------------------------------------------------------------
# Finalize plot
fig.savefig(fig_file, dpi=100, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
