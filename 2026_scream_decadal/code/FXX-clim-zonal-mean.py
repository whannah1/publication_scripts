import os, glob, subprocess as sp, numpy as np, xarray as xr, copy, string, sys, cmocean
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm
import hapy
# from hapy_setres import set_subtitles  # replaced with matplotlib equivalents

#-------------------------------------------------------------------------------
case, opts_list = [], []
def add_case(case_in, **kwargs):
    case.append(case_in)
    case_opts = {}
    for k, val in kwargs.items(): case_opts[k] = val
    opts_list.append(case_opts)

#-------------------------------------------------------------------------------
var_list       = []
eam_var_list   = []
eamxx_var_list = []
unit_list      = []
var_str_list   = []
var_fac_list   = []
def add_var(var_name, eam_var=None, eamxx_var=None, var_str=None, lev=None, fac=None, unit=None):
    var_list.append(var_name)
    eam_var_list.append(eam_var)
    eamxx_var_list.append(eamxx_var)
    unit_list.append(unit)
    var_fac_list.append(fac)
    var_str_list.append(var_name if var_str is None else var_str)
#-------------------------------------------------------------------------------
tmp_path_ne1024 = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal'
tmp_path_ne256  = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
tmp_path_qbo_bm = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu/'

add_case('ERA5', n='ERA5')

add_case('decadal-production-run6',
         n='SCREAM 3-km', comp='eamxx',
         htype='_ANN_199501_200412_climo.nc',
         p=tmp_path_ne1024, s='climo_180x360')

add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',
         n='SCREAM 13-km control', comp='eamxx',
         htype='_ANN_199501_200512_climo',
         p=tmp_path_ne256, s='climo_180x360')

add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1',
         n='SCREAM 13-km tuned', comp='eamxx',
         htype='_ANN_199501_200512_climo',
         p=tmp_path_ne256, s='climo_180x360')

#-------------------------------------------------------------------------------
add_var('Q', eam_var='Q', eamxx_var='qv', fac=1e3, unit='g/kg', var_str='Sp. Humidity')
# add_var('U', eam_var='U', eamxx_var='U',     unit='m/s')
# add_var('T', eam_var='T', eamxx_var='T_mid', unit='K')

#-------------------------------------------------------------------------------
lev = np.array([ 50.,  70., 100., 125., 150., 175., 200., 225., 250., 300., 350., 400.,
                450., 500., 550., 600., 650., 700., 750., 775., 800., 825.,
                850., 875., 900., 925., 950., 975., 1000.])
lev = lev * 1e2   # Pa

#-------------------------------------------------------------------------------
fig_file = f'figs/FXX-clim-zonal-mean'

lat1, lat2 = -15, 15
yr1, yr2   = 1995, 2004

plot_diff = True

recalculate      = True
var_x_case       = False
num_plot_col     = 2
use_common_label_bar = False

#-------------------------------------------------------------------------------
def run_cmd(cmd, verbose=True, indent='  '):
    cmd_str = indent + hapy.tclr.GREEN + cmd + hapy.tclr.END
    if verbose: print('\n' + cmd_str)
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
    (msg, err) = proc.communicate()
    if verbose and msg != '': print(f'  msg: {msg}')
    if err != '' and not verbose: print(cmd_str)
    if err != '': print(f'err: {err}'); exit()
    return msg

#-------------------------------------------------------------------------------
num_case, num_var = len(case), len(var_list)

# Panel layout: [num_case rows, num_var cols] when var_x_case=False
nrows, ncols = (num_var, num_case) if var_x_case else (num_case, num_var)

if 'num_plot_col' in locals():
    nrows, ncols = ( int(np.ceil(num_var*num_case/float(num_plot_col))), num_plot_col )

# Figure size: target ~2048px at 150 dpi → ~13.6 inches per side
dpi     = 150
fig_w   = 10*ncols # max(5 * ncols, 8)
fig_h   =  6*nrows # max(4 * nrows, 6)
fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), dpi=dpi, squeeze=False)
fig.subplots_adjust(wspace=0.15, hspace=0.35)

plots   = [[None] * ncols for _ in range(nrows)]
cfs     = [[None] * ncols for _ in range(nrows)]   # contourf handles for colorbars

#-------------------------------------------------------------------------------
def get_row_col(c,v):
    if 'num_plot_col' in globals():
        ip = v*num_case+c if var_x_case else c*num_var+v
        row = int(np.floor(ip/num_plot_col))
        col = int(ip - row*ncols)
        # print(f'  ip / row / col : {ip} / {row} / {col}')
    else:
        row = c if not var_x_case else v
        col = v if not var_x_case else c
    return (row,col)
#-------------------------------------------------------------------------------
for v in range(num_var):
    print()
    print(f'  var: {hapy.tclr.MAGENTA}{var_list[v]}{hapy.tclr.END}')

    data_list = []
    lat_list_  = []
    lev_list_  = []

    for c in range(num_case):
        case_opts = opts_list[c]
        print('    case: ' + hapy.tclr.GREEN + case[c] + hapy.tclr.END)

        if recalculate:

            if case[c] == 'ERA5':
                obs_root = '/pscratch/sd/w/whannah/Obs/ERA5/e3sm_diags_remap'
                if var_list[v] == 'U': obs_var, input_file_name = 'ua',  f'{obs_root}/ua_197901_201912.remap_180x360.nc'
                if var_list[v] == 'V': obs_var, input_file_name = 'va',  f'{obs_root}/va_197901_201912.remap_180x360.nc'
                if var_list[v] == 'T': obs_var, input_file_name = 'ta',  f'{obs_root}/ta_197901_201912.remap_180x360.nc'
                if var_list[v] == 'Q': obs_var, input_file_name = 'hus', f'{obs_root}/hus_197901_201912.remap_180x360.nc'
                if var_list[v] == 'O3': obs_var, input_file_name = 'tro3', f'{obs_root}/tro3_197901_201912.remap_180x360.nc'

                ds   = xr.open_dataset(input_file_name)
                ds   = ds.where(ds['time.year'] >= yr1, drop=True)
                ds   = ds.where(ds['time.year'] <= yr2, drop=True)
                data = ds[obs_var].mean(dim='time')
                data = data.sel(plev=lev)

            else:
                case_dir = case_opts['p']
                case_sub = case_opts['s']
                htype    = case_opts['htype']

                file_path = f'{case_dir}/{case[c]}/{case_sub}/*{htype}*'
                file_list = sorted(glob.glob(file_path))

                ds = xr.open_mfdataset(file_list)
                if 'time' in ds.dims:
                    if len(ds.time) == 1: ds = ds.isel(time=0, drop=True)

                tvar = None
                if case_opts['comp'] == 'eam':   tvar = eam_var_list[v]
                if case_opts['comp'] == 'eamxx': tvar = eamxx_var_list[v]
                data = ds[tvar]

                data = hapy.vinth2p_simple(data, ds['hyam'], ds['hybm'], lev, ds['ps'],
                                           interp_type='linear', extrapolate=False)
                if len(data.plev) == 1: data = data.isel(plev=0)
        #-----------------------------------------------------------------------
        # yr, mn, dy = ds['time.year'].values, ds['time.month'].values, ds['time.day'].values
        # date1 = f'{yr[ 0]}-{mn[ 0]}-{dy[ 0]}'
        # date2 = f'{yr[-1]}-{mn[-1]}-{dy[-1]}'
        # print(f'\n  date range => {date1} : {date2}')
        # exit()
        #-----------------------------------------------------------------------
        # print(f'\n    latitude range => {data.lat[0].values} : {data.lat[-1].values}\n')
        #-----------------------------------------------------------------------
        if var_fac_list[v] is not None:
            data = data * var_fac_list[v]
        #-----------------------------------------------------------------------
        # Zonal mean
        data = data.mean(dim='lon')
        data_list.append(data.values)
        lat_list_.append(data['lat'].values)
        lev_list_.append(data['plev'].values / 1e2)  # Pa → hPa for display

    #---------------------------------------------------------------------------
    # Contour levels (common across cases for this variable)
    data_min = np.nanmin([np.nanmin(d) for d in data_list])
    data_max = np.nanmax([np.nanmax(d) for d in data_list])
    if data_min == data_max:
        raise ValueError(hapy.tclr.RED + 'WARNING: Difference is zero!' + hapy.tclr.END)

    clev_main = np.linspace(data_min, data_max, 11)
    # norm_main = BoundaryNorm(clev_main, ncolors=256, clip=False)

    # cmap = 'viridis'
    cmap = cmocean.cm.rain
    #---------------------------------------------------------------------------
    if plot_diff and num_case>1:
        for c in range(num_case):
            if c!=0: data_list[c] = data_list[c] - data_list[0]
        # diff_data_min = np.min([np.nanmin(d) for d in data_list[1:]])
        # diff_data_max = np.max([np.nanmax(d) for d in data_list[1:]])
        diff_data_max = np.max([np.nanmax(np.absolute(d)) for d in data_list[1:]])
        diff_data_min = -1*diff_data_max
        clev_diff = np.linspace(diff_data_min, diff_data_max, 12)
        # norm_diff = BoundaryNorm(clev_diff, ncolors=256, clip=False)
    #---------------------------------------------------------------------------
    # Draw panels
    for c in range(num_case):
        case_opts = opts_list[c]
        # row = c if not var_x_case else v
        # col = v if not var_x_case else c
        (row,col) = get_row_col(c,v)
        ax  = axes[row][col]

        img_kwargs = {}
        img_kwargs['cmap']      = cmap
        img_kwargs['levels']    = clev_main
        # img_kwargs['norm']      = norm_main
        # img_kwargs['vmin'] = data_min
        # img_kwargs['vmax'] = data_max

        if plot_diff and c!=0:
            img_kwargs['cmap']     = cmocean.cm.balance
            img_kwargs['levels']   = clev_diff
            # img_kwargs['cmap']     = 'bwr'# match E3SM diags
            # img_kwargs['levels']   = np.array([-1.5,-1.0,-0.5,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.5,1.0,1.5])
            # img_kwargs['norm']      = norm_diff
            # img_kwargs['vmin']   = diff_data_min
            # img_kwargs['vmax']   = diff_data_max

        lats = lat_list_[c]
        levs = lev_list_[c]
        dat  = np.ma.masked_invalid(data_list[c])

        cf = ax.contourf(lats, levs, dat, extend='both', **img_kwargs)
        cfs[row][col] = cf

        ax.set_ylim(levs.max(), levs.min())   # pressure increases downward
        ax.set_xlabel('Latitude',       fontsize=15)
        ax.set_ylabel('Pressure [hPa]', fontsize=15)
        ax.tick_params(labelsize=12)
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter())

        lstr = var_str_list[v]
        if unit_list[v] is not None: lstr = f'{lstr} [{unit_list[v]}]'

        if plot_diff and c!=0: lstr = f'{lstr} (diff)'

        # Subtitles: left = case name, right = variable name (replaces set_subtitles)
        ax.set_title(case_opts['n'], loc='left',  fontsize=20, pad=3)
        ax.set_title(lstr,           loc='right', fontsize=20, pad=3)

        # Panel letter label (a), b), …)
        panel_idx = row * ncols + col
        ax.text(0.02, 0.90, string.ascii_lowercase[panel_idx] + '',
            transform=ax.transAxes, fontsize=20, va='top', ha='left', #fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='black', pad=3))

    # -------------------------------------------------------------------------
    # Colorbar: one per variable column (or per row when var_x_case)
    if use_common_label_bar:
        # Attach colorbar below the bottom panel in this variable's column/row
        ref_ax = axes[-1][v] if not var_x_case else axes[v][-1]
        cbar = fig.colorbar(cfs[nrows - 1 if not var_x_case else v]
                            [v if not var_x_case else ncols - 1],
                            ax=axes[:, v] if not var_x_case else axes[v, :],
                            orientation='horizontal', pad=0.09, fraction=0.04,
                            aspect=30)
        cbar.ax.tick_params(labelsize=7)
        cbar.set_label(f'{var_list[v]} [{unit_list[v]}]', fontsize=8)
    else:
        for c in range(num_case):
            (row,col) = get_row_col(c,v)
            ax = axes[row][col]
            cbar = fig.colorbar(cfs[row][col], ax=ax,
                                orientation='horizontal', pad=0.2, fraction=0.07)
            cbar.ax.tick_params(labelsize=10)

# -------------------------------------------------------------------------
# os.makedirs('figs', exist_ok=True)
fig.savefig(f'{fig_file}.png', dpi=dpi, bbox_inches='tight')
# fig.savefig(f'{fig_file}.png', dpi=dpi)
print(f'\n{fig_file}.png')
# hapy.trim_png(fig_file)
plt.close(fig)