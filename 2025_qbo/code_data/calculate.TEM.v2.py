import os, ngl, copy, glob, xarray as xr, numpy as np
import hapy_common as hc
class tclr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
case,name,case_dir,case_sub,clr,dsh,mrk = [],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n); case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
##------------------------------------------------------------------------------
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01')

# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' )
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' )
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128')

# num_files = 32   # limit the number of files for debugging

#-------------------------------------------------------------------------------

scratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
# scratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'

debug = False

#---------------------------------------------------------------------------------------------------
# constants
H     = 7e3          # m         assumed mean scale heightof the atmosphere
P0    = 101325       # Pa        surface pressure
Rd    = 287.058      # J/kg/K    gas constant for dry air
cp    = 1004.64      # J/kg/K    specific heat for dry air
g     = 9.80665      # m/s       global average of gravity at MSLP
a     = 6.37123e6    # m         Earth's radius
omega = 7.29212e-5   # 1/s       Earth's rotation rate
pi    = 3.14159

#---------------------------------------------------------------------------------------------------
num_case = len(case)
for c in range(num_case):
   print(f'    case: {tclr.CYAN}{case[c]}{tclr.END}')

   idir = f'{scratch}/{case[c]}/data_remap_90x180_prs'
   odir = f'{scratch}/{case[c]}/data_remap_90x180_tem'

   if not os.path.exists(odir): os.mkdir(odir)

   file_path = f'{idir}/*eam.h2.*'
   file_list = sorted(glob.glob(file_path))

   # last file is empty, so remove it
   file_list.remove(file_list[-1]) 

   if 'num_files' in locals(): file_list = file_list[:num_files]

   # print(hc.tcolor.RED+'WARNING - using custom subset to fill in missing data!!!'+hc.tcolor.ENDC)
   # file_list = file_list[8*365:]

   if debug:
      print(hc.tcolor.RED+'WARNING! - only producing a single file for debugging!'+hc.tcolor.ENDC)

   #----------------------------------------------------------------------------
   # loop through files to calculate TEM terms 
   for f in file_list: 

      ds = xr.open_dataset(f)
      nlat = len(ds[ 'lat'].values)
      nlev = len(ds['plev'].values)

      dlat = ds['lat']
      rlat = np.deg2rad(dlat)
      rlat['lat'] = rlat
      cos_lat = np.cos(rlat)
      #-------------------------------------------------------------------------
      # basic zonal means and anomalies

      TH = ds['T'] * np.power( P0/ds['plev'], Rd/cp )
      
      TH_b = TH.mean(dim='lon')
      U_b  = ds['U'].mean(dim='lon')
      V_b  = ds['V'].mean(dim='lon')
      W_b  = ds['OMEGA'].mean(dim='lon')

      TH_p = TH          - TH_b
      U_p  = ds['U']     - U_b
      V_p  = ds['V']     - V_b
      W_p  = ds['OMEGA'] - W_b

      # make sure coordinate data is lat in radians for da.differentiate()
      TH_b['lat'] = rlat
      U_b ['lat'] = rlat
      V_b ['lat'] = rlat
      W_b ['lat'] = rlat

      TH_p['lat'] = rlat
      U_p ['lat'] = rlat
      V_p ['lat'] = rlat
      W_p ['lat'] = rlat

      #-------------------------------------------------------------------------
      # EP flux vectors

      dTHdp = TH_b.differentiate('plev')
      dUdp  =  U_b.differentiate('plev')
      dUdy  = (U_b*cos_lat).differentiate('lat') / (a*cos_lat)

      # eddy stream function
      gamma = (V_p*TH_p).mean(dim='lon') / dTHdp

      fcor = 2*omega*np.sin(rlat)

      ### original version based on Gerber and Manzini
      F_y = a*cos_lat * (       dUdp *gamma - (U_p*V_p).mean(dim='lon') )
      F_z = a*cos_lat * ( (fcor-dUdy)*gamma - (U_p*W_p).mean(dim='lon') )

      F_y = F_y.transpose('time','plev','lat')
      F_z = F_z.transpose('time','plev','lat')

      #-------------------------------------------------------------------------
      # EP flux divergence
      
      dFydy = (F_y*cos_lat).differentiate('lat') / (a*cos_lat)
      dFzdp = F_z.differentiate('plev')

      EP_div = ( dFydy + dFzdp ) / (a*cos_lat)

      #-------------------------------------------------------------------------
      # EP flux divergence transformed to log-pressure

      F_y_lp = F_y * ds['plev']/P0
      F_z_lp = F_z * -H/P0

      dFydy_lp = (F_y_lp*cos_lat).differentiate('lat') / (a*cos_lat)
      dFzdp_lp = F_z_lp.differentiate('plev')

      EP_div_lp = dFydy_lp + dFzdp_lp

      #-------------------------------------------------------------------------
      # TEM meridional and vertical  velocities

      dgamma_dy = (gamma*cos_lat).differentiate('lat') / (a*cos_lat)
      dgamma_dp = gamma.differentiate('plev')

      V_star = V_b - dgamma_dp
      W_star = W_b - dgamma_dy

      #-------------------------------------------------------------------------
      # TEM mass stream function

      dp_int = xr.full_like(ds['plev'],np.nan)
      for k in range(1,nlev-1):
         pint1 = ( ds['plev'][k-1] + ds['plev'][k-0] ) / 2.
         pint2 = ( ds['plev'][k-0] + ds['plev'][k+1] ) / 2.
         dp_int[k] =  pint1 - pint2

      gamma_mass = xr.full_like(U_b,np.nan)
      for k in range(1,nlev-1):
         tmp_integral = xr.full_like(U_b.mean(dim='plev'),0)
         for kk in range(1,k):
            tmp_integral[:,:] = tmp_integral[:,:] + ( V_b[:,kk,:]*dp_int[kk] - gamma[:,kk,] )
         gamma_mass[:,k,:] = ( (2*pi*a*cos_lat[:]/g) * tmp_integral[:,:] ).transpose('time','lat')

      #-------------------------------------------------------------------------
      # TEM northward and upward advection

      dUdt_y = V_star * ( fcor - dUdy )
      dUdt_z = -1 * W_star * dUdp

      #-------------------------------------------------------------------------
      # Create output datset and add variables      
      ds_out = xr.Dataset()

      ds_out['vtem']         = V_star
      ds_out['wtem']         = W_star
      ds_out['psitem']       = gamma_mass 
      ds_out['epfy']         = F_y
      ds_out['epfz']         = F_z
      ds_out['depfydy']      = dFydy 
      ds_out['depfzdz']      = dFzdp 
      # ds_out['depfydy']      = dFydy_lp 
      # ds_out['depfzdz']      = dFzdp_lp 
      ds_out['utendepfd']    = EP_div
      # ds_out['utendepfd_lp'] = EP_div_lp
      ds_out['utendvtem']    = dUdt_y
      ds_out['utendwtem']    = dUdt_z

      ds_out['lat'] = dlat # use latitude in degrees for output

      ds_out['u']            = ds['U' ].mean(dim='lon')
      ds_out['z']            = ds['Z3'].mean(dim='lon')
      ds_out['BUTGWSPEC']    = ds['BUTGWSPEC'].mean(dim='lon') # Beres U tendency - gravity wave spectrum  (convection)
      ds_out['UTGWSPEC' ]    = ds['UTGWSPEC' ].mean(dim='lon') # C&M U tendency - gravity wave spectrum    (frontogenesis)
      ds_out['UTGWORO'  ]    = ds['UTGWORO'  ].mean(dim='lon') # U tendency - orographic gravity wave drag

      #-------------------------------------------------------------------------
      # add variable long names
      ds_out['vtem']        .attrs['long_name'] = 'Transformed Eulerian mean northward wind'
      ds_out['wtem']        .attrs['long_name'] = 'Transformed Eulerian mean upward wind'
      ds_out['psitem']      .attrs['long_name'] = 'Transformed Eulerian mean mass stream function'
      # ds_out['epfy']        .attrs['long_name'] = 'Northward component of the Eliassen-Palm flux'
      # ds_out['epfz']        .attrs['long_name'] = 'Upward component of the Eliassen-Palm flux'
      # ds_out['depfydy']     .attrs['long_name'] = 'Meridional derivative of northward component of the Eliassen-Palm flux'
      # ds_out['depfzdp']     .attrs['long_name'] = 'Vertical derivative of upward component of the Eliassen-Palm flux'
      ds_out['utendepfd']   .attrs['long_name'] = 'Tendency of eastward wind due to TEM Eliassen-Palm flux divergence'
      # ds_out['utendepfd_lp'].attrs['long_name'] = 'Tendency of eastward wind due to TEM Eliassen-Palm flux divergence (log-pressure)'
      ds_out['utendvtem']   .attrs['long_name'] = 'Tendency of eastward wind due to TEM northward wind advection and coriolis'
      ds_out['utendwtem']   .attrs['long_name'] = 'Tendency of eastward wind due to TEM upward wind advection'
      ds_out['u']           .attrs['long_name'] = 'Zonal mean eastward wind'
      ds_out['z']           .attrs['long_name'] = 'Geopotential Height'

      #-------------------------------------------------------------------------
      # add variable units
      ds_out['vtem']        .attrs['units'] = 'm/s'
      ds_out['wtem']        .attrs['units'] = 'm/2'
      ds_out['psitem']      .attrs['units'] = 'kg/s'
      # ds_out['epfy']        .attrs['units'] = 'm3/s2'
      # ds_out['epfz']        .attrs['units'] = 'm3/s2'
      # ds_out['depfydy']     .attrs['units'] = 'm/s2'
      # ds_out['depfzdp']     .attrs['units'] = 'm/s2'
      ds_out['utendepfd']   .attrs['units'] = 'm/s2'
      # ds_out['utendepfd_lp'].attrs['units'] = 'm/s2'
      ds_out['utendvtem']   .attrs['units'] = 'm/s2'
      ds_out['utendwtem']   .attrs['units'] = 'm/s2'
      ds_out['u']           .attrs['units'] = 'm/s'
      ds_out['z']           .attrs['units'] = 'm'

      #-------------------------------------------------------------------------
      # write to file

      f_out = f.replace('.h2.','.h2.tem.').replace(idir,odir)

      print(' '*4+f'writing to file: {f_out}')
      ds_out.to_netcdf(path=f_out,mode='w')

      if debug:
         print(hc.tcolor.RED+'WARNING! - only producing a single file for debugging!'+hc.tcolor.ENDC)
         exit()

      #-------------------------------------------------------------------------

print(); print('done.'); print()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------