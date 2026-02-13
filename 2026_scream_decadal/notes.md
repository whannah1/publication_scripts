
# interactive job commands

```shell
salloc --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account=e3sm
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```

# Data Paths

## NERSC

```shell
# root of analysis code
~/Research/E3SM/pub_figs/2025_scream_decadal
# ne1024 data
/global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6
# ne256 data
/global/cfs/cdirs/e3smdata/simulations/scream-decadal/???
# v3LR data
/global/cfs/cdirs/m3312/whannah/e3smv3_amip/
# DYNANMO MJO hindcasts
/global/cfs/cdirs/e3smdata/simulations/scream-decadal/DYNAMO_MJO_hindcasts

```

## OLCF

```shell
# run directory
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run

```

# SCREAM data on NERSC

## zstash commands

```shell
cd /global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6/run
# zstash ls --hpss=/home/b/bhillma/scream-decadal/decadal-production-run6/run "output.scream.decadal.6hourlyAVG_ne30pg2/*"
# hsi get /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.1995-01-01-00000.nc
# hsi mget /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.1995*-00000.nc
# hsi mget /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.199[6-9]*-00000.nc
# hsi mget /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.200[0-4]*-00000.nc

cd /global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyINST_ne30pg2

# hsi mget /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyINST_ne30pg2/output.scream.decadal.6hourlyINST_ne30pg2.INSTANT.nhours_x6.199[5-9]*.nc
# hsi mget /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyINST_ne30pg2/output.scream.decadal.6hourlyINST_ne30pg2.INSTANT.nhours_x6.200[0-4]*.nc

hsi get /home/b/bhillma/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyINST_ne30pg2/output.scream.decadal.6hourlyINST_ne30pg2.INSTANT.nhours_x6.2001-07-17-21600.nc


cd /global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyAVG_ne30pg2

hsi get /home/b/bhillma/scream-decadal/decadal-production-run6-2001-07-04/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.2001-07-18-00000.nc

hsi get /home/b/bhillma/scream-decadal/decadal-production-run6-2002-12-01/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.2002-12-10-00000.nc

hsi get /home/b/bhillma/scream-decadal/decadal-production-run6-2003-06-04/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.2003-06-28-00000.nc

hsi get /home/b/bhillma/scream-decadal/decadal-production-run6-2004-05-05/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.2004-05-14-00000.nc

```

# v3LR AMIP ensemble data

## zstash commands

```shell
CASE=v3.LR.amip_0101; mkdir /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}
CASE=v3.LR.amip_0151; mkdir /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}
CASE=v3.LR.amip_0201; mkdir /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}

CASE=v3.LR.amip_0101; htype=h3; cd /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}; zstash extract --hpss=/home/w/wlin/E3SMv3/AMIP/${CASE} "archive/atm/hist/*eam.${htype}.199[4-9]*" "archive/atm/hist/*eam.${htype}.200[0-5]*" 
CASE=v3.LR.amip_0151; htype=h3; cd /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}; zstash extract --hpss=/home/w/wlin/E3SMv3/AMIP/${CASE} "archive/atm/hist/*eam.${htype}.199[4-9]*" "archive/atm/hist/*eam.${htype}.200[0-5]*" 
CASE=v3.LR.amip_0201; htype=h2; cd /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}; zstash extract --hpss=/home/w/wlin/E3SMv3/AMIP/${CASE} "archive/atm/hist/*eam.${htype}.199[4-9]*" "archive/atm/hist/*eam.${htype}.200[0-5]*" 

# for reproducing zstash problem
CASE=v3.LR.amip_0101; htype=h3; cd /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}; 

zstash extract --hpss=/home/w/wlin/E3SMv3/AMIP/v3.LR.amip_0101 "archive/atm/hist/v3.LR.amip_0101.eam.h3.1994-01-11-00000.nc"

```

```shell # commands to get restart files
CASE=v3.LR.amip_0201
cd /global/cfs/cdirs/m3312/whannah/e3smv3_amip/${CASE}
zstash ls --hpss=/home/w/wlin/E3SMv3/AMIP/${CASE} "archive/rest/*"
zstash extract --hpss=/home/w/wlin/E3SMv3/AMIP/${CASE} archive/rest/2023-01-01-00000/v3.LR.amip_0201.eam.i.2023-01-01-00000.nc

## v3LR history options from namelist

```

```
 avgflag_pertape = 'A','A','A','A','I','I'
 nhtfrq = 0,-24,-6,-3,-1,0
 mfilt  = 1,30,120,240,720,1

 fincl1 = 'AODALL','AODBC','AODDUST','AODPOM','AODSO4','AODSOA','AODSS','AODVIS',
          'CLDLOW','CLDMED','CLDHGH','CLDTOT',
          'CLDHGH_CAL','CLDLOW_CAL','CLDMED_CAL','CLD_MISR','CLDTOT_CAL',
          'CLMODIS','FISCCP1_COSP','FLDS','FLNS','FLNSC','FLNT','FLUT',
          'FLUTC','FSDS','FSDSC','FSNS','FSNSC','FSNT','FSNTOA','FSNTOAC','FSNTC',
          'ICEFRAC','LANDFRAC','LWCF','OCNFRAC','OMEGA','PRECC','PRECL','PRECSC','PRECSL','PS','PSL','Q',
          'QFLX','QREFHT','RELHUM','SCO','SHFLX','SOLIN','SWCF','T','TAUX','TAUY','TCO',
          'TGCLDLWP','TMQ','TREFHT','TREFMNAV','TREFMXAV','TS','U','U10','V','Z3',
          'dst_a1DDF','dst_a3DDF','dst_c1DDF','dst_c3DDF','dst_a1SFWET','dst_a3SFWET','dst_c1SFWET','dst_c3SFWET',
          'O3','LHFLX',
          'O3_2DTDA_trop','O3_2DTDB_trop','O3_2DTDD_trop','O3_2DTDE_trop','O3_2DTDI_trop','O3_2DTDL_trop',
          'O3_2DTDN_trop','O3_2DTDO_trop','O3_2DTDS_trop','O3_2DTDU_trop','O3_2DTRE_trop','O3_2DTRI_trop',
          'O3_SRF','NO_2DTDS','NO_TDLgt','NO2_2DTDD','NO2_2DTDS','NO2_TDAcf','CO_SRF','TROPE3D_P','TROP_P',
          'CDNUMC','SFDMS','so4_a1_sfgaex1','so4_a2_sfgaex1','so4_a3_sfgaex1','so4_a5_sfgaex1','soa_a1_sfgaex1',
          'soa_a2_sfgaex1','soa_a3_sfgaex1','GS_soa_a1','GS_soa_a2','GS_soa_a3','AQSO4_H2O2','AQSO4_O3',
          'SFSO2','SO2_CLXF','SO2','DF_SO2','AQ_SO2','GS_SO2','WD_SO2','ABURDENSO4_STR','ABURDENSO4_TRO',
          'ABURDENSO4','ABURDENBC','ABURDENDUST','ABURDENMOM','ABURDENPOM','ABURDENSEASALT',
          'ABURDENSOA','AODSO4_STR','AODSO4_TRO',
          'EXTINCT','AODABS','AODABSBC','CLDICE','CLDLIQ','CLD_CAL_TMPLIQ','CLD_CAL_TMPICE','Mass_bc_srf',
          'Mass_dst_srf','Mass_mom_srf','Mass_ncl_srf','Mass_pom_srf','Mass_so4_srf','Mass_soa_srf','Mass_bc_850',
          'Mass_dst_850','Mass_mom_850','Mass_ncl_850','Mass_pom_850','Mass_so4_850','Mass_soa_850','Mass_bc_500',
          'Mass_dst_500','Mass_mom_500','Mass_ncl_500','Mass_pom_500','Mass_so4_500','Mass_soa_500','Mass_bc_330',
          'Mass_dst_330','Mass_mom_330','Mass_ncl_330','Mass_pom_330','Mass_so4_330','Mass_soa_330','Mass_bc_200',
          'Mass_dst_200','Mass_mom_200','Mass_ncl_200','Mass_pom_200','Mass_so4_200','Mass_soa_200',
          'O3_2DTDD','O3_2DCIP','O3_2DCIL','CO_2DTDS','CO_2DTDD','CO_2DCEP','CO_2DCEL','NO_2DTDD',
          'FLNTC','SAODVIS',
          'H2OLNZ',
          'dst_a1SF','dst_a3SF',
          'PHIS','CLOUD','TGCLDIWP','TGCLDCWP','AREL',
          'CLDTOT_ISCCP','MEANCLDALB_ISCCP','MEANPTOP_ISCCP','CLD_CAL',
          'CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN',
          'CLDHGH_CAL_LIQ','CLDHGH_CAL_ICE','CLDHGH_CAL_UN',
          'CLDMED_CAL_LIQ','CLDMED_CAL_ICE','CLDMED_CAL_UN',
          'CLDLOW_CAL_LIQ','CLDLOW_CAL_ICE','CLDLOW_CAL_UN',
          'CLWMODIS','CLIMODIS'

 fincl2 = 'PS', 'FLUT','PRECT','U200','V200','U850','V850',
          'TCO','SCO','TREFHTMN:M','TREFHTMX:X','TREFHT','QREFHT'
 fincl3 = 'PS', 'PSL','PRECT','TUQ','TVQ','UBOT','VBOT','TREFHT','FLUT','OMEGA500','TBOT','U850','V850','U200','V200','T200','T500','Z700'
 fincl4 = 'PRECT'
 fincl5 = 'O3_SRF'
 fincl6 = 'CO_2DMSD','NO2_2DMSD','NO_2DMSD','O3_2DMSD','O3_2DMSD_trop'
```

# SCREAM output variables

## ncdump output

```shell
ncdump -h /global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6/run/output.scream.decadal.3hourlyINST_ne30pg2/output.scream.decadal.3hourlyINST_ne30pg2.INSTANT.nhours_x3.1994-10-04-00000.nc  | grep "float"

  float IceWaterPath(time, ncol) ;
  float LiqWaterPath(time, ncol) ;
  float MeridionalVapFlux(time, ncol) ;
  float PotentialTemperature_at_1000hPa(time, ncol) ;
  float PotentialTemperature_at_700hPa(time, ncol) ;
  float RelativeHumidity_at_300hPa(time, ncol) ;
  float RelativeHumidity_at_500hPa(time, ncol) ;
  float RelativeHumidity_at_700hPa(time, ncol) ;
  float RelativeHumidity_at_850hPa(time, ncol) ;
  float RelativeHumidity_at_925hPa(time, ncol) ;
  float SW_clrsky_flux_up_at_model_top(time, ncol) ;
  float SW_flux_dn_at_model_top(time, ncol) ;
  float SW_flux_up_at_model_top(time, ncol) ;
  float SeaLevelPressure(time, ncol) ;
  float T_2m(time, ncol) ;
  float T_mid_at_300hPa(time, ncol) ;
  float T_mid_at_500hPa(time, ncol) ;
  float T_mid_at_700hPa(time, ncol) ;
  float T_mid_at_850hPa(time, ncol) ;
  float T_mid_at_925hPa(time, ncol) ;
  float U(time, ncol, lev) ;
  float U_at_300hPa(time, ncol) ;
  float U_at_500hPa(time, ncol) ;
  float U_at_700hPa(time, ncol) ;
  float U_at_850hPa(time, ncol) ;
  float U_at_925hPa(time, ncol) ;
  float U_at_model_bot(time, ncol) ;
  float V(time, ncol, lev) ;
  float V_at_300hPa(time, ncol) ;
  float V_at_500hPa(time, ncol) ;
  float V_at_700hPa(time, ncol) ;
  float V_at_850hPa(time, ncol) ;
  float V_at_925hPa(time, ncol) ;
  float V_at_model_bot(time, ncol) ;
  float ZonalVapFlux(time, ncol) ;
  float cldfrac_liq(time, ncol, lev) ;
  float cldfrac_tot_for_analysis(time, ncol, lev) ;
  float omega_at_500hPa(time, ncol) ;
  float omega_at_700hPa(time, ncol) ;
  float omega_at_850hPa(time, ncol) ;
  float precip_ice_surf_mass_flux(time, ncol) ;
  float precip_liq_surf_mass_flux(time, ncol) ;
  float qv_2m(time, ncol) ;
  float wind_speed_10m(time, ncol) ;
  float z_mid_at_300hPa(time, ncol) ;
  float z_mid_at_500hPa(time, ncol) ;
  float z_mid_at_700hPa(time, ncol) ;
  float z_mid_at_850hPa(time, ncol) ;
  float z_mid_at_925hPa(time, ncol) ;
  

ncdump -h /global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6/run/output.scream.decadal.6hourlyAVG_ne30pg2/output.scream.decadal.6hourlyAVG_ne30pg2.AVERAGE.nhours_x6.1995-01-01-00000.nc  | grep "float"

  float LW_flux_dn(time, ncol, ilev) ;
  float LW_flux_dn_at_model_bot(time, ncol) ;
  float LW_flux_up(time, ncol, ilev) ;
  float LW_flux_up_at_model_bot(time, ncol) ;
  float LW_flux_up_at_model_top(time, ncol) ;
  float SW_flux_dn(time, ncol, ilev) ;
  float SW_flux_dn_at_model_bot(time, ncol) ;
  float SW_flux_dn_at_model_top(time, ncol) ;
  float SW_flux_up(time, ncol, ilev) ;
  float SW_flux_up_at_model_bot(time, ncol) ;
  float SW_flux_up_at_model_top(time, ncol) ;
  float T_2m(time, ncol) ;
  float U_at_model_bot(time, ncol) ;
  float V_at_model_bot(time, ncol) ;
  float homme_T_mid_tend(time, ncol, lev) ;
  float homme_qv_tend(time, ncol, lev) ;
  float p3_T_mid_tend(time, ncol, lev) ;
  float p3_qv_tend(time, ncol, lev) ;
  float precip_ice_surf_mass_flux(time, ncol) ;
  float precip_liq_surf_mass_flux(time, ncol) ;
  float ps(time, ncol) ;
  float rrtmgp_T_mid_tend(time, ncol, lev) ;
  float shoc_T_mid_tend(time, ncol, lev) ;
  float shoc_qv_tend(time, ncol, lev) ;
  float surf_evap(time, ncol) ;
  float surf_mom_flux(time, ncol, dim2) ;
  float surf_radiative_T(time, ncol) ;
  float surf_sens_flux(time, ncol) ;
  float surface_upward_latent_heat_flux(time, ncol) ;


```

## yaml files

```shell
> ls -1  /lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/*
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/namelist.nl
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_default_output.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_input.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.15minINST_ARM.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.1dailyAVG_ne1024pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.1dailyMAX_ne1024pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.1dailyMIN_ne1024pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.1hourlyINST_ARM.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.1hourlyINST_ne1024pg2.yaml

/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.3hourlyAVG_ne30pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.3hourlyINST_ne30pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.6hourlyAVG_ne30pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.6hourlyINST_ne30pg2.yaml

/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.6hourlyINST_ne1024pg2.yaml
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.monthlyAVG_ne30pg2.yaml
```

### scream_output.decadal.3hourlyAVG_ne30pg2.yaml

```shell
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.3hourlyAVG_ne30pg2.yaml

    - T_mid
    - qv
    - qc
    - qr
    - qi
    - cldfrac_tot
    - cldfrac_liq
    - omega
    - U
    - V
    - pseudo_density
    - z_mid
    - surf_sens_flux
    - surf_evap
    - surface_upward_latent_heat_flux
    - ps
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - surf_mom_flux
    - surf_radiative_T
    - T_2m
    - sfc_flux_dir_nir
    - sfc_flux_dir_vis
    - sfc_flux_dif_nir
    - sfc_flux_dif_vis
    - sfc_flux_sw_net
    - sfc_flux_lw_dn
    - U_at_model_bot
    - V_at_model_bot
    - T_mid_at_model_bot
    - qv_at_model_bot
    - qc_at_model_bot
    - qi_at_model_bot
    - qr_at_model_bot
    - qm_at_model_bot
    - bm_at_model_bot
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
    - SW_clrsky_flux_up_at_model_top
    - LW_clrsky_flux_up_at_model_top
    - SW_flux_dn_at_model_bot
    - SW_clrsky_flux_dn_at_model_bot
    - SW_flux_up_at_model_bot
    - SW_clrsky_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - LW_clrsky_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - LongwaveCloudForcing
    - ShortwaveCloudForcing
```

### scream_output.decadal.3hourlyINST_ne30pg2.yaml

```shell
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.3hourlyINST_ne30pg2.yaml

    - SeaLevelPressure
    - U_at_model_bot
    - V_at_model_bot
    - wind_speed_10m
    - U_at_925hPa
    - V_at_925hPa
    - U_at_850hPa
    - V_at_850hPa
    - U_at_700hPa
    - V_at_700hPa
    - U_at_500hPa
    - V_at_500hPa
    - U_at_300hPa
    - V_at_300hPa
    - ZonalVapFlux
    - MeridionalVapFlux
    - T_2m
    - qv_2m
    - T_mid_at_925hPa
    - T_mid_at_850hPa
    - T_mid_at_700hPa
    - T_mid_at_500hPa
    - T_mid_at_300hPa
    - RelativeHumidity_at_925hPa
    - RelativeHumidity_at_850hPa
    - RelativeHumidity_at_700hPa
    - RelativeHumidity_at_500hPa
    - RelativeHumidity_at_300hPa
    - z_mid_at_925hPa
    - z_mid_at_850hPa
    - z_mid_at_700hPa
    - z_mid_at_500hPa
    - z_mid_at_300hPa
    - SW_clrsky_flux_up_at_model_top
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LiqWaterPath
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - IceWaterPath
    - cldfrac_tot_for_analysis
    - cldfrac_liq
    - omega_at_500hPa
    - omega_at_700hPa
    - omega_at_850hPa
    - PotentialTemperature_at_700hPa
    - PotentialTemperature_at_1000hPa
    - U
    - V
```

### scream_output.decadal.6hourlyAVG_ne30pg2.yaml

```shell
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.6hourlyAVG_ne30pg2.yaml

    - p3_T_mid_tend
    - shoc_T_mid_tend
    - rrtmgp_T_mid_tend
    - homme_T_mid_tend
    - p3_qv_tend
    - shoc_qv_tend
    - homme_qv_tend
    - SW_flux_dn
    - SW_flux_up
    - LW_flux_dn
    - LW_flux_up
    - surf_sens_flux
    - surf_evap
    - surface_upward_latent_heat_flux
    - ps
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - surf_mom_flux
    - surf_radiative_T
    - T_2m
    - U_at_model_bot
    - V_at_model_bot
    - SW_flux_dn_at_model_bot
    - SW_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
```

### scream_output.decadal.6hourlyINST_ne30pg2.yaml

```shell
/lustre/orion/cli115/proj-shared/brhillman/e3sm_scratch/decadal-production-20240305.ne1024pg2_ne1024pg2.F20TR-SCREAMv1.run1/run/data/scream_output.decadal.6hourlyINST_ne30pg2.yaml

    - T_mid
    - qv
    - RelativeHumidity
    - U
    - V
    - omega
    - qc
    - nc
    - qr
    - qi
    - tke
    - o3_volume_mix_ratio
    - VapWaterPath
    - LiqWaterPath
    - IceWaterPath
    - surf_radiative_T
    - ps
    - qv_2m
    - T_2m
    - ocnfrac
    - landfrac
```


# Mapping

```shell
ncremap -G ttl='Equi-Angular grid 2x2 degree, dimensions 90x180, cell edges on Poles/Equator and Prime Meridian/Date Line'#latlon==90,180#lat_typ=uni#lon_typ=grn_wst -g files_grid/cmip6_90x180_scrip.20250624.nc


SRC_GRID=files_grid/ne30pg2_scrip.nc
DST_GRID=files_grid/cmip6_90x180_scrip.20250624.nc
MAP_FILE=files_maps/map_ne30pg2_to_90x180_traave.20250624.nc
ncremap -a traave --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}

```

## remap IMERG to NOAA OLR grid

```shell
SRC_DATA=/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/IMERG_Daily/PRECT_200101_202012.nc
# DST_DATA=/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/IMERG_Daily/PRECT_200101_202012.remap_73x144.nc
DST_DATA=/global/cfs/cdirs/m3312/whannah/obs_data/IMERG/IMERG_Daily_PRECT_200101_202012.remap_73x144.nc
MAP_FILE=files_maps/map_180x360_to_73x144_traave.20260106.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
```

## remap IMERG to ne30pg2

```shell

SRC_GRID=files_grid/cmip_180x360_scrip.nc
DST_GRID=files_grid/ne30pg2_scrip.nc
MAP_FILE=files_maps/map_180x360_to_ne30pg2_traave.20260112.nc
ncremap -a traave --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}

SRC_DATA=/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/IMERG_Daily/PRECT_200101_202012.nc
DST_DATA=/global/cfs/cdirs/m3312/whannah/obs_data/IMERG/IMERG_Daily_PRECT_200101_202012.remap_ne30pg2.nc
MAP_FILE=files_maps/map_180x360_to_ne30pg2_traave.20260112.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
```

## remap NOAA to ne30pg2

```shell

SRC_GRID=files_grid/73x144_scrip.nc
DST_GRID=files_grid/ne30pg2_scrip.nc
MAP_FILE=files_maps/map_73x144_to_ne30pg2_traave.20260112.nc
ncremap -a traave --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}

SRC_DATA=
DST_DATA=
MAP_FILE=files_maps/map_73x144_to_ne30pg2_traave.20260112.nc
ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}
```
