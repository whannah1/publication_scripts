# E3SMv2 coupled historical ensemble data

--------------------------------------------------------------------------------

## zstash commands to retreive data

```shell
ssh whannah@dtn02.nersc.gov
screen -r
bash
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```

```shell
CASE=v2.LR.historical_0101; zstash ls --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.1950-01*" 
CASE=v2.LR.historical_0151; zstash ls --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.1950-01*" 
CASE=v2.LR.historical_0201; zstash ls --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.1950-01*" 
CASE=v2.LR.historical_0251; zstash ls --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.1950-01*" 
CASE=v2.LR.historical_0301; zstash ls --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.1950-01*" 

CASE=v2.LR.historical_0101; mkdir /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}
CASE=v2.LR.historical_0151; mkdir /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}
CASE=v2.LR.historical_0201; mkdir /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}
CASE=v2.LR.historical_0251; mkdir /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}
CASE=v2.LR.historical_0301; mkdir /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}

CASE=v2.LR.historical_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.*" 
CASE=v2.LR.historical_0151; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.*" 
CASE=v2.LR.historical_0201; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.*" 
CASE=v2.LR.historical_0251; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.*" 
CASE=v2.LR.historical_0301; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h1.*" 

CASE=v2.LR.historical_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.oceanHeatContent.*" 
CASE=v2.LR.historical_0151; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.oceanHeatContent.*" 
CASE=v2.LR.historical_0201; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.oceanHeatContent.*" 
CASE=v2.LR.historical_0251; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.oceanHeatContent.*" 
CASE=v2.LR.historical_0301; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.oceanHeatContent.*" 

CASE=v2.LR.historical_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.19[56789]*" 
CASE=v2.LR.historical_0151; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.19[56789]*" 
CASE=v2.LR.historical_0201; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.19[56789]*" 
CASE=v2.LR.historical_0251; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.19[56789]*" 
CASE=v2.LR.historical_0301; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.19[56789]*"

CASE=v2.LR.historical_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.2*" 
CASE=v2.LR.historical_0151; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.2*" 
CASE=v2.LR.historical_0201; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.2*" 
CASE=v2.LR.historical_0251; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.2*" 
CASE=v2.LR.historical_0301; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*.mpaso.hist.am.meridionalHeatTransport.2*" 

CASE=v2.LR.historical_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "archive/ice/hist/*.mpassi.hist.am.timeSeriesStatsDaily.*"
CASE=v2.LR.historical_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_historical/${CASE}; zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "archive/ice/hist/*.mpassi.hist.am.timeSeriesStatsMonthly.*"

# note that MMF data is on OLCF Kronos - need to use globus => /nl/kronos/olcf/cli115/users/hannah6/hpss_xfer/
CASE=E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1; cd /global/cfs/cdirs/m3312/whannah/2023-CPL/${CASE}; zstash extract --hpss=/home/w/whannah/2023-CPL/${CASE} "archive/ice/hist/*.mpassi.hist.am.timeSeriesStatsMonthly.*"
# CASE=E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1; cd /global/cfs/cdirs/m3312/whannah/2023-CPL/${CASE}; zstash ls --hpss=/home/w/whannah/2023-CPL/${CASE} "archive/ice/hist/*"

CASE=v2.LR.piControl;         mkdir -p /global/cfs/cdirs/m3312/whannah/e3smv2_co2/${CASE} 
CASE=v2.LR.abrupt-4xCO2_0101; mkdir -p /global/cfs/cdirs/m3312/whannah/e3smv2_co2/${CASE} 

CASE=v2.LR.piControl;         cd /global/cfs/cdirs/m3312/whannah/e3smv2_co2/${CASE} ; nohup zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "archive/atm/hist/*eam.h0.*" > zstash_extract.out & 
CASE=v2.LR.abrupt-4xCO2_0101; cd /global/cfs/cdirs/m3312/whannah/e3smv2_co2/${CASE} ; nohup zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "archive/atm/hist/*eam.h0.*" > zstash_extract.out & 

zstash ls --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h*.1950-01-11*"
zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/${CASE} "*eam.h*.1950-01-11*"
```

--------------------------------------------------------------------------------

