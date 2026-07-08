
# NERSC data paths

```shell
# DY1
/global/cfs/cdirs/e3sm/jpaige3/dy1ne1024
# DY2
/global/cfs/cdirs/e3smdata/simulations/ecp-autotune/
```

# SCRIP grid creation

```shell
GRID_ROOT=/pscratch/sd/w/whannah/e3sm_scratch/files_grid
NE=256
GenerateCSMesh --alt --res ${NE} --file ${GRID_ROOT}/ne${NE}.g
GenerateVolumetricMesh --in ${GRID_ROOT}/ne${NE}.g --out ${GRID_ROOT}/ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${GRID_ROOT}/ne${NE}pg2.g --out ${GRID_ROOT}/ne${NE}pg2_scrip.nc
```