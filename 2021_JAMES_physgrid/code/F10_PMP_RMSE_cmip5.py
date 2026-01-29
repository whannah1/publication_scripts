# Script to perfrom comparison with CMIP5 models
# source /global/project/projectdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh
import os,json
import xarray as xr,numpy as np
import matplotlib.pyplot as plt, matplotlib.cbook as cbook
from matplotlib.lines import Line2D

fig_file = 'figs/F10-RMSE-cmip5.png'

verbose = False

top_dir = './'

file_cmip5  = 'data/cmip5_json/cmip5/%s_2.5x2.5_regrid2_regrid2_metrics.json'
file_e3sm_h = 'data/cmip5_json/%s/%s_2.5x2.5_regrid2_regrid2_metrics.json'
file_e3sm_a = 'data/cmip5_json/%s/%s_2.5x2.5_esmf_linear_metrics.json'

seasons = ['ann','djf','mam','jja','son']

variables = [
   {'name':'rt',     'region':'global','title':'Net TOA (W m$^{-2}$)','seasons':seasons},
   {'name':'rstcre', 'region':'global','title':'SW CRE (W m$^{-2}$)','seasons':seasons},
   {'name':'rltcre', 'region':'global','title':'LW CRE (W m$^{-2}$)','seasons':seasons},
   {'name':'pr',     'region':'global','title':'prec (mm day$^{-1}$)','seasons':seasons},
   {'name':'tas',    'region':'land',  'title':'tas (land, K)','seasons':seasons},
   {'name':'tauu',   'region':'ocean', 'title':u'\N{GREEK SMALL LETTER TAU}$_x$ (ocean, Pa)','seasons':seasons},
   {'name':'ua-200', 'region':'global','title':'u200 (m s$^{-1}$)','seasons':seasons},
   {'name':'ua-850', 'region':'global','title':'u850 (m s$^{-1}$)','seasons':seasons},
   {'name':'zg-500', 'region':'global','title':'Zg-500 (m)','seasons':seasons},
   ### {'name':'ta-850', 'region':'global','title':'ta-850','seasons':seasons},
]

# physgrid cases
pg_case = ['E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-cbe53b',
          'E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b',
          'E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-cbe53b',
          'E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-cbe53b',
]

# clr = ['#ff0000','#00ff00','#00aa55','#0055aa','#0000ff']
# clr = ['#ff0000','#00ff00','#00aaaa','#0000ff','#aa00aa']
clr = ['#ff0000','#00ff00','#00ffff','#0000ff','#ff00ff']


scratch  = os.getenv('SCRATCH')

# Dictionaries to hold data
cmip5_models = set()
cmip5 = {}
e3sm_h = {}
e3sm_a = {}
e3sm_np4 = {}
e3sm_pg2 = {}
e3sm_pg3 = {}
e3sm_pg4 = {}

iv = 0
for v in variables:

  # Add entry for variable
  cmip5[v['name']] = {}
  e3sm_h[v['name']] = {}
  e3sm_a[v['name']] = {}
  e3sm_np4[v['name']] = {}
  e3sm_pg2[v['name']] = {}
  e3sm_pg3[v['name']] = {}
  e3sm_pg4[v['name']] = {}

  # CMIP5
  with open(file_cmip5 % (v['name']), 'r') as f:
    data = json.load(f)
  models = sorted(data['RESULTS'].keys())
  cmip5_models |= set(models)
  for m in models:
    for s in v['seasons']:
      try:
        rms = data['RESULTS'][m]['defaultReference']['r1i1p1'][v['region']]['rms_xy'][s]
        if verbose: print('%s,%s,%s,%s' % (v['name'],m,s,rms))
        if s not in cmip5[v['name']]:
           cmip5[v['name']][s] = []
        cmip5[v['name']][s].append(rms)
      except:
        print('%s,%s,%s,%s' % (v['name'],m,s,'missing'))
        pass

  # E3SM historical
  for H in ['H1','H2','H3','H4','H5']:
    with open(file_e3sm_h % (H,v['name']), 'r') as f:
      data = json.load(f)
    for m in sorted(data['RESULTS'].keys()):
      for s in v['seasons']:
        rms = data['RESULTS'][m]['defaultReference']['r1i1p1'][v['region']]['rms_xy_%s' % (s)]
        if verbose: print('%s,%s,%s,%s' % (v['name'],m,s,rms))
        if s not in e3sm_h[v['name']]:
           e3sm_h[v['name']][s] = []
        e3sm_h[v['name']][s].append(rms)

  # E3SM amip
  for H in ['A1',]:
    with open(file_e3sm_a % (H,v['name']), 'r') as f:
      data = json.load(f)
    for m in sorted(data['RESULTS'].keys()):
      for s in v['seasons']:
        rms = data['RESULTS'][m]['defaultReference']['r1i1p1'][v['region']]['rms_xy_%s' % (s)]
        if verbose: print('%s,%s,%s,%s' % (v['name'],m,s,rms))
        if s not in e3sm_a[v['name']]:
           e3sm_a[v['name']][s] = []
        e3sm_a[v['name']][s].append(rms)

  print()

  # file_prefix = 'pmp_rmse'
  file_prefix = 'pmp_rmse_alt' # includes area weighting in RMSE calculation

  ### E3SM physgrid cases

  pg_case = 'E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-cbe53b' 
  pg_rmse_file = f'data/PMP_RMSE/{pg_case}/pmp_rmse/{file_prefix}.{v.get("name")}.nc'
  ds = xr.open_dataset( pg_rmse_file )
  data = ds[v.get("name")]
  print(f'  {v.get("name")}  {pg_case:60}  {data.values}')
  for si,s in enumerate(seasons):
    if s not in e3sm_np4[v['name']]: e3sm_np4[v['name']][s] = []
    e3sm_np4[v['name']][s].append(data[si])

  pg_case = 'E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-cbe53b' 
  pg_rmse_file = f'data/PMP_RMSE/{pg_case}/pmp_rmse/{file_prefix}.{v.get("name")}.nc'
  ds = xr.open_dataset( pg_rmse_file )
  data = ds[v.get("name")]
  print(f'  {v.get("name")}  {pg_case:60}  {data.values}')
  for si,s in enumerate(seasons):
    if s not in e3sm_pg2[v['name']]: e3sm_pg2[v['name']][s] = []
    e3sm_pg2[v['name']][s].append(data[si])

  pg_case = 'E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-cbe53b' 
  pg_rmse_file = f'data/PMP_RMSE/{pg_case}/pmp_rmse/{file_prefix}.{v.get("name")}.nc'
  ds = xr.open_dataset( pg_rmse_file )
  data = ds[v.get("name")]
  print(f'  {v.get("name")}  {pg_case:60}  {data.values}')
  for si,s in enumerate(seasons):
    if s not in e3sm_pg3[v['name']]: e3sm_pg3[v['name']][s] = []
    e3sm_pg3[v['name']][s].append(data[si])

  pg_case = 'E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-cbe53b' 
  pg_rmse_file = f'data/PMP_RMSE/{pg_case}/pmp_rmse/{file_prefix}.{v.get("name")}.nc'
  ds = xr.open_dataset( pg_rmse_file )
  data = ds[v.get("name")]
  print(f'  {v.get("name")}  {pg_case:60}  {data.values}')
  for si,s in enumerate(seasons):
    if s not in e3sm_pg4[v['name']]: e3sm_pg4[v['name']][s] = []
    e3sm_pg4[v['name']][s].append(data[si])


  # Create plot
  
  # CMIP5 data for box and whiskers
  nd = len(cmip5[v['name']][v['seasons'][0]])
  ns = len(v['seasons'])
  data = np.zeros((nd,ns), 'float64')
  i = 0
  labels = []
  for s in v['seasons']:
     data[:,i] = np.array( cmip5[v['name']][s], 'float64' )
     labels.append(s.upper())
     i += 1
  cmip5_stats = cbook.boxplot_stats(data,whis=[0,100],labels=labels)
 
  # E3SM data
  nd = len(e3sm_h[v['name']][v['seasons'][0]])
  ns = len(v['seasons'])
  e3sm_h_data_x = np.zeros((nd,ns), 'float64')
  e3sm_h_data_y = np.zeros((nd,ns), 'float64')
  i = 0
  for s in v['seasons']:
     e3sm_h_data_x[:,i] = i + 1.4
     e3sm_h_data_y[:,i] = np.array( e3sm_h[v['name']][s], 'float64' )
     i += 1

  # E3SM data (amip)
  nd = len(e3sm_a[v['name']][v['seasons'][0]])
  ns = len(v['seasons'])
  e3sm_a_data_x = np.zeros((nd,ns), 'float64')
  e3sm_a_data_y = np.zeros((nd,ns), 'float64')
  i = 0
  for s in v['seasons']:
     e3sm_a_data_x[:,i] = i + 1.1
     e3sm_a_data_y[:,i] = np.array( e3sm_a[v['name']][s], 'float64' )
     i += 1

  # E3SM physgrid data
  nd = len(e3sm_np4[v['name']][v['seasons'][0]])
  ns = len(v['seasons'])
  e3sm_np4_data_x = np.zeros((nd,ns), 'float64')
  e3sm_pg2_data_x = np.zeros((nd,ns), 'float64')
  e3sm_pg3_data_x = np.zeros((nd,ns), 'float64')
  e3sm_pg4_data_x = np.zeros((nd,ns), 'float64')
  e3sm_np4_data_y = np.zeros((nd,ns), 'float64')
  e3sm_pg2_data_y = np.zeros((nd,ns), 'float64')
  e3sm_pg3_data_y = np.zeros((nd,ns), 'float64')
  e3sm_pg4_data_y = np.zeros((nd,ns), 'float64')
  i = 0
  for s in v['seasons']:
    e3sm_np4_data_x[:,i] = i + 1.2
    e3sm_pg2_data_x[:,i] = i + 1.3
    e3sm_pg3_data_x[:,i] = i + 1.4
    e3sm_pg4_data_x[:,i] = i + 1.5
    e3sm_np4_data_y[:,i] = np.array( e3sm_np4[v['name']][s], 'float64' )
    e3sm_pg2_data_y[:,i] = np.array( e3sm_pg2[v['name']][s], 'float64' )
    e3sm_pg3_data_y[:,i] = np.array( e3sm_pg3[v['name']][s], 'float64' )
    e3sm_pg4_data_y[:,i] = np.array( e3sm_pg4[v['name']][s], 'float64' )
    i += 1

  #-----------------------------------------------------------------------------
  # Modified MPL code
  #-----------------------------------------------------------------------------
  if v==variables[0]: 
    fig = plt.figure(figsize=[9,9])
    nsx,nsy = 3,3

  # Plot panel
  ax = plt.subplot(nsy, nsx, iv+1)
  ax.bxp(cmip5_stats,widths=np.array( [0.3]*len(seasons) ))
  s0 = ax.scatter(e3sm_a_data_x,e3sm_a_data_y,    color=clr[0])
  s1 = ax.scatter(e3sm_np4_data_x,e3sm_np4_data_y,color=clr[1])
  s2 = ax.scatter(e3sm_pg2_data_x,e3sm_pg2_data_y,color=clr[2])
  s3 = ax.scatter(e3sm_pg3_data_x,e3sm_pg3_data_y,color=clr[3])
  s4 = ax.scatter(e3sm_pg4_data_x,e3sm_pg4_data_y,color=clr[4])
  ax.set_title('('+chr(97+iv)+')', loc="left")
  ax.set_title(v['title'], loc="right")
  ax.set_xlim([0.4,ns+0.9])

  # if iv==0: cmap_list = [s0.cmap,s1.cmap,s2.cmap,s3.cmap,s4.cmap,]

  iv += 1



#-------------------------------------------------------------------------------
# Create the legend
#-------------------------------------------------------------------------------

legend_elements = [
  Line2D([0],[0],marker='o',color='w',label='G19 AMIP',markerfacecolor=clr[0],markersize=10),
  Line2D([0],[0],marker='o',color='w',label='np4',markerfacecolor=clr[1],markersize=10),
  Line2D([0],[0],marker='o',color='w',label='pg2',markerfacecolor=clr[2],markersize=10),
  Line2D([0],[0],marker='o',color='w',label='pg3',markerfacecolor=clr[3],markersize=10),
  Line2D([0],[0],marker='o',color='w',label='pg4',markerfacecolor=clr[4],markersize=10),
  ]
fig.legend(handles=legend_elements, loc='lower center', ncol=len(legend_elements))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
fig.subplots_adjust(wspace=0.3,hspace=0.3)
fig.savefig(fig_file,bbox_inches='tight', dpi=300)
#fig.savefig("test1.pdf")
   
# print('CMIP5 Models used for plot:')
# for i,c5mod in enumerate(cmip5_models): print(f'  {i}  {c5mod}')

print(f'\n  {fig_file} \n')
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
