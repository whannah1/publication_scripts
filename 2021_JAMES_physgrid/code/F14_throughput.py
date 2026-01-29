#!/usr/bin/env python
import sys, os, fileinput, re, subprocess as sp, glob
import hapy_common as hc, hapy_setres as hs
import ngl, numpy as np, xarray as xr
# Set up terminal colors
class bcolor:
   ENDC    = '\033[0m';  BLACK  = '\033[30m'; RED   = '\033[31m'  
   GREEN   = '\033[32m'; YELLOW = '\033[33m'; BLUE  = '\033[34m'
   MAGENTA = '\033[35m'; CYAN   = '\033[36m'; WHITE = '\033[37m'
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


### physgrid validation
name,cases,clr,git_hash = [],[],[],'cbe53b'
cases.append(f'E3SM.PGVAL.ne30_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')         ; clr.append('black')
cases.append(f'E3SM.PGVAL.ne30pg2_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')      ; clr.append('black')
cases.append(f'E3SM.PGVAL.ne30pg3_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')      ; clr.append('black')
cases.append(f'E3SM.PGVAL.ne30pg4_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')      ; clr.append('black')
cases.append(f'E3SM.PGVAL.conusx4v1_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}')    ; clr.append('blue')
cases.append(f'E3SM.PGVAL.conusx4v1pg2_r05_oECv3.F2010SC5-CMIP6.master-{git_hash}') ; clr.append('blue')
for c in cases:
   if 'E3SM.PGVAL.ne30_r05_oECv3'         in c: name.append('ne30np4')
   if 'E3SM.PGVAL.ne30pg2_r05_oECv3'      in c: name.append('ne30pg2')
   if 'E3SM.PGVAL.ne30pg3_r05_oECv3'      in c: name.append('ne30pg3')
   if 'E3SM.PGVAL.ne30pg4_r05_oECv3'      in c: name.append('ne30pg4')
   if 'E3SM.PGVAL.conusx4v1_r05_oECv3'    in c: name.append('RRM np4')
   if 'E3SM.PGVAL.conusx4v1pg2_r05_oECv3' in c: name.append('RRM pg2')



param_list = []
param_list.append('Throughput')
# param_list.append('physics_fraction')
param_list.append('dynamics_cost')
param_list.append('physics_cost')
# param_list.append('Run Time    :')
# param_list.append('ATM Run Time')
# param_list.append('\"a:crm\"')
# param_list.append('a:CAM_run1')
# param_list.append('a:CAM_run2')
# param_list.append('a:CAM_run3')
# param_list.append('a:dyn_run')
# param_list.append('a:ac_physics')
# param_list.append('a:bc_physics')

print_data = False

min_num_days = 4

#-------------------------------------------------------------------------------
# setup plot stuff
#-------------------------------------------------------------------------------
# num_param = len(param)

fig_file = 'figs/F14-throughput'

wks = ngl.open_wks('png',fig_file)
plot = []
res = hs.res_xy()
res.vpHeightF = 0.3
res.xyMarkLineMode = 'Markers'


lres = hs.res_xy()
lres.xyLineThicknessF      = 1
lres.xyDashPattern         = 0
lres.xyLineColor           = 'gray'


# pmres = hs.res_default()
pmres  = ngl.Resources()
pmres.gsMarkerSizeF       = 5.0 
pmres.gsMarkerColor       = 'red'
# pmres.gsMarkerIndex       = 15 # circle with X
pmres.gsMarkerIndex       = 16
pmres.gsMarkerThicknessF  = 2

#-------------------------------------------------------------------------------
# make sure all logs are unzipped
#-------------------------------------------------------------------------------
max_case_len = 0
for case in cases:
   # timing_dir = os.getenv('HOME')+f'/E3SM/Cases/{case}/timing'
   timing_dir = f'timing_data/{case}/'
   timing_stat_gz_path = f'{timing_dir}/*.gz'
   if len(glob.glob(timing_stat_gz_path))>0: os.system(f'gunzip {timing_stat_gz_path} ')
   # find max char len for case name
   if len(case)>max_case_len: max_case_len = len(case)

#-------------------------------------------------------------------------------
# Loop through parameters and cases
#-------------------------------------------------------------------------------
davg_list_list = []
for param in param_list :
   if print_data: print()
   data_list,dpos_list,davg_list = [],[],[]
   for c,case in enumerate(cases):
      if print_data: print()
      timing_dir = os.getenv('HOME')+f'/E3SM/Cases/{case}/timing'
      
      case_name_fmt = bcolor.CYAN+f'{case:{max_case_len}}'+bcolor.ENDC

      # check that the timing files exist
      timing_file_path = f'{timing_dir}/*'
      if len(glob.glob(timing_file_path))==0: 
         print(f'{case_name_fmt}  No files!')
         continue

      tparam = param
      if param=='physics_fraction': tparam = 'a:CAM_run'
      # if param=='physics_fraction': tparam = 'a:[p,d][h,y].*_run'
      # if param=='physics_cost':     tparam = 'a:CAM_run'
      # if param=='dynamics_cost':    tparam = 'a:CAM_run3'
      # if param=='physics_cost':     tparam = 'a:phys_run'
      if param=='physics_cost':     tparam = 'a:[a,b]c_physics'
      if param=='dynamics_cost':    tparam = 'a:dyn_run'
      #-------------------------------------------------------------------------
      # grep for the appendropriate line in the stat files
      cmd = 'grep  \''+tparam+f'\'  {timing_dir}/*'
      proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True, universal_newlines=True)
      (msg, err) = proc.communicate()
      msg = msg.split('\n')

      ### check what grep found
      # for m in msg: print(m)
      # exit()

      stat_file = msg[0].split(':')[0]

      #-------------------------------------------------------------------------
      # Clean up message but don't print yet
      for m in range(len(msg)): 
         line = msg[m]
         line = line.replace(timing_dir+'/','')
         line = line.replace(f'e3sm_timing.{case}.','e3sm_timing        ')
         line = line.replace('e3sm_timing_stats.'  ,'e3sm_timing_stats  ')
         # Add case name
         if line!='': line = case_name_fmt+' '+line
         line = line.replace(f'e3sm_timing      ','')
         msg[m] = line
      #-------------------------------------------------------------------------
      # print stat file header with indentation
      if 'a:' in tparam : 
         head = sp.check_output(['head',stat_file],universal_newlines=True).split('\n')
         for line in head: 
            if 'walltotal' in line:
               indent = len(msg[0].split(':')[0])+1
               line = ' '*indent+line
               hline = line
               # Get rid of some dead space
               line = line.replace('name        ','name')
               if print_data: print(hline)
               break
      #-------------------------------------------------------------------------
      # set up character indices for color
      if 'a:' in tparam:
         
         # use for avg time per task
         n0 = hline.find('on')           +len('on')
         n1 = hline.find('processes')    +len('processes')
         n2 = hline.find('count')        +len('count')
         n3 = hline.find('walltotal')    +len('walltotal')
         n4 = hline.find('wallmax')      +len('wallmax')
         n5 = hline.find('wallmax')      +len('wallmax (proc   thrd  )')
         n6 = hline.find('wallmin')      +len('wallmin')
         
         # use for wallmax
         # n1 = hline.find('walltotal')    +len('walltotal')
         # n2 = hline.find('wallmax')      +len('wallmax')
         # n3 = hline.find('wallmax')      +len('wallmax (proc   thrd  )')
         # n4 = hline.find('wallmin')      +len('wallmin')
      elif 'ATM Run Time' in tparam:
         line = msg[0]
         n1 = line.replace(':','', 1).find(':')+2
         num_in_list = re.findall(r'\d+\.\d+', line[n1:])
         n1 = line.find(num_in_list[2])
         n2 = line.find(num_in_list[2])+len(num_in_list[2])
      else:
         line = msg[0]
         n1 = line.replace(':','', 1).find(':')+2
         num_in_list = re.findall(r'\d+\.\d+', line[n1:])
         n2 = line.find(num_in_list[0])+len(num_in_list[0])

      #-------------------------------------------------------------------------
      # figure out run length for each file
      # if param=='physics_cost':
      num_sim_days = []

      cmd = f'grep  \'run length  :\'  {timing_dir}/*'
      proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True, universal_newlines=True)
      (tmp_msg, err) = proc.communicate()
      tmp_msg = tmp_msg.split('\n')

      # Clean up grep message
      for m in range(len(tmp_msg)): 
         line = tmp_msg[m]
         line = line.replace(f'{timing_dir}/e3sm_timing.{case}.','')
         tmp_msg[m] = line

      for line in tmp_msg[-len(tmp_msg):] : 
         if line!='':
            # num_in_list = re.findall(r'\d+\d+', line)
            # num_sim_days.append( float(num_in_list[3]) )
            num_in_list = re.findall(r'\d+\.\d+', line)
            num_sim_days.append( float(num_in_list[-1]) )
      
      # print(num_sim_days)

      # # print number of days
      # for num_days in num_sim_days: 
      #    print(f'  num_days: {num_days}')

      # # compare the counts of num_days to timing data
      # mcnt = 0
      # for m in msg[-len(msg):]: 
      #    if m=='': continue
      #    mcnt += 1
      # print(len(num_sim_days))
      # print(mcnt)
      # exit()
      #-------------------------------------------------------------------------
      # print the timing data
      num_file = -len(msg)
      data,dpos = [],[]
      frac_cnt,phys_cost,tot_cost = 0,0,0
      nday_cnt = 0
      dcnt = 0
      for l,line in enumerate(msg[num_file:]) : 
         if line=='': continue
         data_tmp = None
         if param=='physics_fraction': 
            if frac_cnt==0 and print_data: print()
            if any([s in line for s in ['a:CAM_run1','a:CAM_run2','a:CAM_run3']]):
            # if any([s in line for s in ['a:phys_run1','a:phys_run2','a:dyn_run']]):
               # tot_cost += float(line[n1:n2])
               processes = float(line[n0:n1])
               walltime = float(line[n2:n3])  # wall total
               tot_cost += walltime / processes
               # tot_cost += float(line[n3:n4])  # wall max
            if any([s in line for s in ['a:CAM_run1','a:CAM_run2']]):
            # if any([s in line for s in ['a:phys_run1','a:phys_run2']]):
               # phys_cost += float(line[n1:n2])
               processes = float(line[n0:n1])
               walltime = float(line[n2:n3])  # wall total
               phys_cost += walltime / processes
               # phys_cost += float(line[n3:n4])  # wall max
            frac_cnt += 1
            if frac_cnt==4:
            # if frac_cnt==3:
               data_tmp = phys_cost / tot_cost
               frac_cnt,phys_cost,tot_cost = 0,0,0
         elif param=='physics_cost': 
            if frac_cnt==0 and print_data: print()
            # if any([s in line for s in ['a:CAM_run1','a:CAM_run2']]):
            # if any([s in line for s in ['a:phys_run1','a:phys_run2']]):
            if any([s in line for s in ['a:ac_physics','a:bc_physics']]):
               processes = float(line[n0:n1])
               walltotal = float(line[n2:n3])
               phys_cost += walltotal / processes
            frac_cnt += 1
            # if frac_cnt==4:
            if frac_cnt==2:
               # convert to sypd
               phys_cost = num_sim_days[dcnt]/phys_cost*(3600*24)/365
               data_tmp = phys_cost
               frac_cnt,phys_cost,tot_cost = 0,0,0
         elif param=='dynamics_cost': 
            # convert to sypd
            processes = float(line[n0:n1])
            walltotal = float(line[n2:n3])
            data_tmp = num_sim_days[dcnt]/(walltotal/processes)*(3600*24)/365
         elif 'a:' in param:
            processes = float(line[n0:n1])
            walltotal = float(line[n2:n3])
            data_tmp =  walltotal / processes
         else:
            data_tmp = float(line[n1:n2])

         if data_tmp is not None:
            if num_sim_days[dcnt]>min_num_days:
               data.append( data_tmp )
               dpos.append(float(c))

            dcnt += 1

         #----------------------------------------------------------------------
         # format and print the data
         if print_data: 
            if 'a:' in tparam : 
               line = line[:n1] \
                    +bcolor.CYAN  +line[n1:n2]+bcolor.CYAN +line[n2:n3] \
                    +bcolor.GREEN +line[n3:n4]+bcolor.ENDC +line[n4:]
               # Get rid of some dead space
               line = line.replace('        ','')
            else:
               line = line[:n1] \
                    +bcolor.GREEN +line[n1:n2]+bcolor.ENDC \
                    +line[n2:]
               # Print conversion to min aand hours for specific params
               offset = len(bcolor.GREEN)
               if tparam=='Run Time    :' and line[n1+offset:n2+offset]!='' :
                  sec = float( line[n1+offset:n2+offset] )
                  line = line+'  ('+bcolor.GREEN+f'{sec/60:.2f}'     +bcolor.ENDC+' min)'
                  line = line+'  ('+bcolor.GREEN+f'{sec/60/60:.2f}'  +bcolor.ENDC+' hour)'
            # print the line
            print(line)
         #----------------------------------------------------------------------

      data = np.array(data)
      avg = np.average(data)
      data_list.append(data)
      dpos_list.append(np.array(dpos))
      davg_list.append(avg)
      if print_data: print(f'average: {bcolor.GREEN}{avg:6.3}{bcolor.ENDC}')

   #-------------------------------------------------------------------------------

   apx = np.arange(0,len(davg_list))
   apy = np.array(davg_list)
   
   #-------------------------------------------------------------------------------
   # plot timing data

   # res.tiXAxisString = xvar
   # res.tiYAxisString = yvar

   if param in ['Throughput','physics_cost','dynamics_cost']:
      res.tiYAxisString = 'Throughput [sypd]'
   elif param in ['physics_fraction']:
      res.tiYAxisString = ''
   else:
      res.tiYAxisString = ''

   # limit the number of data points - make sure number of points is consistent
   max_list_length = 100
   for l,tlist in enumerate(dpos_list):
      max_list_length = np.minimum(max_list_length,len(tlist))
   for l,tlist in enumerate(dpos_list):
      if len(dpos_list[l])>max_list_length:
         dpos_list[l] = dpos_list[l][:max_list_length]
         data_list[l] = data_list[l][:max_list_length]

   px = xr.DataArray(np.array(dpos_list)).stack().values
   py = xr.DataArray(np.array(data_list)).stack().values

   res.trXMinF = np.min(px)-0.25
   res.trXMaxF = np.max(px)+0.25

   dy = np.max(py) - np.min(py)
   
   res.trYMinF = 0
   # res.trYMinF = np.min(py)-0.1*dy
   res.trYMaxF = np.max(py)+0.1*dy

   res.tmXBMode = 'Explicit'
   res.tmXBValues = np.arange(0,len(name))
   res.tmXBLabels = name

   # res.tiYAxisString = '[sypd]'

   # # plot individual timing estimates
   # res.xyMarker = 16
   # res.xyMarkerSizeF = 0.005
   # res.xyMarkerColor = 'gray'
   # plot.append( ngl.xy(wks, px, py, res)  )


   if param=='Throughput'   : lmax,dl1,dl2 =  6,1,10
   if param=='dynamics_cost': lmax,dl1,dl2 = 36,5,10
   if param=='physics_cost' : lmax,dl1,dl2 = 20,2,10
   res.tmYLMode = 'Explicit'
   if param=='Throughput'   : res.tmYLValues = np.arange(0, 6,dl1)
   if param=='dynamics_cost': res.tmYLValues = np.arange(0,36,dl1)
   if param=='physics_cost' : res.tmYLValues = np.arange(0,20,dl1)
   res.tmYLLabels = res.tmYLValues
   res.tmYLLabelFontHeightF = 0.010

   # plot average value
   res.xyMarker = 16 # large filled dot
   # res.xyMarker = 4  # open circle
   res.xyMarkerSizeF = 0.012
   res.xyMarkerColor = 'black'
   # ngl.overlay( plot[len(plot)-1], ngl.xy(wks,px,py,res) )
   plot.append( ngl.xy(wks, apx, apy, res)  )


   # Add horizontal lines
   xx = np.array([-1e2,1e2])
   lres.xyLineThicknessF = 2
   for hline_y in range(0,50,dl2):
      yy = np.array([float(hline_y),float(hline_y)])
      ngl.overlay( plot[len(plot)-1], ngl.xy(wks, xx, yy, lres) )

   # Add more horizontal lines
   xx = np.array([-1e2,1e2])
   lres.xyLineThicknessF = 1
   for hline_y in range(0,50,dl1):
      yy = np.array([float(hline_y),float(hline_y)])
      ngl.overlay( plot[len(plot)-1], ngl.xy(wks, xx, yy, lres) )


   for c,case in enumerate(cases):
      res.xyMarkerColor = clr[c]
      ngl.overlay( plot[len(plot)-1], ngl.xy(wks, [apx[c],apx[c]], [apy[c],apy[c]], res) )

   # overlay individual timing estimates
   res.xyMarker = 16
   res.xyMarkerSizeF = 0.005
   res.xyMarkerColor = 'gray'
   ngl.overlay( plot[len(plot)-1], ngl.xy(wks,px,py,res) )

   # # Add horizontal lines
   # ngl.overlay( plot[len(plot)-1], ngl.xy(wks,[0,0],[-1e8,1e8],lres) )
   # ngl.overlay( plot[len(plot)-1], ngl.xy(wks,[-1e8,1e8],[0,0],lres) )
   
   param_name = ''
   if param=='Throughput':       param_name = 'Overall Model Throughput'
   if param=='ATM Run Time':     param_name = 'Atmosphere Model Throughput [sypd]'
   if param=='physics_fraction': param_name = 'Proportional Cost of Physics'
   # if param=='physics_cost':     param_name = 'Throughput of Physics (excl. I/O and mapping)'
   if param=='physics_cost':     param_name = 'Physics Throughput'
   if param=='dynamics_cost':    param_name = 'Dynamics Throughput'
   
   subtitle_font_height = 0.025
   if len(param_list)==2: subtitle_font_height = 0.015
   if len(param_list)==3: subtitle_font_height = 0.010
   hs.set_subtitles(wks, plot[len(plot)-1], '', param_name, '', font_height=subtitle_font_height)

   # create list for summary table
   davg_list_list.append(davg_list)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
print()
max_case_len = 0
for c,case in enumerate(cases):
   if len(case)>max_case_len: max_case_len = len(case)

for p,param in enumerate(param_list) :
   davg_list = davg_list_list[p]
   print(f'{param}  (averages)')
   for c,case in enumerate(cases):
      avg = davg_list[c]
      print(f'  {case:{max_case_len}} : {bcolor.GREEN}{avg:6.3}{bcolor.ENDC}')

      

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
pres = hs.setres_panel()
pres.nglPanelFigureStrings            = ['a','b','c','d','e','f','g','h']
pres.nglPanelFigureStringsJust        = "TopLeft"
if len(param_list)==1: pres.nglPanelFigureStringsFontHeightF = 0.015
if len(param_list)==2: pres.nglPanelFigureStringsFontHeightF = 0.015
if len(param_list)>=3: pres.nglPanelFigureStringsFontHeightF = 0.01
# pres.nglPanelYWhiteSpacePercent = 5
# pres.nglPanelXWhiteSpacePercent = 5

# layout = [1,len(plot)]
layout = [len(plot),1]
if len(plot)==4: layout = [2,2]
ngl.panel(wks,plot[0:len(plot)],layout,pres)
ngl.end()

hc.trim_png(fig_file)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
