import os, ngl, sys, numba, copy, xarray as xr, numpy as np, itertools
np.seterr(divide='ignore', invalid='ignore')
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

deg_to_rad = np.pi/180.
rad_to_deg = 180./np.pi

num_neighbors = 8

all_possible_sets = xr.DataArray(np.array(
      [[0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,1], [0,0,0,0,0,0,1,1], [0,0,0,0,0,1,0,1],
       [0,0,0,0,0,1,1,1], [0,0,0,0,1,0,0,1], [0,0,0,0,1,0,1,1], [0,0,0,0,1,1,0,1],
       [0,0,0,0,1,1,1,1], [0,0,0,1,0,0,0,1], [0,0,0,1,0,0,1,1], [0,0,0,1,0,1,0,1],
       [0,0,0,1,0,1,1,1], [0,0,0,1,1,0,0,1], [0,0,0,1,1,0,1,1], [0,0,0,1,1,1,0,1],
       [0,0,0,1,1,1,1,1], [0,0,1,0,0,1,0,1], [0,0,1,0,0,1,1,1], [0,0,1,0,1,0,1,1],
       [0,0,1,0,1,1,0,1], [0,0,1,0,1,1,1,1], [0,0,1,1,0,0,1,1], [0,0,1,1,0,1,0,1],
       [0,0,1,1,0,1,1,1], [0,0,1,1,1,0,1,1], [0,0,1,1,1,1,0,1], [0,0,1,1,1,1,1,1],
       [0,1,0,1,0,1,0,1], [0,1,0,1,0,1,1,1], [0,1,0,1,1,0,1,1], [0,1,0,1,1,1,1,1],
       [0,1,1,0,1,1,1,1], [0,1,1,1,0,1,1,1], [0,1,1,1,1,1,1,1], [1,1,1,1,1,1,1,1]])
      ,dims=['set','neighbors'] )

chx_only_sets = xr.DataArray(np.array(
      [[0,1,0,1,0,1,0,1], [1,0,1,0,1,0,1,0]])
      ,dims=['set','neighbors'] )

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def get_set_labels(set_array):
   (num_set,set_len) = set_array.shape
   set_labels = []
   for s in range(num_set):
      label = f'{set_array[s,:].values}'
      set_labels.append(label)
   return set_labels

#---------------------------------------------------------------------------------------------------
# routine for reorganizing SCRIP file coordinates with grid_rank==2
#---------------------------------------------------------------------------------------------------
@numba.njit
def reorder_coords(lat,lon,ilat,ilon):
   for j in range(nlat):
      for i in range(nlon):
         n = j*nlon + i
         lat[n],lon[n] = ilat[j],ilon[i]
#---------------------------------------------------------------------------------------------------
# single routine to return sorted neighbors and bearings using methods below
#---------------------------------------------------------------------------------------------------
def get_neighbors_and_bearings(scrip_file):
   
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   lat = scrip_file['grid_center_lat']
   lon = scrip_file['grid_center_lon']
   grid_rank = len(scrip_file['grid_dims'])

   # SCRIP data is always 1D? Otherwise we need the block below
   ncol = len(lat)

   # if grid_rank==1:
   #    ncol = len(lat)
   # if grid_rank==2:
   #    nlat,nlon = len(lat), len(lon)
   #    ncol = len(lat) * len(lon)
   #    ilat,ilon = lat,lon
   #    lat,lon = np.zeros(ncol), np.zeros(ncol)
   #    reorder_coords(lat,lon,ilat,ilon)
   #    lat = xr.DataArray(lat,coords=[np.arange(ncol)],dims='ncol')
   #    lon = xr.DataArray(lon,coords=[np.arange(ncol)],dims='ncol')

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   corner_lon = scrip_file['grid_corner_lon'][:,:].values
   corner_lat = scrip_file['grid_corner_lat'][:,:].values

   neighbors = np.full([ncol,num_neighbors+1],-1)
   bearings  = np.full([ncol,num_neighbors+1],-1)
   edge_flag = np.full([ncol,num_neighbors+1],-1)

   # Find adjacent neighbors
   find_neighbors(lat.values, lon.values, corner_lat, corner_lon, neighbors, bearings, edge_flag)

   neighbors = neighbors[:,1:num_neighbors+1]
   bearings  = bearings [:,1:num_neighbors+1]
   edge_flag = edge_flag[:,1:num_neighbors+1]
   
   #----------------------------------------------------------------------------
   # sort neighbors by bearing
   #----------------------------------------------------------------------------
   bear = bearings[:,:]
   bear = np.where(bear<0,bear+360,bear)
   neighbors_sorted = sort_neighbors( neighbors[:,:], bear )
   bearings_sorted  = sort_neighbors( bearings[:,:],  bear )
   edge_flag_sorted = sort_neighbors( edge_flag[:,:], bear )

   # make sure the northernmost neighbor (first in list) is an edge neighbor, not a corner neighbor
   adjust_sorted_neighbors( ncol, neighbors_sorted, bearings_sorted, edge_flag_sorted )

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   return ( neighbors_sorted.astype(int), bearings_sorted )

#---------------------------------------------------------------------------------------------------
# routines for finding neighbors using SCRIP corner data
#---------------------------------------------------------------------------------------------------
# @numba.njit
def calc_great_circle_distance(lat1,lat2,lon1,lon2):
   '''
   input should be in degrees
   '''
   dlon = lon2 - lon1
   cos_dist = np.sin(lat1*deg_to_rad)*np.sin(lat2*deg_to_rad) + \
              np.cos(lat1*deg_to_rad)*np.cos(lat2*deg_to_rad)*np.cos(dlon*deg_to_rad)
   # print( str(cos_dist.min()) +"   "+ str(cos_dist.max()) )
   cos_dist = np.where(cos_dist> 1.0, 1.0,cos_dist)
   cos_dist = np.where(cos_dist<-1.0,-1.0,cos_dist)
   dist = np.arccos( cos_dist )
   return dist

@numba.njit
def calc_great_circle_bearing(lat1_in,lat2_in,lon1_in,lon2_in):
   deg_to_rad = np.pi/180.
   rad_to_deg = 180./np.pi
   lat1 = lat1_in * deg_to_rad
   lat2 = lat2_in * deg_to_rad
   lon1 = lon1_in * deg_to_rad
   lon2 = lon2_in * deg_to_rad

   dlon = lon1 - lon2

   atan_tmp1 = np.sin(lon2-lon1)*np.cos(lat2)
   atan_tmp2 = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
   bearing = np.arctan2( atan_tmp1, atan_tmp2 )

   bearing = bearing * rad_to_deg
   return bearing

@numba.njit()
def find_neighbors(center_lat, center_lon, corner_lat, corner_lon, neighbors, bearings, edge_flag):
   '''
   Use FV cell information from a SCRIP file to identify adjacent cells
   center_lat[num_col,num_corner]    - input latitude at cell centers
   center_lon[num_col,num_corner]    - input longitude at cell centers
   corner_lat[num_col,num_corner]    - input latitude at cell corners
   corner_lon[num_col,num_corner]    - input longitude at cell corners
   neighbors [num_col,num_neighbors] - output column indices of neighbors                       
   bearings  [num_col,num_neighbors] - output great circle bearing of each neighbor (for sorting)
   edge_flag [num_col,num_neighbors] - output flag indicating if neighbor shares edge or corner
   '''
   ncol = len(center_lat)

   ### Exact equality doesn't always work so we need to round the coords a little
   for n in range(ncol):
      for c in range(4):
         corner_lat[n,c] = np.round(corner_lat[n,c],num_neighbors)
         corner_lon[n,c] = np.round(corner_lon[n,c],num_neighbors)
         # also make sure all prime meridian points have lon=0 rather than 360
         if corner_lon[n,c]==360: corner_lon[n,c] = 0
   
   # loop through all points that will define the center of the neighborhood
   for n in range(ncol):
      neighbors[n,0] = n
      cnt = 1
      # loop over all 4 corners of the FV cell
      for c in range(4):
         # loop over all points again to find points that share each corner
         for nn in range(ncol):
            # make sure this potential neighbor has not been found already
            if nn not in neighbors[n,0:cnt] :
               # check if corner locations match 
               if  corner_lat[n,c] in corner_lat[nn,:] \
               and corner_lon[n,c] in corner_lon[nn,:] :
                  neighbors[n,cnt] = nn
                  bearings[n,cnt] = calc_great_circle_bearing( center_lat[n], center_lat[nn],\
                                                               center_lon[n], center_lon[nn] )
                  # count how many cell corners are shared
                  common_corner_cnt = 0
                  for cc in range(4):
                     if (corner_lat[n,cc] in corner_lat[nn,:]) and \
                        (corner_lon[n,cc] in corner_lon[nn,:]):
                        common_corner_cnt += 1
                  # if 2 corners are shared then this is an edge neighbor
                  edge_flag[n,cnt] = 1 if common_corner_cnt==2 else 0
                  # increment the neighbor count
                  cnt += 1
         # stop after finding 8 neighbors
         if cnt==9: break

            # if nn==ncol: print('neighbor not found?!?!')
   return

#---------------------------------------------------------------------------------------------------
# routines for sorting neighbor sets
#---------------------------------------------------------------------------------------------------
# @numba.njit() - argsort doesnt work with numba?
def sort_neighbors( neighbors, bearing ):
   num_col = neighbors.shape[0]
   neighbors_sorted = np.zeros(neighbors.shape)
   for i in range(num_col): neighbors_sorted[i,:] = neighbors[i, np.argsort(bearing[i,:]) ]
   return neighbors_sorted

# make sure the northernmost neighbor (first in list) is an edge neighbor, not a corner neighbor
@numba.njit()
def adjust_sorted_neighbors(ncol, neighbors_sorted, bearings, edge_flag):
   for n in range(ncol):
      neighbors_out = neighbors_sorted
      if np.any(neighbors_sorted[n,:]==-1): continue
      if edge_flag[n,0]==0 :
         b1 = np.absolute( bearings[n,7] )
         b2 = np.absolute( bearings[n,1] )
         if b1<b2:           neighbors_out[n,:] = np.concatenate( ( neighbors_sorted[n,7:8], neighbors_sorted[n,0:7] ) )
         if b1>b2 or b1==b2: neighbors_out[n,:] = np.concatenate( ( neighbors_sorted[n,1:8], neighbors_sorted[n,0:1] ) )
   neighbors_sorted = neighbors_out
   return

#---------------------------------------------------------------------------------------------------
# routines for counting the number of possible ordered binary sets (0's and 1's)
#---------------------------------------------------------------------------------------------------
def get_factors(x):
   """
   return list of all multiplicative factors of x
   """
   factors = []
   for i in range(1, x + 1):
      if x % i == 0: factors.append(i)
   return factors[::-1]

def count_unique_sets_no_repeat(n):
   """
   count unique binary sets of length n with no repeating patterns
   subtract sets with repeating patterns smaller than length n
   """
   count = np.power(2.,n)
   if n>1:
      for f in get_factors(n)[1:]: count -= count_unique_sets_no_repeat(f)
   return count

def count_rotated_sets(n):
   """
   count unique binary sets for each factor f of n 
   divide each by f to avoid double counting cyclicly identical patterns
   """
   count = 0
   for f in get_factors(n): count += count_unique_sets_no_repeat(f)/f
   return int(count)

#---------------------------------------------------------------------------------------------------
# routines for removing the mean of the local neighborhood
#---------------------------------------------------------------------------------------------------
@numba.njit()
def remove_mean_numba(nn_values,center_value):
   """
   remove the mean of the nearest neighborhood
   """
   (num_time,num_neighbors) = nn_values.shape
   anomalies = np.zeros((num_time,num_neighbors))
   for t in range(num_time):
      mean = ( np.sum(nn_values[t,:]) + center_value[t] ) / ( num_neighbors+1 )
      for n in range(num_neighbors):
         anomalies[t,n] = nn_values[t,n] - mean
   return anomalies

@numba.njit()
def remove_neighbor_mean_numba( data, neighbors, bearings, num_neighbors=8, sort_by_bearing=False ):
   (num_time,num_col) = data.shape
   nn_anomalies = np.zeros( ( num_time, num_col, num_neighbors ) )
   nn_values = np.zeros(num_neighbors)
   for t in range(num_time):
      for i in range(num_col) :
         # get nearest neighbors
         nn = neighbors[i,:]
         if sort_by_bearing:
            sort_ind = np.argsort( bearings[i,:] )
            nn = nn[sort_ind]
         # calculate nearest neighbor data as anomalies from center point
         for n in range(num_neighbors) : nn_values[n] = data[t,nn[n]]
         nn_anomalies[t,i,:] = nn_values[:] - data[t,i]
   return nn_anomalies

#---------------------------------------------------------------------------------------------------
# routines for counting the occurrence of a collection of sets
#---------------------------------------------------------------------------------------------------
@numba.njit()
def count_sets_numba(states,sets,valid_flag=None,rotate_sets=False):
   """
   count the occurrences of each set in states
   input states should be single column and have dims (time,neighbor)
   """
   (num_set,set_len) = sets.shape
   (num_state,state_len) = states.shape
   if state_len!=set_len: raise ValueError('Inputs must have identical second dimension!')
   count = np.zeros(num_set)
   for n in range(num_state):
      if valid_flag is not None:
         if valid_flag[n]!=1: continue
      found = False
      for s in range(num_set):
         if found: break
         if rotate_sets:
            for l in range(set_len):
               if found: break
               s_rot = np.concatenate( ( sets[s,l:], sets[s,:l] ) )
               # s_rot = xr.concat( ( sets[s,l:], sets[s,:l] ), dim='neighbors' )
               ### check if rotation works correctly
               # print(f'  {states[n,:]}  {sets[s,:]}   {s_rot}   {np.all(states[n,:] == s_rot)}')            
               if np.all(states[n,:] == s_rot[:]):
                  count[s] += 1
                  found = True
         else:
            if np.all(states[n,:] == sets[s,:]):
               count[s] += 1
               found = True
   return count

def count_sets(states,sets):
   """
   count the occurrences of each set in states
   """
   (num_set,set_len) = sets.shape
   count = xr.DataArray( np.zeros(num_set) )
   for state in states:
      found = False
      for s in range(num_set):
         if found: break
         for l in range(set_len):
            if found: break
            # s_rot = np.concatenate( ( sets[s,l:], sets[s,:l] ) )
            s_rot = xr.concat( ( sets[s,l:], sets[s,:l] ), dim='neighbors' )
            # print(f'    {state.values}  {sets[s,:].values}   {s_rot.values}   {np.all(state.values == s_rot.values)}')            
            if np.all(state.values == s_rot.values):
               count[s] += 1
               found = True
   return count

@numba.njit()
def count_sets_per_col_numba(nn_states,sets,ntime,ncol,rotate_sets,valid_flag=None):
   cnt_out = np.zeros((ncol,len(sets)))
   # time_ind = np.arange(ntime)
   if valid_flag is None: valid_flag = np.ones((ntime,ncol))
   for i in range(ncol) :
      ### count the occurrence of sets and calculate frequency
      count = count_sets_numba( nn_states[:,i,:], sets, valid_flag[:,i], rotate_sets=rotate_sets )
      cnt_out[i,:] = count
   return cnt_out

@numba.njit()
def count_sets_per_col_keep_time_numba(nn_states,sets,ntime,ncol,rotate_sets,valid_flag=None):
   '''
   Note that nn_states and valid_flag include a dummy dimension (axis=0)
   '''
   (num_set,set_len) = sets.shape
   # (ntime,ncol,state_len) = nn_states.shape
   cnt_out = np.zeros((ntime,ncol,len(sets)))
   # time_ind = np.arange(ntime)
   if valid_flag is None: valid_flag = np.ones((ntime,ncol))
   for i in range(ncol) :
         # count = count_sets_numba( nn_states[:,t,i,:], sets, 
         #                           valid_flag[:,t,i] , rotate_sets=rotate_sets )
         # cnt_out[t,i,:] = count
         # def count_sets_numba(states,sets,valid_flag=None,rotate_sets=False):
         # (num_set,set_len) = sets.shape
         # (num_state,state_len) = nn_states.shape
         # if state_len!=set_len: raise ValueError('Inputs must have identical second dimension!')
         # count = np.zeros(num_set)
         for t in range(ntime):
            if valid_flag[t,i]!=1: continue
            found = False
            for s in range(num_set):
               if found: break
               if rotate_sets:
                  for l in range(set_len):
                     if found: break
                     s_rot = np.concatenate( ( sets[s,l:], sets[s,:l] ) )
                     if np.all(nn_states[t,i,:] == s_rot[:]):
                        cnt_out[t,i,s] += 1
                        found = True
               else:
                  if np.all(nn_states[t,i,:] == sets[s,:]):
                     cnt_out[t,i,s] += 1
                     found = True
   return cnt_out

#---------------------------------------------------------------------------------------------------
# routines for checking stuff
#---------------------------------------------------------------------------------------------------
def is_partial_checkerboard(set_in,subset_length=6):
   """
   check for partial checkerboard given a minimum alternating sequence length
   """
   found = False
   for l in range(len(set_in)):
      if found: break
      s_rot = np.concatenate( ( set_in[l:], set_in[:l] ) ) # rotate the set
      subset1 = [0,1]*int(subset_length/2)
      subset2 = [1,0]*int(subset_length/2)
      if np.all(s_rot[:subset_length]==subset1): found = True
      if np.all(s_rot[:subset_length]==subset2): found = True
   return found 

def set_in_rotated_sets(test_set,sets):
   """
   brute force check if test_set is contained within sets
   also check n-1 rotations of each set in sets
   """
   (num_set,set_len) = sets.shape
   found = False
   for s in range(num_set):
      if all(sets[s,:]<0): continue
      for l in range(set_len):
         s_rot = np.concatenate( ( sets[s,l:], sets[s,:l] ) )
         if (test_set == s_rot).all():
            found = True
         if found: break
      if found: break
   return found

#---------------------------------------------------------------------------------------------------
# routine for brute force generation of all possible sets
#---------------------------------------------------------------------------------------------------
def generate_sets(n):
   # generate all possible perturbations
   permutations = [p for p in itertools.product([0,1], repeat=n)]
   # initialize output array to hold sets
   sets = xr.DataArray( [ np.array(permutations[0]) ], dims=['set','neighbors'] )
   # brute force check to throw out sets that are identical to other rotated ones
   for p in permutations[1:]:
      test_set = xr.DataArray( [ np.array(p) ], dims=['set','neighbors'] )
      if not set_in_rotated_sets(test_set,sets): 
         sets = xr.concat( [sets, test_set], dim='set' )
   # check the final product and assign coordinate
   num_sets = count_rotated_sets(n)
   num_set_found = sets.shape[0]
   if num_set_found != num_sets: 
      raise ValueError(f'number of sets is incorrect {num_set_found} vs {num_sets}')
   set_coord,nn_coord = np.arange(num_sets),np.arange(n)
   sets.assign_coords(coords={'set':set_coord,'neighbors':nn_coord})
   return(sets)

#---------------------------------------------------------------------------------------------------
# generate flag to indicate valid points
#---------------------------------------------------------------------------------------------------
@numba.njit()
def get_valid_flag_numba(nn_anomalies,neighbors,ntime,ncol,num_neighbor):
   valid_flag = np.ones((ntime,ncol))
   for i in range(ncol):
      for t in range(ntime):
         tmp_anom = nn_anomalies[t,i,:]
         if np.any(np.isnan(tmp_anom)): valid_flag[t,i] = 0
         if np.any(np.isinf(tmp_anom)): valid_flag[t,i] = 0
         if np.any(neighbors[i,:]==-1): valid_flag[t,i] = 0
         # GPM data has a lot of missing data as 9.99999961690316e+35 so...
         if np.any(tmp_anom> 1e35): valid_flag[t,i] = 0
         if np.any(tmp_anom<-1e35): valid_flag[t,i] = 0
   return valid_flag
#---------------------------------------------------------------------------------------------------
# top level routine to encapsulate the pattern detection algorithm
#---------------------------------------------------------------------------------------------------
def find_neighbors_and_count_sets( data, scripfile_path, sets, 
                                   rotate_sets, keep_time=False, 
                                   verbose=False, indent=' '*4 ):
   """
   top level routine to encapsulate the pattern detection algorithm
   """

   # use scrip data to locate adjacent neighbor
   if verbose: print(indent+'Loading scrip file data and finding adjacent neighbors')
   scripfile = xr.open_dataset(scripfile_path)
   ( neighbors, bearings ) = get_neighbors_and_bearings(scripfile)

   # calculate anomalies relative to local neighbor
   if verbose: print(indent+'Calculating anomalies from local neighborhood mean')
   nn_anomalies = remove_neighbor_mean_numba( data.values, neighbors, bearings, sort_by_bearing=False)

   # define array of neighbor states
   nn_states = np.sign( nn_anomalies )
   nn_states = xr.where(nn_states<=0,0,nn_states)
   nn_states = nn_states.astype(type(sets.values[0,0]))

   # flag valid points considering the neighbor states - also deal with cube sphere corners
   valid_flag = get_valid_flag_numba( nn_anomalies, neighbors, 
                                      len(data['time']), len(data['ncol']), len(sets))

   # Count sets
   if verbose: print(indent+'Counting neighbor state sets')
   if keep_time:
      # add dummy dimension as a hacky workaround 
      cnt_list = count_sets_per_col_keep_time_numba( nn_states, sets.values, 
                                                     len(data['time']), len(data['ncol']), 
                                                     rotate_sets, valid_flag )
   else:
      cnt_list = count_sets_per_col_numba( nn_states, sets.values, 
                                           len(data['time']), len(data['ncol']), 
                                           rotate_sets, valid_flag )

   num_valid = xr.DataArray(valid_flag,dims=('time','ncol')).sum(dim='time')

   # put set count data into an xarray dataset
   if keep_time:
      if verbose: print(indent+'Creating final dataset')
      cnt_ds = xr.Dataset()
      cnt_ds['cnt'] = ( ('time','ncol','set'), cnt_list )
      cnt_ds.coords['time'] = ('time',data['time'].values)
      cnt_ds.coords['ncol'] = ('ncol',data['ncol'].values)
      cnt_ds.coords['set'] = ('set',sets['set'])
      cnt_ds['num_time'] = len(data['time'].values)
      cnt_ds['num_valid'] = ( ('ncol'), num_valid )
   else:
      if verbose: print(indent+'Creating final dataset')
      cnt_ds = xr.Dataset()
      cnt_ds['cnt'] = ( ('ncol','set'), cnt_list )
      cnt_ds.coords['ncol'] = ('ncol',data['ncol'].values)
      cnt_ds.coords['set'] = ('set',sets['set'])
      # add the neighbor, comes in handy for debugging
      cnt_ds['neighbors'] = ( ('ncol','neighbor'), neighbors )
      cnt_ds.coords['neighbor'] = ('neighbor',np.arange(8))
      # add valid time length for later normalization
      cnt_ds['num_time'] = len(data['time'].values)
      cnt_ds['num_valid'] = ( ('ncol'), num_valid )

   return cnt_ds
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

