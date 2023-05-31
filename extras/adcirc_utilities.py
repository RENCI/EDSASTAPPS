# Derived from git@github-work:BrianOBlanton/py_utils

import os
import re
import yaml
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy import spatial as sp

from datetime import date, datetime
import netCDF4 as nc4
from netCDF4 import num2date, date2num, date2index

import matplotlib.tri as Tri
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import time as tm

# import cartopy as ctpi
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# dict for correspondance between internal netCDF variable names and netCDF file names.
# these should be determined by querying the standard names ...
varnamedict = {
  "zeta_max": "maxele.63.nc",
  "vel_max":  "maxvel.63.nc",
  "inun_max": "maxinundepth.63.nc",
  "zeta":     "fort.63.nc",
  "swan_HS":  "swan_HS.63.nc"
}
print("utilities:netCDF4 Version = {}".format(nc4.__version__))

def load_config(area='DEFAULT', yaml_file='config.yml'):
    if not os.path.exists(yaml_file):
        raise IOError("Failed to load yaml config file {}".format(yaml_file))
    with open(yaml_file, 'r') as stream:
        config = yaml.safe_load(stream)
    return config[area]

def validate_url(url):  # TODO
    """
    """
    pass

def get_adcirc_grid(nc):
    """
    """
    agdict = {}
    agdict['lon'] = nc.variables['x'][:]
    agdict['lat'] = nc.variables['y'][:]
    agdict['ele'] = nc.variables['element'][:,:] - 1
    agdict['depth'] = nc.variables['depth'][:]
    agdict['latmin'] = np.mean(nc.variables['y'][:])  # needed for scaling lon/lat plots
    return agdict

def get_adcirc_grid_from_ds(ds):
    """
    """
    agdict = {}
    agdict['lon'] = ds['x'][:]
    agdict['lat'] = ds['y'][:]
    agdict['ele'] = ds['element'][:,:] - 1
    agdict['depth'] = ds['depth'][:]
    agdict['latmin'] = np.mean(ds['y'][:])  # needed for scaling lon/lat plots
    return agdict

def attach_element_areas(agdict):
    """
    """
    x=agdict['lon'].values
    y=agdict['lat'].values
    e=agdict['ele'].values
    
    # COMPUTE GLOBAL DX,DY, Len, angles
    i1=e[:,0]
    i2=e[:,1]
    i3=e[:,2]

    x1=x[i1];x2=x[i2];x3=x[i3];
    y1=y[i1];y2=y[i2];y3=y[i3];

    # coordinate deltas
    dx23=x2-x3
    dx31=x3-x1
    dx12=x1-x2
    dy23=y2-y3
    dy31=y3-y1
    dy12=y1-y2

    # lengths of sides
    a = np.sqrt(dx12*dx12 + dy12*dy12)
    b = np.sqrt(dx31*dx31 + dy31*dy31)
    c = np.sqrt(dx23*dx23 + dy23*dy23)
    
    agdict['areas'] = ( x1*dy23 + x2*dy31 + x3*dy12 )/2.
    agdict['edge_lengths']=[a, b, c];
    agdict['dl']=np.mean(agdict['edge_lengths'],axis=0)

    return agdict
    
def basis2d_withinElement(phi):
    """
    """
    interior_status = np.all(phi[:]<=1,axis=1) & np.all(phi[:]>=0,axis=1)
    return interior_status

def basis2d(agdict,xylist,j):
    """
    """
    # check length of j and xylist
    # check for needed arrays in agdict
    phi=[]
    #nodes for the elements in j
    n3=agdict['ele'][j]
    x=agdict['lon'][n3].values     
    x1=x[:,0];x2=x[:,1];x3=x[:,2];
    y=agdict['lat'][n3].values
    y1=y[:,0];y2=y[:,1];y3=y[:,2];  

    areaj=agdict['areas'][j]
    xp=xylist[:,0]
    yp=xylist[:,1]

    # Basis function 1
    a=(x2*y3-x3*y2)
    b=(y2-y3)
    c=-(x2-x3)
    phi0=(a+b*xp+c*yp)/(2.0*areaj)
    # Basis function 2
    a=(x3*y1-x1*y3)
    b=(y3-y1)
    c=-(x3-x1)
    phi1=(a+b*xp+c*yp)/(2.0*areaj)
    # Basis function 3
    a=(x1*y2-x2*y1)
    b=(y1-y2)
    c=-(x1-x2)
    phi2=(a+b*xp+c*yp)/(2.0*areaj)
    
    return np.array([phi0, phi1, phi2]).T

    
def get_adcirc_trifxn(agdict):
    """
    """
    return Tri.Triangulation(agdict['lon'],agdict['lat'], triangles=agdict['ele'])
    
def get_adcirc_time(nc):
    """
    """
    atdict = {}
    temp = nc.variables['time']
    dates = num2date(temp[:], temp.units) # , only_use_python_datetimes=False)
    #dates = [date.strftime('%Y%m%d%H%M') for date in dates]
    dates = [pd.to_datetime(date.strftime('%Y%m%d%H%M')) for date in dates]
    atdict['time'] = dates    
    return atdict

def get_adcirc_time_from_ds(ds):
    """
    """
    return {'time': ds['time']}

def get_adcirc_slice(nc,v,it=None):
    """
    """
    advardict = {}
    var = nc.variables[v]
    if re.search('max', v):
        var_d = var[:] # the actual data
    else:
        if var.dimensions[0] == 'node':
            #print('transposed')
            var_d = var[it,:].T
        elif var.dimensions[0] == 'time':
            var_d = var[:,it] # the actual data
        else:
            print(f'netCDF4: Unexpected leading variable name {nc[v].dimensions}: Abort')
            sys.exit(1)
    var_d[var_d.mask] = np.nan
    advardict['var'] = var_d.data
    return advardict

def get_adcirc_slice_from_ds(ds,v,it=0):
    """
    """
    advardict = {}
    var = ds.variables[v]
    if re.search('max', v) or re.search('depth', v):
        var_d = var[:] # the actual data
    else:
        if ds.variables[v].dims[0] == 'node':
            #print('transposed')
            var_d = var[it,:].T # the actual data
        elif ds.variables[v].dims[0] == 'time':
            var_d = var[:,it] # the actual data
        else:
            print(f'Unexpected leading variable name {ds.variables[v].dims}: Abort')
            sys.exit(1)
    #var_d[var_d.mask] = np.nan
    advardict['var'] = var_d.data
    return advardict

def f63_to_xr(url):
    """
    """
    dropvars=['neta', 'nvel',  'max_nvdll', 'max_nvell']
    return xr.open_dataset(url,drop_variables=dropvars)

def f61_to_xr(url):
    """
    """
    ds=xr.open_dataset(url)
    ids=[]
    for i,s in enumerate(ds['station_name'].values):
        s=s.decode().strip()
        ids.append(s)
#     temp={'station_name': ids}
#     print(temp)
    ds=ds.assign_coords(station_name=ids)
    ids=[]
    for k,v in enumerate(ds.station_name):
        stationID=np.array2string(v.values).split()[0][1:]
        ids.append(stationID)
    
    return ds.assign_coords(station_ids=ids)
    
def f61_to_dict(url):
    nc=nc4.Dataset(url)
    var_name='zeta'
    time=nc.variables['time']
    dtime=nc4.num2date(time[:],time.units)
    tstart=dtime[0].strftime('%Y-%b-%d')
    var=nc.variables[var_name]
    data=var[:,:]
    sn=nc.variables['station_name']
    lstr=[]
    for i in range(sn.shape[0]) :
        temp=str(nc4.chartostring(sn[i,:])).strip() 
        lstr.append(temp)

    return {'nc': nc,
            'time': dtime,
            'tstart': tstart,
            'data': data,
            'station_name': lstr}

def discrete_cmap(N, base_cmap=None):
    """
    Create an N-bin discrete colormap from the specified input map
    """

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

# def set_gridLines(ax, xlabels_top,xlabels_bottom,ylabels_left,ylabels_right, dl=5):
#     """
#     """
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.25, color='black', 
#                       alpha=0.5, linestyle='--', draw_labels=True)
#     gl.x_inline = False
#     gl.y_inline = False
#     gl.top_labels = xlabels_top
#     gl.bottom_labels=xlabels_bottom
#     gl.left_labels = ylabels_left
#     gl.right_labels=ylabels_right
#     gl.xlines = True
#     gl.ylines = True
#     gl.xlocator = mticker.FixedLocator(np.arange(-180,180,dl))
#     gl.ylocator = mticker.FixedLocator(np.arange(-90,90,dl))
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlabel_style = {'color': 'k', 'size': 14} #, 'weight': 'bold'}    
#     gl.ylabel_style = {'color': 'k', 'size': 14} #, 'weight': 'bold'}    

# def addFeatures(ax):
#     """
#     """
#     kwargs = {"scale": '50m', "facecolor": "none",  "edgecolor": "k", "linewidth": 1.5} 
#     ax.add_feature(cfeature.NaturalEarthFeature(category='physical', name='coastline', **kwargs))
#     kwargs = {"scale": '110m', "facecolor": "none",  "edgecolor": "k", "linewidth": .5} 
#     ax.add_feature(cfeature.NaturalEarthFeature(category='physical', name='lakes', **kwargs))
#     ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', **kwargs))
#     ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', **kwargs))

def fft_lowpass(signal, low, high):
    """
    Performs a low pass filer on the series.
    low and high specifies the boundary of the filter.

    >>> from oceans.filters import fft_lowpass
    >>> import matplotlib.pyplot as plt
    >>> t = np.arange(500)  # Time in hours.
    >>> x = 2.5 * np.sin(2 * np.pi * t / 12.42)
    >>> x += 1.5 * np.sin(2 * np.pi * t / 12.0)
    >>> x += 0.3 * np.random.randn(len(t))
    >>> filtered = fft_lowpass(x, low=1/30, high=1/40)
    >>> fig, ax = plt.subplots()
    >>> l1, = ax.plot(t, x, label='original')
    >>> l2, = ax.plot(t, filtered, label='filtered')
    >>> legend = ax.legend()

    """

    if len(signal) % 2:
        result = np.fft.rfft(signal, len(signal))
    else:
        result = np.fft.rfft(signal)

    freq = np.fft.fftfreq(len(signal))[:len(signal) // 2 + 1]
    factor = np.ones_like(freq)
    factor[freq > low] = 0.0
    sl = np.logical_and(high < freq, freq < low)
    a = factor[sl]

    # Create float array of required length and reverse.
    a = np.arange(len(a) + 2).astype(float)[::-1]

    # Ramp from 1 to 0 exclusive.
    a = (a / a[0])[1:-1]

    # Insert ramp into factor.
    factor[sl] = a

    result = result * factor

    return np.fft.irfft(result, len(signal))

def lanc(numwt, haf):
    """
    Generates a numwt + 1 + numwt lanczos cosine low pass filter with -6dB
    (1/4 power, 1/2 amplitude) point at haf

    Parameters
    ----------
    numwt : int
            number of points
    haf : float
          frequency (in 'cpi' of -6dB point, 'cpi' is cycles per interval.
          For hourly data cpi is cph,

    Examples
    --------
    >>> from oceans.filters import lanc
    >>> import matplotlib.pyplot as plt
    >>> t = np.arange(500)  # Time in hours.
    >>> h = 2.5 * np.sin(2 * np.pi * t / 12.42)
    >>> h += 1.5 * np.sin(2 * np.pi * t / 12.0)
    >>> h += 0.3 * np.random.randn(len(t))
    >>> wt = lanc(96+1+96, 1./40)
    >>> low = np.convolve(wt, h, mode='same')
    >>> high = h - low
    >>> fig, (ax0, ax1) = plt.subplots(nrows=2)
    >>> _ = ax0.plot(high, label='high')
    >>> _ = ax1.plot(low, label='low')
    >>> _ = ax0.legend(numpoints=1)
    >>> _ = ax1.legend(numpoints=1)

    """
    summ = 0
    numwt += 1
    wt = np.zeros(numwt)

    # Filter weights.
    ii = np.arange(numwt)
    wt = 0.5 * (1.0 + np.cos(np.pi * ii * 1. / numwt))
    ii = np.arange(1, numwt)
    xx = np.pi * 2 * haf * ii
    wt[1:numwt + 1] = wt[1:numwt + 1] * np.sin(xx) / xx
    summ = wt[1:numwt + 1].sum()
    xx = wt.sum() + summ
    wt /= xx
    return np.r_[wt[::-1], wt[1:numwt + 1]]

def crosscorr(datax, datay, lag=0):
    """ 
    Lag-N cross correlation.
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    return datax.corr(datay.shift(lag))

def AdcircCppForward(lon, lat, lon0=-70, lat0=32):
    """
    Platte Carre Forward projection
    """
    r = 6378206.4
    x = r*(lon-lon0)*np.pi/180.*np.cos(lat0*np.pi/180)
    y = lat*np.pi/180*r
    return x, y 

def AdcircCppInverse(x, y, lon0=-70, lat0=32):
    """
    Platte Carre Inverse projection
    """    
    r = 6378206.4
    alpha = np.cos(lat0*np.pi/180)
    lam = lon0 + 180 / np.pi*(np.divide(x,r*alpha))
    phi = 180/np.pi*y/r
    return lam, phi

def ComputeTree(agdict):
    """
    Given the pre-loaded lon,lat,ele entries to agdict compute the means and 
    genrate the full ADCIRC grid KDTree
    """
    t0=tm.time()
    try:
        x=agdict['lon'].values.ravel() # ravel; not needed
        y=agdict['lat'].values.ravel()
        e=agdict['ele'].values
    except Exception as e:
        print('Failed in finding lon,lat,ele data on agdict. Perhaps check for preprocesing using attach_element_areas method. {}'.format(e))
        sys.exit(1)
    xe=np.mean(x[e],axis=1)
    ye=np.mean(y[e],axis=1)
    agdict['tree']=tree = sp.KDTree(np.c_[xe,ye])
    print(f'Build KDTree time is {tm.time()-t0}s')
    return agdict

def ComputeQuery(xylist, agdict, kmax=10):
    """
    Generate the kmax-set of nearest neighbors to each lon,lat pair xylst.
    Each test point (each lon/lat pair) gets associaed with a distance (dd) and element (j) object 
    At this stage it is possible that some test points are not interior to the nearest element. We will
    subsequently check that.

    dd: num points by neighbors
    j: num points by neighbors
    """
    t0=tm.time()
    agresults=dict()
    dd, j = agdict['tree'].query(xylist, k=kmax)
    if kmax==1:
        dd=dd.reshape(-1,1)
        j=j.reshape(-1,1)
    agresults['distance']=dd
    agresults['elements']=j
    agresults['number_neighbors']=kmax
    agresults['geopoints']=xylist # We shall use this later
    print(f'KDTree query {kmax} took {tm.time()-t0}s')
    return agresults

def ComputeBasisRepresentation(xylist, agdict, agresults):
    """
    For each test point with kmax number_neighbors, compute basis representations for
    each neighbor. Then, check which, if any, element the test point actually resides within.
    If not, then the returned basis function (weights) set to nans. Doing it this way
    will retain the order of test points, overall
    """
    # First build all the basis weights and determine if it was interior or not
    t0=tm.time()
    kmax = agresults['number_neighbors']
    j = agresults['elements']
    phival_list=list()
    within_interior=list()
    # Let's try to keep the FIRST True value of the lists
    for k_value in range(0,kmax):
        phival=basis2d(agdict,xylist,j[:,k_value])
        phival_list.append(phival)
        within_interior.append(basis2d_withinElement(phival))
    # Second only keep the "interior" results or nans if none
    final_weights= np.full( (phival_list[0].shape[0],phival_list[0].shape[1]),np.nan)
    final_jvals = np.full( j.T[0].shape[0],-99999)
    final_status = np.full( within_interior[0].shape[0],False)
    for pvals,jvals,testvals in zip(phival_list, j.T, within_interior): 
        final_weights[testvals] = pvals[testvals]
        final_jvals[testvals]=jvals[testvals]
        final_status[testvals] = testvals[testvals]
    agresults['final_weights']=final_weights
    agresults['final_jvals']=final_jvals
    agresults['final_status']=final_status
    print(f'Identify best basis weights took {tm.time()-t0}s')
    outside_elements = np.argwhere(np.isnan(final_weights).all(axis=1)).ravel()
    agresults['outside_elements']=outside_elements
    return agresults

def WaterLevelReductions(t, data_list, final_weights):
    """
    Each data_list is a df for a single test point containing 3 columns. 
    These columns are reduced using the final_weights previously calculated
    
    A final df is returned with index=time and a single column for each of the
    input test points (some of which may be partially or completely nan)
    """
    final_list = list()
    for index,dataseries,weights in zip(range(0,len(data_list)), data_list,final_weights):
        reduced_data = np.matmul(dataseries.values, weights.T)
        df = pd.DataFrame(reduced_data, index=t, columns=[f'P{index+1}'])
        final_list.append(df)
    df_final_data = pd.concat(final_list, axis=1)
    return df_final_data

# An unused alternatve method specific to water levels
#def WaterLevelFromDs(ds, elements):
#    """
#    A simple reading of ADCIRC data but also checking for transposed or not
#    No need to worry about the times nor the column names. They get applied at the end
#    """
#    zeta = ds.zeta
#    if zeta.dims[0] == 'node':
#        nnodes, ntimes = zeta.shape
#        print('file is transposed')
#        print(elements)
#        df=pd.DataFrame(zeta[elements].values)
#    elif zeta.dims[0] == 'time':
#        nnodes, ntimes = zeta.shape
#        print(eleements)
#        df=pd.DataFrame(zeta[:,elements].values)
#    else:
#        print(f'Unexpected leading variable name {zeta.dims[0]}: Abort')
#        sys.exit(1)
#    print(f"Nnodes,Ntimes = {nnodes, ntimes}")
#    return df.T

def GenerateMetadata(agresults):
    """
    """
    df_meta=pd.DataFrame(agresults['geopoints'], columns=['LON','LAT'])
    df_meta['GEOPOINT']=df_meta.index+1
    df_meta.set_index('GEOPOINT', inplace=True)
    df_meta.rename('P{}'.format, inplace=True)
    return df_meta
## Can get the equivelent data using the calls:
# advardict = get_adcirc_slice(nc_data,'zeta',it=e[vstation])
# df = pd.DataFrame(advardict['var'])
## Which is slightly faster than calling the ds data object

def ConstructReducedWaterLevelData_from_ds(ds, agdict, agresults, variable_name=None): # Maybe sed in "ds" instead of url
    """
    This method acquires ADCIRC water levels (three for the current grids) for each specified element of the grid
    The resulting time series are reduced to a single time series using a (basis2d) weighted sum. For a non-nan value to 
    result in the final data product the data must:
    1) Be non-nan for each time series' at the specified time tick
    2) The test point muct be interior to the specified element
    """
    if variable_name is None:
        print('User MUST supply the correct variable name')
        sys.exit(1)
    print(f'Variable name is {variable_name}')
    t0 = tm.time()
    data_list=list()
    final_weights = agresults['final_weights']
    final_jvals = agresults['final_jvals']
    #nc_data = nc4.Dataset(url)
    #acdict = get_adcirc_time(nc_data)
    acdict=get_adcirc_time_from_ds(ds)
    t=acdict['time'].values
    e = agdict['ele'].values
    # Need to call the internal reader
    for vstation in final_jvals:
        ## Not used anymore df = WaterLevelFromDs(ds, e[vstation].T)
        advardict = get_adcirc_slice_from_ds(ds,variable_name,it=e[vstation])
        ## Can access using netCDF4 - This is slightly faster
        ##advardict = get_adcirc_slice(nc_data,variable_name,it=e[vstation])
        df = pd.DataFrame(advardict['var'])
        data_list.append(df)
    print(f'Time to fetch all test station (triplets) was {tm.time()-t0}s')
    df_final=WaterLevelReductions(t, data_list, final_weights)
    df_meta=GenerateMetadata(agresults) # This is here mostly for future considerations
    print(f'Time to reduce {len(final_jvals)} test stations is {tm.time()-t0}s')
    agresults['final_reduced_data']=df_final
    agresults['final_meta_data']=df_meta
    return agresults

def load_ds(url):
    """
    """
    ds=f63_to_xr(url)
    if 'fort.63' in url:
        variable_name='zeta'
    elif 'swan_HS' in url:
        variable_name='swan_HS'
    else:
        print(f'load_ds error: cannot deterine variable_name {url}')
        sys.exit(1)
    return ds, variable_name

def Combined_pipeline(url, geopoints, nearest_neighbors=10):
    """
    Simply run the series of steps as part of this method to shield users from some
    of the details of the processing
    """
    ds, variable_name = load_ds(url)
    agdict=get_adcirc_grid_from_ds(ds)
    agdict=attach_element_areas(agdict)

    print('Start KDTree pipeline')
    agdict=ComputeTree(agdict)
    agresults=ComputeQuery(geopoints, agdict, kmax=nearest_neighbors)
    agresults=ComputeBasisRepresentation(geopoints, agdict, agresults)
    agresults=ConstructReducedWaterLevelData_from_ds(ds, agdict, agresults, variable_name=variable_name)

    print(f'List of {len(agresults["outside_elements"])} stations not assigned to any grid element follows for kmax {nearest_neighbors}')
    print(geopoints.iloc[agresults['outside_elements']])

    df_product_data=agresults['final_reduced_data']
    df_product_metadata=agresults['final_meta_data']
    print('Finished')
    return df_product_data, df_product_metadata


