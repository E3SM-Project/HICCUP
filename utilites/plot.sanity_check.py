#!/usr/bin/env python
# ------------------------------------------------------------------------------
# This script will plot the data output from HICCUP as a sanity check (uses matplotlib)
# ------------------------------------------------------------------------------
import xarray as xr, numpy as np, os
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from optparse import OptionParser

# Coastlines are drawn when cartopy is available, but it is not required
try:
   import cartopy.crs as ccrs
   HAVE_CARTOPY = True
except ImportError:
   HAVE_CARTOPY = False
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def main(fig_file='sanity_check',fig_type='png',ifile=None,gfile=None,var_list=None):

   if ifile is None: ifile = '../data/HICCUP_TEST.output.atm.nc'
   if gfile is None: gfile = '../files_grid/scrip_ne30np4.nc'

   # Specify list of variables to plot
   if var_list is None:
      var = ['PS']
      # var = ['PS','TS','T','Q','U','CLDLIQ']
      # var = ['T','Q','U','V']
      # var = ['t','q','u','v']
   else:
      var = var_list

   # specify level to use for data with "lev" dimension
   # starts at TOA, but negative values can be used to start from surface
   klev = -20

   # Make plot subtitle font size vary with number of plot panels
   title_fs = max(7, 14/np.sqrt(len(var)))

   #----------------------------------------------------------------------------
   # Create dataset objects
   ds = xr.open_dataset(ifile)
   scrip_ds = xr.open_dataset(gfile)

   # grid cell corner/center coordinates used to draw the unstructured cells
   corner_lon = scrip_ds['grid_corner_lon'].values
   corner_lat = scrip_ds['grid_corner_lat'].values
   center_lon = scrip_ds['grid_center_lon'].values

   # print(); print(ds[var[0]])
   # print(); print(scrip_ds)
   # exit()

   #----------------------------------------------------------------------------
   # Set up plot layout
   if len(var)<4 :
      ncol = 1
   else:
      ncol = 2                          # 2-column layout
   nrow = int(np.ceil(len(var)/ncol))

   fig = plt.figure(figsize=(7*ncol,3.5*nrow))
   proj_kw = {'projection':ccrs.PlateCarree(central_longitude=180)} if HAVE_CARTOPY else {}

   #----------------------------------------------------------------------------
   # load data and create plot
   for v in range(len(var)):

      data = ds[var[v]].isel(time=0)

      lev_str = ''
      lev_name = None
      if 'lev' in data.dims : lev_name = 'lev'
      if 'plev' in data.dims : lev_name = 'plev'
      if 'level' in data.dims : lev_name = 'level'
      if lev_name is not None:
         data = data.isel({lev_name:klev})
         plev = ds[lev_name].isel({lev_name:klev}).values
         lev_str = f'{plev:6.2f} hPa'

      # Print some statistics of the data
      print(f'\nvar: {var[v]}')
      print(f'  min : {data.min().values} ')
      print(f'  avg : {data.mean().values} ')
      print(f'  max : {data.max().values} ')

      # exit()

      # Create map plot
      ax = fig.add_subplot(nrow,ncol,v+1,**proj_kw)
      pc = add_cell_fill(ax, corner_lon, corner_lat, center_lon, data.values)

      if HAVE_CARTOPY:
         ax.coastlines(linewidth=0.5)
         ax.set_global()
      else:
         ax.set_xlim(0,360); ax.set_ylim(-90,90)
         ax.set_aspect('equal')

      cbar = fig.colorbar(pc, ax=ax, orientation='horizontal', pad=0.05, shrink=0.9)
      cbar.ax.tick_params(labelsize=title_fs)

      left_str = data.attrs['long_name'] if 'long_name' in data.attrs else var[v]
      ax.set_title(left_str, loc='left',  fontsize=title_fs)
      ax.set_title(lev_str,  loc='right', fontsize=title_fs)

   #----------------------------------------------------------------------------
   # Combine plot panels and save
   fig.tight_layout()
   fig.savefig(f'{fig_file}.{fig_type}',dpi=150,bbox_inches='tight')
   plt.close(fig)
   print(f'\n{fig_file}.{fig_type}\n')

#---------------------------------------------------------------------------------------------------
# Draw unstructured-grid data as filled polygons using the SCRIP cell corner bounds
#---------------------------------------------------------------------------------------------------
def add_cell_fill(ax, corner_lon, corner_lat, center_lon, values, cmap='viridis'):
   """ Draw each grid cell as a filled polygon (equivalent to NGL CellFill) """
   # unwrap each cell's corner longitudes relative to its center so cells that
   # straddle the dateline do not stretch all the way across the plot
   lon = corner_lon - 360.0*np.round( (corner_lon - center_lon[:,None])/360.0 )
   verts = np.stack([lon, corner_lat], axis=-1)   # (ncell, ncorner, 2)
   values = np.asarray(values)
   kw = {'transform':ccrs.PlateCarree()} if HAVE_CARTOPY else {}
   pc = PolyCollection(verts, array=values, cmap=cmap, edgecolors='face', linewidths=0.0, **kw)
   pc.set_clim(np.nanmin(values), np.nanmax(values))
   ax.add_collection(pc)
   return pc
#---------------------------------------------------------------------------------------------------
# Binning routine for calculating zonal mean on unstructured grid
#---------------------------------------------------------------------------------------------------
def bin_YbyX (Vy,Vx,bins=[],bin_min=0,bin_max=1,bin_spc=1,wgt=[],keep_time=False):
   """ Average Vy into bins of Vx values. """
   #----------------------------------------------------------------------------
   # use min, max, and spc (i.e. stride) to define bins
   nbin    = np.round( ( bin_max - bin_min + bin_spc )/bin_spc ).astype(int)
   bins    = np.linspace(bin_min,bin_max,nbin)
   bin_coord = xr.DataArray( bins )
   #----------------------------------------------------------------------------
   # create output data arrays
   nlev  = len(Vy['lev'])  if 'lev'  in Vy.dims else 1
   ntime = len(Vy['time']) if 'time' in Vy.dims else 1
   if ntime==1 and keep_time==True : keep_time = False

   shape,dims,coord = (nbin,),'bin',[('bin', bin_coord)]
   if nlev >1 and keep_time==False : shape,coord,dims = (nbin,nlev), [('bin',bin_coord),('lev',Vy['lev'])], ['bin','lev']
   if nlev==1 and keep_time==False : shape,dims,coord = (nbin,),'bin',[('bin',bin_coord)]

   mval = np.nan
   bin_val = xr.DataArray( np.full(shape,mval,dtype=Vy.dtype), coords=coord, dims=dims )
   bin_std = xr.DataArray( np.full(shape,mval,dtype=Vy.dtype), coords=coord, dims=dims )
   bin_cnt = xr.DataArray( np.zeros(shape,    dtype=Vy.dtype), coords=coord, dims=dims )
   #----------------------------------------------------------------------------
   levchk = False
   if 'lev' in Vy.dims and len(Vy.lev)>1 : levchk = True

   if levchk :
      avg_dims = ['ncol']
      if 'time' in Vy.dims : avg_dims = ['time','ncol']
      avg_dims_wgt = ['ncol']

   val_chk = np.isfinite(Vx.values)
   #----------------------------------------------------------------------------
   # Loop through bins
   for b in range(nbin):
      bin_bot = bin_min - bin_spc/2. + bin_spc*(b  )
      bin_top = bin_min - bin_spc/2. + bin_spc*(b+1)

      condition = xr.DataArray( np.full(Vx.shape,False,dtype=bool), coords=Vx.coords )
      condition.values = ( np.where(val_chk,Vx.values,bin_bot-1e3) >=bin_bot ) \
                        &( np.where(val_chk,Vx.values,bin_bot-1e3)  <bin_top )

      if np.sum(condition)>0 :
         if levchk :
            if len(wgt)==0 :
               bin_val[b,:] = Vy.where(condition,drop=True).mean( dim=avg_dims, skipna=True )
            else:
               if wgt.dims != Vy.dims :
                  wgt, *__ = xr.broadcast(wgt, Vy)
                  if 'time' in Vy.dims :
                     wgt = wgt.transpose('time','lev','ncol')
                  else :
                     wgt = wgt.transpose('lev','ncol')
               if 'time' in Vy.dims :
                  bin_val[b,:] = ( (Vy*wgt).where(condition,drop=True).sum( dim='ncol', skipna=True ) \
                                      / wgt.where(condition,drop=True).sum( dim='ncol', skipna=True ) ).mean(dim='time', skipna=True )
               else:
                  bin_val[b,:] = ( (Vy*wgt).where(condition,drop=True).sum( dim='ncol', skipna=True ) \
                                      / wgt.where(condition,drop=True).sum( dim='ncol', skipna=True ) )
            bin_std[b,:] = Vy.where(condition,drop=True).std(  dim=avg_dims, skipna=True )
            bin_cnt[b,:] = Vy.where(condition,drop=True).count(dim=avg_dims)
         else:
            bin_val[b] = Vy.where(condition).mean(skipna=True)
            bin_std[b] = Vy.where(condition).std()
            bin_cnt[b] = np.sum( condition )
   #----------------------------------------------------------------------------
   # use a dataset to hold all the output
   dims = ('bins','lev') if levchk else ('bins',)
   bin_ds = xr.Dataset()
   bin_ds['bin_val'] = (dims, bin_val )
   bin_ds['bin_std'] = (dims, bin_std )
   bin_ds['bin_cnt'] = (dims, bin_cnt )
   bin_ds['bin_pct'] = (dims, bin_cnt/bin_cnt.sum()*1e2 )
   bin_ds.coords['bins'] = ('bins',bin_coord)
   if levchk : bin_ds.coords['lev'] = ( 'lev', xr.DataArray(Vy['lev']) )
   #----------------------------------------------------------------------------
   return bin_ds
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
   # Parse the command line options
   help_header = 'usage: ./%prog [file] [file] ...\n'
   help_header += '\nThis script will plot the data output from HICCUP as a sanity check (uses matplotlib)'
   parser = OptionParser(usage=help_header)
   parser.add_option('-i',dest='ifile',default=None,help='input file name')
   parser.add_option('--grid_file',dest='gfile',default=None,help='grid file name')
   parser.add_option('--vars',dest='vars',default=None,help='comma separated list of variables to plot')
   (opts, args) = parser.parse_args()

   var_list = None
   if opts.vars is not None:
      var_list = opts.vars.split(',')

   main(ifile=opts.ifile,gfile=opts.gfile,var_list=var_list)
