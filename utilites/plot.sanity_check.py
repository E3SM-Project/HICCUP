#!/usr/bin/env python
# ------------------------------------------------------------------------------
# This scripts will plot the data output from HICCUP as a sanity check (requires PyNGL)
# ------------------------------------------------------------------------------
import xarray as xr, numpy as np, ngl, os
from optparse import OptionParser
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
   font_height = 0.015/np.sqrt(len(var))

   #----------------------------------------------------------------------------
   # Create dataset objects
   ds = xr.open_dataset(ifile)
   scrip_ds = xr.open_dataset(gfile)

   # print(); print(ds[var[0]])
   # print(); print(scrip_ds)
   # exit()
   
   #----------------------------------------------------------------------------
   # Set up plot stuff
   plot = []
   wks = ngl.open_wks(fig_type,fig_file)
   res = get_resources(map_plot=True)
   res.cnFillMode    = 'CellFill'
   # res.sfXCellBounds = scrip_ds['grid_corner_lon'].values
   # res.sfYCellBounds = scrip_ds['grid_corner_lat'].values
   res.sfXArray         = scrip_ds.variables['grid_center_lon'].values
   res.sfXCellBounds    = scrip_ds.variables['grid_corner_lon'].values
   res.sfYArray         = scrip_ds.variables['grid_center_lat'].values
   res.sfYCellBounds    = scrip_ds.variables['grid_corner_lat'].values

   # separate resources for zonal mean plot
   res2 = get_resources()
   res2.vpHeightF = 0.4
   res2.trYReverse = True
   res2.tiXAxisString = 'Latitude'
   res2.tiYAxisString = 'Pressure [hPa]'

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
      plot.append( ngl.contour_map(wks,data.values,res) )
      if 'long_name' in data.attrs:
         left_str = data.attrs['long_name']
      else:
         left_str = var[v]
      set_subtitles(wks, plot[len(plot)-1], left_str, '', lev_str, font_height=font_height )
      
      # if 'lev' in ds[var[v]].dims :
      #    # Create zonal mean plot (height vs lat) using area-weighted averaging of columns
      #    bin_ds = bin_YbyX( ds[var[v]].isel(time=0), ds['lat'], bin_min=-88, bin_max=88, bin_spc=2, wgt=ds['area'] )

      #    # Use sin(lat) axis to minimize distortion
      #    sin_lat_bins = np.sin( bin_ds['bins'].values*np.pi/180. )
      #    lat_tick = np.array([-90,-60,-30,0,30,60,90])
      #    res2.tmXBMode, res2.tmXBValues, res2.tmXBLabels = "Explicit", np.sin( lat_tick*3.14159/180. ), lat_tick

      #    res2.sfXCStartV, res2.sfXCEndV = -1.0, 1.0
      #    res2.sfYCStartV, res2.sfYCEndV = min( bin_ds['lev'].values ), max( bin_ds['lev'].values )

      #    plot.append( ngl.contour(wks, bin_ds['bin_val'].transpose().values, res2) )
      #    set_subtitles(wks, plot[len(plot)-1], data.attrs['long_name'], '', 'Zonal Mean', font_height=font_height )

      # else: 

      #    # Create map plot
      #    plot.append( ngl.contour_map(wks,data.values,res) )
      #    set_subtitles(wks, plot[len(plot)-1], data.attrs['long_name'], '', lev_str, font_height=font_height )

   #----------------------------------------------------------------------------
   # Combine plots panels
   if len(plot)<4 :
      layout = [len(plot),1]
   else:
      # layout = [ np.ceil(np.sqrt(len(plot))), np.ceil(np.sqrt(len(plot))) ]    # square layout
      layout = [np.ceil(len(plot)),2]     # 2-column layout

   pres = ngl.Resources()
   pres.nglPanelXWhiteSpacePercent = 10
   pres.nglPanelYWhiteSpacePercent = 10
   ngl.panel(wks,plot,layout,pres)
   ngl.end()

   #----------------------------------------------------------------------------
   # trim white space from image using imagemagik
   if fig_type == 'png' :
      os.system(f'convert -trim +repage {fig_file}.png {fig_file}.png')
      print(f'\n{fig_file}.png\n')

#---------------------------------------------------------------------------------------------------
# define function to add subtitles to the top of plot
#---------------------------------------------------------------------------------------------------
def set_subtitles(wks, plot, left_string='', center_string='', right_string='', font_height=0.01):
   ttres         = ngl.Resources()
   ttres.nglDraw = False

   ### Use plot extent to call ngl.text(), otherwise you will see this error: 
   ### GKS ERROR NUMBER   51 ISSUED FROM SUBROUTINE GSVP  : --RECTANGLE DEFINITION IS INVALID
   strx = ngl.get_float(plot,'trXMinF')
   stry = ngl.get_float(plot,'trYMinF')
   ttres.txFontHeightF = font_height

   # Set annotation resources to describe how close text is to be attached to plot
   amres = ngl.Resources()
   amres.amOrthogonalPosF = -0.52   # Top of plot plus a little extra to stay off the border
   if hasattr(ttres,'amOrthogonalPosF'): amres.amOrthogonalPosF = ttres.amOrthogonalPosF   

   # Add left string
   amres.amJust,amres.amParallelPosF = 'BottomLeft', -0.5   # Left-justified
   tx_id_l   = ngl.text(wks, plot, left_string, strx, stry, ttres)
   anno_id_l = ngl.add_annotation(plot, tx_id_l, amres)
   # Add center string
   amres.amJust,amres.amParallelPosF = 'BottomCenter', 0.0   # Centered
   tx_id_c   = ngl.text(wks, plot, center_string, strx, stry, ttres)
   anno_id_c = ngl.add_annotation(plot, tx_id_c, amres)
   # Add right string
   amres.amJust,amres.amParallelPosF = 'BottomRight', 0.5   # Right-justified
   tx_id_r   = ngl.text(wks, plot, right_string, strx, stry, ttres)
   anno_id_r = ngl.add_annotation(plot, tx_id_r, amres)

   return
#---------------------------------------------------------------------------------------------------
# Define function for setting the plot resources
#---------------------------------------------------------------------------------------------------
def get_resources(map_plot=False):
   res = ngl.Resources()
   res.nglDraw                  = False
   res.nglFrame                 = False
   res.tmXTOn                   = False
   res.tmXBMajorOutwardLengthF  = 0.
   res.tmXBMinorOutwardLengthF  = 0.
   res.tmYLMajorOutwardLengthF  = 0.
   res.tmYLMinorOutwardLengthF  = 0.
   res.tmYLLabelFontHeightF     = 0.015
   res.tmXBLabelFontHeightF     = 0.015
   res.tiXAxisFontHeightF       = 0.015
   res.tiYAxisFontHeightF       = 0.015
   res.tmXBMinorOn              = False
   res.tmYLMinorOn              = False
   res.cnFillPalette            = "MPL_viridis"
   res.cnFillOn                 = True
   res.cnLinesOn                = False
   res.cnLineLabelsOn           = False
   res.cnInfoLabelOn            = False
   res.lbLabelFontHeightF       = 0.015
   if map_plot==True :
      res.lbOrientation            = "Horizontal"
      res.mpGridAndLimbOn          = False
      res.mpCenterLonF             = 180
      res.mpLimitMode              = "LatLon" 
   return res
#---------------------------------------------------------------------------------------------------
# Binning routine for calculating zonal mean on unstructured grid
#---------------------------------------------------------------------------------------------------
def bin_YbyX (Vy,Vx,bins=[],bin_min=0,bin_max=1,bin_spc=1,wgt=[],keep_time=False):
   """ Average Vy into bins of Vx values. """
   #----------------------------------------------------------------------------
   # use min, max, and spc (i.e. stride) to define bins   
   nbin    = np.round( ( bin_max - bin_min + bin_spc )/bin_spc ).astype(np.int)
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
   help_header += '\nThis scripts will plot the data output from HICCUP as a sanity check (requires PyNGL)'
   parser = OptionParser(usage=help_header)
   parser.add_option('-i',dest='ifile',default=None,help='input file name')
   parser.add_option('--grid_file',dest='gfile',default=None,help='grid file name')
   parser.add_option('--vars',dest='vars',default=None,help='comma separated list of variables to plot')
   (opts, args) = parser.parse_args()

   var_list = None
   if opts.vars is not None:
      var_list = opts.vars.split(',')

   main(ifile=opts.ifile,gfile=opts.gfile,var_list=var_list)
