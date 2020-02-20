#!/usr/bin/env python
import xarray as xr
import numpy as np
import ngl
import os

def main(fig_file='sanity_check',fig_type='png'):

   var = ['PS','T','Q','U']

   klev = -5

   ds = xr.open_dataset('HICCUP_TEST.output.nc' )
   scrip_ds = xr.open_dataset('scrip_ne30pg2.nc')

   plot = []
   wks = ngl.open_wks(fig_type,fig_file)
   res = get_resources()
   res.cnFillMode    = 'CellFill'
   res.sfXCellBounds = scrip_ds['grid_corner_lon'].values
   res.sfYCellBounds = scrip_ds['grid_corner_lat'].values

   for v in range(len(var)):
      data = ds[var[v]].isel(time=0)
      if 'lev' in data.dims: 
         data = data.isel(lev=klev)
         plev = ds['lev'].isel(lev=klev).values
         lev_str = f'{plev:6.2f} hPa'
      else:
         lev_str = ''

      print(f'\nvar: {var[v]}')
      print(f'  min : {data.min().values} ')
      print(f'  avg : {data.mean().values} ')
      print(f'  max : {data.max().values} ')

      plot.append( ngl.contour_map(wks,data.values,res) )
      set_subtitles(wks, plot[len(plot)-1], data.attrs['long_name'], '', lev_str )

   if len(plot)<4 :
      layout = [len(plot),1]
   else:
      pdim = np.ceil(len(plot)/2)
      layout = [pdim,pdim]

   ngl.panel(wks,plot,layout,False)
   ngl.end()

   ### trim white space from image using imagemagik
   if fig_type == 'png' :
      os.system(f'convert -trim +repage {fig_file}.png {fig_file}.png')
      print(f'\n{fig_file}.png\n')

#---------------------------------------------------------------------------------------------------
# define function to add subtitles to the top of plot
#---------------------------------------------------------------------------------------------------
def set_subtitles(wks, plot, left_string='', center_string='', right_string='', font_height=0.015):
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
def get_resources():
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
   res.lbOrientation            = "Horizontal"
   res.lbLabelFontHeightF       = 0.015
   res.mpGridAndLimbOn          = False
   res.mpCenterLonF             = 180
   res.mpLimitMode              = "LatLon" 
   return res
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__': main()