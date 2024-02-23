#!/usr/bin/env python3
import xarray, numpy, dask, sys
from subprocess import call

'''
usage: convert_eami_to_screami.py [-h] inputfile topofile gasfile outputfile

topofile=/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne256np4_16xdel2_consistentSGH_20220429.nc
gasfile=/global/cfs/cdirs/e3sm/whannah/init_data/gas_constituents_ne256np4L128_20230222.nc
inputfile=/global/cfs/cdirs/m2637/jsgoodni/HICCUP.atm_era5.2006-08-09.Saomai_2006_ne128x8_lon130E_lat25N.L128.nc
outputfile=/global/cfs/cdirs/m2637/whannah/HICCUP.atm_era5.2006-08-09.Saomai_2006_ne128x8_lon130E_lat25N.L128.converted.nc
python3 convert_eami_to_screami.py $inputfile $topofile $gasfile $outputfile
'''

def print_message(message):
    print(message); sys.stdout.flush()


def main(inputfile, gasfile, outputfile):

    # Stuff we need to get from cami file
    screami_from_cami_coordvars = {  # These we want to be doubles
        'time': 'time',
        'hyam': 'hyam',
        'hybm': 'hybm',
        'hyai': 'hyai',
        'hybi': 'hybi',
        'pref_mid': 'lev', #hybm
        'lat' : 'lat',
        'lon' : 'lon',
        'P0'  : 'P0',
    }
    screami_from_cami_datavars = {  # These can be floats
        'qv'  : 'Q',
        'ps'  : 'PS',
        'T_mid': 'T',
        'qc': 'CLDLIQ',
        'qi': 'CLDICE',
        'qr': 'RAINQM',
        'nc': 'NUMLIQ',
        'ni': 'NUMICE',
        'nr': 'NUMRAI',
    }
    screami_from_cami = {**screami_from_cami_coordvars, **screami_from_cami_datavars}

    # Other things we need:
    #   gas concentrations for trace gas constituents (only o3 now)
    #   surface geopotential
    gas_names = ('o3',)

    with xarray.open_dataset(inputfile, decode_cf=False) as ds_in:

        ds_out = xarray.Dataset()

        ntime = ds_in.sizes['time']
        ncol = ds_in.sizes['ncol']
        nlev = ds_in.sizes['lev']

        # Fields we can copy directly
        print_message('Copy fields from cami')
        for scream_name, eam_name in screami_from_cami.items():
            print_message(f'  {eam_name} -> {scream_name}')
            if eam_name in ds_in.variables.keys():
                ds_out[scream_name] = ds_in[eam_name]
            else:
                print_message(f'WARNING: {eam_name} not found in inputfile; hope you did not need that one...')

        # Make sure ps field contains a time dimension
        if 'time' not in ds_out['ps'].dims:
            ds_out['ps'] = ds_out['ps'].expand_dims('time', axis=0)

        # Handle U,V to horiz_winds(time,ncol,dim2,lev)
        print_message('Map U,V to horiz_winds')
        ds_out['horiz_winds'] = xarray.concat([ds_in['U'], ds_in['V']], 'dim2')

        # Grab trace gas concentrations from file
        print_message('Get trace gases')
        with xarray.open_dataset(gasfile, decode_cf=False).isel(time=0).drop('time') as ds_gas:
            for v in gas_names:
                print_message(f'  get {v} from gas file...')
                ds_out[v + '_volume_mix_ratio'] = xarray.DataArray(
                    numpy.zeros([ncol, nlev], dtype=float), 
                    dims=('ncol', 'lev'),
                    attrs=ds_gas[v].attrs
                )
                ds_out[v + '_volume_mix_ratio'].values = ds_gas[v].transpose('ncol', 'lev').values
                ds_out[v + '_volume_mix_ratio'], *__ = xarray.broadcast(ds_out[v + '_volume_mix_ratio'], ds_out['ps'])
                ds_out[v + '_volume_mix_ratio'].attrs['units'] = 'mol/mol'
                ds_out[v + '_volume_mix_ratio'].attrs['long_name'] = v + ' volume mixing ratio'

        # Convert datavars to single precision float
        for v in ds_out.data_vars.keys():
            if v not in screami_from_cami_coordvars.keys():
                print_message(f'convert {v} to float...')
                ds_out[v] = ds_out[v].astype(numpy.float32)

        # Permute dimensions
        print_message('Permute dimensions')
        try:
            if 'nbnd' in ds_out.dims:
                ds_out = ds_out.transpose('time','ncol','dim2','lev','ilev','nbnd')
            else:
                ds_out = ds_out.transpose('time','ncol','dim2','lev','ilev')
        except:
            print_message('  permute dimensions failed, but continuing anyways.')

        # Apply clipping
        for v in ('qv', 'qc', 'qi', 'qr', 'nc', 'ni', 'nr'):
            if v in ds_out.variables.keys():
                ds_out[v].data = numpy.where(ds_out[v] < 0, 0, ds_out[v])

        # Check values
        for v in ds_out.variables:
            numnan = ds_out[v].isnull().sum().values
            print_message(f'{v} numnan = {numnan}')

        print_message('Write to netcdf')
        ds_out.to_netcdf(outputfile) #format="NETCDF3_64BIT")

        # Convert to cdf5
        print_message('Convert to cdf5')
        call(f'ncks -5 {outputfile} {outputfile}.cdf5'.split(' '))
        call(f'mv {outputfile}.cdf5 {outputfile}'.split(' '))

        print_message('Done.')

        print_message('')
        print_message(f'  outputfile: {outputfile}')
        print_message('')

if __name__ == '__main__':
    import plac; plac.call(main)