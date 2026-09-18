#!/usr/bin/env python
# --------------------------------------------------------------------------------------------------
# Renormalize an offline map so that each destination row's weights sum to 1.
#
# A map produced by GenerateTransposeMap is the flux adjoint of the forward map: it is
# conservative (area-weighted column sums are preserved) but NOT consistent (row sums are
# not 1), so a constant source field does not map to a constant destination field. That is
# the correct property for remapping fluxes and the wrong one for remapping state variables
# (T, ps, qv, winds), which is what initial-condition generation does.
#
# This rescales each row by 1/rowsum, which is the standard conversion from a conservative
# map to a consistent one. The trade is explicit: the result is consistent but no longer
# exactly conservative. Keep the original file for any flux application.
#
# Only rows that actually deviate (|rowsum-1| > tol) are touched, so every well-formed row
# stays bit-identical and a diff against the original isolates exactly what changed.
# --------------------------------------------------------------------------------------------------
import os
import sys
import shutil
import argparse
import numpy as np
import netCDF4 as nc

CHUNK = 50_000_000

# --------------------------------------------------------------------------------------------------
# weight accumulation
# --------------------------------------------------------------------------------------------------
def compute_row_sums(ds, n_s, n_b):
  """accumulate the sum of weights for each destination row, in chunks so a map with
  hundreds of millions of entries never has to be resident all at once"""
  rowsum = np.zeros(n_b, dtype='float64')
  for i0 in range(0, n_s, CHUNK):
    i1 = min(i0+CHUNK, n_s)
    row = ds.variables['row'][i0:i1].astype('int64') - 1
    S   = ds.variables['S'][i0:i1].astype('float64')
    rowsum += np.bincount(row, weights=S, minlength=n_b)
    print(f'    row sums: {i1}/{n_s}', flush=True)
  return rowsum

# --------------------------------------------------------------------------------------------------
def main():
  p = argparse.ArgumentParser(description='renormalize an offline map to be row-consistent')
  p.add_argument('--in_map',  required=True, help='input map file (not modified)')
  p.add_argument('--out_map', required=True, help='output map file (created)')
  p.add_argument('--tol', type=float, default=1e-12,
                 help='rows with |rowsum-1| <= tol are left bit-identical (default 1e-12)')
  p.add_argument('--dry_run', action='store_true', help='report what would change, write nothing')
  args = p.parse_args()

  if os.path.abspath(args.in_map) == os.path.abspath(args.out_map):
    raise ValueError('in_map and out_map must differ - the original is needed for flux use')

  # ----------------------------------------------------------------------------
  # pass 1 - measure
  # ----------------------------------------------------------------------------
  print(f'reading {args.in_map}', flush=True)
  with nc.Dataset(args.in_map, 'r') as ds:
    n_s = len(ds.dimensions['n_s'])
    n_b = len(ds.dimensions['n_b'])
    print(f'  n_s={n_s}  n_b={n_b}', flush=True)
    rowsum = compute_row_sums(ds, n_s, n_b)

  if not np.all(np.isfinite(rowsum)):
    n_invalid = int((~np.isfinite(rowsum)).sum())
    raise ValueError(f'{n_invalid} rows have a non-finite weight sum and cannot be rescaled')
  dev = np.abs(rowsum - 1.0)
  bad = dev > args.tol
  print()
  print(f'  row sums : min={rowsum.min():.8f} max={rowsum.max():.8f}')
  print(f'  rows to rescale (|rowsum-1|>{args.tol:g}): {int(bad.sum())} ({100*bad.mean():.4f}%)')
  if rowsum.min() <= 0:
    n_zero = int((rowsum <= 0).sum())
    raise ValueError(f'{n_zero} rows have a non-positive weight sum and cannot be rescaled; '
                     f'these indicate a genuine mesh/overlap failure, not a normalization issue')
  if args.dry_run:
    print('  dry run - nothing written'); return

  # scale factor of exactly 1.0 for every row we are leaving alone, so those weights
  # come through the multiply bit-identical
  scale = np.ones(n_b, dtype='float64')
  scale[bad] = 1.0/rowsum[bad]

  # ----------------------------------------------------------------------------
  # pass 2 - copy, then rescale the weights in place
  # ----------------------------------------------------------------------------
  print(f'\ncopying to {args.out_map}', flush=True)
  shutil.copyfile(args.in_map, args.out_map)

  with nc.Dataset(args.out_map, 'a') as ds:
    for i0 in range(0, n_s, CHUNK):
      i1 = min(i0+CHUNK, n_s)
      row = ds.variables['row'][i0:i1].astype('int64') - 1
      S   = ds.variables['S'][i0:i1].astype('float64')
      ds.variables['S'][i0:i1] = S*scale[row]
      print(f'    rescaling: {i1}/{n_s}', flush=True)

    # frac_b records the fraction of each destination cell covered by the source grid,
    # which for a consistent map is unity wherever the row was rescaled
    if 'frac_b' in ds.variables:
      frac_b = ds.variables['frac_b'][:]
      frac_b[bad] = 1.0
      ds.variables['frac_b'][:] = frac_b
      print('    frac_b updated for rescaled rows', flush=True)

    hist = f'renormalized row sums to 1 ({int(bad.sum())} rows) for state-variable remapping; ' \
           f'no longer exactly conservative - see {os.path.basename(args.in_map)} for the flux form'
    ds.setncattr('history', hist + '\n' + getattr(ds, 'history', ''))
    ds.setncattr('renormalized_from', os.path.abspath(args.in_map))

  # ----------------------------------------------------------------------------
  # verify
  # ----------------------------------------------------------------------------
  print('\nverifying', flush=True)
  with nc.Dataset(args.out_map, 'r') as ds:
    rowsum2 = compute_row_sums(ds, n_s, n_b)
  dev2 = np.abs(rowsum2 - 1.0)
  print()
  print(f'  row sums after: min={rowsum2.min():.10f} max={rowsum2.max():.10f}')
  print(f'  rows still off 1.0 by >1e-10: {int((dev2>1e-10).sum())}')
  print('  PASS' if (dev2 > 1e-10).sum() == 0 else '  FAIL - inspect before using')

# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  main()
