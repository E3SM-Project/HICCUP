#!/usr/bin/env python
# ==================================================================================================
# HICCUP - Hindcast Initial Condition Creation Utility/Processor
# Create an ensemble of perturbed initial condition files from a single final
# IC file. Each member is a copy of the input file with small random
# perturbations added to the state variables; all other variables and file
# metadata are preserved byte-for-byte.
#
# The perturbations can be spatially correlated (low-pass filtered so they look
# synoptic-scale and coherent) rather than grid-point noise, which is useful for
# hindcast ensembles - e.g. of heat waves. Set spatially_correlated=False to fall
# back to the original independent grid-point perturbations.
# ==================================================================================================
import os
from hiccup import hiccup
# ------------------------------------------------------------------------------
# User settings

# final IC file to perturb (output of a normal HICCUP workflow)
input_file = os.getenv('SCRATCH')+'/HICCUP/data/HICCUP.eam_i_ne30np4_L80.nc'

# directory and file name pattern for the ensemble members
output_dir     = os.getenv('SCRATCH')+'/HICCUP/data/ensemble'
output_pattern = 'HICCUP.eam_i_ne30np4_L80.pert{member:03d}.nc'

# number of ensemble members to create
num_members = 100

# variables to perturb (only those present in the file are perturbed)
var_list = ['T','PS','U','V']

# spatial coherence of the perturbations
spatially_correlated = True   # False => original grid-point (IID) noise
corr_length_km       = 1000.  # approx. 1/e spatial correlation length [km]

# ------------------------------------------------------------------------------
# Create the ensemble

if not os.path.exists(output_dir): os.makedirs(output_dir)

# Reuse the smoothing operator across members so the (relatively expensive)
# KD-tree is only built once. create_perturbed_file returns the smoother, which
# we pass back in on each subsequent call. seed=member gives reproducible,
# distinct members.
smoother = None
for member in range(num_members):
    output_file = f'{output_dir}/'+output_pattern.format(member=member)

    smoother = hiccup.create_perturbed_file(input_file, output_file,
                                            var_list=var_list,
                                            spatially_correlated=spatially_correlated,
                                            corr_length_km=corr_length_km,
                                            seed=member,
                                            smoother=smoother,
                                            clobber=True,
                                            verbose=True)

# ------------------------------------------------------------------------------
print(f'\ncreated {num_members} perturbed members in {output_dir}\n')
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
