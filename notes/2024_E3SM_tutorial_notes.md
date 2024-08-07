
# Let’s Create an Initial Condition!

Clone HICCUP repo
`git clone https://github.com/E3SM-Project/HICCUP.git `


Move into HICCUP directory
```shell
cd HICCUP/
```

1. Activate the E3SM unified conda environment
```shell
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```

1. Install HICCUP as a package
```shell
pip install ./
```

1. Copy the template HICCUP script to the user_scripts directory
```shell
cp  template_hiccup_scripts/E3SM_tutorial_create_IC_from_ERA5-NOAA.py  user_scripts/
```

1. Edit the script to modify output_root 
```shell
vim user_scripts/E3SM_tutorial_create_IC_from_ERA5-NOAA.py
```

1. Run HICCUP to Generate files
```shell
python user_scripts/E3SM_tutorial_create_IC_from_ERA5-NOAA.py
```


# Now Let’s Run E3SM!

1. Deactivate the conda environment - or create a new session - Otherwise you will get errors!!!
```shell
conda deactivate
```

1. Make a copy the template run script
```shell
cp  template_run_scripts/run_E3SM.2024_tutorial_hindcast.nersc.py  user_scripts/
```

1. Edit the file to update allocation code and source code path
```shell
vim user_scripts/run_E3SM.2024_tutorial_hindcast.nersc.py
```
  - `acct        = ntrain6`
  - `reservation = e3sm_day4`
  - `src_dir     = ????`


1. Execute the run script
```shell
user_scripts/run_E3SM.2024_tutorial_hindcast.nersc.py  
```

1. Use this command to monitor the queue:
  ```shell
  squeue --user=${USER}
  ```
  Or if you want to get fancy:
  ```shell
  squeue --user=${USER} --format="%.18i  %.8u  %.70j  %.9P  %.8T  %.10M  %.12l  %.6D  %V "
  ```

