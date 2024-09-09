# constant parameter values for HICCUP calculations

tk_zero = 273.15 # value for converting between celsius and Kelvin

std_lapse = 0.0065           # std. atmosphere lapse rate              ~ -6.5 K/km
gravit    = 9.80616          # acceleration of gravity                 ~ m/s^2
boltz     = 1.38065e-23      # boltzmann's constant                    ~ J/k/molecule
avogad    = 6.02214e26       # avogadro's number                       ~ molecules/kmole
MW_dryair = 28.966           # molecular weight of dry air             ~ g/mol
MW_ozone  = 47.998           # molecular weight of ozone               ~ g/mol
MW_vapor  = 18.0             # molecular weight of dry air             ~ g/mol
Rgas      = avogad*boltz     # universal gas constant                  ~ J/k/kmole
Rdair     = Rgas/MW_dryair   # gas constant for dry air                ~ J/k/kg
Rvapor    = Rgas/MW_vapor    # gas constant for water vapor            ~ J/(kg/kg)
P0        = 1e5              # reference pressure

