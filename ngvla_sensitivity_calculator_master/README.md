ngVLA sensitivity calculator written in Python by Viviana Rosero.


This version uses  Tsys and aperture efficiency values that are frequency dependent (from Wes Grammer's front-end cascade model 2021.10.20). Also, it calculates the line sensitivity and brightness temperature at the desired frequency instead of the band center. At band 6 it tries to center the 20 GHz of instantaneous bandwidth on the desired frequency instead of using the band average.


Requires:   scipy version >= 1.0.0

It should work in Python versions 2 and 3 provided the above
requirement is satisfied.  It will not work in CASA version 5 because
the scipy version is too old, but will work in CASA version 6.

It can be run from the command line or the individual functions
can be imported and used directly.

See the command line help for further instructions, i.e.,

./ngVLA_sensitivity_calculator.py -h
