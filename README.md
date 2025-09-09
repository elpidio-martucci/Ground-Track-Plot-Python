# Ground-Track-Plot-Python
Python tool for computing and visualizing satellite ground track from orbital parameters

This is my first Python project. The main goal was to get comfortable with Python for scientific computing, while applying concepts I already know from aerospace engineering: orbital mechanics, reference frame transformations, and vector operations (previously practiced in MATLAB).

The script takes as input the orbital parameters of a satellite (sma, inc, raan, ecc, aop, ta) and computes the evolution of its ground track. Initial time (hours) and step resolution can also be adjusted.

Required libraries: numpy, matplotlib, cartopy.

Once the parameters are inserted, the logic order is:
 - computing 'r' for each true anomaly step
 - reference frame 3D rotation from orbital plane to ECI
 - reference frame 3D rotation from ECI to ECEF
 - coordinates transformation from ECEF to Geodetic, accounting for Earth ellipticity (WGS84 model)
 - plotting all the lat7long position over Earth map

(NOTE: This project was developed by me, with help from ChatGPT for debugging plotting functions and tidying up the code.
Future improvements may include support for multiple satellites, TLE inputs, and interactive Jupyter notebooks.)
