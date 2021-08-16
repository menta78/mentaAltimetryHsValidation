This is a "micro library" for validation of significant wave height (Hs) as reproduced by a model with Hs observed by satellite atlimeters.

It is developed in python, and requires netCDF4, numpy, scipy, matplotlib.

The subdirectory altHsValidation contain the library functions, while the subdirectory scripts contain a sample script (validateWwmUOST_ICE.py) which does the validation.

In validateWwmUOST_ICE.py 3 functions are invoked:

- coarsenSatData: the raw altimetry estimations of Hs come with high frequency (one every few seconds) and noise. This routine averages them on segments along the track with latitudinal extent
 given by the parameter latdelta (0.4 deg in the example). The results are saved in the crsSatDataDir in the numpy format. The result of this operation does not depend on the model, and can be
 performed once, unless the user wants to change the parameter latdelta, or the time extent of the sat data.
 The input data for this routine are provided by the GlobWave database
 (http://globwave.ifremer.fr/).<br />
!!!!!!!!!!!!!!!!!!!!!!!!!<br />
 FOR ARON: for our stuff you can skip this step, as you can find its output here:
   https://www.dropbox.com/sh/y8e364if6m9bs9f/AAAJySeZRLuGXixHcEOuo186a?dl=0 <br />
!!!!!!!!!!!!!!!!!!!!!!!!!

- interpolateModelToCoarsenedSatData_schismWWM: this routine, located in the module interpolateModelToCoarsenedSatData, interpolates the model results at the observation points, along the satellite tracks. This function does the job with the results of the schismWWM model, while similar functions (interpolateModelToCoarsenedSatData_WW3 and interpolateModelToCoarsenedSatData_WWMExperimental) do the job for WW3 and the standalone WWM. The results (the couple of obsrved and modelled Hs along the tracks) are saved in the numpy binary format.

- function elaborateMeasures in module computeStats: elaborates the statistics and saves them in a text format in the directory statsDir. The statistics are generated on a global regular mesh, latitudinally limited by the latlims parameter, and with x,y steps dx, dy
    These files are saved by elaborateMeasures:
    - lons.csv, lats.csv: coordinates of the cells of the statistics grid.
    - dtcount.csv, grid with the data count for each cell.
    - bias.csv, grid of normalized bias.
    - absbias.csv, grid of normalized absolute bias.
    - nrmse.csv, grid of nrmse.
    - hh.csv, grid of hh as defined by Mentaschi et al. 2013.
    - modelYearMaxMean.csv, grid of model mean annual maxima (useful for multi-year validations).
    - obsYearMaxMean.csv, grid of observation mean annual maxima (useful for multi-year validations).
    - modelTotMax.csv, grid of model total maxima.
    - obsTotMax.csv, grid of observation total maxima.
    - totalIndicators.csv, list with global statistical indicators.
 

