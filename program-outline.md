# Program outline

this wiki explains ideas on how the program should be organized.

## Programs

1. computePaths
  * in: params
  * out: list stats
2. minimizeDist
  * in: final stats
  * out: params
3. calcStats
  * in: multiple out stats
  * out: final stats

## data structures

### out-stats
This class will need to hold lists of data that is output from the simulation. Also the number of loops will need to be saved as to normalize the lists with other lists correctly

* ev{0,4} : lists of length Nb
  * list of floats (length Nb)
  * holds the KL distance calculations
* xbarPath
  * list of floats (length Nb)
  * holds the mean path for the run
* xxbarPath
  * list of floats (length Nb)
  * holds the average position squared
  * used to compute the standard deviation of the path (along with the xbarPath)
* loops
  * integer
  * holds the number of loops that were calculated 
  * used for normailization

### params
This class acts more like a struct and holds the parameters for the run. This info will be passed in at the beginning of the run and minimizeDist routine will calculate the next runs parameters.

* mean
  * mplus
  * gammam
  * tm
  * gammamH
  * mHlist (list of Hermite coefficients) 
* A
  * Aplus
  * Amid
  * gammaA
  * tA
  * gammaAH
  * AHlist (list of Hermite coefficients) 
