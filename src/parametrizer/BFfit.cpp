/*

This is version 1.4 of the C++ BusyFunc library that fits the Busy Function to data.
See Westmeier, Jurek, Obreschkow, Koribalski & Staveley-Smith (2013) for more details about the implementation
and the Busy Function.

v1.4 updates
------------
1. Removed the (>= 0.0) constraint on the power-law scaling.
2. Changed the 2 <= n <= 8 power-law index constraint to 1 <= n <= 8.
3. BUG FIX: The partial derivative of the BF with respect to the power-law exponent, n, was previously over-estimated 
in magnitude by 33%. 
4. BUG FIX: The function ApproxObsCovar was using an incorrect range when applying/removing the power-law exponent's 
internal re-mapping.

v1.3 updates
------------
None.

v1.2 updates
------------
1. Created a version of the v1.1 C/C++ library that doesn't use OMP.
2. Python setup.py script uses the files, PyBFfit.h and libPyBFfit.a. These are symbolic links to either the 
OMP or non-OMP enabled C libraries.

v1.1 updates
------------
1. Added the ability to generate random variants of a Busy Function fit from a fit's covariance matrix.
2. The number of LVM seeds and maximum iterations per seed *must* now be specified. Previously, 1000 seeds and 30 
iterations/seed were used.
3. A single, 8 parameter BF curve can be used as the sole LVM seed. This is achieved by populating the fit_params array 
and requesting -1 LVM seeds.
4. BUG FIX: In v1.0, the best fit of previous BF variants were changing the values of x_min and x_max when 
FitBusyFunc_engine was called for a given BF variant. This has been corrected. x_min and x_max are now calculated from 
the values mid and amp, which are calculated only once at the beginning of FitBusyFunc.
5. FitBusyFunc now returns an integer value corresponding to the BF variant that returned the best fit. This change was 
made, because the fit dimensionality is not a unique identifier. The fit dimensionality is now returned in a variable
that is passed by reference, best_NOp.

Email: Russell.Jurek@gmail.com

Version 1.4 created September 23rd, 2014 by Russell J. Jurek.
Version 1.3 created September 15th, 2014 by Russell J. Jurek.
Version 1.2 created December 26th, 2013 by Russell J. Jurek.
Version 1.1 created November 9th, 2013 by Russell J. Jurek.
Version 1.0 created June 30th, 2013 by Russell J. Jurek.

*/

#include "BFfit.h"

// macro-like inline functions used by NR-based functions --- called by the other functions in various places

inline float MAX(const double &a, const float &b)
        {return ((b > a) ? (b) : float(a));}

inline float MAX(const float &a, const double &b)
        {return ((b > a) ? float(b) : (a));}

inline float MIN(const double &a, const float &b)
        {return ((b < a) ? (b) : float(a));}

inline float MIN(const float &a, const double &b)
        {return ((b < a) ? float(b) : (a));}

inline float SIGN(const float &a, const double &b)
        {return ((b >= 0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

inline float SIGN(const double &a, const float &b)
        {return (float)((b >= 0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

// chain of functions that fit the Busy Function

void BusyFunc(double x, double a[], double *y, double dyda[], int fit_mode, double mid, double amp){

  double erf_1, erf_2, power_law;
  const double pi = 3.1415926535;
  double map_1, map_2, map_3, map_4, map_5, map_6, map_7, map_8;
  int i;

  //calculate the y value and partial derivatives for the different Busy Function models
  switch(fit_mode){
  case 1:
    map_1 = exp(a[0]);
    if(std::isinf(map_1)){ map_1 = 9E30; }
    if((std::isnan(map_1)) || (map_1 < 0.0)){ map_1 = 1E-10; }
    map_2 = 0.0;
    if(a[2] >= (2.0 * pi)){ a[2]-=(2.0 * pi); }
    if(a[2] <= (-2.0 * pi)){ a[2]+=(2.0 * pi); }
    map_3 = mid + (amp * (sin(a[2])));
    if(std::isnan(map_3)){ map_3 = 0.0; }
    if(std::isinf(map_3)){ map_3 = (map_3 > 0.0) ? 9E30 : -9E30; }
    if(a[3] >= (2.0 * pi)){ a[3]-=(2.0 * pi); }
    if(a[3] <= (-2.0 * pi)){ a[3]+=(2.0 * pi); }
    map_4 = mid + (amp * (sin(a[3])));
    if(std::isnan(map_4)){ map_4 = 0.0; }
    if(std::isinf(map_4)){ map_4 = (map_4 > 0.0) ? 9E30 : -9E30; }
    map_7 = exp(a[1]);
    if(std::isinf(map_7)){ map_7 = 9E30; }
    if((std::isnan(map_7)) || (map_7 < 0.0)){ map_7 = 1E-10; }
    erf_1 = 1.0 + (erff((map_7 * (x - map_3))));
    erf_2 = 1.0 + (erff((map_7 * (map_4 - x))));
    power_law = 1.0;
    *y = 0.25 * map_1 * erf_1 * erf_2 * power_law;
    dyda[0] = map_1 * 0.25 * erf_1 * erf_2 * power_law;
    dyda[1] = map_7 * 0.25 * map_1 * power_law * 2.0 * ((erf_2 * exp((-1.0 * map_7 * map_7 * (x - map_3) * (x - map_3))) * (x - map_3)) + (erf_1 * exp((-1.0 * map_7 * map_7 * (map_4 - x) * (map_4 - x))) * (map_4 - x))) / (sqrt(pi));
    dyda[2] = amp * (cos(a[2])) * 0.25 * map_1 * erf_2 * power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (x - map_3) * (x - map_3))) * -1.0 * map_7 / (sqrt(pi));
    dyda[3] = amp * (cos(a[3])) * 0.25 * map_1 * erf_1 * power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (map_4 - x) * (map_4 - x))) * map_7 / (sqrt(pi));
    for(i = 0; i < 4; i++){
      if(std::isinf(a[i])){ a[i] = (a[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(a[i])){ a[i] = 0.0; }
      if(std::isinf(dyda[i])){ dyda[i] = (dyda[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(dyda[i])){ dyda[i] = 0.0; }
    }
    a[4] = a[5] = a[6] = a[7] = 0.0;
    break;
  case 2:
    map_1 = exp(a[0]);
    if(std::isinf(map_1)){ map_1 = 9E30; }
    if((std::isnan(map_1)) || (map_1 < 0.0)){ map_1 = 1E-10; }
    map_2 = 0.0;
    if(a[2] >= (2.0 * pi)){ a[2]-=(2.0 * pi); }
    if(a[2] <= (-2.0 * pi)){ a[2]+=(2.0 * pi); }
    map_3 = mid + (amp * (sin(a[2])));
    if(std::isnan(map_3)){ map_3 = 0.0; }
    if(std::isinf(map_3)){ map_3 = (map_3 > 0.0) ? 9E30 : -9E30; }
    if(a[4] >= (2.0 * pi)){ a[4]-=(2.0 * pi); }
    if(a[4] <= (-2.0 * pi)){ a[4]+=(2.0 * pi); }
    map_4 = mid + (amp * (sin(a[4])));
    if(std::isnan(map_4)){ map_4 = 0.0; }
    if(std::isinf(map_4)){ map_4 = (map_4 > 0.0) ? 9E30 : -9E30; }
    map_7 = exp(a[1]);
    if(std::isinf(map_7)){ map_7 = 9E30; }
    if((std::isnan(map_7)) || (map_7 < 0.0)){ map_7 = 1E-10; }
    map_8 = exp(a[3]);
    if(std::isinf(map_8)){ map_8 = 9E30; }
    if((std::isnan(map_8)) || (map_8 < 0.0)){ map_8 = 1E-10; }
    erf_1 = 1.0 + (erff((map_7 * (x - map_3))));
    erf_2 = 1.0 + (erff((map_8 * (map_4 - x))));
    power_law = 1.0;
    *y = 0.25 * map_1 * erf_1 * erf_2 * power_law;
    dyda[0] = map_1 * 0.25 * erf_1 * erf_2 * power_law;
    dyda[1] = map_7 * 0.25 * map_1 * erf_2 * power_law * 2.0 * (exp((-1.0 * map_7 * map_7 * (x - map_3) * (x - map_3)))) * (x - map_3) / (sqrt(pi));
    dyda[2] = amp * (cos(a[2])) * 0.25 * map_1 * erf_2 * (power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (x - map_3) * (x - map_3))) * -1.0 * map_7 / (sqrt(pi)));
    dyda[3] = map_8 * 0.25 * map_1 * erf_1 * power_law * 2.0 * (exp((-1.0 * map_8 * map_8 * (map_4 - x) * (map_4 - x)))) * (map_4 - x) / (sqrt(pi));
    dyda[4] = amp * (cos(a[4])) * 0.25 * map_1 * erf_1 * (power_law * 2.0 * exp((-1.0 * map_8 * map_8 * (map_4 - x) * (map_4 - x))) * map_8 / (sqrt(pi)));
    for(i = 0; i < 5; i++){
      if(std::isinf(a[i])){ a[i] = (a[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(a[i])){ a[i] = 0.0; }
      if(std::isinf(dyda[i])){ dyda[i] = (dyda[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(dyda[i])){ dyda[i] = 0.0; }
    }
    a[5] = a[6] = a[7] = 0.0;
    break;
  case 3:
    map_1 = exp(a[0]);
    if(std::isinf(map_1)){ map_1 = 9E30; }
    if((std::isnan(map_1)) || (map_1 < 0.0)){ map_1 = 1E-10; }
    //map_2 = exp(a[4]);
    map_2 = a[4];
    if(std::isinf(map_2)){ map_2 = 9E30; }
    if((std::isnan(map_2)) || (map_2 < 0.0)){ map_2 = 1E-10; }    
    if(a[2] >= (2.0 * pi)){ a[2]-=(2.0 * pi); }
    if(a[2] <= (-2.0 * pi)){ a[2]+=(2.0 * pi); }
    map_3 = mid + (amp * (sin(a[2])));
    if(std::isnan(map_3)){ map_3 = 0.0; }
    if(std::isinf(map_3)){ map_3 = (map_3 > 0.0) ? 9E30 : -9E30; }
    if(a[3] >= (2.0 * pi)){ a[3]-=(2.0 * pi); }
    if(a[3] <= (-2.0 * pi)){ a[3]+=(2.0 * pi); }
    map_4 = mid + (amp * (sin(a[3])));
    if(std::isnan(map_4)){ map_4 = 0.0; }
    if(std::isinf(map_4)){ map_4 = (map_4 > 0.0) ? 9E30 : -9E30; }
    map_7 = exp(a[1]);
    if(std::isinf(map_7)){ map_7 = 9E30; }
    if((std::isnan(map_7)) || (map_7 < 0.0)){ map_7 = 1E-10; }
    erf_1 = 1.0 + (erff((map_7 * (x - map_3))));
    erf_2 = 1.0 + (erff((map_7 * (map_4 - x))));
    power_law = 1.0 + (map_2 * (pow((fabs(x - (0.5 * (map_3 + map_4)))),4.0)));
    if(std::isinf(power_law)){ power_law = (power_law > 0.0) ? 9E30 : -9E30; }
    if((power_law < 0.0) || (std::isnan(power_law))){ power_law = 0.0; }
    *y = 0.25 * map_1 * erf_1 * erf_2 * power_law;
    dyda[0] = map_1 * 0.25 * erf_1 * erf_2 * power_law;
    dyda[1] = map_7 * 0.25 * map_1 * power_law * 2.0 * ((erf_2 * exp((-1.0 * map_7 * map_7 * (x - map_3) * (x - map_3))) * (x - map_3)) + (erf_1 * exp((-1.0 * map_7 * map_7 * (map_4 - x) * (map_4 - x))) * (map_4 - x))) / (sqrt(pi));
    dyda[2] = amp * (cos(a[2])) * 0.25 * map_1 * erf_2 * ((power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (x - map_3) * (x - map_3))) * -1.0 * map_7 / (sqrt(pi))) + (erf_1 * -2.0 * map_2 * (pow((fabs(x - (0.5 * (map_3 + map_4)))),3.0)) * (x - (0.5 * (map_3 + map_4))) / (fabs((x - (0.5 * (map_3 + map_4)))))));
    dyda[3] = amp * (cos(a[3])) * 0.25 * map_1 * erf_1 * ((power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (map_4 - x) * (map_4 - x))) * map_7 / (sqrt(pi))) + (erf_2 * -2.0 * map_2 * (pow((fabs(x - (0.5 * (map_3 + map_4)))),3.0)) * (x - (0.5 * (map_3 + map_4))) / (fabs((x - (0.5 * (map_3 + map_4)))))));
    //dyda[4] = map_2 * 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - (0.5 * (map_3 + map_4)))),4.0));
    dyda[4] = 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - (0.5 * (map_3 + map_4)))),4.0));
    for(i = 0; i < 5; i++){
      if(std::isinf(a[i])){ a[i] = (a[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(a[i])){ a[i] = 0.0; }
      if(std::isinf(dyda[i])){ dyda[i] = (dyda[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(dyda[i])){ dyda[i] = 0.0; }
    }
    a[5] = a[6] = a[7] = 0.0;
    break;
  case 4:
    map_1 = exp(a[0]);
    if(std::isinf(map_1)){ map_1 = 9E30; }
    if((std::isnan(map_1)) || (map_1 < 0.0)){ map_1 = 1E-10; }
    //map_2 = exp(a[4]);
    map_2 = a[4];
    if(std::isinf(map_2)){ map_2 = 9E30; }
    if((std::isnan(map_2)) || (map_2 < 0.0)){ map_2 = 1E-10; }
    if(a[5] >= (2.0 * pi)){ a[5]-=(2.0 * pi); }
    if(a[5] <= (-2.0 * pi)){ a[5]+=(2.0 * pi); }
    //map_3 = 5.0 + (3.0 * (sin(a[5])));
    map_3 = 4.5 + (3.5 * (sin(a[5])));
    if(std::isnan(map_3)){ map_3 = 0.0; }
    if(std::isinf(map_3)){ map_3 = (map_3 > 0.0) ? 9E30 : -9E30; }
    if(a[2] >= (2.0 * pi)){ a[2]-=(2.0 * pi); }
    if(a[2] <= (-2.0 * pi)){ a[2]+=(2.0 * pi); }
    map_4 = mid + (amp * (sin(a[2])));
    if(std::isnan(map_4)){ map_4 = 0.0; }
    if(std::isinf(map_4)){ map_4 = (map_4 > 0.0) ? 9E30 : -9E30; }
    if(a[3] >= (2.0 * pi)){ a[3]-=(2.0 * pi); }
    if(a[3] <= (-2.0 * pi)){ a[3]+=(2.0 * pi); }
    map_5 = mid + (amp * (sin(a[3])));
    if(std::isnan(map_5)){ map_5 = 0.0; }
    if(std::isinf(map_5)){ map_5 = (map_5 > 0.0) ? 9E30 : -9E30; }
    map_7 = exp(a[1]);
    if(std::isinf(map_7)){ map_7 = 9E30; }
    if((std::isnan(map_7)) || (map_7 < 0.0)){ map_7 = 1E-10; }
    erf_1 = 1.0 + (erff((map_7 * (x - map_4))));
    erf_2 = 1.0 + (erff((map_7 * (map_5 - x))));
    power_law = 1.0 + (map_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),map_3)));
    if(std::isinf(power_law)){ power_law = (power_law > 0.0) ? 9E30 : -9E30; }
    if((std::isnan(power_law)) || (power_law < 0.0)){ power_law = 0.0; }
    *y = 0.25 * map_1 * erf_1 * erf_2 * power_law;
    dyda[0] = map_1 * 0.25 * erf_1 * erf_2 * power_law;
    dyda[1] = map_7 * 0.25 * map_1 * power_law * 2.0 * ((erf_2 * exp((-1.0 * map_7 * map_7 * (x - map_4) * (x - map_4))) * (x - map_4)) + (erf_1 * exp((-1.0 * map_7 * map_7 * (map_5 - x) * (map_5 - x))) * (map_5 - x))) / (sqrt(pi));
    dyda[2] = amp * (cos(a[2])) * 0.25 * map_1 * erf_2 * ((power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (x - map_4) * (x - map_4))) * -1.0 * map_7 / (sqrt(pi))) + (erf_1 * -0.5 * map_3 * map_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),(map_3 - 1.0))) * (x - (0.5 * (map_4 + map_5))) / (fabs((x - (0.5 * (map_4 + map_5)))))));
    dyda[3] = amp * (cos(a[3])) * 0.25 * map_1 * erf_1 * ((power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (map_5 - x) * (map_5 - x))) * map_7 / (sqrt(pi))) + (erf_2 * -0.5 * map_3 * map_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),(map_3 - 1.0))) * (x - (0.5 * (map_4 + map_5))) / (fabs((x - (0.5 * (map_4 + map_5)))))));
    //dyda[4] = map_2 * 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),map_3));
    dyda[4] = 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),map_3));
    //dyda[5] = 4.0 * (cos(a[5])) * 0.25 * map_1 * erf_1 * erf_2 * map_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),map_3)) * (log((fabs((x - (0.5 * (map_4 + map_5)))))));
    dyda[5] = 3.5 * (cos(a[5])) * 0.25 * map_1 * erf_1 * erf_2 * map_2 * (pow((fabs(x - (0.5 * (map_4 + map_5)))),map_3)) * (log((fabs((x - (0.5 * (map_4 + map_5)))))));
    for(i = 0; i < 6; i++){
      if(std::isinf(a[i])){ a[i] = (a[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(a[i])){ a[i] = 0.0; }
      if(std::isinf(dyda[i])){ dyda[i] = (dyda[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(dyda[i])){ dyda[i] = 0.0; }
    }
    a[6] = a[7] = 0.0;
    break;
  case 5:
    map_1 = exp(a[0]);
    if(std::isinf(map_1)){ map_1 = 9E30; }
    if((std::isnan(map_1)) || (map_1 < 0.0)){ map_1 = 1E-10; }
    //map_2 = exp(a[4]);
    map_2 = a[4];
    if(std::isinf(map_2)){ map_2 = 9E30; }
    if((std::isnan(map_2)) || (map_2 < 0.0)){ map_2 = 1E-10; }
    if(a[6] >= (2.0 * pi)){ a[6]-=(2.0 * pi); }
    if(a[6] <= (-2.0 * pi)){ a[6]+=(2.0 * pi); }
    //map_3 = 5.0 + (3.0 * (sin(a[6])));
    map_3 = 4.5 + (3.5 * (sin(a[6])));
    if(std::isnan(map_3)){ map_3 = 0.0; }
    if(std::isinf(map_3)){ map_3 = (map_3 > 0.0) ? 9E30 : -9E30; }
    if(a[2] >= (2.0 * pi)){ a[2]-=(2.0 * pi); }
    if(a[2] <= (-2.0 * pi)){ a[2]+=(2.0 * pi); }
    map_4 = mid + (amp * (sin(a[2])));
    if(std::isnan(map_4)){ map_4 = 0.0; }
    if(std::isinf(map_4)){ map_4 = (map_4 > 0.0) ? 9E30 : -9E30; }
    if(a[3] >= (2.0 * pi)){ a[3]-=(2.0 * pi); }
    if(a[3] <= (-2.0 * pi)){ a[3]+=(2.0 * pi); }
    map_5 = mid + (amp * (sin(a[3])));
    if(std::isnan(map_5)){ map_5 = 0.0; }
    if(std::isinf(map_5)){ map_5 = (map_5 > 0.0) ? 9E30 : -9E30; }
    if(a[5] >= (2.0 * pi)){ a[5]-=(2.0 * pi); }
    if(a[5] <= (-2.0 * pi)){ a[5]+=(2.0 * pi); }
    map_6 = mid + (amp * (sin(a[5])));
    if(std::isnan(map_6)){ map_6 = 0.0; }
    if(std::isinf(map_6)){ map_6 = (map_6 > 0.0) ? 9E30 : -9E30; }
    map_7 = exp(a[1]);
    if(std::isinf(map_7)){ map_7 = 9E30; }
    if((std::isnan(map_7)) || (map_7 < 0.0)){ map_7 = 1E-10; }
    erf_1 = 1.0 + (erff((map_7 * (x - map_4))));
    erf_2 = 1.0 + (erff((map_7 * (map_5 - x))));
    power_law = 1.0 + (map_2 * (pow((fabs(x - map_6)),map_3)));
    if(std::isinf(power_law)){ power_law = (power_law > 0.0) ? 9E30 : -9E30; }
    if((std::isnan(power_law)) || (power_law < 0.0)){ power_law = 0.0; }
    *y = 0.25 * map_1 * erf_1 * erf_2 * power_law;
    dyda[0] = map_1 * 0.25 * erf_1 * erf_2 * power_law;
    dyda[1] = map_7 * 0.25 * map_1 * power_law * 2.0 * ((erf_2 * exp((-1.0 * map_7 * map_7 * (x - map_4) * (x - map_4))) * (x - map_4)) + (erf_1 * exp((-1.0 * map_7 * map_7 * (map_5 - x) * (map_5 - x))) * (map_5 - x))) / (sqrt(pi));
    dyda[2] = amp * (cos(a[2])) * 0.25 * map_1 * erf_2 * power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (x - map_4) * (x - map_4))) * -1.0 * map_7 / (sqrt(pi));
    dyda[3] = amp * (cos(a[3])) * 0.25 * map_1 * erf_1 * power_law * 2.0 * exp((-1.0 * map_7 * map_7 * (map_5 - x) * (map_5 - x))) * map_7 / (sqrt(pi));
    //dyda[4] = map_2 * 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - map_6)),map_3));
    dyda[4] = 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - map_6)),map_3));
    dyda[5] = amp * (cos(a[5])) * 0.25 * map_1 * erf_1 * erf_2 * map_3 * map_2 * (pow((fabs(x - map_6)),(map_3 - 1.0))) * -1.0 * (x - map_6) / (fabs((x - map_6)));
    //dyda[6] = 4.0 * (cos(a[6])) * 0.25 * map_1 * erf_1 * erf_2 * map_2 * (pow((fabs(x - map_6)),map_3)) * (log((fabs((x - map_6)))));
    dyda[6] = 3.5 * (cos(a[6])) * 0.25 * map_1 * erf_1 * erf_2 * map_2 * (pow((fabs(x - map_6)),map_3)) * (log((fabs((x - map_6)))));
    for(i = 0; i < 7; i++){
      if(std::isinf(a[i])){ a[i] = (a[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(a[i])){ a[i] = 0.0; }
      if(std::isinf(dyda[i])){ dyda[i] = (dyda[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(dyda[i])){ dyda[i] = 0.0; }
    }
    a[7] = 0.0;
    break;
  default:
    map_1 = exp(a[0]);
    if(std::isinf(map_1)){ map_1 = 9E30; }
    if((std::isnan(map_1)) || (map_1 < 0.0)){ map_1 = 1E-10; }
    //map_2 = exp(a[5]);
    map_2 = a[5];
    if(std::isinf(map_2)){ map_2 = 9E30; }
    if((std::isnan(map_2)) || (map_2 < 0.0)){ map_2 = 1E-10; }
    if(a[7] >= (2.0 * pi)){ a[7]-=(2.0 * pi); }
    if(a[7] <= (-2.0 * pi)){ a[7]+=(2.0 * pi); }
    //map_3 = 5.0 + (3.0 * (sin(a[7])));
    map_3 = 4.5 + (3.5 * (sin(a[7])));
    if(std::isnan(map_3)){ map_3 = 0.0; }
    if(std::isinf(map_3)){ map_3 = (map_3 > 0.0) ? 9E30 : -9E30; }
    if(a[2] >= (2.0 * pi)){ a[2]-=(2.0 * pi); }
    if(a[2] <= (-2.0 * pi)){ a[2]+=(2.0 * pi); }
    map_4 = mid + (amp * (sin(a[2])));
    if(std::isnan(map_4)){ map_4 = 0.0; }
    if(std::isinf(map_4)){ map_4 = (map_4 > 0.0) ? 9E30 : -9E30; }
    if(a[4] >= (2.0 * pi)){ a[4]-=(2.0 * pi); }
    if(a[4] <= (-2.0 * pi)){ a[4]+=(2.0 * pi); }
    map_5 = mid + (amp * (sinf(a[4])));
    if(std::isnan(map_5)){ map_5 = 0.0; }
    if(std::isinf(map_5)){ map_5 = (map_5 > 0.0) ? 9E30 : -9E30; }
    if(a[6] >= (2.0 * pi)){ a[6]-=(2.0 * pi); }
    if(a[6] <= (-2.0 * pi)){ a[6]+=(2.0 * pi); }
    map_6 = mid + (amp * (sin(a[6])));
    if(std::isnan(map_6)){ map_6 = 0.0; }
    if(std::isinf(map_6)){ map_6 = (map_6 > 0.0) ? 9E30 : -9E30; }
    map_7 = exp(a[1]);
    if(std::isinf(map_7)){ map_7 = 9E30; }
    if((std::isnan(map_7)) || (map_7 < 0.0)){ map_7 = 1E-10; }
    map_8 = exp(a[3]);
    if(std::isinf(map_8)){ map_8 = 9E30; }
    if((std::isnan(map_8)) || (map_8 < 0.0)){ map_8 = 1E-10; }
    erf_1 = 1.0 + (erff((map_7 * (x - map_4))));
    erf_2 = 1.0 + (erff((map_8 * (map_5 - x))));
    power_law = 1.0 + (map_2 * (pow((fabs(x - map_6)),map_3)));
    if(std::isinf(power_law)){ power_law = (power_law > 0.0) ? 9E30 : -9E30; }
    if((std::isnan(power_law)) || (power_law < 0.0)){ power_law = 0.0; }
    *y = 0.25 * map_1 * erf_1 * erf_2 * power_law;
    dyda[0] = map_1 * 0.25 * erf_1 * erf_2 * power_law;
    dyda[1] = map_7 * 0.25 * map_1 * erf_2 * power_law * 2.0 * (exp((-1.0 * map_7 * map_7 * (x - map_4) * (x - map_4)))) * (x - map_4) / (sqrt(pi));
    dyda[2] = amp * (cos(a[2])) * 0.25 * map_1 * erf_2 * power_law * 2.0 * (exp((-1.0 * map_7 * map_7 * (x - map_4) * (x - map_4)))) * -1.0 * map_7 / (sqrt(pi));
    dyda[3] = map_8 * 0.25 * map_1 * erf_1 * power_law * 2.0 * (exp((-1.0 * map_8 * map_8 * (map_5 - x) * (map_5 - x)))) * (map_5 - x) / (sqrt(pi));
    dyda[4] = amp * (cos(a[4])) * 0.25 * map_1 * erf_1 * power_law * 2.0 * (exp((-1.0 * map_8 * map_8 * (map_5 - x) * (map_5 - x)))) * map_8 / (sqrt(pi));
    //dyda[5] = map_2 * 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - map_6)),map_3));
    dyda[5] = 0.25 * map_1 * erf_1 * erf_2 * (pow((fabs(x - map_6)),map_3));
    dyda[6] = amp * (cos(a[6])) * 0.25 * map_1 * erf_1 * erf_2 * map_3 * map_2 * (pow((fabs(x - map_6)),(map_3 - 1.0))) * -1.0 * (x - map_6) / (fabs((x - map_6)));
    //dyda[7] = 4.0 * (cos(a[7])) * 0.25 * map_1 * erf_1 * erf_2 * map_2 * (log(fabs((x - map_6)))) * (pow((fabs(x - map_6)),map_3));
    dyda[7] = 3.5 * (cos(a[7])) * 0.25 * map_1 * erf_1 * erf_2 * map_2 * (log(fabs((x - map_6)))) * (pow((fabs(x - map_6)),map_3));
    for(i = 0; i < 8; i++){
      if(std::isinf(a[i])){ a[i] = (a[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(a[i])){ a[i] = 0.0; }
      if(std::isinf(dyda[i])){ dyda[i] = (dyda[i] > 0.0) ? 9E30 : -9E30; } 
      if(std::isnan(dyda[i])){ dyda[i] = 0.0; }
    }
    break;
  }

  // test range of values
  if(std::isinf(*y)){ *y = (*y > 0.0) ? 9E30 : -9E30; }
  if(std::isnan(*y)){ *y = 0.0; }

}

void covsrt(double **covar, int ma, int ia[], int mfit){

  int i,j,k;

  for (i=mfit;i<ma;i++)
    for (j=0;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit - 1;
  for (j=ma-1;j>=0;j--) {
    if (ia[j]) {
      for (i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j]);
      for (i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i]);
      k--;
    }
  }

}

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]){
  
  int jj,j,i;
  double s,*tmp;
  
  tmp = new double[n];
  for (j=0;j<n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=0;i<m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<n;j++) {
    s=0.0;
    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  delete [] tmp;
  tmp = NULL;

}

void svdcmp(double **a, int m, int n, double w[], double **v){
  
  double pythag(double a, double b);
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1 = new double[n];
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for(k=i;k<m;k++){ scale+=(fabs(a[k][i])); }
      if (scale != 0.0) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s+=(a[k][i]*a[k][i]);
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l-1;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++){ s+=(a[k][i]*a[k][j]); }
	  f=s/h;
	  for(k=i;k<m;k++){ a[k][j]+=(f*a[k][i]); }
	}
	for(k=i;k<m;k++){ a[k][i]*=scale; }
      }
    }
    w[i]=scale * g;
    g=s=scale=0.0;
    if((i < m) && (i != (n - 1))){
      for(k=l-1;k<n;k++){ scale+=(fabs(a[i][k])); }
      if(scale != 0.0){
	for (k=l-1;k<n;k++) {
	  a[i][k] /= scale;
	  s+=(a[i][k]*a[i][k]);
	}
	f=a[i][(l-1)];
	g = -SIGN(sqrt(s),f);
	h=(f*g)-s;
	a[i][(l - 1)]=f-g;
	for(k=l-1;k<n;k++){ rv1[k]=a[i][k]/h; }
	for(j=l-1;j<m;j++){
	  for (s=0.0,k=l-1;k<n;k++){ s+=(a[j][k]*a[i][k]); }
	  for (k=l-1;k<n;k++){ a[j][k]+=(s*rv1[k]); }
	}
	for (k=l-1;k<n;k++){ a[i][k]*=scale; }
      }
    }
    anorm=MAX(anorm,((fabs(w[i])+fabs(rv1[i]))));
  }
  for (i=n-1;i>=0;i--) {
    if (i < (n - 1)) {
      if (g != 0.0) {
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++){ s+=(a[i][k]*v[k][j]); }
	  for (k=l;k<n;k++){ v[k][j]+=(s*v[k][i]); }
	}
      }
      for (j=l;j<n;j++){ v[i][j]=v[j][i]=0.0; }
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for(i=(MIN(m,n) - 1);i>=0;i--){
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++){ a[i][j]=0.0; }
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++){ s+=(a[k][i]*a[k][j]); }
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++){ a[k][j]+=(f*a[k][i]); }
      }
      for (j=i;j<m;j++){ a[j][i]*=g; }
    } else for (j=i;j<m;j++){ a[j][i]=0.0; }
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<100;its++) {
      flag=1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	if ((l == 0) || ((double)(fabs(rv1[l])+anorm) == anorm)) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=(y*c)+(z*s);
	    a[j][i]=(z*c)-(y*s);
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if(its == 200){ 

	std::cout << "No convergence in 200 svdcmp iterations. Exiting." << std::endl; 
	delete [] rv1;
	rv1 = NULL; 
	return; 

      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  delete [] rv1;
  rv1 = NULL;

}
 
double ran2(long &idum){

  #define IM1 2147483563
  #define IM2 2147483399
  #define AM (1.0/IM1)
  #define IMM1 (IM1-1)
  #define IA1 40014
  #define IA2 40692
  #define IQ1 53668
  #define IQ2 52774
  #define IR1 12211
  #define IR2 3791
  #define NTAB 32
  #define NDIV (1+IMM1/NTAB)
  #define EPS 1.2e-307
  #define RNMX (1.0-EPS)

  long j,k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (idum <= 0) {
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    idum2=(idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(idum)/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=(idum)/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;

  #undef IM1 
  #undef IM2 
  #undef AM 
  #undef IMM1 
  #undef IA1 
  #undef IA2 
  #undef IQ1 
  #undef IQ2 
  #undef IR1 
  #undef IR2 
  #undef NTAB 
  #undef NDIV 
  #undef EPS 
  #undef RNMX 

}

double pythag(double a, double b){

  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if(absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));

}

// functions used to generate an arbitrary number of correlated Busy Function fit parameters,
// and derive uncertainties in observational properties

void choldc(double ** a, int n){

  int i,j,k;
  double sum;

  for(i=0;i<n;i++){
    if(abs(a[i][i]) <= 1.0E-12){ continue; }
    for(j=i;j<n;j++){
      for(sum=a[i][j],k=i-1;k>=0;k--){ sum-=(a[i][k]*a[j][k]); }
      if(i == j){
	if(sum <= 0.0){
	  a[i][i] = 0.0;
	} else {
	  a[i][i] = sqrt(sum);
	} 
      } else {
	a[j][i]=sum/a[i][i];
      }
    }
    for(j=i-1;j>=0;j--){ a[j][i] = 0.0; }
  }

  // range test resultant matrix
  for(j = 0; j < n; j++){
    for(i = 0; i < n; i++){
      if(std::isnan(a[j][i])){ a[j][i] = 0.0; }
      if(std::isinf(a[j][i])){ a[j][i] = (a[j][i] > 0.0) ? 9E30 : -9E30; } 
    }
  }

}

double inverf(double p){

  int j;
  double x,err,t,pp;

  p = 1.0 - p;
  if(p >= 2.0){ return -100.0; }
  if(p <= 0.0){ return 100.0; }
  pp = (p < 1.0) ? p : (2.0 - p);
  t = sqrt(-2.0*log(pp/2.0));
  x = -0.70711*((2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481)) - t);
  for(j=0;j<2;j++){
    err=erfc(x) - pp;
    x+=(err/(1.12837916709551257*exp(-SQR(x))-x*err));
  }
  return  (p < 1.0 ? x : (-1.0*x));

}





