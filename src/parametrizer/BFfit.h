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

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>

#ifndef _BFfit_
#define _BFfit_

// templated macro-like inline functions used by NR-based functions --- called by the other functions in various places

template<class T>
inline T SQR(const T a) {return ((a == 0) ? 0.0 : (a*a));}

template<class T>
inline const T &MAX(const T &a, const T &b)
{return ((b > a) ? (b) : (a));}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return ((b < a) ? (b) : (a));}

template<class T>
inline T SIGN(const T &a, const T &b)
{return ((b >= 0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

// Chain of functions that carry out Busy Function fitting --- templated at the highest level

template <class T_count, class T_xvals, class T_data>
void mrqcof(T_xvals x[], T_data y[], T_data sig[], T_count ndata, double a[], int ia[], int ma, double **alpha, double beta[], double *chisq, void (*funcs) (double, double [], double *, double [], int, double, double), int fit_mode, double mid, double amp){

  int i,j,k,l,m,mfit=0;
  double ymod,wt,sig2i,dy,*dyda;

  dyda = new double[ma];
  for(j=0;j<ma;j++){ if(ia[j]){ mfit++; } }
  for(j=0;j<mfit;j++){
    for(k=0;k<=j;k++){ alpha[j][k]=0.0; }
    beta[j]=0.0;
  }
  *chisq=0.0;
  for(i=0;i<ndata;i++){
    (*funcs)(((double) x[i]),a,&ymod,dyda,fit_mode,mid,amp);
    sig2i=1.0/((double)(sig[i]*sig[i]));
    dy=((double) y[i])-ymod;
    for (j=0,l=0;l<ma;l++) {
      if (ia[l]) {
	wt=dyda[l]*sig2i;
	for (k=0,m=0;m<=l;m++){
	  if(ia[m]){ alpha[j][k++]+=(wt*dyda[m]); }
	}
	beta[j++]+=(dy*wt);
      }
    }
    *chisq+=(dy*dy*sig2i);
  }
  if((std::isinf(*chisq)) || (std::isnan(*chisq))){ *chisq = 9E30; }
  for (j=1;j<mfit;j++){
    for (k=0;k<j;k++){ alpha[k][j]=alpha[j][k]; }
  }
 
  // test for inf and nan values in alpha, beta and dyda
  for(j = 0; j < mfit; j++){

    if(std::isinf(beta[j])){ beta[j] = (beta[j] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(beta[j])){ beta[j] = 0.0; }

    for(k = 0; k < mfit; k++){

      if(std::isinf(alpha[j][k])){ alpha[j][k] = (alpha[j][k] > 0.0) ? 9E30 : -9E30; }
      if(std::isnan(alpha[j][k])){ alpha[j][k] = 0.0; }

    }

  }

  delete [] dyda;
  dyda = NULL;
  
}

template <class T_count, class T_xvals, class T_data>
int mrqmin(T_xvals x[], T_data y[], T_data sig[], T_count ndata, double a[], int ia[], int ma, double **covar, double **alpha, double *chisq, void (*funcs)(double, double [], double *, double [], int, double, double), double *alamda, int fit_mode, double mid, double amp, int &mfit, double &ochisq, double * atry, double * beta, double * da, double ** oneda){

  void covsrt(double **covar, int ma, int ia[], int mfit);
  void svdcmp(double **a, int m, int n, double w[], double **v);
  void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
  int i,j,k,l,state_flag;
  double ** temp_covar, * w, ** v, * soln, * b, wmax, wmin;

  // set the value of state_flag 
  state_flag = 0;
  if((*alamda > -2.0) && (*alamda < 0.0)){ 

    *alamda = 0.0; 
    state_flag = 1; 

  } else if(*alamda < -2.0){
    
    state_flag = -1;

  }

  // test range of *alamda to ensure it's not ludicrously small or large during mid-stages
  if(state_flag >= 0){

    if(std::isnan(*alamda)){ *alamda = 1E-10; }
    if(std::isinf(*alamda)){ *alamda = (*alamda > 0.0) ? 1E20 : 1E-10; }
    if(*alamda < 1E-10){ *alamda = 1E-10; }
    if(*alamda > 1E20){ *alamda = 1E20; }
    
  }

  temp_covar = new double * [ma];
  for(j = 0; j < ma; j++){ temp_covar[j] = new double[ma]; }
  w = new double[ma];
  v = new double * [ma];
  for(j = 0; j < ma; j++){ v[j] = new double[ma]; }
  soln = new double[ma];
  b = new double[ma];

  if(state_flag < 0){
    for (mfit=0,j=0;j<ma;j++){ if(ia[j] > 0){ mfit++; } }
    *alamda=0.001;
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs,fit_mode,mid,amp);
    ochisq=(*chisq);
    for(j=0;j<ma;j++){ atry[j]=a[j]; }
  }
  for(j=0;j<mfit;j++){
    for(k=0;k<mfit;k++){ covar[j][k]=alpha[j][k]; }
    covar[j][j]=alpha[j][j]*(1.0+(*alamda));
    oneda[j][0]=beta[j];
  }

  // calculate inverse of covar[] using SVD 
  for(j=0;j<mfit;j++){
    for(k=0;k<mfit;k++){
      temp_covar[j][k]=covar[j][k];
    }
  }

  svdcmp(temp_covar,mfit,mfit,w,v);

  for(j=0;j<mfit;j++){
    if(std::isnan(w[j])){ w[j] = 0.0; }
    if(std::isinf(w[j])){ w[j] = (w[j] > 0.0) ? 9E30 : -9E30; }
    for(k=0;k<mfit;k++){
      if(std::isnan(temp_covar[j][k])){ temp_covar[j][k] = 0.0; }
      if(std::isinf(temp_covar[j][k])){ temp_covar[j][k] = (temp_covar[j][k] > 0.0) ? 9E30 : -9E30; }
      if(std::isnan(v[j][k])){ v[j][k] = 0.0; }
      if(std::isinf(v[j][k])){ v[j][k] = (v[j][k] > 0.0) ? 9E30 : -9E30; }
    }
  }

  // threshold useless weights --- uncomment the for loop to apply thresholding 
  wmax=0.0;
  //for(j=0;j<mfit;j++){ if(w[j] > wmax){ wmax = w[j]; } }
  //wmin=wmax * 1E-6;
  //wmin=wmax * 1E-9;
  wmin=wmax * 1E-12;
  for(j=0;j<mfit;j++){ if(w[j] < wmin){ w[j] = 0.0; } }

  for(j=0;j<mfit;j++){
    for(k=0;k<mfit;k++){ b[k] = soln[k] = 0.0; }
    b[j]=1.0;
    svbksb(temp_covar,w,v,mfit,mfit,b,soln);
    for(k=0;k<mfit;k++){ covar[k][j] = soln[k]; }

  }
  
  // ensure oneda[j][1] contains the solution, x, to A * x = b
  // using svbksb 
  for(j=0;j<mfit;j++){ b[j]=oneda[j][0]; soln[j] = 0.0; }
  svbksb(temp_covar,w,v,mfit,mfit,b,soln);
  for(j=0;j<mfit;j++){ da[j]=soln[j]; }

  if(state_flag > 0) {
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
    for(i = 0; i < ma; i++){ delete [] temp_covar[i]; temp_covar[i] = NULL; }
    delete [] temp_covar;
    temp_covar = NULL;
    delete [] w;
    w = NULL;
    for(i = 0; i < ma; i++){ delete [] v[i]; v[i] = NULL; }
    delete [] v;
    v = NULL;
    delete [] soln;
    soln = NULL;
    delete [] b;
    b = NULL;
    return 1;
  }

  for(j=-1,l=0;l<ma;l++){ if(ia[l]){ atry[l]=a[l]+da[++j]; } }

  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs,fit_mode,mid,amp);

  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq=(*chisq);
    for (j=0;j<mfit;j++) {
      for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
      beta[j]=da[j];
    }
    for (l=0;l<ma;l++) a[l]=atry[l];
  } else {
    *alamda *= 10.0;
    *chisq=ochisq;
  }

  for(i = 0; i < ma; i++){ delete [] temp_covar[i]; temp_covar[i] = NULL; }
  delete [] temp_covar;
  temp_covar = NULL;
  delete [] w;
  w = NULL;
  for(i = 0; i < ma; i++){ delete [] v[i]; v[i] = NULL; }
  delete [] v;
  v = NULL;
  delete [] soln;
  soln = NULL;
  delete [] b;
  b = NULL;
  
  return 1;

}

void BusyFunc(double x, double a[], double *y, double dyda[], int fit_mode, double mid, double amp);

template <class T_count, class T_xvals, class T_data>
  void FitBusyFunc_engine(T_count NOvals, T_xvals * x_vals, T_data * y_vals, T_data * n_vals, double * model_params, int fit_mode, int NOs, double ** start_vals, double mid, double amp, double ** fit_covar, int iter_max, int vb_flag){

  void BusyFunc(double x, double a[], double *y, double dyda[], int fit_mode, double mid, double amp);
  double ran2(long &idum);
  double chi2_best, chi2_sum, prev_chi2_val, chi2_val, ** best_params, *** best_fit_covar, * model_vals, model_max, y_max;
  // c arrays keep track of the following
  // c1 -- global scaling ---> re-mapped to exp(global scaling)
  // c2 -- first erff slope ---> re-mapped to exp(first erff slope)
  // c3 -- first erff position ---> re-mapped to mid + amp * sin()
  // c4 -- second erff slope ---> re-mapped to exp(second erff slope)
  // c5 -- second erff position ---> re-mapped to mid + amp * sin()
  // c6 -- power-law scaling ---> re-mapped to exp(power-law scaling) prior to v1.4, no longer re-mapped
  // c7 -- power-law position ---> re-mapped to mid + amp * sin()
  // c8 -- power-law exponent ---> re-mapped to 5 + 3 * sin()
  int i, j, c[8], NOc[8], SVD_flag, btr_cnt, thrd, NOthrd;
  double c_val[8], c_min[8], c_max[8], c_step[8];    
  // variables used to do LVM optimisation
  double * alamda, ** a;
  long int seed;
  int ma = 8, iter, ** ia, s, s_done;
  double *** covar, *** alpha;
  void (*funcs)(double x, double a[], double *y, double dyda[], int fit_mode, double mid, double amp);
  funcs = &BusyFunc;
  // variables used to show start seed progress
  double progress;
  // variables used to implement MRQMIN using multi-threading
  int mfit;
  double ochisq, ** atry, ** beta, ** da, *** oneda;

  // 0. test the range of model_params
  for(i = 0; i < 17; i++){

    if(std::isinf(model_params[i])){ model_params[i] = (model_params[i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[i])){ model_params[i] = 0.1; }

  }

  // 1. set the number of values to be used for each variable
  NOc[0] = 8; // over-all scaling
  NOc[1] = NOc[3] = 13; // roll-off 
  NOc[2] = NOc[4] = 9; // positions of roll-offs
  NOc[5] = 14; // power-law scaling --- minimum is 11, starts at 1E-8 --> powers to 1E3 is 13, maximal is 17 @ 1E7
  NOc[6] = 7; // position of power-law
  NOc[7] = 3; // power-law exponent

  // 2. set values according to input data
  c_min[0] = c_step[0] = (1.0 / ((double) NOc[0])); c_max[0] = 1.0;    
  c_min[1] = c_min[3] = 0.1; c_max[1] = c_max[3] = 10.0; c_step[1] = c_step[3] = 0.5; // this information isn't used in the current implementation --- retained for completeness
  c_min[2] = c_min[4] = c_min[6] = ((mid - amp) >= (double) x_vals[0]) ? (mid - amp) : (double) x_vals[0]; 
  c_max[2] = c_max[4] = c_max[6] = ((mid + amp) <= (double) x_vals[(NOvals - 1)]) ? (mid + amp) : (double) x_vals[(NOvals - 1)]; 
  c_step[2] = c_step[4] = c_step[6] = (c_max[2] - c_min[2]) / ((double) NOc[2] - 1.0);
  c_min[5] = -9.0; c_max[5] = 7.0; c_step[5] = 1.0; 
  c_min[7] = 2.0; c_max[7] = 8.0; c_step[7] = 2.0;  

  // 3. If the fit_mode is > 0, use brute force chi^2 sampling to find a new starting point.
  // Otherwise, use the previous starting point in start_vals[].
  if(fit_mode > 0){
    
    // 3.0 Create array model_vals[] used to generate and store curves
    model_vals = new double[NOvals];

    // 3.1. Find the maximum value of the input data
    y_max = -9E20;
    for(i = 0; i < NOvals; i++){ if((double) y_vals[i] >= y_max){ y_max = (double) y_vals[i]; } }

    // 3.2 seed the random number generator
    srand (time(NULL));
    seed = (long) (-1*(rand()));

    // 3.3 randomly generate NOs LVM seed values
    if(vb_flag > 0){

      std::cout << "Generating start seeds . . . " << std::endl;
      std::cout << "0         25        50        75        100%" << std::endl;
      std::cout << "| | | | | : | | | | : | | | | : | | | | :" << std::endl;

    }
    progress = 0.0;
    for(s = 0; s < NOs; s++){

      // randomly pick all the c_val[], then construct a curve and re-scale
      // c_val[0] accordingly
      for(i = 0; i < 8; i++){ 

	c[i] = (int) floor((ran2(seed) * ((double) NOc[i]))); 
	c_val[i] = c_min[i] + (((double) c[i]) * c_step[i]);
	if(std::isinf(c_val[i])){ c_val[i] = (c_val[i] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(c_val[i])){ c_val[i] = 0.0; }
	
      }
      c[4] = c[2] + (int) floor((ran2(seed) * ((double) (NOc[4] - c[2])))); 
      c_val[4] = c_min[4] + (((double) c[4]) * c_step[4]);
      if(c[5] > 0){
	
	c_val[5] = pow(10.0,c_val[5]);
	c_val[5] = (ran2(seed) < 0.5) ? c_val[5] : -c_val[5];
      
      } else { 

	c_val[5] = 0.0; 

      }
      if(c[2] == c[4]){ c[6] = 0; }
      c_val[6] = c_val[2] + ((double) c[6] * ((c_val[4] - c_val[2]) / (double)(NOc[6] - 1)));
      for(i = 1; i < 4; i+=2){

	switch (c[i]) {
	case 0:
	  c_val[i] = 0.1;
	  break;
	case 1:
	  c_val[i] = 0.2;
	  break;
	case 2:
	  c_val[i] = 1.6;
	  break;
	case 3:
	  c_val[i] = 6.4;
	  break;
	case 4:
	  c_val[i] = 0.4;
	  break;
	case 5:
	  c_val[i] = 0.8;
	  break;
	case 6:
	  c_val[i] = 3.2;
	  break;
	case 7:
	  c_val[i] = 0.747806 / 1.25;
	  break;
	case 8:
	  c_val[i] = 0.747806 / 2.5;
	  break;
	case 9:
	  c_val[i] = 0.747806 / 5.0;
	  break;
	case 10:
	  c_val[i] = 0.747806 / 10.0;
	  break;
	case 11:
	  c_val[i] = 0.747806 / 15.0;
	  break;
	case 12:
	  c_val[i] = 0.747806 / 20.0;
	default:
	  c_val[i] = 10.0;
	  break;
	}

      }

      // test range of c_val[]
      for(i = 0; i < 8; i++){ 

	if(std::isinf(c_val[i])){ c_val[i] = (c_val[i] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(c_val[i])){ c_val[i] = 0.0; }
	
      }      

      // calculate curve
      model_max = -9E20;
      for(i = 0; i < NOvals; i++){
	
	model_vals[i] = 0.25 * 1E-3 * (1.0 + erff((c_val[1] * ((double) x_vals[i] - c_val[2])))) * (1.0 + erff((c_val[3] * (c_val[4] - (double) x_vals[i])))) * (1.0 + (c_val[5] * (pow((fabs((double) x_vals[i] - c_val[6])),c_val[7]))));
	if(std::isnan(model_vals[i])){ model_vals[i] = 0.0; }
	if(std::isinf(model_vals[i])){ model_vals[i] = (model_vals[i] > 0.0) ? 9E30 : -9E30; } 
	if(model_vals[i] >= model_max){ model_max = model_vals[i]; }
	
	// for(i = 0; i < NOvals; i++)
      }
		      
      // calculate chi^2
      chi2_sum = 0.0;
      for(i = 0; i < NOvals; i++){
	
	chi2_val = ((double) y_vals[i] - (c_val[0] * y_max * model_vals[i] / model_max)) / (double) n_vals[i];
	chi2_val = chi2_val * chi2_val;
	if((!(std::isnan(chi2_val))) && (!(std::isinf(chi2_val)))){
	  
	  chi2_sum+=chi2_val;
	  
	} else if(std::isinf(chi2_val)){ chi2_sum+=1E20; }
	
	// for(i = 0; i < NOvals; i++)
      }
     
      if(std::isinf(chi2_sum)){ chi2_sum = (chi2_sum > 0.0) ? 9E30 : -9E30; }

      // re-scale c_val[0] and c_val[5]
      c_val[0]*=(y_max * 1E-3 / model_max);
      c_val[5]*=(y_max * 1E-3 / model_max);

      // update start_vals
      for(i = 0; i < 8; i++){ start_vals[s][i] = c_val[i]; }
      start_vals[s][8] = chi2_sum;

      // display progress
      if(vb_flag > 0){ while((((float) (s + 1)) / ((float) NOs)) >= progress){ std::cout << "*" << std::flush; progress+=0.025; } }

      // for(s = 0; s < NOs; s++)
    }
    if(vb_flag > 0){ std::cout << " done." << std::endl; }
  
    // 3.4 free up memory
    delete [] model_vals;
    model_vals = NULL;
    
    // if(fit_mode > 0)
  } else {

    if(vb_flag > 0){ std::cout << "Using existing LVM starting seeds . . . " << std::endl; }
    
  }
 
  // 4. Try to optimise each starting position using the Levenberg-Marquardt algorithm 
  // --- implemented using SVD --- and retain the global chi^2 minium
  NOthrd = 1;
  ia = new int * [NOthrd];
  alamda = new double[NOthrd];
  a = new double * [NOthrd];
  covar = new double ** [NOthrd];
  alpha = new double ** [NOthrd];
  best_params = new double * [NOthrd];
  best_fit_covar = new double ** [NOthrd];
  atry = new double * [NOthrd];
  beta = new double * [NOthrd];
  da = new double * [NOthrd];
  oneda = new double ** [NOthrd];
  for(thrd = 0; thrd < NOthrd; thrd++){
    a[thrd] = new double[8];
    covar[thrd] = new double * [8];
    alpha[thrd] = new double * [8];
    ia[thrd] = new int[8];
    best_params[thrd] = new double[17];
    best_fit_covar[thrd] = new double * [8];
    atry[thrd] = new double[8];
    beta[thrd] = new double[8];
    da[thrd] = new double[8];
    oneda[thrd] = new double * [8];
    for(i = 0; i < 8; i++){
      
      covar[thrd][i] = new double[8];
      alpha[thrd][i] = new double[8];
      best_fit_covar[thrd][i] = new double[8];
      oneda[thrd][i] = new double;

    }

  }
  for(thrd = 0; thrd < NOthrd; thrd++){
    for(i = 0; i < 17; i++){ best_params[thrd][i] = -99.0; }
    for(j = 0; j < 8; j++){
      for(i = 0; i < 8; i++){
	best_fit_covar[thrd][j][i] = -99.0;
      }
    }
  }
  if(vb_flag > 0){

    std::cout << "Applying LVM to start seeds . . . " << std::endl;
    std::cout << "0         25        50        75        100%" << std::endl;
    std::cout << "| | | | | : | | | | : | | | | : | | | | :" << std::endl;

  }
  progress = 0.0;
  SVD_flag = -1;
  s_done = 0;
  chi2_best = -99.0;
  fit_mode = abs(fit_mode);
  for(s = 0; s < NOs; s++){

    // get thread number
    thrd = 0;
    
    // 4a. use fit_mode to set the number of free parameters and initialise a[] array
    for(i = 0; i < 8; i++){ ia[thrd][i] = 1; }
    switch(fit_mode){
    case 1:    
      ma = 4;
      if(start_vals[s][0] > 0.0){ a[thrd][0] = log(start_vals[s][0]); } else { a[thrd][0] = -9E10; }
      if(start_vals[s][1] > 0.0){ a[thrd][1] = log(start_vals[s][1]); } else { a[thrd][1] = -9E10; }
      a[thrd][2] = (start_vals[s][2] - mid)/amp;
      if(a[thrd][2] > 1.0){ a[thrd][2] = 1.0; }
      if(a[thrd][2] < -1.0){ a[thrd][2] = -1.0; }
      a[thrd][2] = asin(a[thrd][2]);
      a[thrd][3] = (start_vals[s][4] - mid)/amp;
      if(a[thrd][3] > 1.0){ a[thrd][3] = 1.0; }
      if(a[thrd][3] < -1.0){ a[thrd][3] = -1.0; }
      a[thrd][3] = asin(a[thrd][3]);
      a[thrd][4] = a[thrd][5] = a[thrd][6] = a[thrd][7] = 0.0;
      break;
    case 2:    
      ma = 5;
      if(start_vals[s][0] > 0.0){ a[thrd][0] = log(start_vals[s][0]); } else { a[thrd][0] = -9E10; }
      if(start_vals[s][1] > 0.0){ a[thrd][1] = log(start_vals[s][1]); } else { a[thrd][1] = -9E10; }
      a[thrd][2] = (start_vals[s][2] - mid)/amp;
      if(a[thrd][2] > 1.0){ a[thrd][2] = 1.0; }
      if(a[thrd][2] < -1.0){ a[thrd][2] = -1.0; }
      a[thrd][2] = asin(a[thrd][2]);
      if(start_vals[s][3] > 0.0){ a[thrd][3] = log(start_vals[s][3]); } else { a[thrd][3] = -9E10; }
      a[thrd][4] = (start_vals[s][4] - mid)/amp;
      if(a[thrd][4] > 1.0){ a[thrd][4] = 1.0; }
      if(a[thrd][4] < -1.0){ a[thrd][4] = -1.0; }
      a[thrd][4] = asin(a[thrd][4]);
      a[thrd][5] = a[thrd][6] = a[thrd][7] = 0.0;
      break;
    case 3:    
      ma = 5;
      if(start_vals[s][0] > 0.0){ a[thrd][0] = log(start_vals[s][0]); } else { a[thrd][0] = -9E10; }
      if(start_vals[s][1] > 0.0){ a[thrd][1] = log(start_vals[s][1]); } else { a[thrd][1] = -9E10; }
      a[thrd][2] = (start_vals[s][2] - mid)/amp;
      if(a[thrd][2] > 1.0){ a[thrd][2] = 1.0; }
      if(a[thrd][2] < -1.0){ a[thrd][2] = -1.0; }
      a[thrd][2] = asin(a[thrd][2]);
      a[thrd][3] = (start_vals[s][4] - mid)/amp;
      if(a[thrd][3] > 1.0){ a[thrd][3] = 1.0; }
      if(a[thrd][3] < -1.0){ a[thrd][3] = -1.0; }
      a[thrd][3] = asin(a[thrd][3]);
      //if(start_vals[s][5] > 0.0){ a[thrd][4] = log(start_vals[s][5]); } else { a[thrd][4] = -9E10; }
      a[thrd][4] = start_vals[s][5];
      a[thrd][5] = a[thrd][6] = a[thrd][7] = 0.0;
      break;
    case 4:
      ma = 6;
      if(start_vals[s][0] > 0.0){ a[thrd][0] = log(start_vals[s][0]); } else { a[thrd][0] = -9E10; }
      if(start_vals[s][1] > 0.0){ a[thrd][1] = log(start_vals[s][1]); } else { a[thrd][1] = -9E10; }
      a[thrd][2] = (start_vals[s][2] - mid)/amp;
      if(a[thrd][2] > 1.0){ a[thrd][2] = 1.0; }
      if(a[thrd][2] < -1.0){ a[thrd][2] = -1.0; }
      a[thrd][2] = asin(a[thrd][2]);
      a[thrd][3] = (start_vals[s][4] - mid)/amp;
      if(a[thrd][3] > 1.0){ a[thrd][3] = 1.0; }
      if(a[thrd][3] < -1.0){ a[thrd][3] = -1.0; }
      a[thrd][3] = asin(a[thrd][3]);
      //if(start_vals[s][5] > 0.0){ a[thrd][4] = log(start_vals[s][5]); } else { a[thrd][4] = -9E10; }
      a[thrd][4] = start_vals[s][5];
      //a[thrd][5] = (start_vals[s][7] - 5.0)/3.0;
      a[thrd][5] = (start_vals[s][7] - 4.5)/3.5;
      if(a[thrd][5] > 1.0){ a[thrd][5] = 1.0; }
      if(a[thrd][5] < -1.0){ a[thrd][5] = -1.0; }
      a[thrd][5] = asin(a[thrd][5]);
      a[thrd][6] = a[thrd][7] = 0.0;
      break;
    case 5:
      ma = 7;
      if(start_vals[s][0] > 0.0){ a[thrd][0] = log(start_vals[s][0]); } else { a[thrd][0] = -9E10; }
      if(start_vals[s][1] > 0.0){ a[thrd][1] = log(start_vals[s][1]); } else { a[thrd][1] = -9E10; }
      a[thrd][2] = (start_vals[s][2] - mid)/amp;
      if(a[thrd][2] > 1.0){ a[thrd][2] = 1.0; }
      if(a[thrd][2] < -1.0){ a[thrd][2] = -1.0; }
      a[thrd][2] = asin(a[thrd][2]);
      a[thrd][3] = (start_vals[s][4] - mid)/amp;
      if(a[thrd][3] > 1.0){ a[thrd][3] = 1.0; }
      if(a[thrd][3] < -1.0){ a[thrd][3] = -1.0; }
      a[thrd][3] = asin(a[thrd][3]);
      //if(start_vals[s][5] > 0.0){ a[thrd][4] = log(start_vals[s][5]); } else { a[thrd][4] = -9E10; }
      a[thrd][4] = start_vals[s][5];
      a[thrd][5] = (start_vals[s][6] - mid)/amp;
      if(a[thrd][5] > 1.0){ a[thrd][5] = 1.0; }
      if(a[thrd][5] < -1.0){ a[thrd][5] = -1.0; }
      a[thrd][5] = asin(a[thrd][5]);
      //a[thrd][6] = (start_vals[s][7] - 5.0)/3.0;
      a[thrd][6] = (start_vals[s][7] - 4.5)/3.5;
      if(a[thrd][6] > 1.0){ a[thrd][6] = 1.0; }
      if(a[thrd][6] < -1.0){ a[thrd][6] = -1.0; }
      a[thrd][6] = asin(a[thrd][6]);
      a[thrd][7] = 0.0;
      break;
    default:
      ma = 8;
      if(start_vals[s][0] > 0.0){ a[thrd][0] = log(start_vals[s][0]); } else { a[thrd][0] = -9E10; }
      if(start_vals[s][1] > 0.0){ a[thrd][1] = log(start_vals[s][1]); } else { a[thrd][1] = -9E10; }
      a[thrd][2] = (start_vals[s][2] - mid)/amp;
      if(a[thrd][2] > 1.0){ a[thrd][2] = 1.0; }
      if(a[thrd][2] < -1.0){ a[thrd][2] = -1.0; }
      a[thrd][2] = asin(a[thrd][2]);
      if(start_vals[s][3] > 0.0){ a[thrd][3] = log(start_vals[s][3]); } else { a[thrd][3] = -9E10; }
      a[thrd][4] = (start_vals[s][4] - mid)/amp;
      if(a[thrd][4] > 1.0){ a[thrd][4] = 1.0; }
      if(a[thrd][4] < -1.0){ a[thrd][4] = -1.0; }
      a[thrd][4] = asin(a[thrd][4]);
      //if(start_vals[s][5] > 0.0){ a[thrd][5] = log(start_vals[s][5]); } else { a[thrd][5] = -9E10; }
      a[thrd][5] = start_vals[s][5];
      a[thrd][6] = (start_vals[s][6] - mid)/amp;
      if(a[thrd][6] > 1.0){ a[thrd][6] = 1.0; }
      if(a[thrd][6] < -1.0){ a[thrd][6] = -1.0; }
      a[thrd][6] = asin(a[thrd][6]);
      //a[thrd][7] = (start_vals[s][7] - 5.0)/3.0;
      a[thrd][7] = (start_vals[s][7] - 4.5)/3.5;
      if(a[thrd][7] > 1.0){ a[thrd][7] = 1.0; }
      if(a[thrd][7] < -1.0){ a[thrd][7] = -1.0; }
      a[thrd][7] = asin(a[thrd][7]);
      break;
    }
    
    // 4b. iterate until convergence is achieved for this starting position
    alamda[thrd] = prev_chi2_val = chi2_val = -99.0;
    iter = btr_cnt = 0;
    while((prev_chi2_val < 0.0) || (chi2_val >= prev_chi2_val) || ((prev_chi2_val - chi2_val) > (prev_chi2_val * 1E-4)) || (chi2_val > 8.999E30) || (btr_cnt < 5)){
      
      prev_chi2_val = chi2_val;
      SVD_flag = mrqmin(x_vals,y_vals,n_vals,NOvals,a[thrd],ia[thrd],ma,covar[thrd],alpha[thrd],&chi2_val,funcs,&alamda[thrd],fit_mode,mid,amp,mfit,ochisq,atry[thrd],beta[thrd],da[thrd],oneda[thrd]);
      if(std::isinf(chi2_val)){ chi2_val = (chi2_val > 0.0) ? 9E30 : prev_chi2_val + 1.0; }
      if(std::isnan(chi2_val)){ chi2_val = prev_chi2_val + 1.0; }
      if(chi2_val <= prev_chi2_val){ btr_cnt++; } else { btr_cnt = 0; }
      iter++;
      if(iter >= iter_max){ break; }
      
    }
      
    // 4c. calculate final covariance matrix
    alamda[thrd] = -1.0;
    SVD_flag = mrqmin(x_vals,y_vals,n_vals,NOvals,a[thrd],ia[thrd],ma,covar[thrd],alpha[thrd],&chi2_val,funcs,&alamda[thrd],fit_mode,mid,amp,mfit,ochisq,atry[thrd],beta[thrd],da[thrd],oneda[thrd]);
    if(std::isinf(chi2_val)){ chi2_val = (chi2_val > 0.0) ? 9E30 : prev_chi2_val + 1.0; }
    if(std::isnan(chi2_val)){ chi2_val = prev_chi2_val + 1.0; }
    
    if((chi2_val > 0.0) && ((chi2_val < chi2_best) || (chi2_best < 0.0))){
            
      chi2_best = chi2_val;
           
      for(i = 0; i < 8; i++){
	for(j = 0; j < 8; j++){  
	  best_fit_covar[thrd][i][j] = 0.0;
	}
      }
           
      best_params[thrd][0] = exp(a[thrd][0]);
      if(std::isinf(best_params[thrd][0])){ best_params[thrd][0] = (best_params[thrd][0] > 0.0) ? 9E30 : -9E30; }
      if(std::isnan(best_params[thrd][0])){ best_params[thrd][0] = 0.0; }
      best_params[thrd][1] = best_params[thrd][0] * sqrt(covar[thrd][0][0]);
      switch(fit_mode){
      case 1:
	best_params[thrd][2] = exp(a[thrd][1]);
	if(std::isinf(best_params[thrd][2])){ best_params[thrd][2] = (best_params[thrd][2] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][2])){ best_params[thrd][2] = 0.0; }
	best_params[thrd][3] = best_params[thrd][2] * sqrt(covar[thrd][1][1]);
	best_params[thrd][4] = mid + (amp * (sin(a[thrd][2])));
	best_params[thrd][5] = amp * fabs(cos(a[thrd][2])) * (sqrt(covar[thrd][2][2]));
	best_params[thrd][6] = best_params[thrd][2];
	best_params[thrd][7] = best_params[thrd][3];
	best_params[thrd][8] = mid + (amp * (sin(a[thrd][3])));
	best_params[thrd][9] = amp * fabs(cos(a[thrd][3])) * (sqrt(covar[thrd][3][3]));
	best_params[thrd][10] = 0.0;
	best_params[thrd][11] = 0.0;
	best_params[thrd][12] = 0.5 * (best_params[thrd][4] + best_params[thrd][8]);
	best_params[thrd][13] = 0.0;
	best_params[thrd][14] = 4.0;
	best_params[thrd][15] = 0.0;
	best_params[thrd][16] = chi2_val;
	for(i = 0; i < 4; i++){
	  for(j = 0; j < 4; j++){
	    best_fit_covar[thrd][i][j] = covar[thrd][i][j];
	  }
	}
	for(j = 7; j > 3; j--){
	  best_fit_covar[thrd][j][j] = best_fit_covar[thrd][(j - 1)][(j - 1)];
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][i][j] = best_fit_covar[thrd][i][(j - 1)]; }
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][j][i] = best_fit_covar[thrd][(j - 1)][i]; }
	}
	for(i = 0; i < 8; i++){ best_fit_covar[thrd][i][3] = best_fit_covar[thrd][3][i] = best_fit_covar[thrd][i][5] = best_fit_covar[thrd][5][i] = best_fit_covar[thrd][i][6] = best_fit_covar[thrd][6][i] = best_fit_covar[thrd][i][7] = best_fit_covar[thrd][7][i] = 0.0; }
	break;
      case 2:
	best_params[thrd][2] = exp(a[thrd][1]);
	if(std::isinf(best_params[thrd][2])){ best_params[thrd][2] = (best_params[thrd][2] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][2])){ best_params[thrd][2] = 0.0; }
	best_params[thrd][3] = best_params[thrd][2] * sqrt(covar[thrd][1][1]);
	best_params[thrd][4] = mid + (amp * (sin(a[thrd][2])));
	best_params[thrd][5] = amp * fabs(cos(a[thrd][2])) * (sqrt(covar[thrd][2][2]));
	best_params[thrd][6] = exp(a[thrd][3]);
	if(std::isinf(best_params[thrd][6])){ best_params[thrd][6] = (best_params[thrd][6] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][6])){ best_params[thrd][6] = 0.0; }
	best_params[thrd][7] = best_params[thrd][6] * sqrt(covar[thrd][3][3]);
	best_params[thrd][8] = mid + (amp * (sin(a[thrd][4])));
	best_params[thrd][9] = amp * fabs(cos(a[thrd][4])) * (sqrt(covar[thrd][4][4]));
	best_params[thrd][10] = 0.0;
	best_params[thrd][11] = 0.0;
	best_params[thrd][12] = 0.5 * (best_params[thrd][4] + best_params[thrd][8]);
	best_params[thrd][13] = 0.0;
	best_params[thrd][14] = 4.0;
	best_params[thrd][15] = 0.0;
	best_params[thrd][16] = chi2_val;
	for(i = 0; i < 5; i++){
	  for(j = 0; j < 5; j++){
	    best_fit_covar[thrd][i][j] = covar[thrd][i][j];
	  }
	}
	for(i = 0; i < 8; i++){ best_fit_covar[thrd][i][5] = best_fit_covar[thrd][5][i] = best_fit_covar[thrd][i][6] = best_fit_covar[thrd][6][i] = best_fit_covar[thrd][i][7] = best_fit_covar[thrd][7][i] = 0.0; }
	break;
      case 3:
	best_params[thrd][2] = exp(a[thrd][1]);
	if(std::isinf(best_params[thrd][2])){ best_params[thrd][2] = (best_params[thrd][2] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][2])){ best_params[thrd][2] = 0.0; }
	best_params[thrd][3] = best_params[thrd][2] * sqrt(covar[thrd][1][1]);
	best_params[thrd][4] = mid + (amp * (sin(a[thrd][2])));
	best_params[thrd][5] = amp * fabs(cos(a[thrd][2])) * (sqrt(covar[thrd][2][2]));
	best_params[thrd][6] = best_params[thrd][2];
	best_params[thrd][7] = best_params[thrd][3];
	best_params[thrd][8] = mid + (amp * (sin(a[thrd][3])));
	best_params[thrd][9] = amp * fabs(cos(a[thrd][3])) * (sqrt(covar[thrd][3][3]));
	//best_params[thrd][10] = exp(a[thrd][4]);
	best_params[thrd][10] = a[thrd][4];
	if(std::isinf(best_params[thrd][10])){ best_params[thrd][10] = (best_params[thrd][10] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][10])){ best_params[thrd][10] = 0.0; }
	//best_params[thrd][11] = best_params[thrd][10] * sqrt(covar[thrd][4][4]);
	best_params[thrd][11] = sqrt(covar[thrd][4][4]);
	best_params[thrd][12] = 0.5 * (best_params[thrd][4] + best_params[thrd][8]);
	best_params[thrd][13] = 0.0;
	best_params[thrd][14] = 4.0;
	best_params[thrd][15] = 0.0;
	best_params[thrd][16] = chi2_val;
	for(i = 0; i < 5; i++){
	  for(j = 0; j < 5; j++){
	    best_fit_covar[thrd][i][j] = covar[thrd][i][j];
	  }
	}
	for(j = 7; j > 3; j--){
	  best_fit_covar[thrd][j][j] = best_fit_covar[thrd][(j - 1)][(j - 1)];
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][i][j] = best_fit_covar[thrd][i][(j - 1)]; }
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][j][i] = best_fit_covar[thrd][(j - 1)][i]; }
	}
	for(i = 0; i < 8; i++){ best_fit_covar[thrd][i][3] = best_fit_covar[thrd][3][i] = best_fit_covar[thrd][i][6] = best_fit_covar[thrd][6][i] = best_fit_covar[thrd][i][7] = best_fit_covar[thrd][7][i] = 0.0; }
	break;
      case 4:
	best_params[thrd][2] = exp(a[thrd][1]);
	if(std::isinf(best_params[thrd][2])){ best_params[thrd][2] = (best_params[thrd][2] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][2])){ best_params[thrd][2] = 0.0; }
	best_params[thrd][3] = best_params[thrd][2] * sqrt(covar[thrd][1][1]);
	best_params[thrd][4] = mid + (amp * (sin(a[thrd][2])));
	best_params[thrd][5] = amp * fabs(cos(a[thrd][2])) * (sqrt(covar[thrd][2][2]));
	best_params[thrd][6] = best_params[thrd][2];
	best_params[thrd][7] = best_params[thrd][3];
	best_params[thrd][8] = mid + (amp * (sin(a[thrd][3])));
	best_params[thrd][9] = amp * fabs(cos(a[thrd][3])) * (sqrt(covar[thrd][3][3]));
	//best_params[thrd][10] = exp(a[thrd][4]);
	best_params[thrd][10] = a[thrd][4];
	if(std::isinf(best_params[thrd][10])){ best_params[thrd][10] = (best_params[thrd][10] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][10])){ best_params[thrd][10] = 0.0; }
	//best_params[thrd][11] = best_params[thrd][10] * sqrt(covar[thrd][4][4]);
	best_params[thrd][11] = sqrt(covar[thrd][4][4]);
	best_params[thrd][12] = 0.5 * (best_params[thrd][4] + best_params[thrd][8]);
	best_params[thrd][13] = 0.0;
	//best_params[thrd][14] = 5.0 + (4.0 * (sin(a[thrd][5])));
	best_params[thrd][14] = 4.5 + (3.5 * (sin(a[thrd][5])));
	//best_params[thrd][15] = amp * fabs(cos(a[thrd][5])) * (sqrt(covar[thrd][5][5]));
	best_params[thrd][15] = 3.5 * fabs(cos(a[thrd][5])) * (sqrt(covar[thrd][5][5]));
	best_params[thrd][16] = chi2_val;
	for(i = 0; i < 6; i++){
	  for(j = 0; j < 6; j++){
	    best_fit_covar[thrd][i][j] = covar[thrd][i][j];
	  }
	}
	for(j = 7; j > 3; j--){
	  best_fit_covar[thrd][j][j] = best_fit_covar[thrd][(j - 1)][(j - 1)];
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][i][j] = best_fit_covar[thrd][i][(j - 1)]; }
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][j][i] = best_fit_covar[thrd][(j - 1)][i]; }
	}
	for(j = 7; j > 5; j--){
	  best_fit_covar[thrd][j][j] = best_fit_covar[thrd][(j - 1)][(j - 1)];
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][i][j] = best_fit_covar[thrd][i][(j - 1)]; }
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][j][i] = best_fit_covar[thrd][(j - 1)][i]; }
	}
	for(i = 0; i < 8; i++){ best_fit_covar[thrd][i][3] = best_fit_covar[thrd][3][i] = best_fit_covar[thrd][i][6] = best_fit_covar[thrd][6][i] = 0.0; }
	break;
      case 5:
	best_params[thrd][2] = exp(a[thrd][1]);
	if(std::isinf(best_params[thrd][2])){ best_params[thrd][2] = (best_params[thrd][2] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][2])){ best_params[thrd][2] = 0.0; }
	best_params[thrd][3] = best_params[thrd][2] * sqrt(covar[thrd][1][1]);
	best_params[thrd][4] = mid + (amp * (sin(a[thrd][2])));
	best_params[thrd][5] = amp * fabs(cos(a[thrd][2])) * (sqrt(covar[thrd][2][2]));
	best_params[thrd][6] = best_params[thrd][2];
	best_params[thrd][7] = best_params[thrd][3];
	best_params[thrd][8] = mid + (amp * (sin(a[thrd][3])));
	best_params[thrd][9] = amp * fabs(cos(a[thrd][3])) * (sqrt(covar[thrd][3][3]));
	//best_params[thrd][10] = exp(a[thrd][4]);
	best_params[thrd][10] = a[thrd][4];
	if(std::isinf(best_params[thrd][10])){ best_params[thrd][10] = (best_params[thrd][10] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][10])){ best_params[thrd][10] = 0.0; }
	//best_params[thrd][11] = best_params[thrd][10] * sqrt(covar[thrd][4][4]);
	best_params[thrd][11] = sqrt(covar[thrd][4][4]);
	best_params[thrd][12] = mid + (amp * (sin(a[thrd][5])));
	best_params[thrd][13] = amp * fabs(cos(a[thrd][5])) * (sqrt(covar[thrd][5][5]));
	//best_params[thrd][14] = 5.0 + (3.0 * (sin(a[thrd][6])));
	best_params[thrd][14] = 4.5 + (3.5 * (sin(a[thrd][6])));
	//best_params[thrd][15] = amp * fabs(cos(a[thrd][6])) * (sqrt(covar[thrd][6][6]));
	best_params[thrd][15] = 3.5 * fabs(cos(a[thrd][6])) * (sqrt(covar[thrd][6][6]));
	best_params[thrd][16] = chi2_val;
	for(i = 0; i < 7; i++){
	  for(j = 0; j < 7; j++){
	    best_fit_covar[thrd][i][j] = covar[thrd][i][j];
	  }
	}
	for(j = 7; j > 3; j--){
	  best_fit_covar[thrd][j][j] = best_fit_covar[thrd][(j - 1)][(j - 1)];
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][i][j] = best_fit_covar[thrd][i][(j - 1)]; }
	  for(i = 0; i < j; i++){ best_fit_covar[thrd][j][i] = best_fit_covar[thrd][(j - 1)][i]; }
	}
	for(i = 0; i < 8; i++){ best_fit_covar[thrd][i][3] = best_fit_covar[thrd][3][i] = 0.0; }
	break;
      default:
	best_params[thrd][2] = exp(a[thrd][1]);
	if(std::isinf(best_params[thrd][2])){ best_params[thrd][2] = (best_params[thrd][2] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][2])){ best_params[thrd][2] = 0.0; }
	best_params[thrd][3] = best_params[thrd][2] * sqrt(covar[thrd][1][1]);
	best_params[thrd][4] = mid + (amp * (sin(a[thrd][2])));
	best_params[thrd][5] = amp * fabs(cos(a[thrd][2])) * (sqrt(covar[thrd][2][2]));
	best_params[thrd][6] = exp(a[thrd][3]);
	if(std::isinf(best_params[thrd][6])){ best_params[thrd][6] = (best_params[thrd][6] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][6])){ best_params[thrd][6] = 0.0; }
	best_params[thrd][7] = best_params[thrd][6] * sqrt(covar[thrd][3][3]);
	best_params[thrd][8] = mid + (amp * (sin(a[thrd][4])));
	best_params[thrd][9] = amp * fabs(cos(a[thrd][4])) * (sqrt(covar[thrd][4][4]));
	//best_params[thrd][10] = exp(a[thrd][5]);
	best_params[thrd][10] = a[thrd][5];
	if(std::isinf(best_params[thrd][10])){ best_params[thrd][10] = (best_params[thrd][10] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(best_params[thrd][10])){ best_params[thrd][10] = 0.0; }
	//best_params[thrd][11] = best_params[thrd][10] * sqrt(covar[thrd][5][5]);
	best_params[thrd][11] = sqrt(covar[thrd][5][5]);
	best_params[thrd][12] = mid + (amp * (sin(a[thrd][6])));
	best_params[thrd][13] = amp * fabs(cos(a[thrd][6])) * (sqrt(covar[thrd][6][6]));
	//best_params[thrd][14] = 5.0 + (3.0 * (sin(a[thrd][7])));
	best_params[thrd][14] = 4.5 + (3.5 * (sin(a[thrd][7])));
	//best_params[thrd][15] = amp * fabs(cos(a[thrd][7])) * (sqrt(covar[thrd][7][7]));
	best_params[thrd][15] = 3.5 * fabs(cos(a[thrd][7])) * (sqrt(covar[thrd][7][7]));
	best_params[thrd][16] = chi2_val;
	for(i = 0; i < 8; i++){
	  for(j = 0; j < 8; j++){
	    best_fit_covar[thrd][i][j] = covar[thrd][i][j];
	  }
	}
	break;
      }
            
      // if(chi2_val < chi2_best)
    }
    
    // display progress
    if(vb_flag > 0){
      s_done++;
      while((((float) (s_done + 1)) / ((float) NOs)) >= progress){ std::cout << "*" << std::flush; progress+=0.025; }
    }    
    
    // for(s = 0; s < NOs; s++)
  }
  if(vb_flag > 0){ std::cout << " done." << std::endl; }
  
  // 5. If any starting position succeeded, return the LVM results. Notify the user if no start position converged.
  if(SVD_flag < 0){ std::cout << "WARNING!!! None of the LVM-J starting positions converged. Refinement STRONGLY recommended." << std::endl; }

  // update model_params array
  for(s = 0, thrd = 1; thrd < NOthrd; thrd++){ if(((best_params[s][16] <= 0.0) || (best_params[thrd][16] < best_params[s][16])) && (best_params[thrd][16] > 0.0)){ s = thrd; } }
  for(i = 0; i < 17; i++){ model_params[i] = best_params[s][i]; }
  for(i = 0; i < 8; i++){ model_params[(1 + (2*i))] = fabs(model_params[(1 + (2*i))]); }
  for(i = 0; i < 8; i++){
    for(j = 0; j < 8; j++){
      fit_covar[i][j] = best_fit_covar[s][i][j];
    }
  }
  
  delete [] alamda;
  alamda = NULL;
  for(thrd = 0; thrd < NOthrd; thrd++){
    
    for(i = 0; i < 8; i++){

      delete [] covar[thrd][i];
      covar[thrd][i] = NULL;
      delete [] alpha[thrd][i];
      alpha[thrd][i] = NULL;
      delete [] oneda[thrd][i];
      oneda[thrd][i] = NULL;

    }

    delete [] covar[thrd];
    covar[thrd] = NULL;
    delete [] alpha[thrd];
    alpha[thrd] = NULL;
    delete [] a[thrd];
    a[thrd] = NULL;
    delete [] ia[thrd];
    ia[thrd] = NULL;

    delete [] best_params[thrd];
    best_params[thrd] = NULL;
    delete [] best_fit_covar[thrd];
    best_fit_covar[thrd] = NULL;
    delete [] oneda[thrd];
    oneda[thrd] = NULL;
    delete [] da[thrd];
    da[thrd] = NULL;
    delete [] beta[thrd];
    beta[thrd] = NULL;
    delete [] atry[thrd];
    atry[thrd] = NULL;

  }
  delete [] ia;
  ia = NULL;
  delete [] a;
  a = NULL;
  delete [] covar;
  covar = NULL;
  delete [] alpha;
  alpha = NULL;
  delete [] best_params;
  best_params = NULL;
  delete [] best_fit_covar;
  best_fit_covar = NULL;
  delete [] oneda;
  oneda = NULL;
  delete [] da;
  da = NULL;
  delete [] atry;
  atry = NULL;
  delete [] beta;
  beta = NULL;

  return;

}

template <class T_count, class T_xvals, class T_data, class T_result>
  int FitBusyFunc(T_count NOvals, T_xvals * x_vals, T_data * y_vals, T_data * n_vals, T_result * fit_params, T_result ** fit_covar, int & best_NOp, int NOs, int iter_max, int vb_flag){

  double model_params[2][17], y_max, y_scale = 1.0, ** start_vals, ** temp_covar, mid, amp;
  //double chi2_val, model_val;
  int i, j, best_BFtype, y_scale_flag;

  // range test the values of NOs and iter_max
  if(iter_max <= 0){ 
    std::cout << "ERROR: iter_max must be greater than zero. Exiting." << std::endl;
    return 0;
  }
  if((NOs < 0) && (NOs != -1)){
    std::cout << "ERROR: NOs is not greater than 0 or equal to -1 (flag value). Exiting." << std::endl;
    return 0;
  }

  if(vb_flag >= 0){ std::cout << "Fitting Busy Function . . . " << std::endl; }

  // create arrays 
  temp_covar = new double * [8];
  for(i = 0; i < 8; i++){ temp_covar[i] = new double[8]; }
  start_vals = new double * [(int) abs(NOs)];
  for(i = 0; i < (int) abs(NOs); i++){ start_vals[i] = new double[9]; }

  // initialise the start_vals array if NOs == -1
  if(NOs == -1){
 
    for(i = 0; i < 8; i++){ start_vals[0][i] = fit_params[(2 * i)]; }
    start_vals[0][8] = 9E10;
    
  }

  // initialise the mid-point and amplitude used to limit the roll-off and power-law positions during fitting using variable re-mapping
  if(NOs > 0){
    fit_params[4] = (fit_params[4] >= x_vals[0]) ? fit_params[4] : x_vals[0];
    fit_params[8] = (fit_params[8] <= x_vals[(NOvals - 1)]) ? fit_params[8] : x_vals[(NOvals - 1)];
    mid = 0.5 * (double) (fit_params[4] + fit_params[8]);
    amp = 0.5 * (double) (fit_params[8] - fit_params[4]);
  } else {
    mid = 0.5 * (double) (fit_params[4] + fit_params[8]);
    amp = (double) fit_params[4];
    fit_params[4] = mid - (double) (fit_params[8] - fit_params[4]);
    fit_params[8] = mid + (double) fit_params[8] - amp;
    fit_params[4] = (fit_params[4] >= x_vals[0]) ? fit_params[4] : x_vals[0];
    fit_params[8] = (fit_params[8] <= x_vals[(NOvals - 1)]) ? fit_params[8] : x_vals[(NOvals - 1)];
    mid = 0.5 * (double) (fit_params[4] + fit_params[8]);
    amp = 0.5 * (double) (fit_params[8] - fit_params[4]);
  }

  // initialise the fit_params array
  if(vb_flag >= 0){ std::cout << "Initialising fit params to: " << std::flush; }
  for(i = 0; i < 17; i++){ 

    model_params[0][i] = model_params[1][i] = (double) fit_params[i]; 
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
    if(std::isinf(model_params[1][i])){ model_params[1][i] = (model_params[1][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[1][i])){ model_params[1][i] = 0.0; }
    if(vb_flag >= 0){ std::cout << model_params[0][i] << " " << std::flush; }
    
  }
  if(vb_flag >= 0){ std::cout << std::endl; }

  // determine the y-axis range of the input data
  y_max = -9E30;
  for(i = 0; i < NOvals; i++){
    if((fabs((double) y_vals[i])) >= y_max){ y_max = fabs((double) y_vals[i]); }
  }

  // determine if the y-axis needs to be re-scaled and do so 
  // (includes re-scaling model_params[] and fit_params[] (+start_vals[] for NOs == -1))
  y_scale_flag = -1;
  if(y_max >= 1E7){

    y_scale = y_max / 1E7;
    if(vb_flag >= 0){ std::cout << "Re-scaling with " << y_scale << " . . . " << std::endl; }
    for(i = 0; i < NOvals; i++){ 

      y_vals[i]/=((T_data) y_scale); 
      n_vals[i]/=((T_data) y_scale);

    }
    model_params[0][0]/=y_scale;
    model_params[1][0]/=y_scale;
    if(NOs == -1){ start_vals[0][0]/=y_scale; }
    y_scale_flag = 1;

  }
 
  // fit_mode == 1 -- 4 parameters; maximum symmetry and no power-law
  if(NOs > 0){
    if(vb_flag >= 0){ std::cout << "Trying: fit mode 1 . . .  " << std::endl; }
    FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],1,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  } else {
    if(vb_flag >= 0){ std::cout << "Trying: fit mode -1 . . .  " << std::endl; }
    NOs = 1;
    FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],-1,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  }
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
    model_params[1][i] = model_params[0][i];
    if(vb_flag >= 0){ std::cout << model_params[1][i] << " " << std::flush; }
  }
  if(vb_flag >= 0){ std::cout << ": " << model_params[1][16] << std::endl; }
  for(i = 0; i < 8; i++){
    for(j = 0; j < 8; j++){
      fit_covar[i][j] = (T_result) temp_covar[i][j];
      if(std::isinf(fit_covar[i][j])){ fit_covar[i][j] = (fit_covar[i][j] > 0.0) ? 9E30 : -9E30; }
      if(std::isnan(fit_covar[i][j])){ fit_covar[i][j] = 0.0; }
    }
  }
  best_NOp = 4;
  best_BFtype = 1;
  if(vb_flag >= 0){ std::cout << "best_NOp = " << best_NOp << ", best_BFtype = " << best_BFtype << std::endl; }

  // fit_mode == 2 -- 5 parameters; no forced symmetry, but no power-law
  if(vb_flag >= 0){ std::cout << "Trying: fit mode -2 . . .  " << std::endl; }
  FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],-2,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
  }
  if((model_params[0][16] + (2 * (5 - best_NOp))) < model_params[1][16]){
    for(i = 0; i < 17; i++){
      model_params[1][i] = model_params[0][i];
      if(vb_flag >= 0){ std::cout << model_params[1][i] << " " << std::flush; }
    }
    if(vb_flag >= 0){ std::cout << ": " << model_params[1][16] << " ---> " << (model_params[0][16] + (2 * (5 - best_NOp))) << std::endl; }
    for(i = 0; i < 8; i++){
      for(j = 0; j < 8; j++){
	fit_covar[i][j] = (T_result) temp_covar[i][j];
	if(std::isinf(fit_covar[i][j])){ fit_covar[i][j] = (fit_covar[i][j] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(fit_covar[i][j])){ fit_covar[i][j] = 0.0; }
      }
    }
    best_NOp = 5;
    best_BFtype = 2;
  } else {
    if(vb_flag >= 0){ 
      std::cout << "Too many parameters! " << model_params[0][16] << " + " << (2 * (5 - best_NOp)) << " >= " << model_params[1][16] << " for model: " << std::flush; 
      for(i = 0; i < 17; i++){
	std::cout << model_params[0][i] << " (" << model_params[1][i] << ") " << std::flush;
      }
      std::cout << std::endl;
    }
  }
  if(vb_flag >= 0){ std::cout << "best_NOp = " << best_NOp << ", best_BFtype = " << best_BFtype << std::endl; }

  // fit_mode == 3 -- 5 parameters; maximum symmetry and power-law = 4
  if(vb_flag >= 0){ std::cout << "Trying: fit mode -3 . . .  " << std::endl; }
  FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],-3,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
  }
  if((model_params[0][16] + (2 * (5 - best_NOp))) < model_params[1][16]){
    for(i = 0; i < 17; i++){
      model_params[1][i] = model_params[0][i];
      if(vb_flag >= 0){ std::cout << model_params[1][i] << " " << std::flush; }
    }
    if(vb_flag >= 0){ std::cout << ": " << model_params[1][16] << " ---> " << (model_params[0][16] + (2 * (5 - best_NOp))) << std::endl; }
    for(i = 0; i < 8; i++){
      for(j = 0; j < 8; j++){
	fit_covar[i][j] = (T_result) temp_covar[i][j];
	if(std::isinf(fit_covar[i][j])){ fit_covar[i][j] = (fit_covar[i][j] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(fit_covar[i][j])){ fit_covar[i][j] = 0.0; }
      }
    }
    best_NOp = 5;
    best_BFtype = 3;
  } else {
    if(vb_flag >= 0){
      std::cout << "Too many parameters! " << model_params[0][16] << " + " << (2 * (5 - best_NOp)) << " >= " << model_params[1][16] << " for model: " << std::flush;
      for(i = 0; i < 17; i++){
	std::cout << model_params[0][i] << " (" << model_params[1][i] << ") " << std::flush;
      }
      std::cout << std::endl;
    }
  }
  if(vb_flag >= 0){ std::cout << "best_NOp = " << best_NOp << ", best_BFtype = " << best_BFtype << std::endl; }

  // fit_mode == 4 -- 6 parameters; maximum symmetry and power-law = ? ==> use previous fitting results
  if(vb_flag >= 0){ std::cout << "Trying: fit mode -4 . . .  " << std::endl; }
  FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],-4,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
  }
  if((model_params[0][16] + (2 * (6 - best_NOp))) < model_params[1][16]){
    for(i = 0; i < 17; i++){
      model_params[1][i] = model_params[0][i];
      if(vb_flag >= 0){ std::cout << model_params[1][i] << " " << std::flush; }
    }
    if(vb_flag >= 0){ std::cout << ": " << model_params[1][16] << " ---> " << (model_params[0][16] + (2 * (6 - best_NOp))) << std::endl; }
    for(i = 0; i < 8; i++){
      for(j = 0; j < 8; j++){
	fit_covar[i][j] = (T_result) temp_covar[i][j];
	if(std::isinf(fit_covar[i][j])){ fit_covar[i][j] = (fit_covar[i][j] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(fit_covar[i][j])){ fit_covar[i][j] = 0.0; }
      }
    }
    best_NOp = 6;
    best_BFtype = 4;
  } else {
    if(vb_flag >= 0){
      std::cout << "Too many parameters! " << model_params[0][16] << " + " << (2 * (6 - best_NOp)) << " >= " << model_params[1][16] << " for model: " << std::flush;
      for(i = 0; i < 17; i++){
	std::cout << model_params[0][i] << " (" << model_params[1][i] << ") " << std::flush;
      }
      std::cout << std::endl;
    }
  }
  if(vb_flag >= 0){ std::cout << "best_NOp = " << best_NOp << ", best_BFtype = " << best_BFtype << std::endl; }
   
  // fit_mode == 5 -- 7 parameters; intermediate symmetry (symmetric error functions, but power-law is off-centre) and power-law = ? ==> use previous fitting results
  if(vb_flag >= 0){ std::cout << "Trying: fit mode -5 . . .  " << std::endl; }
  FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],-5,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
  }
  if((model_params[0][16] + (2 * (7 - best_NOp))) < model_params[1][16]){
    for(i = 0; i < 17; i++){
      model_params[1][i] = model_params[0][i];
      if(vb_flag >= 0){ std::cout << model_params[1][i] << " " << std::flush; }
    }
    if(vb_flag >= 0){ std::cout << ": " << model_params[1][16] << " ---> " << (model_params[0][16] + (2 * (7 - best_NOp))) << std::endl; }
    for(i = 0; i < 8; i++){
      for(j = 0; j < 8; j++){
	fit_covar[i][j] = (T_result) temp_covar[i][j];
	if(std::isinf(fit_covar[i][j])){ fit_covar[i][j] = (fit_covar[i][j] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(fit_covar[i][j])){ fit_covar[i][j] = 0.0; }
      }
    }
    best_NOp = 7;
    best_BFtype = 5;
  } else {
    if(vb_flag >= 0){ 
      std::cout << "Too many parameters! " << model_params[0][16] << " + " << (2 * (7 - best_NOp)) << " >= " << model_params[1][16] << " for model: " << std::flush;
      for(i = 0; i < 17; i++){
	std::cout << model_params[0][i] << " (" << model_params[1][i] << ") " << std::flush;
      }
      std::cout << std::endl;
    }
  }
  if(vb_flag >= 0){ std::cout << "best_NOp = " << best_NOp << ", best_BFtype = " << best_BFtype << std::endl; }
  
  // fit_mode == 6 -- 8 parameters; no forced symmetry and power-law = ? ==> use previous fitting results
  if(vb_flag >= 0){ std::cout << "Trying: fit mode -6 . . .  " << std::endl; }
  FitBusyFunc_engine(NOvals,x_vals,y_vals,n_vals,model_params[0],-6,NOs,start_vals,mid,amp,temp_covar,iter_max,vb_flag);
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[0][i])){ model_params[0][i] = (model_params[0][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[0][i])){ model_params[0][i] = 0.0; }
  }
  if((model_params[0][16] + (2 * (8 - best_NOp))) < model_params[1][16]){
    for(i = 0; i < 17; i++){
      model_params[1][i] = model_params[0][i];
      if(vb_flag >= 0){ std::cout << model_params[1][i] << " " << std::flush; }
    }
    if(vb_flag >= 0){ std::cout << ": " << model_params[1][16] << " ---> " << (model_params[0][16] + (2 * (8 - best_NOp))) << std::endl; }
    for(i = 0; i < 8; i++){
      for(j = 0; j < 8; j++){
	fit_covar[i][j] = (T_result) temp_covar[i][j];
	if(std::isinf(fit_covar[i][j])){ fit_covar[i][j] = (fit_covar[i][j] > 0.0) ? 9E30 : -9E30; }
	if(std::isnan(fit_covar[i][j])){ fit_covar[i][j] = 0.0; }
      }
    }
    best_NOp = 8;
    best_BFtype = 6;
  } else {
    if(vb_flag >= 0){
      std::cout << "Too many parameters! " << model_params[0][16] << " + " << (2 * (8 - best_NOp)) << " >= " << model_params[1][16] << " for model: " << std::flush;
      for(i = 0; i < 17; i++){
	std::cout << model_params[0][i] << " (" << model_params[1][i] << ") " << std::flush;
      }
      std::cout << std::endl;
    }
  }

  if(vb_flag >= 0){ 

    std::cout << "best_NOp = " << best_NOp << ", best_BFtype = " << best_BFtype << std::endl; 
    std::cout << "Best fitting model is fit_mode = " << best_BFtype << std::endl;

  }

  // free up memory 
  for(i = 0; i < NOs; i++){

    delete [] start_vals[i];
    start_vals[i] = NULL;

  }
  delete [] start_vals;
  start_vals = NULL;
  for(i = 0; i < 8; i++){

    delete [] temp_covar[i];
    temp_covar[i] = NULL;

  }
  delete [] temp_covar;
  temp_covar = NULL;

  // if required, de-scale the input y_vals[] and the model parameters
  if(y_scale_flag > 0){

    if(vb_flag >= 0){ std::cout << "De-scaling with " << y_scale << " . . . " << std::endl; }
    for(i = 0; i < NOvals; i++){ 

      y_vals[i]*=((T_data) y_scale);
      n_vals[i]*=((T_data) y_scale);

    }
    model_params[1][0]*=y_scale;
    model_params[1][1]*=y_scale;

  }

  // return best fitting model (incorporating AIC penalties) and the number of model parameters
  for(i = 0; i < 17; i++){
    if(std::isinf(model_params[1][i])){ model_params[1][i] = (model_params[1][i] > 0.0) ? 9E30 : -9E30; }
    if(std::isnan(model_params[1][i])){ model_params[1][i] = 0.0; }
    fit_params[i] = (T_result) model_params[1][i];
  }
  return best_BFtype;
  
}

void covsrt(double **covar, int ma, int ia[], int mfit);

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

void svdcmp(double **a, int m, int n, double w[], double **v);
 
double ran2(long &idum);

double pythag(double a, double b);

// functions used to generate an arbitrary number of correlated Busy Function fit parameters,
// and derive uncertainties in observational properties

void choldc(double ** a, int n);

double inverf(double p);

template <class T_fit, class T_range, class T_obs>
void CreateRandFits(int NOr, T_obs ** rand_fits, int fit_type, T_fit * BF_fit, T_fit ** BF_covar, T_range x_min, T_range x_max, int vb_flag){
  
  int i,j,r,r_done;
  float progress;
  double ** BF_corr, BF_errs[8], rand_val, z_val, cfd_min[8], cfd_max[8], rand_set[8], mid, amp, temp_cos[4];
  long int seed;

  void choldc(double ** a, int n);
  double ran2(long &seed);

  if(vb_flag >= 0){ std::cout << "Generating random variants of Busy Function fit . . . " << std::endl; }

  if(vb_flag >= 0){ std::cout << "Identified Busy Function type = " << fit_type << std::endl; }

  // calculate mid and amp --- used in variable mapping
  mid = 0.5 * (double) (x_min + x_max);
  amp = 0.5 * (double) (x_max - x_min);

  // calculate temp_cos values --- used in variable mappping
  temp_cos[0] = cos(asin(((double) BF_fit[4] - mid) / amp));
  if(std::isnan(temp_cos[0])){ temp_cos[0] = 0.0; }
  if(std::isinf(temp_cos[0])){ temp_cos[0] = (temp_cos[0] > 0.0) ? 9E30 : -9E30; }
  if(temp_cos[0] > 1.0){ temp_cos[0] = 1.0; }
  if(temp_cos[0] < -1.0){ temp_cos[0] = -1.0; }
  temp_cos[1] = cos(asin(((double) BF_fit[8] - mid) / amp));
  if(std::isnan(temp_cos[1])){ temp_cos[1] = 0.0; }
  if(std::isinf(temp_cos[1])){ temp_cos[1] = (temp_cos[1] > 0.0) ? 9E30 : -9E30; }
  if(temp_cos[1] > 1.0){ temp_cos[1] = 1.0; }
  if(temp_cos[1] < -1.0){ temp_cos[1] = -1.0; }
  temp_cos[2] = cos(asin(((double) BF_fit[12] - mid) / amp));
  if(std::isnan(temp_cos[2])){ temp_cos[2] = 0.0; }
  if(std::isinf(temp_cos[2])){ temp_cos[2] = (temp_cos[2] > 0.0) ? 9E30 : -9E30; }
  if(temp_cos[2] > 1.0){ temp_cos[2] = 1.0; }
  if(temp_cos[2] < -1.0){ temp_cos[2] = -1.0; }
  //temp_cos[3] = cos(asin(((double) BF_fit[14] - 5.0) / 3.0));
  temp_cos[3] = cos(asin(((double) BF_fit[14] - 4.5) / 3.5));
  if(std::isnan(temp_cos[3])){ temp_cos[3] = 0.0; }
  if(std::isinf(temp_cos[3])){ temp_cos[3] = (temp_cos[3] > 0.0) ? 9E30 : -9E30; }
  if(temp_cos[3] > 1.0){ temp_cos[3] = 1.0; }
  if(temp_cos[3] < -1.0){ temp_cos[3] = -1.0; }

  // calculate uncertainties in mapped BF variables
  for(i = 0; i < 8; i++){

    BF_errs[i] = sqrt((double) BF_covar[i][i]);
    if(std::isnan(BF_errs[i])){ BF_errs[i] = 0.0; }
    if(std::isinf(BF_errs[i])){ BF_errs[i] = (BF_errs[i] > 0.0) ? 9E30 : -9E30; }

  }

  // calculate minimum and maximum cfd values for each parameter
  if(BF_errs[0] > 0.0){
    cfd_min[0] = 0.5 * (1.0 + erf((double) (-1.0 / (BF_errs[0] * sqrt(2.0)))));
    if(std::isnan(cfd_min[0]) || (cfd_min[0] < 0.0)){ cfd_min[0] = 0.0; }
    if(cfd_min[0] > 1.0){ cfd_min[0] = 1.0; }
    cfd_max[0] = 1.0;
  } else {
    cfd_min[0] = cfd_max[0] = 0.5;
  }

  if(BF_errs[1] > 0.0){
    cfd_min[1] = 0.5 * (1.0 + erf((double) (-1.0 / (BF_errs[1] * sqrt(2.0)))));
    if(std::isnan(cfd_min[1]) || (cfd_min[1] < 0.0)){ cfd_min[1] = 0.0; }
    if(cfd_min[1] > 1.0){ cfd_min[1] = 1.0; }
    cfd_max[1] = 1.0;
  } else {
    cfd_max[1] = cfd_min[1] = 0.5;
  }

  if(BF_errs[2] > 0.0){
    cfd_min[2] = 0.5 * (1.0 + erf((double) ((x_min - BF_fit[4]) / (amp * BF_errs[2] * temp_cos[0] * sqrt(2.0)))));
    if(std::isnan(cfd_min[2]) || (cfd_min[2] < 0.0)){ cfd_min[2] = 0.0; }
    if(cfd_min[2] > 1.0){ cfd_min[2] = 1.0; }
    cfd_max[2] = 0.5 * (1.0 + erf((double) ((x_max - BF_fit[4]) / (amp * BF_errs[2] * temp_cos[0] * sqrt(2.0)))));
    if(std::isnan(cfd_max[2]) || (cfd_max[2] < 0.0)){ cfd_max[2] = 0.0; }
    if(cfd_max[2] > 1.0){ cfd_max[2] = 1.0; }
  } else {
    cfd_min[2] = cfd_max[2] = 0.5;
  }

  if(BF_errs[3] > 0.0){
    cfd_min[3] = 0.5 * (1.0 + erf((double) (-1.0 / (BF_errs[3] * sqrt(2.0)))));
    if((std::isnan(cfd_min[3])) || (cfd_min[3] < 0.0)){ cfd_min[3] = 0.0; }
    if(cfd_min[3] > 1.0){ cfd_min[3] = 1.0; }
    cfd_max[3] = 1.0;
  } else {
    cfd_min[3] = cfd_max[3] = 0.5;
  }

  if(BF_errs[4] > 0.0){
    cfd_min[4] = 0.5 * (1.0 + erf((double) ((x_min - BF_fit[8]) / (amp * BF_errs[4] * temp_cos[1] * sqrt(2.0)))));
    if(std::isnan(cfd_min[4]) || (cfd_min[4] < 0.0)){ cfd_min[4] = 0.0; }
    if(cfd_min[4] > 1.0){ cfd_min[4] = 1.0; }
    cfd_max[4] = 0.5 * (1.0 + erf((double) ((x_max - BF_fit[8]) / (amp * BF_errs[4] * temp_cos[1] * sqrt(2.0)))));
    if(std::isnan(cfd_max[4]) || (cfd_max[4] < 0.0)){ cfd_max[4] = 0.0; }
    if(cfd_max[4] > 1.0){ cfd_max[4] = 1.0; }
  } else {
    cfd_min[4] = cfd_max[4] = 0.5;
  }

  if(BF_errs[5] > 0.0){
    //cfd_min[5] = 0.5 * (1.0 + erf((double) (-1.0 / (BF_errs[5] * sqrt(2.0)))));
    //if((std::isnan(cfd_min[5])) || (cfd_min[5] < 0.0)){ cfd_min[5] = 0.0; }
    //if(cfd_min[5] > 1.0){ cfd_min[5] = 1.0; }
    cfd_min[5] = 0.0;
    cfd_max[5] = 1.0;
  } else {
    cfd_min[5] = cfd_max[5] = 0.5;
  }

  if(BF_errs[6] > 0.0){
    cfd_min[6] = 0.5 * (1.0 + erf((double) ((x_min - BF_fit[12]) / (amp * BF_errs[6] * temp_cos[2] * sqrt(2.0)))));
    if(std::isnan(cfd_min[6]) || (cfd_min[6] < 0.0)){ cfd_min[6] = 0.0; }
    if(cfd_min[6] > 1.0){ cfd_min[6] = 1.0; }
    cfd_max[6] = 0.5 * (1.0 + erf((double) ((x_max - BF_fit[12]) / (amp * BF_errs[6] * temp_cos[2] * sqrt(2.0)))));
    if(std::isnan(cfd_max[6]) || (cfd_max[6] < 0.0)){ cfd_max[6] = 0.0; }
    if(cfd_max[6] > 1.0){ cfd_max[6] = 1.0; }
  } else {
    cfd_min[6] = cfd_max[6] = 0.5;
  }

  if(BF_errs[7] > 0.0){
    //cfd_min[7] = 0.5 * (1.0 + erf((double) ((2.0 - BF_fit[14]) / (3.0 * BF_errs[7] * temp_cos[3] * sqrt(2.0)))));
    cfd_min[7] = 0.5 * (1.0 + erf((double) ((1.0 - BF_fit[14]) / (3.5 * BF_errs[7] * temp_cos[3] * sqrt(2.0)))));
    if(std::isnan(cfd_min[7]) || (cfd_min[7] < 0.0)){ cfd_min[7] = 0.0; }
    if(cfd_min[7] > 1.0){ cfd_min[7] = 1.0; }
    //cfd_max[7] = 0.5 * (1.0 + erf((double) ((8.0 - BF_fit[14]) / (3.0 * BF_errs[7] * temp_cos[3] * sqrt(2.0)))));
    cfd_max[7] = 0.5 * (1.0 + erf((double) ((8.0 - BF_fit[14]) / (3.5 * BF_errs[7] * temp_cos[3] * sqrt(2.0)))));
    if(std::isnan(cfd_max[7]) || (cfd_max[7] < 0.0)){ cfd_max[7] = 0.0; }
    if(cfd_max[7] > 1.0){ cfd_max[7] = 1.0; }
  } else {
    cfd_min[7] = cfd_max[7] = 0.5;
  }

  // initialise random number generation
  srand (time(NULL));
  seed = (long) (-1*(rand()));
  ran2(seed);

  // generate a correlation matrix from the covariance matrix --- set correlation to zero when they
  // involve parameters with uncertainties of zero, because they are most likely fixed values
  BF_corr = new double * [8];
  for(i = 0; i < 8; i++){ BF_corr[i] = new double[8]; }
  for(r = 0; r < 8; r++){
    for(i = 0; i < 8; i++){
      if(i == r){
	BF_corr[r][i] = 1.0;
      } else {
	if((BF_covar[r][r] > 0.0) && (BF_covar[i][i] > 0.0)){
	  BF_corr[r][i] = (double) BF_covar[r][i] / (BF_errs[i] * BF_errs[r]);
	  if(std::isnan(BF_corr[r][i])){ BF_corr[r][i] = 0.0; }
	  if(std::isinf(BF_corr[r][i])){ BF_corr[r][i] = (BF_corr[r][i] > 0.0) ? 9E30 : -9E30; }
	} else { BF_corr[r][i] = 0.0; }
      }
    }
  }
  
  // adjust the correlation matrix for the fit type
  switch(fit_type){
  case 1:
    for(i = 0; i < 8; i++){
      if(i != 3){
	BF_corr[3][i] = BF_corr[1][i];
	BF_corr[i][3] = BF_corr[i][1];
      } 
    }
    for(i = 0; i < 8; i++){
      if(i != 5){
	BF_corr[5][i] = BF_corr[i][5] = 0.0;
      }
    }
    for(i = 0; i < 8; i++){
      if(i != 6){
	BF_corr[6][i] = 0.5 * (BF_corr[2][i] + BF_corr[4][i]);
	BF_corr[i][6] = 0.5 * (BF_corr[i][2] + BF_corr[i][4]);
      }
    }
    for(i = 0; i < 7; i++){ BF_corr[7][i] = BF_corr[i][7] = 0.0; }    
    break;
  case 2:
    for(i = 0; i < 8; i++){
      if(i != 5){
	BF_corr[5][i] = BF_corr[i][5] = 0.0;
      }
    }
    for(i = 0; i < 8; i++){
      if(i != 6){
	BF_corr[6][i] = 0.5 * (BF_corr[2][i] + BF_corr[4][i]);
	BF_corr[i][6] = 0.5 * (BF_corr[i][2] + BF_corr[i][4]);
      }
    }
    for(i = 0; i < 7; i++){ BF_corr[7][i] = BF_corr[i][7] = 0.0; }    
    break;
  case 3:
    for(i = 0; i < 8; i++){
      if(i != 3){
	BF_corr[3][i] = BF_corr[1][i];
	BF_corr[i][3] = BF_corr[i][1];
      } 
    }
    for(i = 0; i < 8; i++){
      if(i != 6){
	BF_corr[6][i] = 0.5 * (BF_corr[2][i] + BF_corr[4][i]);
	BF_corr[i][6] = 0.5 * (BF_corr[i][2] + BF_corr[i][4]);
      }
    }
    for(i = 0; i < 7; i++){ BF_corr[7][i] = BF_corr[i][7] = 0.0; }    
    break;
  case 4:
    for(i = 0; i < 8; i++){
      if(i != 3){
	BF_corr[3][i] = BF_corr[1][i];
	BF_corr[i][3] = BF_corr[i][1];
      } 
    }
    for(i = 0; i < 8; i++){
      if(i != 6){
	BF_corr[6][i] = 0.5 * (BF_corr[2][i] + BF_corr[4][i]);
	BF_corr[i][6] = 0.5 * (BF_corr[i][2] + BF_corr[i][4]);
      }
    }
    break;
  case 5:
    for(i = 0; i < 8; i++){
      if(i != 3){
	BF_corr[3][i] = BF_corr[1][i];
	BF_corr[i][3] = BF_corr[i][1];
      } 
    }
    break;
  default:
    break;
  }
  
  // calculate the Cholesky decomposition of the correlation matrix
  choldc(BF_corr,8);

  // generate NOr correlated, random variations of the Busy Function fit 
  if(vb_flag > 0){

    std::cout << "0         25        50        75        100%" << std::endl;
    std::cout << "| | | | | : | | | | : | | | | : | | | | :" << std::endl;

  }
  r_done = 0;
  progress = 0.0;
  for(r = 0; r < NOr; r++){
    
    // generate a random array of 8 values following a Gaussian distribution with zero mean and standard deviation of 1 ---
    // use binary search and Donald Knuth's cfd. mapping trick
    for(i = 0; i < 8; i++){
      
      // it's crucial that the same seed value isn't shared between threads, otherwise the threads will use
      // identical sequences of pseudo-random numbers
      rand_val = ran2(seed);

      // convert random uniform deviate value into a z-value for a gaussian distribution
      // with a mean of 0 and a standard deviation of 1
      rand_val = cfd_min[i] + ((cfd_max[i] - cfd_min[i])*rand_val);
      z_val = inverf(((2.0 * (double) rand_val) - 1.0)) * sqrt(2.0);
      if(std::isnan(z_val)){ z_val = 0.0; }
      if(std::isinf(z_val)){ z_val = (z_val > 0.0) ? 4.5 : -4.5; }
      rand_set[i] = z_val;

      // for(i = 0; i < 8; i++)
    }
      
    // apply the fit type constraints
    switch(fit_type){
    case 1:
      rand_set[3] = rand_set[1];
      rand_set[5] = 0.0;
      rand_set[6] = 0.5 * (rand_set[2] + rand_set[4]);
      rand_set[7] = 4.0;
      break;
    case 2:
      rand_set[5] = 0.0;
      rand_set[6] = 0.5 * (rand_set[2] + rand_set[4]);
      rand_set[7] = 4.0;
      break;
    case 3:
      rand_set[3] = rand_set[1];
      rand_set[6] = 0.5 * (rand_set[2] + rand_set[4]);
      rand_set[7] = 4.0;
      break;
    case 4:
      rand_set[3] = rand_set[1];
      rand_set[6] = 0.5 * (rand_set[2] + rand_set[4]);
      break;
    case 5:
      rand_set[3] = rand_set[1];
      break;
    default:
      break;
    }
  
    // multiply by the correlation matrix Cholesky decomposition
    for(j = 0; j < 8; j++){
      rand_val = 0.0;
      for(i = 0; i < 8; i++){
	rand_val+=(BF_corr[j][i]*rand_set[i]); 
      }
      rand_fits[r][j] = (T_obs) rand_val;
    }
    
    // apply random offsets to BF fits to generate random BF variants --- incorporate
    // variable mappings at the same time
    rand_fits[r][0] = (T_obs) (BF_fit[0] * (1.0 + (rand_fits[r][0] * BF_errs[0])));
    if(std::isnan((double) rand_fits[r][0])){ rand_fits[r][0] = 0.0; }
    if(std::isinf((double) rand_fits[r][0])){ rand_fits[r][0] = (rand_fits[r][0] > 0.0) ? 9E30 : -9E30; }
    if(rand_fits[r][0] < 0.0){ rand_fits[r][0] = 0.0; }

    rand_fits[r][1] = (T_obs) (BF_fit[2] * (1.0 + (rand_fits[r][1] * BF_errs[1])));
    if(std::isnan((double) rand_fits[r][1])){ rand_fits[r][1] = 0.0; }
    if(std::isinf((double) rand_fits[r][1])){ rand_fits[r][1] = (rand_fits[r][1] > 0.0) ? 9E30 : -9E30; }
    if(rand_fits[r][1] < 0.0){ rand_fits[r][1] = 0.0; }

    rand_fits[r][2] = (T_obs) (BF_fit[4] + (amp * temp_cos[0] * rand_fits[r][2] * BF_errs[2]));
    if(std::isnan((double) rand_fits[r][2])){ rand_fits[r][2] = 0.0; }
    if(std::isinf((double) rand_fits[r][2])){ rand_fits[r][2] = (rand_fits[r][2] > 0.0) ? 9E30 : -9E30; }
    if(rand_fits[r][2] < x_min){ rand_fits[r][2] = x_min; }
    if(rand_fits[r][2] > x_max){ rand_fits[r][2] = x_max; }

    rand_fits[r][3] = (T_obs) (BF_fit[6] * (1.0 + (rand_fits[r][3] * BF_errs[3])));
    if(std::isnan((double) rand_fits[r][3])){ rand_fits[r][3] = 0.0; }
    if(std::isinf((double) rand_fits[r][3])){ rand_fits[r][3] = (rand_fits[r][3] > 0.0) ? 9E30 : -9E30; }
    if(rand_fits[r][3] < 0.0){ rand_fits[r][3] = 0.0; }

    rand_fits[r][4] = (T_obs) (BF_fit[8] + (amp * temp_cos[1] * rand_fits[r][4] * BF_errs[4]));
    if(std::isnan((double) rand_fits[r][4])){ rand_fits[r][4] = 0.0; }
    if(std::isinf((double) rand_fits[r][4])){ rand_fits[r][4] = (rand_fits[r][4] > 0.0) ? 9E30 : -9E30; }
    if(rand_fits[r][4] < x_min){ rand_fits[r][4] = x_min; }
    if(rand_fits[r][4] > x_max){ rand_fits[r][4] = x_max; }    

    rand_fits[r][5] = (T_obs) (BF_fit[10] * (1.0 + (rand_fits[r][5] * BF_errs[5])));
    if(std::isnan((double) rand_fits[r][5])){ rand_fits[r][5] = 0.0; }
    if(std::isinf((double) rand_fits[r][5])){ rand_fits[r][5] = (rand_fits[r][5] > 0.0) ? 9E30 : -9E30; }
    //if(rand_fits[r][5] < 0.0){ rand_fits[r][5] = 0.0; }

    rand_fits[r][6] = (T_obs) (BF_fit[12] + (amp * temp_cos[2] * rand_fits[r][6] * BF_errs[6]));
    if(std::isnan((double) rand_fits[r][6])){ rand_fits[r][6] = 0.0; }
    if(std::isinf((double) rand_fits[r][6])){ rand_fits[r][6] = (rand_fits[r][6] > 0.0) ? 9E30 : -9E30; }
    if(rand_fits[r][6] < x_min){ rand_fits[r][6] = x_min; }
    if(rand_fits[r][6] > x_max){ rand_fits[r][6] = x_max; }

    //rand_fits[r][7] = (T_obs) (BF_fit[14] + (3.0 * temp_cos[3] * rand_fits[r][7] * BF_errs[7]));    
    rand_fits[r][7] = (T_obs) (BF_fit[14] + (3.5 * temp_cos[3] * rand_fits[r][7] * BF_errs[7]));    
    if(std::isnan((double) rand_fits[r][7])){ rand_fits[r][7] = 0.0; }
    if(std::isinf((double) rand_fits[r][7])){ rand_fits[r][7] = (rand_fits[r][7] > 0.0) ? 9E30 : -9E30; }
    //if(rand_fits[r][7] < 2.0){ rand_fits[r][7] = 2.0; }
    if(rand_fits[r][7] < 1.0){ rand_fits[r][7] = 1.0; }
    if(rand_fits[r][7] > 8.0){ rand_fits[r][7] = 8.0; }

    // apply constraints of BF fit type
    switch(fit_type){
    case 1:
      rand_fits[r][3] = rand_fits[r][1];
      rand_fits[r][5] = 0.0;
      rand_fits[r][6] = 0.5 * (rand_fits[r][2] + rand_fits[r][4]);
      rand_fits[r][7] = 4.0;
      break;
    case 2:
      rand_fits[r][5] = 0.0;
      rand_fits[r][6] = 0.5 * (rand_fits[r][2] + rand_fits[r][4]);
      rand_fits[r][7] = 4.0;
      break;
    case 3:
      rand_fits[r][3] = rand_fits[r][1];
      rand_fits[r][6] = 0.5 * (rand_fits[r][2] + rand_fits[r][4]);
      rand_fits[r][7] = 4.0;
      break;
    case 4:
      rand_fits[r][3] = rand_fits[r][1];
      rand_fits[r][6] = 0.5 * (rand_fits[r][2] + rand_fits[r][4]);
      break;
    case 5:
      rand_fits[r][3] = rand_fits[r][1];
      break;
    default:
      break;
    }

    // display progress
    if(vb_flag > 0){
      while((((float) (r_done + 1)) / ((float) NOr)) >= progress){ std::cout << "*" << std::flush; progress+=0.025; }
      r_done++;
    }    

    // for(r = 0; r < NOr; r++)
  }
  if(vb_flag > 0){ std::cout << std::endl; }

  // free up memory
  for(i = 0; i < 8; i++){ delete [] BF_corr[i]; }
  delete [] BF_corr;

}

template <class T_x, class T_obs, class T_calc>
  void CalcObsParams(int NOr, T_obs ** rand_fits, int NOvals, T_x * x_vals, T_calc ** obs_vals, int vb_flag){

  int r,i,j,k,r_done;
  float progress;
  double erf_1,erf_2,plaw,BF_val,w20_min,w20_max,w50_min,w50_max, ** y_vals;

  // allocate memory
  y_vals = new double * [NOr];
  for(r = 0; r < NOr; r++){ y_vals[r] = new double[NOvals]; }

  // for each random variation (correctly correlated) of a Busy Function fit, calculate the total flux,
  // peak flux, w_50 and w_20 values (conventional velocity definitions ala Duchamp)
  if(vb_flag >= 0){   std::cout << "Calculating observational properties of BF fits . . . " << std::endl;  }
  if(vb_flag > 0){

    std::cout << "0         25        50        75        100%" << std::endl;
    std::cout << "| | | | | : | | | | : | | | | : | | | | :" << std::endl;

  }
  r_done = 0;
  progress = 0.0;
  for(r = 0; r < NOr; r++){
    
    // initialise flux values
    obs_vals[r][0] = 0.0; // total flux
    obs_vals[r][1] = 0.0; // peak flux height
    obs_vals[r][2] = (T_calc) x_vals[0]; // peak flux position
   
    // calculate total flux and find peak flux --- loop over all x values
    for(i = 0; i < NOvals; i++){

      // calculate Busy Function value for this x value
      erf_1 = 1.0 + erf((double) (rand_fits[r][1] * (x_vals[i] - rand_fits[r][2])));
      if(std::isnan(erf_1)){ erf_1 = 0.0; }
      if(std::isinf(erf_1)){ erf_1 = (erf_1 > 0.0) ? 1.0 : 0.0; }
      erf_2 = 1.0 + erf((double) (rand_fits[r][3] * (rand_fits[r][4] - x_vals[i])));
      if(std::isnan(erf_2)){ erf_2 = 0.0; }
      if(std::isinf(erf_2)){ erf_2 = (erf_2 > 0.0) ? 1.0 : 0.0; }
      plaw = 1.0 + (rand_fits[r][5] * pow(abs((rand_fits[r][6] - x_vals[i])),rand_fits[r][7]));
      if(std::isnan(plaw)){ plaw = 0.0; } 
      if(std::isinf(plaw)){ plaw = (plaw > 0.0) ? 9E30 : -9E30; }
      BF_val = 0.25 * rand_fits[r][0] * erf_1 * erf_2 * plaw;
      if(std::isnan(BF_val)){ BF_val = 0.0; } 
      if(std::isinf(BF_val)){ BF_val = (BF_val > 0.0) ? 9E30 : -9E30; }

      // store the BF val for this x
      y_vals[r][i] = BF_val;

      // update total flux and peak flux
      obs_vals[r][0]+=((T_calc) BF_val);
      if(BF_val >= (double) obs_vals[r][1]){ 

	obs_vals[r][1] = (T_calc) BF_val;
	obs_vals[r][2] = (T_calc) x_vals[i];
	
      }
	
    }
   
    // calculate w_20 and w_50 in conventional manner --- loop over all input x values
    w20_min = 0.0;
    j = 0;
    for(i = 0; i < NOvals; i++){ 
      
      if((double) y_vals[r][i] >= (0.2 * (double) obs_vals[r][1])){ 
	
	if(i > 0){
	  
	  w20_min = (double) x_vals[(i - 1)] + ((double) (x_vals[i] - x_vals[(i - 1)]) * ((0.2 * (double) obs_vals[r][1]) - (double) y_vals[r][(i - 1)]) / (double) (y_vals[r][i] - y_vals[r][(i - 1)])); 
	  
	} else {
	  
	  w20_min = (double) x_vals[i]; 
	  
	}
	
	break; 
	
      } else {

	j = i;

      }
      
    }
    w20_max = (double) x_vals[(NOvals - 1)];
    k = NOvals - 1;
    for(i = NOvals - 1; i >= 0; i--){ 
	
      if(y_vals[r][i] >= (0.2 * (double) obs_vals[r][1])){ 
	
	if(i < (NOvals - 1)){
	  
	  w20_max = (double) x_vals[(i + 1)] - ((double) (x_vals[(i + 1)] - x_vals[i]) * ((0.2 * (double) obs_vals[r][1]) - (double) y_vals[r][(i + 1)]) / (double) (y_vals[r][i] - y_vals[r][(i + 1)]));
	  
	} else {
	  
	  w20_max = (double) x_vals[i]; 
	  
	}
	
	break; 
	
      } else {

	k = i;
	
      }
      
    }
    w50_min = 0.0;
    for(i = j; i < NOvals; i++){ 
      
      if(y_vals[r][i] >= (0.5 * (double) obs_vals[r][1])){ 
	
	if(i > 0){
	  
	  w50_min = (double) x_vals[(i - 1)] + ((double) (x_vals[(i + 1)] - x_vals[i]) * ((0.5 * (double) obs_vals[r][1]) - (double) y_vals[r][(i - 1)]) / (double) (y_vals[r][i] - y_vals[r][(i - 1)])); 
	  
	} else {
	  
	  w50_min = (double) x_vals[i]; 
	  
	}	
	
	break; 
	
      } 
      
    }
    w50_max = (double) x_vals[(NOvals - 1)];
    for(i = k; i >= 0; i--){ 
      
      if(y_vals[r][i] >= (double) (0.5 * obs_vals[r][1])){ 
	
	if(i < (NOvals - 1)){
	  
	  w50_max = (double) x_vals[(i + 1)] - ((double) (x_vals[(i + 1)] - x_vals[i]) * ((0.5 * (double) obs_vals[r][1]) - (double) y_vals[r][(i + 1)]) / (double) (y_vals[r][i] - y_vals[r][(i + 1)]));
	  
	} else {
	  
	  w50_max = (double) x_vals[i]; 
	  
	}
	
	break; 
	
      } 
      
    }
    obs_vals[r][3] = (T_calc) (w50_max - w50_min);
    obs_vals[r][4] = (T_calc) (0.5 * (w50_min + w50_max));
    obs_vals[r][5] = (T_calc) (w20_max - w20_min);
    obs_vals[r][6] = (T_calc) (0.5 * (w20_min + w20_max));

    // range test values
    for(i = 0; i < 7; i++){
      if(std::isnan(obs_vals[r][i])){ obs_vals[r][i] = 0.0; }
      if(std::isinf(obs_vals[r][i])){ obs_vals[r][i] = (obs_vals[r][i] > 0.0) ? 9E30 : -9E30; }
    }

    // display progress
    if(vb_flag > 0){
      while((((float) (r_done + 1)) / ((float) NOr)) >= progress){ std::cout << "*" << std::flush; progress+=0.025; }
      r_done++;
    }    

    // for(r = 0; r < NOr; r++)
  }
  if(vb_flag > 0){ std::cout << std::endl; }

  // free up memory
  for(r = 0; r < NOr; r++){ delete [] y_vals[r]; }
  delete [] y_vals;
  
}

template <class T_obs>
void heapsort(int n, T_obs ra[]){

  // This function is called by analyse_obs_vals. It sorts an array of n values into ascending
  // order using the heapsort algorithm. Note that it uses base-1 array indexing. This is 
  // most easily dealt with by passing the input array address minus 1 e.g. ra - 1.

  int i,ir,j,l;
  T_obs rra;

  if (n<2) return;
  l = (n >> 1) + 1;
  ir = n;
  for(;;){

    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1]=rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	i=j;
	j <<= 1;
      } else break;
    }
    ra[i]=rra;

  }

}

template <class T_obs>
void AnalyseObsVals(int NOr, int NOprops, T_obs ** obs_vals, T_obs ** obs_stats, int vb_flag){

  T_obs ** temp_vals;
  int i,j,p;
  float progress;

  // calculate the median and standard deviation (from IQR) of each property 
  //
  // for a given property, the 0th, 1st, 2nd, 3rd, 4th, 5th and 6th array elements are the median, std.dev. (from IQR), minimum, maximum,
  //                       mean, std. dev., skewness and kurtosis

  // 0. Exit if NOr <= 1. The sample isn't big enough to process.
  if(NOr <= 1){

    std::cout << "WARNING: NOr is too small. Exiting AnalyseObsVals function." << std::endl;
    return;

  }

  // 1. Create array to store obs_vals in opposite row-column order. This improves memory access patterns for the rest
  // of the function.
  if(vb_flag >= 0){ std::cout << "Creating temporary working arrays . . . " << std::endl; }
  temp_vals = new T_obs * [NOprops];
  for(p = 0; p < NOprops; p++){ 

    temp_vals[p] = new T_obs[NOr]; 
    for(i = 0; i < NOr; i++){ temp_vals[p][i] = obs_vals[i][p]; }

  }

  // 2. Process each property. Sort the values, then calculate the min, max, median and IQR --> convert to standard deviation.
  if(vb_flag >= 0){ std::cout << "Analysing observational properties of BF fits . . . " << std::endl; }
  if(vb_flag > 0){

    std::cout << "0         25        50        75        100%" << std::endl;
    std::cout << "| | | | | : | | | | : | | | | : | | | | :" << std::endl;
  
  }
  i = 0;
  progress = 0.0;
  for(p = 0; p < NOprops; p++){

    // sort values
    heapsort(NOr,temp_vals[p] - 1);

    // calculate median
    if((NOr % 2) == 0){

      obs_stats[p][0] = 0.5 * (temp_vals[p][(int) floorf((float) NOr / 2.0)] + temp_vals[p][(int) floorf(((float) NOr / 2.0) - 1.0)]);

    } else {

      obs_stats[p][0] = temp_vals[p][(int) floorf(((float) NOr - 1.0)/2.0)];

    }

    // calculate standard deviation from IQR
    obs_stats[p][1] = (temp_vals[p][(int)(floorf((0.5 + (0.75 * (float) NOr))))] - temp_vals[p][(int)(floorf((0.5 + (0.25 * (float) NOr))))]) / 1.349;

    // assign minimum and maximum
    obs_stats[p][2] = temp_vals[p][0];
    obs_stats[p][3] = temp_vals[p][(NOr - 1)];

    // calcuate mean and conventional standard deviation
    obs_stats[p][4] = obs_stats[p][5] = 0.0;
    for(j = 0; j < NOr; j++){

      obs_stats[p][4]+=temp_vals[p][j];
      obs_stats[p][5]+=(temp_vals[p][j]*temp_vals[p][j]);

    }
    obs_stats[p][4]/=((T_obs) NOr);
    obs_stats[p][5] = (T_obs) sqrt((double) ((obs_stats[p][5] / ((T_obs) NOr)) - (obs_stats[p][4] * obs_stats[p][4])));

    // calculate skewness and kurtosis
    obs_stats[p][6] = obs_stats[p][7] = 0.0;
    for(j = 0; j < NOr; j++){

      obs_stats[p][6]+=((temp_vals[p][j] - obs_stats[p][4])*(temp_vals[p][j] - obs_stats[p][4])*(temp_vals[p][j] - obs_stats[p][4]));
      obs_stats[p][7]+=((temp_vals[p][j] - obs_stats[p][4])*(temp_vals[p][j] - obs_stats[p][4])*(temp_vals[p][j] - obs_stats[p][4])*(temp_vals[p][j] - obs_stats[p][4]));

    }
    obs_stats[p][6]/=(obs_stats[p][5]*obs_stats[p][5]*obs_stats[p][5] * ((T_obs) NOr));
    obs_stats[p][7]/=(obs_stats[p][5]*obs_stats[p][5]*obs_stats[p][5]*obs_stats[p][5] * ((T_obs) NOr));

    // range test stats
    for(j = 0; j < 8; j++){
      if(std::isnan(obs_stats[p][i])){ obs_stats[p][i] = 0.0; }
      if(std::isinf(obs_stats[p][i])){ obs_stats[p][i] = (obs_stats[p][i] > 0.0) ? 9E30 : -9E30; }
    }

    // display progress
    if(vb_flag > 0){
      while((((float) (i + 1)) / ((float) NOprops)) >= progress){ std::cout << "*" << std::flush; progress+=0.025; }
      i++;
    }    

    // for(p = 0; p < NOprops; p++)
  }
  if(vb_flag > 0){ std::cout << std::endl; }

  // 3. Free up memory assigned to temp_vals.
  for(p = 0; p < NOprops; p++){ delete [] temp_vals[p]; }
  delete [] temp_vals;

}

// function to calculate a linear approximation of the observable parameter covariance matrix
template <class T_x, class T_obs, class T_approx, class T_range>
  void ApproxObsCovar(int fit_type, int NOvals, T_x * x_vals, T_obs * fit_params, T_obs ** fit_covar, T_approx ** obs_covar, T_range x_min, T_range x_max){
  
  double ** model_vals, ** obs_vals,temp_matrix[7][8],scaling[8],mid,amp,jacobian[8][7], use_fit_covar[8][8];
  int i,j,k;

  // set mid and amp
  mid = 0.5 * (x_min + x_max);
  amp = 0.5 * (x_max - x_min);

  // create arrays
  obs_vals = new double * [9];
  model_vals = new double * [9];
  for(i = 0; i < 9; i++){ 
    obs_vals[i] = new double[7]; 
    model_vals[i] = new double[8];
  }

  // convert BF fit covar from internal, mapped variables to external variables
  for(j = 0; j < 8; j++){
    for(i = 0; i < 8; i++){
      use_fit_covar[j][i] = fit_covar[j][i];
    }
  }
  scaling[0] = fit_params[0];
  scaling[1] = fit_params[2];
  scaling[2] = amp * cos(asin(((double) fit_params[4] - mid) / amp));
  scaling[3] = fit_params[6];
  scaling[4] = amp * cos(asin(((double) fit_params[8] - mid) / amp));
  scaling[5] = fit_params[10];
  scaling[6] = amp * cos(asin(((double) fit_params[12] - mid) / amp));
  //scaling[7] = amp * cos(asin(((double) fit_params[14] - mid) / amp));
  scaling[7] = 3.5 * cos(asin(((double) fit_params[14] - 4.5) / 3.5));
  for(i = 0; i < 8; i++){
    if(std::isnan(scaling[i])){ scaling[i] = 0.0; }
    if(std::isinf(scaling[i])){ scaling[i] = (scaling[i] > 0.0) ? 9E30 : -9E30; }
  }
  for(j = 0; j < 8; j++){
    for(i = 0; i < 8; i++){
    use_fit_covar[j][i]*=(scaling[j] * scaling[i]);
    }
  }

  // generate 9 BF fits --- the first is the original, the remaining 8 have an offset applied to 
  // one of the parameters
  for(i = 0; i < 8; i++){ model_vals[0][i] = (double) fit_params[(2 * i)]; }  
  switch(fit_type){
  case 1:
    for(j = 1; j < 9; j++){
      for(i = 0; i < 8; i++){ model_vals[j][i] = (double) fit_params[(2 * i)]; }
    }
    model_vals[1][0]+=(1.0E-5 * fabs((double) fit_params[0]));
    model_vals[2][1]+=(1.0E-5 * fabs((double) fit_params[2]));
    model_vals[2][3] = model_vals[2][1];
    model_vals[3][2]+=(1.0E-5 * fabs((double) fit_params[4]));
    model_vals[4][3]+=(1.0E-5 * fabs((double) fit_params[6]));
    model_vals[4][1] = model_vals[4][3];
    model_vals[5][4]+=(1.0E-5 * fabs((double) fit_params[8]));
    for(j = 0; j < 9; j++){
      for(i = 0; i < 8; i++){ 
	if(std::isnan(model_vals[j][i])){ model_vals[j][i] = 0.0; }
	if(std::isinf(model_vals[j][i])){ model_vals[j][i] = (model_vals[j][i] > 0.0) ? 9E30 : -9E30; }
      }    
    }
    break;
  case 2:
    for(j = 1; j < 9; j++){
      for(i = 0; i < 8; i++){ model_vals[j][i] = (double) fit_params[(2 * i)]; }
    }
    model_vals[1][0]+=(1.0E-5 * fabs((double) fit_params[0]));
    model_vals[2][1]+=(1.0E-5 * fabs((double) fit_params[2]));
    model_vals[3][2]+=(1.0E-5 * fabs((double) fit_params[4]));
    model_vals[4][3]+=(1.0E-5 * fabs((double) fit_params[6]));
    model_vals[5][4]+=(1.0E-5 * fabs((double) fit_params[8]));
    for(j = 0; j < 9; j++){
      for(i = 0; i < 8; i++){ 
	if(std::isnan(model_vals[j][i])){ model_vals[j][i] = 0.0; }
	if(std::isinf(model_vals[j][i])){ model_vals[j][i] = (model_vals[j][i] > 0.0) ? 9E30 : -9E30; }
      }    
    }
    break;
  case 3:
    for(j = 1; j < 9; j++){
      for(i = 0; i < 8; i++){ model_vals[j][i] = (double) fit_params[(2 * i)]; }
    }
    model_vals[1][0]+=(1.0E-5 * fabs((double) fit_params[0]));
    model_vals[2][1]+=(1.0E-5 * fabs((double) fit_params[2]));
    model_vals[2][3] = model_vals[2][1];
    model_vals[3][2]+=(1.0E-5 * fabs((double) fit_params[4]));
    model_vals[3][6]+=(0.5 * 1.0E-5 * fabs((double) fit_params[4]));
    model_vals[4][3]+=(1.0E-5 * fabs((double) fit_params[6]));
    model_vals[4][1] = model_vals[4][3];
    model_vals[5][4]+=(1.0E-5 * fabs((double) fit_params[8]));
    model_vals[5][6]+=(0.5 * 1.0E-5 * fabs((double) fit_params[8]));
    model_vals[6][5]+=(1.0E-5 * fabs((double) fit_params[10]));
    for(j = 0; j < 9; j++){
      for(i = 0; i < 8; i++){ 
	if(std::isnan(model_vals[j][i])){ model_vals[j][i] = 0.0; }
	if(std::isinf(model_vals[j][i])){ model_vals[j][i] = (model_vals[j][i] > 0.0) ? 9E30 : -9E30; }
      }    
    }
    break;
  case 4:
    for(j = 1; j < 9; j++){
      for(i = 0; i < 8; i++){ model_vals[j][i] = (double) fit_params[(2 * i)]; }
    }
    model_vals[1][0]+=(1.0E-5 * fabs((double) fit_params[0]));
    model_vals[2][1]+=(1.0E-5 * fabs((double) fit_params[2]));
    model_vals[2][3] = model_vals[2][1];
    model_vals[3][2]+=(1.0E-5 * fabs((double) fit_params[4]));
    model_vals[3][6]+=(0.5 * 1.0E-5 * fabs((double) fit_params[4]));
    model_vals[4][3]+=(1.0E-5 * fabs((double) fit_params[6]));
    model_vals[4][1] = model_vals[4][3];
    model_vals[5][4]+=(1.0E-5 * fabs((double) fit_params[8]));
    model_vals[5][6]+=(0.5 * 1.0E-5 * fabs((double) fit_params[8]));
    model_vals[6][5]+=(1.0E-5 * fabs((double) fit_params[10]));
    model_vals[8][7]+=(1.0E-5 * fabs((double) fit_params[14]));
    for(j = 0; j < 9; j++){
      for(i = 0; i < 8; i++){ 
	if(std::isnan(model_vals[j][i])){ model_vals[j][i] = 0.0; }
	if(std::isinf(model_vals[j][i])){ model_vals[j][i] = (model_vals[j][i] > 0.0) ? 9E30 : -9E30; }
      }    
    }
    break;
  case 5:
    for(j = 1; j < 9; j++){
      for(i = 0; i < 8; i++){ model_vals[j][i] = (double) fit_params[(2 * i)]; }
    }
    model_vals[1][0]+=(1.0E-5 * fabs((double) fit_params[0]));
    model_vals[2][1]+=(1.0E-5 * fabs((double) fit_params[2]));
    model_vals[3][2]+=(1.0E-5 * fabs((double) fit_params[4]));
    model_vals[2][3] = model_vals[2][1];
    model_vals[4][3]+=(1.0E-5 * fabs((double) fit_params[6]));
    model_vals[4][1] = model_vals[4][3];
    model_vals[5][4]+=(1.0E-5 * fabs((double) fit_params[8]));
    model_vals[6][5]+=(1.0E-5 * fabs((double) fit_params[10]));
    model_vals[7][6]+=(1.0E-5 * fabs((double) fit_params[12]));
    model_vals[8][7]+=(1.0E-5 * fabs((double) fit_params[14]));
    for(j = 0; j < 9; j++){
      for(i = 0; i < 8; i++){ 
	if(std::isnan(model_vals[j][i])){ model_vals[j][i] = 0.0; }
	if(std::isinf(model_vals[j][i])){ model_vals[j][i] = (model_vals[j][i] > 0.0) ? 9E30 : -9E30; }
      }    
    }
    break;
  default:
    for(j = 1; j < 9; j++){
      for(i = 0; i < 8; i++){ model_vals[j][i] = (double) fit_params[(2 * i)]; }
      model_vals[j][(j - 1)]+=(1.0E-5 * fabs((double) fit_params[(2 * (j - 1))]));
      if(std::isnan(model_vals[j][(j - 1)])){ model_vals[j][(j - 1)] = 0.0; }
      if(std::isinf(model_vals[j][(j - 1)])){ model_vals[j][(j - 1)] = (model_vals[j][(j - 1)] > 0.0) ? 9E30 : -9E30; }
    }    
    for(j = 0; j < 9; j++){
      for(i = 0; i < 8; i++){ 
	if(std::isnan(model_vals[j][i])){ model_vals[j][i] = 0.0; }
	if(std::isinf(model_vals[j][i])){ model_vals[j][i] = (model_vals[j][i] > 0.0) ? 9E30 : -9E30; }
      }    
    }
    break;
  }

  // calculate observational parameters for 9 BF fits
  CalcObsParams(9,model_vals,NOvals,x_vals,obs_vals,-1);

  // calculate Jacobian matrix
  for(j = 0; j < 8; j++){
    for(i = 0; i < 7; i++){
      jacobian[j][i] = (obs_vals[(j + 1)][i] - obs_vals[0][i])/(1.0E-5 * fabs((double) fit_params[(2 * j)]));
      if(std::isnan(jacobian[j][i])){ jacobian[j][i] = 0.0; }
      if(std::isinf(jacobian[j][i])){ jacobian[j][i] = (jacobian[j][i] > 0.0) ? 9E30 : -9E30; }
    }
  }

  // apply BF fit type constraints
  switch(fit_type){
  case 1:
    for(i = 0; i < 7; i++){ jacobian[3][i] = jacobian[1][i]; }
    for(i = 0; i < 7; i++){ jacobian[5][i] = 0.0; }
    for(i = 0; i < 7; i++){ jacobian[6][i] = 0.5 * (jacobian[2][i] + jacobian[4][i]); }
    for(i = 0; i < 7; i++){ jacobian[7][i] = 0.0; }
    break;
  case 2:
    for(i = 0; i < 7; i++){ jacobian[5][i] = 0.0; }
    for(i = 0; i < 7; i++){ jacobian[6][i] = 0.5 * (jacobian[2][i] + jacobian[4][i]); }
    for(i = 0; i < 7; i++){ jacobian[7][i] = 0.0; }
    break;
  case 3:
    for(i = 0; i < 7; i++){ jacobian[3][i] = jacobian[1][i]; }
    for(i = 0; i < 7; i++){ jacobian[6][i] = 0.5 * (jacobian[2][i] + jacobian[4][i]); }
    for(i = 0; i < 7; i++){ jacobian[7][i] = 0.0; }
    break;
  case 4:
    for(i = 0; i < 7; i++){ jacobian[3][i] = jacobian[1][i]; }
    for(i = 0; i < 7; i++){ jacobian[6][i] = 0.5 * (jacobian[2][i] + jacobian[4][i]); }
    break;
  case 5:
    for(i = 0; i < 7; i++){ jacobian[3][i] = jacobian[1][i]; }
    break;
  default:
    break;
  }
  
  // create obs_covar approximation: jacobian transpose * use_fit_covar * jacobian
  for(j = 0; j < 7; j++){
    for(i = 0; i < 8; i++){
      temp_matrix[j][i] = 0.0;
      for(k = 0; k < 8; k++){ temp_matrix[j][i]+=(jacobian[k][j] * use_fit_covar[k][i]); }
      if(std::isnan(temp_matrix[j][i])){ temp_matrix[j][i] = 0.0; }
      if(std::isinf(temp_matrix[j][i])){ temp_matrix[j][i] = (temp_matrix[j][i] > 0.0) ? 9E30 : -9E30; }      
    }
  }
  for(j = 0; j < 7; j++){
    for(i = 0; i < 7; i++){
      obs_covar[j][i] = 0.0;
      for(k = 0; k < 8; k++){ obs_covar[j][i]+=(T_approx)(temp_matrix[j][k] * jacobian[k][i]); }
      if(std::isnan(obs_covar[j][i])){ obs_covar[j][i] = 0.0; }
      if(std::isinf(obs_covar[j][i])){ obs_covar[j][i] = (obs_covar[j][i] > 0.0) ? 9E30 : -9E30; }            
    }
  }
  
  // ensure diagonal elements are positive definite
  for(i = 0; i < 7; i++){ obs_covar[i][i] = (T_obs) fabs((double) obs_covar[i][i]); }

  // free up memory
  for(i = 0; i < 9; i++){ 
    delete [] obs_vals[i]; 
    delete [] model_vals[i];
  }
  delete [] obs_vals;
  delete [] model_vals;

}

#endif




