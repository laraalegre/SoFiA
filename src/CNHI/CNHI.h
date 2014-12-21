#include<iostream>
#include<cmath>

#ifndef _CNHIsfinder_
#define _CNHIsfinder_

// definitions used by probkuiper
#define EPS1 0.00001
#define EPS2 1.0e-08

using namespace std;

template <typename dtype>
dtype probkuiper(dtype alam){
  
  int j;
  dtype pa2,na2,fac=2.0,sum=0.0,term,termbf=0.0;

  if (alam < 0.4){ return 1.0; }

  pa2 = 4.0*alam*alam;
  na2 = -2.0*alam*alam;
  for(j=1;j<=100;j++){
    term=fac*((j*j*pa2) - 1)*exp(na2*j*j);
    sum += term;
    if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
    termbf=fabs(term);
  }
  return 1.0;
}

template <typename dtype>
void CNHI_heapsort(int n, dtype ra[], int indices[]){

  int i,ir,j,l;
  dtype rra;
  int rraindex;

  if (n<2) return;
  l = (n >> 1) + 1;
  ir = n;
  for(;;){

    if (l > 1) {
      rra = ra[indices[--l]];
      rraindex = indices[l];
    } else {
      rra=ra[indices[ir]];
      rraindex = indices[ir];
      indices[ir] = indices[1];
      if (--ir == 1) {
	indices[1] = rraindex;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[indices[j]] < ra[indices[j+1]]) j++;
      if (rra < ra[indices[j]]) {
	indices[i] = indices[j];
	i = j;
	j <<= 1;
      } else break;
    }
    indices[i] = rraindex;
  }
}

template<typename dtype, typename mtype, typename rtype>
  void CNHI_find_sources(dtype * data, mtype * mask, int verbose, int metric[9], rtype p_req, int tr_min, int tr_max, mtype flag_val, bool median_test, rtype q_req){
 
  int * map_XtoI_pos, * map_ItoX_pos, * tr_order, i, j, k, p, tr_start, tr_size, proc_min, proc_max;
  int * p_scale, above, below, NOx, NOy, NOz, median_rem_pos;
  dtype * data_LoS;
  double D_min, D_max, dist, alam, * p_val, signif, median_tr, * median_rem, progress;
  bool rem_flag, run_flag;

  // assign values from the metric
  NOx = metric[6];
  NOy = metric[7];
  NOz = metric[8];

  // Create memory.
  data_LoS = new dtype[NOz];
  map_XtoI_pos = new int[NOz];
  map_ItoX_pos = new int[NOz];
  tr_order = new int[NOz];
  p_val = new double[NOz];
  p_scale = new int[NOz];
  median_rem = new double[NOz];

  // set test region size range
  if(tr_min < 1){ tr_min = 1; }
  if(tr_max <= 0){ tr_max = (int) ceilf((float) NOz / 2.0); }
  if(tr_max > (NOz - tr_min)){ tr_max = NOz - tr_min; }

  // Calculate the minimum test region size that delivers the required Q value
  while((tr_min * (NOz - tr_min) / NOz) < q_req){ ++tr_min; }    

  // Check tr_max against tr_min
  if(tr_max < tr_min){ tr_max = tr_min; }

  if(verbose > 1){

    cout << "Using tr_min = " << tr_min << " to achieve Q_req = " << q_req << endl;
    cout << "Testing scales up to tr_max = " << tr_max << endl;
    cout << "Requiring sources to have a probability <= " << p_req << endl;
    if(median_test){ cout << "Requiring sources to have a median >= that of the remaining data." << endl; }
    cout << "Updating mask with flag value = " << flag_val << endl;
    cout << "Using metric: ";
    for(i = 0; i < 9; ++i){ cout << metric[i] << " "; }
    cout << endl;

  }

  // Initialise the progress variable to 0.
  progress = 0.0;

  // Process each line of sight through the datacube individually.
  if(verbose > 0){
   
    cout << "Processing each line of sight through the datacube . . . " << endl;
    cout << "0   :   : | :   :   |   :   : | :   :   100% done." << endl;
  
  }
  for(j = 0; j < NOy; ++j){

    for(i = 0; i < NOx; ++i){
      
      // 1. Initilise arrays for this line-of-sight.
      for(k = 0; k < NOz; ++k){ data_LoS[k] = data[i + metric[3] + ((j + metric[4]) * metric[0]) + ((k + metric[5]) * metric[0] * metric[1])]; }
      for(k = 0; k < NOz; ++k){ map_ItoX_pos[k] = k; }
      for(k = 0; k < NOz; ++k){ p_val[k] = -99.0; }
      for(k = 0; k < NOz; ++k){ p_scale[k] = 0; }
            
      // 2. Populate sorted index array. --> Maps intensity order to physical location.
      CNHI_heapsort(NOz,data_LoS,map_ItoX_pos - 1);
      
      // 3. Generate inversion array for sorted index array. --> Maps physical location to intensity order.
      for(k = 0; k < NOz; ++k){ map_XtoI_pos[map_ItoX_pos[k]] = k; }

      // 4. For each test region size:
      for(tr_size = tr_min; tr_size <= tr_max; ++tr_size){

	// 4. Initialise the position of the `left edge/start' of the test region
	tr_start = 0;

	// b. Build initial test region.
	for(k = 0; k < NOz && k < tr_size; ++k){ tr_order[k] = k; }

	// c. Sort initial test region.
	CNHI_heapsort(tr_size,data_LoS,tr_order - 1);

	if(median_test){

	  // d. Calculate the tr_size+1 possible medians of NOz-tr_size data points
	  if(((NOz - tr_size) % 2) == 0){

	    for(k = 0; k < tr_size+1; ++k){

	      median_rem[k] = 0.5 * (data_LoS[map_ItoX_pos[k + ((NOz - tr_size)/2) - 1]] + data_LoS[map_ItoX_pos[k + ((NOz - tr_size)/2)]]);

	    }

	  } else {

	    for(k = 0; k < tr_size+1; ++k){

	      median_rem[k] = data_LoS[map_ItoX_pos[k + ((NOz - tr_size - 1)/2)]];
	      
	    }

	  }
	
	  // e. Initialise the median_rem_pos value that indexes the median_rem array
	  median_rem_pos = 0;
	  while(data_LoS[tr_order[median_rem_pos]] <= median_rem[median_rem_pos]){ if(median_rem_pos < tr_size){ ++median_rem_pos; } else { break; } }

	} else { median_rem_pos = 0; }

	// f. Apply median test to test region if requested. Must be greater than or equal to median of remainder.
	run_flag = true;
	if(median_test){
	  
	  // calculate median of test region
	  if((tr_size % 2) == 0){

	    median_tr = 0.5 * (data_LoS[tr_order[0] + tr_order[(tr_size/2) - 1]] + data_LoS[tr_order[0] + tr_order[(tr_size/2)]]);

	  } else {
	    
	    median_tr = data_LoS[tr_order[0] + tr_order[(tr_size - 1)/2]];

	  }

	  // update run_flag
	  if(median_tr < median_rem[median_rem_pos]){ run_flag = false; }
	  
	}

	if(run_flag){

	  // g. Initialise "to process" test range: 0 to size of test region.
	  proc_min = 0;
	  proc_max = tr_size - 1;

	  // h. For each element in "test range", test adjacent c.f.d. points and update D values.
	  D_min = 9.0E9;
	  D_max = -9.0E9;
	  for(p = proc_min; p <= proc_max; ++p){
	    
	    // Notes on mapping between physical and intensity positions
	    //
	    // physical location of intensity sorted elements in test region = tr_order[p]
	    // overall intensity position of this physical location = map_XtoI_pos[tr_order[p]]
	    // physical location of intensity adjacent elements = map_ItoX_pos[map_XtoI_pos[tr_order[p]] +/- 1]
	    // 
	    // below and above store the physical location of intensity adjacent elements
	    // overall intensity position of these physical locations = map_XtoI_pos[above or below]
	    
	    // I) find "adjacent" non-identical value below the test point, which is not in the test region
	    k = 0;
	    if(map_XtoI_pos[tr_order[p]] > 0){
	      
	      below = map_XtoI_pos[tr_order[p]] - 1;
	      while((data_LoS[map_ItoX_pos[below]] >= data_LoS[tr_order[p]]) || ((map_ItoX_pos[below] >= tr_start) && (map_ItoX_pos[below] < (tr_start + tr_size)))){ 

		if((map_ItoX_pos[below] >= tr_start) && (map_ItoX_pos[below] < (tr_start + tr_size))){ ++k; }
		--below; if(below < 0){ break; } 

	      }
	      
	    } else {
	      
	      below = -1;
	      
	    }
	    
	    // II) calculate c.f.d. distance to adjacent-below point
	    if(below == -1){
	      
	      dist = (double) (1 + p) / (double) tr_size;
	      
	    } else {
	      
	      dist = ((double) p / (double) tr_size) - ((double) (below + k - p) / (double) (NOz - tr_size));
	      
	    }
	    
	    // III) update D_min and D_max for the c.f.d. distance
	    D_min = (dist <= D_min) ? dist : D_min;
	    D_max = (dist >= D_max) ? dist : D_max;
	    
	    // IV) find "adjacent" non-identical value above the test point, which is not in the test region
	    k = 0;
	    if(map_XtoI_pos[tr_order[p]] < (NOz - 1)){
	      
	      above = map_XtoI_pos[tr_order[p]] + 1;
	      while((data_LoS[map_ItoX_pos[above]] <= data_LoS[tr_order[p]]) || ((map_ItoX_pos[above] >= tr_start) && (map_ItoX_pos[above] < (tr_start + tr_size)))){

		if((map_ItoX_pos[above] >= tr_start) && (map_ItoX_pos[above] < (tr_start + tr_size))){ ++k; }
		++above; if(above >= NOz){ break; } }
	      
	    } else {
	      
	      above = NOz;
	      
	    }
	    
	    // V) calculate c.f.d. distance to adjacent-above point
	    if(above == NOz){
	      
	      dist = 1.0 - ((double) (p + 1) / (double) tr_size);
	      
	    } else {
	      
	      dist = ((double) (p + 1) / (double) tr_size) - ((double) (above - p - k - 1) / (double) (NOz - tr_size));
	      
	    }
	    
	    // VI) update D_min and D_max for the c.f.d. distance
	    D_min = (dist <= D_min) ? dist : D_min;
	    D_max = (dist >= D_max) ? dist : D_max;
	    
	    // for(p = proc_min; p <= proc_max; ++p){
	  }	  
	  
	  // i. If p <= p_req, then update results array [factoring in previous results].
	  dist = 0.0;
	  if(!(isnan(D_min)) && !(isinf(D_min)) && (D_min <= 0.0)){ dist+=(-1.0*D_min); }
	  if(!(isnan(D_max)) && !(isinf(D_max)) && (D_max >= 0.0)){ dist+=D_max; }	
	  alam = sqrtf(((tr_size * (NOz - tr_size)) / NOz));
	  alam = (alam+0.155+(0.24/alam))*dist;
	  signif = probkuiper(alam);
	  alam = tr_size * (NOz - tr_size) / NOz;
	  if((alam >= q_req) && (signif <= p_req) && ((signif <= p_val[tr_start]) || (p_val[tr_start] < 0.0) || (p_scale[tr_start] <= 0))){ 
	    
	    p_val[tr_start] = signif; 
	    p_scale[tr_start] = tr_size;
	    
	  }

	  // if(run_flag)
	}

	// j. While new points are available:
	while((tr_start + tr_size) < NOz){

	  // I) Replace "oldest/closest to origin" point with new point, and initialise "to process" to location of replaced point.
	  // II) Initialise "to process" to location of replaced point.
	  for(p = 0; p < tr_size; ++p){
	      
	    if(tr_order[p] == tr_start){
		
	      tr_order[p] = tr_start + tr_size;
	      proc_min = proc_max = p;
	      break;
	      
	    }
	    
	    // for(p = 0; p < tr_size; ++p)
	  }
	  
	  // III) Use insertion sort to move new point to correct position.
	  // IV) Update "to process" for final position of new point.
	  if(proc_max < (tr_size - 1)){
	    
	    while((map_XtoI_pos[tr_order[proc_max]] > map_XtoI_pos[tr_order[proc_max + 1]])){
	      
	    
	      p = tr_order[proc_max + 1];
	      tr_order[proc_max + 1] = tr_order[proc_max];
	      tr_order[proc_max] = p;
	      ++proc_max;
	      if(proc_max == (tr_size - 1) || proc_max == (NOz - 1)){ break; }
	      
	    }
	  
	  }
	  if(proc_min > 0){
	  
	    while((map_XtoI_pos[tr_order[proc_min]] < map_XtoI_pos[tr_order[proc_min - 1]])){
		
	      p = tr_order[proc_min - 1];
	      tr_order[proc_min - 1] = tr_order[proc_min];
	      tr_order[proc_min] = p;
	      --proc_min;
	      if(proc_min == 0){ break; }
	      
	    }
	  
	  }
	  
	  // V) Apply  median test to test region if requested. Must be greater than or equal to median of remainder.
	  run_flag = true;
	  if(median_test){
	  
	    // calculate median of test region
	    if((tr_size % 2) == 0){

	      median_tr = 0.5 * (data_LoS[tr_order[0] + tr_order[(tr_size/2) - 1]] + data_LoS[tr_order[0] + tr_order[(tr_size/2)]]);
	      
	    } else {
	      
	      median_tr = data_LoS[tr_order[0] + tr_order[(tr_size - 1)/2]];
	      
	    }
	    
	    // update median_rem_pos 
	    if((data_LoS[tr_order[proc_min]] <= median_rem[median_rem_pos]) && (data_LoS[tr_order[proc_max]] > median_rem[median_rem_pos])){ --median_rem_pos; }
	    if((data_LoS[tr_order[proc_min]] > median_rem[median_rem_pos]) && (data_LoS[tr_order[proc_max]] <= median_rem[median_rem_pos])){ ++median_rem_pos; }

	    // update run_flag
	    if(median_tr < median_rem[median_rem_pos]){ run_flag = false; }
	    
	  }

	  if(run_flag){

	    // VI) For each element of test range in "to process" interval, test adjacent c.f.d. points and update D values.
	    D_min = 9.0E30;
	    D_max = -9.0E30;
	    for(p = proc_min; p <= proc_max; ++p){

	      // Notes on mapping between physical and intensity positions
	      //
	      // physical location of intensity sorted elements in test region = tr_order[p]
	      // overall intensity position of this physical location = map_XtoI_pos[tr_order[p]]
	      // physical location of intensity adjacent elements = map_ItoX_pos[map_XtoI_pos[tr_order[p]] +/- 1]
	      // 
	      // below and above store the physical location of intensity adjacent elements
	      // overall intensity position of these physical locations = map_XtoI_pos[above or below]
	      
	      // find "adjacent" non-identical value below the test point
	      k = 0;
	      if(map_XtoI_pos[tr_order[p]] > 0){
	      
		below = map_XtoI_pos[tr_order[p]] - 1;
		while((data_LoS[map_ItoX_pos[below]] >= data_LoS[tr_order[p]]) || ((map_ItoX_pos[below] >= tr_start) && (map_ItoX_pos[below] < (tr_start + tr_size)))){ 
		  
		  if((map_ItoX_pos[below] >= tr_start) && (map_ItoX_pos[below] < (tr_start + tr_size))){ ++k; }
		  --below; if(below < 0){ break; } 
		  
		}
		
	      } else {
		
		below = -1;
		
	      }

	      // calculate c.f.d. distance to adjacent-below point
	      if(below == -1){
	      
		dist = (double) (1 + p) / (double) tr_size;
	      
	      } else {
	      
		dist = ((double) p / (double) tr_size) - ((double) (below + k - p) / (double) (NOz - tr_size));
		
	      }

	      // update D_min and D_max for the c.f.d. distance
	      D_min = (dist <= D_min) ? dist : D_min;
	      D_max = (dist >= D_max) ? dist : D_max;
	      
	      // find "adjacent" non-identical value above the test point, that is not in the 
	      k = 0;
	      if(map_XtoI_pos[tr_order[p]] < (NOz - 1)){
	      
		above = map_XtoI_pos[tr_order[p]] + 1;
		while((data_LoS[map_ItoX_pos[above]] <= data_LoS[tr_order[p]]) || ((map_ItoX_pos[above] >= tr_start) && (map_ItoX_pos[above] < (tr_start + tr_size)))){

		  if((map_ItoX_pos[above] >= tr_start) && (map_ItoX_pos[above] < (tr_start + tr_size))){ ++k; }
		  ++above; if(above >= NOz){ break; } }
		
	      } else {
				
		above = NOz;
		
	      }

	      // calculate c.f.d. distance to adjacent-above point
	      if(above == NOz){
	      
		dist = 1.0 - ((double) (p + 1) / (double) tr_size);
	      
	      } else {
	      
		dist = ((double) (p + 1) / (double) tr_size) - ((double) (above - p - k - 1) / (double) (NOz - tr_size));
	      
	      }
	      
	      // update D_min and D_max for the c.f.d. distance
	      D_min = (dist <= D_min) ? dist : D_min;
	      D_max = (dist >= D_max) ? dist : D_max;
	      
	      // for(p = proc_min; p <= proc_max; ++p){
	    }	  
	    
	    // VII) If p <= p_req, then update results array [factoring in previous results].
	    dist = 0.0;
	    if(!(isnan(D_min)) && !(isinf(D_min)) && (D_min <= 0.0)){ dist+=(-1.0*D_min); }
	    if(!(isnan(D_max)) && !(isinf(D_max)) && (D_max >= 0.0)){ dist+=D_max; }	
	    alam = sqrtf(((tr_size * (NOz - tr_size)) / NOz));
	    alam = (alam+0.155+(0.24/alam))*dist;
	    signif = probkuiper(alam);
	    alam = tr_size * (NOz - tr_size) / NOz;
	    if((alam >= q_req) && (signif <= p_req) && ((signif <= p_val[tr_start]) || (p_val[tr_start] < 0.0) || (p_scale[tr_start] <= 0))){ 
	      
	      p_val[tr_start] = signif; 
	      p_scale[tr_start] = tr_size;
	      
	    }

	    // if(run_flag)
	  }

	  // increment tr_start
	  ++tr_start;

	  // while((tr_start + tr_size) < NOz)
	}
	
	//for(tr_size = tr_min; tr_size <= tr_max; ++tr_size)
      }

      // Collapse results into a unique set of non-overlapping detections
      for(k = 0; k < NOz; ++k){

	// test if a detection starts at this position
	if(p_scale[k] > 0){
	  
	  // initialise remove flag
	  rem_flag = false;
	  
	  // if any other detections start in this detection's range, set the rem_flag to true if the other detection is 'centered' within this detection and is more significant
	  // remove the other detection if it is less significant
	  for(p = k + 1; p < (k + p_scale[k]); ++p){

	    // process another detection if it starts at point p within the current detection's range
	    if(p_scale[p] > 0){

	      if((((double) p) + 0.5 * ((double) p_scale[p])) <= (double) k){

		if(p_val[p] < p_val[k]){

		  rem_flag = true;

		} else {

		  p_scale[k] = 0;
		  p_val[k] = -99.0;

		}

		// if((((float) p) + 0.5 * ((float) p_scale[p])) <= (float) k)
	      }

	      // if(p_scale[p] > 0)
	    }

	    // for(p = k + 1; p < (k + p_scale[k]); ++p){
	  }
	  
	  // remove this detection if rem_flag is true
	  if(rem_flag){

	    p_val[k] = -99.0;
	    p_scale[k] = 0;

	  }

	  // if(p_scale[k] > 0)
	}

	// for(k = 0; k < NOz; ++k)
      }

      // Update the mask array
      for(k = 0; k < NOz; ++k){

	// test if a detection starts at this position
	if(p_scale[k] > 0){

	  for(p = k; p < (k + p_scale[k]); ++p){

	    mask[i + metric[3] + ((j + metric[4]) * metric[0]) + ((p + metric[5]) * metric[0] * metric[1])] = flag_val;

	  }

	}

	// for(k = 0; k < NOz; ++k)
      }
      
      // update progress
      if(verbose > 0){

	while(progress <= ((double)(1 + i + (j * NOx)) / (double)(NOx * NOy))){ 

	  cout << "*"; 
	  cout << flush; 
	  progress+=0.025; 
	  
	}

      }

      // for(i = 0; i < NOx; ++i)
    }

    // for (j = 0; j < NOy; ++j)
  }

  // update progress
  if(verbose > 0){

    while(progress <= 1.0){ 
    
      cout << "*"; 
      cout << flush; 
      progress+=0.025; 
      
    }    
    cout << endl;

  }

  // 5. Free up memory.
  delete [] data_LoS;
  delete [] map_XtoI_pos;
  delete [] map_ItoX_pos;
  delete [] tr_order;
  delete [] p_val;
  delete [] p_scale;
  delete [] median_rem;

  // 6. End of function.
  return;

}

#endif
