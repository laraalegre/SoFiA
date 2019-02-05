#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

// functions using floats

void ThresholdObjs(vector< object_props *> & detections, int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int max_x_size, int max_y_size, int max_z_size, int min_v_size, int max_v_size, int min_los, int max_los, float min_fill, float max_fill, float intens_thresh_min, float intens_thresh_max){
  
  float progress;
  int i, j, k, obj_batch;

  progress = 0.0;
  std::cout << "Thresholding objects . . . " << std::endl;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  for(k = 0; k < NOobj; ++k){
	
    // calculate the obj_batch value for the existing object
    obj_batch = (int) floorf(((float) k / (float) obj_limit));
	
    // move to the next object if this one has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 

      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 

    }
	
    // count the number of LOSs through this object that contain an object section
    i = 0;
    for(j = 0; j < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){

      if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid((j + 1)) > detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid(j)){ ++i; }

    }

    // calculate the size of the detection's bounding box
    j = (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1);

    // apply the thresholds, and if it fails re-initialise the object
    if(((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) < min_x_size)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) < min_y_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) < min_z_size) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < min_v_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) > max_x_size && max_x_size >= 1)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) > max_y_size && max_y_size >= 1) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > max_z_size && max_z_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > max_v_size && max_v_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < intens_thresh_min) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() > intens_thresh_max) || (i < min_los) || (i > max_los && max_los >= 1) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) < min_fill) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) > max_fill && max_fill >= 1)){

      // re-initialise object
      if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_srep();
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_mini();
	
      }
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_size();
      detections[obj_batch][(k - (obj_batch * obj_limit))].Set_srep_update(0);
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit();

    }
    
    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
    // for(i = 0; i < NO_check_obj_ids; ++i)
  }
  std::cout << "* done." << std::endl;

}

void ThresholdObjs(vector< object_props *> & detections, long int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int max_x_size, int max_y_size, int max_z_size, int min_v_size, int max_v_size, int min_los, int max_los, float min_fill, float max_fill, float intens_thresh_min, float intens_thresh_max){
  
  float progress;
  int i, j;
  long int k, obj_batch;

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  for(k = 0; k < NOobj; ++k){
	
    // calculate the obj_batch value for the existing object
    obj_batch = (long int) floor(((double) k / (double) obj_limit));
	
    // move to the next object if this one has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 

      while(progress <= (((double) (k + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 

    }
	
    // count the number of LOSs through this object that contain an object section
    i = 0;
    for(j = 0; j < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){

      if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid((j + 1)) > detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid(j)){ ++i; }

    }

    // calculate the size of the detection's bounding box
    j = (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1);

    // apply the thresholds, and if it fails re-initialise the object
    if(((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) < min_x_size)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) < min_y_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) < min_z_size) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < min_v_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) > max_x_size && max_x_size >= 1)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) > max_y_size && max_y_size >= 1) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > max_z_size && max_z_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > max_v_size && max_v_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < intens_thresh_min) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() > intens_thresh_max) || (i < min_los) || (i > max_los && max_los >= 1) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) < min_fill) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) > max_fill && max_fill >= 1)){

      // re-initialise object
      if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_srep();
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_mini();
	
      }
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_size();
      detections[obj_batch][(k - (obj_batch * obj_limit))].Set_srep_update(0);
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit();

    }
    
    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
    // for(k = 0; k < NOobj; ++k)
  }
  std::cout << "* done." << std::endl;

}

// functions using doubles

void ThresholdObjs(vector< object_props_dbl *> & detections, int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int max_x_size, int max_y_size, int max_z_size, int min_v_size, int max_v_size, int min_los, int max_los, float min_fill, float max_fill, double intens_thresh_min, double intens_thresh_max){
  
  float progress;
  int i, j, k, obj_batch;

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  for(k = 0; k < NOobj; ++k){
	
    // calculate the obj_batch value for the existing object
    obj_batch = (int) floorf(((float) k / (float) obj_limit));
	
    // move to the next object if this one has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 

      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 

    }
	
    // count the number of LOSs through this object that contain an object section
    i = 0;
    for(j = 0; j < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){

      if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid((j + 1)) > detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid(j)){ ++i; }

    }

    // calculate the size of the detection's bounding box
    j = (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1);

    // apply the thresholds, and if it fails re-initialise the object
    if(((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) < min_x_size)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) < min_y_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) < min_z_size) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < min_v_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) > max_x_size && max_x_size >= 1)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) > max_y_size && max_y_size >= 1) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > max_z_size && max_z_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > max_v_size && max_v_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < intens_thresh_min) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() > intens_thresh_max) || (i < min_los) || (i > max_los && max_los >= 1) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) < min_fill) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) > max_fill && max_fill >= 1)){

      // re-initialise object
      if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_srep();
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_mini();
	
      }
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_size();
      detections[obj_batch][(k - (obj_batch * obj_limit))].Set_srep_update(0);
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit();

    }
    
    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
    // for(i = 0; i < NO_check_obj_ids; ++i)
  }
  std::cout << "* done." << std::endl;

}

void ThresholdObjs(vector< object_props_dbl *> & detections, long int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int max_x_size, int max_y_size, int max_z_size, int min_v_size, int max_v_size, int min_los, int max_los, float min_fill, float max_fill, double intens_thresh_min, double intens_thresh_max){
  
  float progress;
  int i, j;
  long int k, obj_batch;

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  for(k = 0; k < NOobj; ++k){
	
    // calculate the obj_batch value for the existing object
    obj_batch = (long int) floor(((double) k / (double) obj_limit));
	
    // move to the next object if this one has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 

      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 

    }
	
    // count the number of LOSs through this object that contain an object section
    i = 0;
    for(j = 0; j < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){

      if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid((j + 1)) > detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid(j)){ ++i; }

    }

    // calculate the size of the detection's bounding box
    j = (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1);

    // apply the thresholds, and if it fails re-initialise the object
    if(((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) < min_x_size)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) < min_y_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) < min_z_size) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < min_v_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) > max_x_size && max_x_size >= 1)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) > max_y_size && max_y_size >= 1) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > max_z_size && max_z_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > max_v_size && max_v_size >= 1) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < intens_thresh_min) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() > intens_thresh_max) || (i < min_los) || (i > max_los && max_los >= 1) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) < min_fill) || ((((float) detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels())/((float) j)) > max_fill && max_fill >= 1)){

      // re-initialise object
      if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_srep();
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_mini();
	
      }
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_size();
      detections[obj_batch][(k - (obj_batch * obj_limit))].Set_srep_update(0);
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit();

    }
    
    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
    // for(k = 0; k < NOobj; ++k)
  }
  std::cout << "* done." << std::endl;

}

