/**
 * @file real_example.c
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au>
 *
 * @section DESCRIPTION
 *
 * This file provides an example of running the planar CDI 
 * reconstruction on real data (from the file test_dat.tif)
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Projection.h"
#include "Config.h"
//#include "FourierT.h"
//#include <google/profiler.h>
#include "shrinkwrap.h"

using namespace std;


int main(void){

  string data_file_name = "../examples/real_example_iter_120.ppm";
  //string data_file_name = "test.ppm";
  int size_x,size_y;
  double **data;
  int status = read_ppm(data_file_name, &size_x, &size_y, &data);  
  

  //write out the image before and after the crop and threashold to see 
  //what they look. "true" is used to indicate we want it ouput on log scale.

  cout << "writing before (log scale)"<<endl;
  write_ppm("before.ppm", size_x, size_y, data, false);
  apply_shrinkwrap(size_x,size_y,&data,1.5,0.05);

  //apply_threshold(size_x,size_y,&data,0.2);
  //convolve(size_x,size_y,&data,3);
  cout << "writing after (log scale)" <<endl;
  write_ppm("after.ppm", size_x, size_y, data, false);


  return 0;

}

