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
//#include "FourierT.h"
#include <google/profiler.h>

using namespace std;


int main(void){

  ProfilerStart("profiler.prof");

  //define some constants which will be used in the code:

  //the data file name
  const static char * data_file_name = "image_files/test_dat.tif";

  //the approx. level of background noise in the image
  const double noise_level = 30;

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  const static char * support_file_name = "image_files/support_2.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 200;
  
  //number of error reduction iterations to perform after the HIO.
  const int er_iterations = 0;

  //output the current image ever "output_iterations"
  const int output_iterations = 20;


  /*******  get the diffraction data from file and read into an array *****/
  int size_x, size_y;
  double ** data;
  int status = read_tiff(data_file_name, &size_x, &size_y, &data);  
  
  //check that the file could be opened okay
  if(!status){
    cout << "failed to get data from "<< data_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }

  /***** set up a 2d array which whill hold the magnitude information ****/
  //The image is a bit big which would slow down the fourier transforms.
  //So lets crop the edges; the image will go from 2048x2048 to 1024x1024.

  int nx = size_x/2;
  int ny = size_y/2;

  //declare the array, allocate some memory for it and fill it.
  double ** intensity = new double*[nx];

  //loop over the pixels in the image
  for(int i=0; i < nx; i++){
    intensity[i]= new double[ny];    
    for(int j=0; j<ny; j++){
      //copy to the new array
      intensity[i][j] = data[i+size_x/4][j+size_y/4];

      //apply a threashold to the data to remove background
      intensity[i][j]-=noise_level;
      if(intensity[i][j]<0)
	intensity[i][j]=0;
    }
  }

  //write out the image before and after the crop and threashold to see 
  //what they look. "true" is used to indicate we want it ouput on log scale.
  write_ppm("data_before.ppm", size_x, size_y, data, true);
  write_ppm("data_after.ppm", nx, ny, intensity, true);


  /******* get the support from file and read it into an array *****/

  double ** support;
  int nx_s, ny_s;
  status = read_tiff(support_file_name, &nx_s, &ny_s, &support);  
  if(!status){
    cout << "failed to get data from "<< support_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }
  if( nx_s != nx || ny_s != ny){
    cout << "dimensions of the support to not match ... exiting"  << endl;
    return(1);
  }

  /*******  set up the reconstuction *********************/

  //create the projection object which will be used to
  //perform the reconstuction.
  Complex_2D object_estimate(nx,ny);
  Projection proj(&object_estimate);
 
  //set the support and intensity
  proj.set_support(support);
  proj.set_intensity(intensity);
  //set the algorithm to hybrid input-output
  proj.set_algorithm(HIO);


  //Initialise the current object ESW with a random numbers
  proj.initialise_estimate(0);
  
  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.
  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];


  /******* now get the autocorrelation ************/

  double ** autoc = proj.get_intensity_autocorrelation();
  write_ppm("test_autocorrelation.ppm", nx, ny, autoc, true);

  /*** run the reconstruction ************/

  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
    proj.iterate(); 
    
    if(i%output_iterations==0){
      //output the current estimate of the object
      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,&result);
      temp_str << "real_example_iter_" << i << ".ppm";
      write_ppm(temp_str.str(), nx, ny, result);
      temp_str.clear();

      //uncomment to output the estimated 
      /**Complex_2D * temp = object_estimate.clone();
      fft.perform_forward_fft(temp);
      temp->get_2d(MAG_SQ,&result);
      temp_str << "diffraction.ppm";
      write_ppm(temp_str.str(), nx, ny, result, true); 
      delete temp;
      **/
    }

  }
  
  //now change to the error reduction algorithm 
  proj.set_algorithm(ER);

  for(int i=hio_iterations; i<(hio_iterations+er_iterations+1); i++){

    cout << "iteration " << i << endl;
    
    proj.iterate(); 
    
    if(i%output_iterations==0){
      //output the current estimate of the object
      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,&result);
      temp_str << "real_example_iter_" << i << ".ppm";
      write_ppm(temp_str.str(), nx, ny, result);
      temp_str.clear();
    }
    
  }

  //clean up
  for(int i=0; i< nx; i++){
    delete[] intensity[i];
    delete[] result[i];
    delete[] support[i];
    delete[] autoc[i];
  }

  for(int i=0; i< size_x; i++){
    delete[] data[i];
  }

  delete[] intensity;
  delete[] result;
  delete[] support;
  delete[] autoc;
  delete[] data;

  ProfilerStop();

  return 0;
}

