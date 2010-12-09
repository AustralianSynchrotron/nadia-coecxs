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
#include "FresnelCDI_WF.h"
#include "Config.h"
#include "shrinkwrap.h"
//#include "FourierT.h"
#include <google/profiler.h>

using namespace std;


int main(void){

  //the data file name
  //string data_file_name = "image_files/nadia_wf_data.tif";
  string data_file_name = "../tools/wf_A_real.ppm";


  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = "image_files/wf_A_support.tiff";

  /*******  get the diffraction data from file and read into an array *****/
  int nx, ny;
  double ** data;
  int status = read_ppm(data_file_name, &nx, &ny, &data);  
   
  //check that the file could be opened okay
  if(!status){
    cout << "failed to get data from "<< data_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }


  //write out the image before and after the crop and threashold to see 
  //what they look. "true" is used to indicate we want it ouput on log scale.
  //write_ppm("data_before.ppm", size_x, size_y, data, true);
  //write_ppm("data_after.ppm", nx, ny, intensity, true);


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
  Complex_2D zone_estimate(nx,ny);

  //dimensions in micron
  /** double wavelength = 1.240/(2.54*1000); //0.000488
  cout << "wavelength " << wavelength << endl;

  double focallength = 2*80*0.05/wavelength; //16,400
  cout << "focallength " << focallength << endl;

  FresnelCDI_WF proj(&zone_estimate,
		       wavelength,
		       focallength,
		       800000,
		       13.5); **/
 
  
  FresnelCDI_WF proj(&zone_estimate,
		     4.892e-10,
		     16.353e-3,
		     0.909513-16.353e-3,
		     13.5e-3);

  //set the support and intensity
  proj.set_support(support);
  proj.set_intensity(data);
  //set the algorithm to hybrid input-output
  //proj.set_algorithm(HIO);

  //Initialise the current object ESW with a random numbers
  proj.initialise_estimate(0);
  
  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.
  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];

  ProfilerStart("my_prof.prof");

  /*** run the reconstruction ************/
  for(int i=0; i<100; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
    if(i%5==0){

      //output the current estimate of the detector
      ostringstream temp_str ( ostringstream::out ) ;
      zone_estimate.get_2d(PHASE,&result);
      temp_str << "detector_phase_iter_" << i << ".ppm";
      write_ppm(temp_str.str(), nx, ny, result);
      temp_str.clear();

      //uncomment to output the estimated 
      //Complex_2D * temp = zone_estimate.clone();
      //fft.perform_forward_fft(temp);
      //temp->get_2d(MAG_SQ,&result);
      //temp_str << "diffraction.ppm";
      //write_ppm(temp_str.str(), nx, ny, result, true); 
      //delete temp;**/
      
      
      //apply the shrinkwrap algorithm
      //apply_shrinkwrap(nx,ny,&result,2,0.1);
      //proj.set_support(result);
      //write_ppm("shrink.ppm", nx, ny, result);
    }
    proj.iterate(); 
    

  }

  ProfilerStop();

  //now change to the error reduction algorithm 

  //clean up
  for(int i=0; i< nx; i++){
    delete[] data[i];
    delete[] result[i];
    delete[] support[i];
  }

  delete[] data;
  delete[] result;
  delete[] support;

  //ProfilerStop();

  return 0;
  }

