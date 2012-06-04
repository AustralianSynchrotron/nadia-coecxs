// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.
//
// date last modified: 17/01/2011

/**
 * @file PartialCDI_simulation_example.c
 *
 * \a PartialCDI_simulation_example.c The object from the
 * PartialCDI_example (in fact I used the reconstructed image as the
 * object) is used to simulate a diffraction pattern for partial
 * CDI. The diffraction pattern is thresholded to make it more
 * realistic and then CDI reconstruction is performed.
 *
 */

#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "PartialCDI.h"
#include <sstream>

using namespace std;


/**************************************/
int main(void){

  //define some constants which will be used in the code:

  //the data file name
  const static char * data_file_name = "image_files/object.tiff";

  //the approx. level of background noise in the image
  const double noise_level = 15;

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  const static char * support_file_name = "image_files/planar_support.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 0;

  //the number of error-reduction iterations to perform.
  const int er_iterations = 0;

  //output the current image ever "output_iterations"
  const int output_iterations = 50;

  int nleg=7;
  int nmodes=7;


  /****** get the object from an image file ****************/

  //get the data from file
  
  Double_2D data;

  //read the data into an array
  int status = read_tiff(data_file_name, data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  int n_x = data.get_size_x();
  int n_y = data.get_size_y();
  
  //fill the complex no. with image data
  Complex_2D input(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      input.set_value(i,j,REAL, 1.0/sqrt(2.0)*(data.get(i,j)));
      input.set_value(i,j,IMAG, 1.0/sqrt(2.0)*(data.get(i,j)));
    }
  }

  /**** create the projection/reconstruction object *****/

  Complex_2D pattern(n_x,n_y);
  PartialCDI my_partial(pattern, 0.9, 1.0, 100000.0, 8.0, 8.0, 4, 0);

  my_partial.set_threshold(+1.0e-11);


  my_partial.initialise_matrices(nleg, nmodes);

  Double_2D result(n_x, n_y);

  ostringstream temp_str0 ( ostringstream::out ) ;
  pattern.get_2d(MAG,result);
  temp_str0 << "field.ppm";
  write_image(temp_str0.str(), result, false);

  ostringstream temp_str1 ( ostringstream::out ) ;
  temp_str1 << "fieldl.ppm";
  write_image(temp_str1.str(), result, true);


  for(int i=0; i< nmodes*nmodes; i++){
    ostringstream temp_str0 ( ostringstream::out ) ;
    my_partial.get_mode(i).get_2d(REAL,result);
    temp_str0 << "mode"<<i<<".ppm";
    write_image(temp_str0.str(), result);
  }


  //propagate to the detector plane
  my_partial.set_transmission(input);

  Double_2D intensity=my_partial.propagate_modes_to_detector();

  //write the fourier transform to file.
  //Double_2D intensity(n_x,n_y);
  //intensity=my_partial.get_intensity();

  //apply a threshold to make the simulation a bit more realistic
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      intensity.set(i,j,intensity.get(i,j)-noise_level);
      if(intensity.get(i,j)<0)
	intensity.set(i,j,0);
    }
  }


  //write the output to file (use log scale)
  write_ppm("part_sim_intensity.ppm",intensity,true);
  write_tiff("part_sim_intensity.tiff",intensity,false);

  //Double_2D tmp(n_x, n_y);
  //input.get_2d(PHASE, tmp);
  //write_ppm("part_sim_intensity.ppm",tmp,true);

  /******** get the support from file ****************************/

  Double_2D support(n_x,n_y);
  status = read_tiff(support_file_name, support);

  /*************** do the reconstruction *******************/

  //create a project object and set the options.
  my_partial.set_support(support);
  my_partial.set_intensity(intensity);
  my_partial.set_algorithm(HIO);

  //set the inital guess to be random inside the support
  //and zero outside. Note that this must be called
  //after "my_partial.set_support()"
  //my_partial.initialise_estimate(0);

  //make a temporary arrary
  //Double_2D result(n_x,n_y);

  //apply the projection operators
  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    my_partial.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      pattern.get_2d(MAG,result);
      temp_str << "part_sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      my_partial.apply_shrinkwrap();

    }
  }

  my_partial.set_algorithm(ER);

  for(int i=hio_iterations; i<(er_iterations+hio_iterations+1); i++){

    cout << "iteration " << i << endl;

    my_partial.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      pattern.get_2d(MAG,result);
      temp_str << "part_sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      my_partial.apply_shrinkwrap();

    }
  }




  return 0;
}

