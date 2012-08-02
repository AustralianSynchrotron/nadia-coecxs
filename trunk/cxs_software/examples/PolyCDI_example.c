// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.

/**
 * @file PolyCDI_example.c
 *
 * \a PolyCDI_example.c This example reconstructs some poly
 * diffraction data (Lachie's data). The shrinkwrap algorithm is used
 * to improve the reconstruction. A combination of HIO and the
 * error-reduction algorithm are used.
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
#include "PolyCDI.h"
#include "Double_2D.h"

//#include "google/profiler.h"

using namespace std;

int main(void){

  //Define some constants which will be used in the code.

  //the data file name
  string data_file_name = 
    "poly_sim_intensity.tiff";
    //"lowcoherence.dbin";
  //"/home/tdjempire/Desktop/polychromatic_03667.ppm";//"real_sim_intensity.tiff";//"image_files/planar_data.tif";
  //"part_sim_intensity.tiff";
//  "image_files/planar_data.tif";

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = //"image_files/planar_support2.tiff";
  //"/home/tdjempire/Desktop/support.tiff";
  "image_files/planar_support.tiff";

  string data_spectrum_name = "tdj.txt";

  const int cycles=5;
  //number of error reduction iterations to perform before the HIO.
  const int er_iterations1 = 50;

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 100;

  //number of error reduction iterations to perform after the HIO.
  const int er_iterations2 = 50;

  //output the current image every "output_iterations"
  int output_iterations = 10;

  //apply the shrinkwrap algorithm every "shrinkwrap iterations"
  int shrinkwrap_iterations = 200;

  //the number of pixels in x and y
  int nx = 1024;

  int ny = 1024;

  /**** get the diffraction data from file and read into an array *****/

  Double_2D data;
  read_image(data_file_name, data, nx, ny);  

  /****** get the support from a file and read it into an array *****/

  //Every pixel with a zero value is intepreted as being outside
  //the support

  Double_2D support;
  read_image(support_file_name, support, 50, 2);

  //check that the support image has the same dimensions
  //as the data.
  if(support.get_size_x()!=nx || support.get_size_y()!=ny){
    cout << "support file has the wrong dimensions"  << endl;
    return(1);
  }

  /*******  set up the reconstuction *********************/

  //Create a complex 2D field which will hold the result of
  //the reconstruction.
  Complex_2D object_estimate(nx,ny);

  Double_2D spectrum;
  read_txt(data_spectrum_name, spectrum);


  /*double n=4000.0;

    Double_2D spectrum(n, 2);

    double sigma, mean, scale, del;

    sigma=0.035;
    mean=1.4;
    scale=1.0/(sqrt(2.0*3.14159));
    del=(1.6-1.2)/n;

  /*  for(double i=0; i<n; i++){


  spectrum.set(int(i), 0, 1.0/(mean-del*(n/2.0+i)));
  spectrum.set(int(i), 1, scale*exp(-del*((-n/2.0)+i)*del*(-n/2.0+i)/2.0/sigma/sigma));

  std::cout<<spectrum.get(i, 0)<<" "<<spectrum.get(i, 1)<<"\n";

  }

  //Create the spectrum

  /*Double_2D spectrum(1, 1);

  spectrum.set(0, 0, 1.0/1.4);
  spectrum.set(0, 1, 1.0);
   */
  /*  Double_2D spectrum(7, 7);

      spectrum.set(0, 0, 1.0/1.3);
      spectrum.set(0, 1, 0.3);
      spectrum.set(1, 0, 1.0/1.35);
      spectrum.set(1, 1, 0.0);
      spectrum.set(2, 0, 1.0/1.40);
      spectrum.set(2, 1, 0.4);
      spectrum.set(3, 0, 1.0/1.45);
      spectrum.set(3, 1, 0.0);
      spectrum.set(4, 0, 1.0/1.5);
      spectrum.set(4, 1, 0.3);
   */
  //create the poly CDI object which will be used to
  //perform the reconstuction.
  PolyCDI poly(object_estimate, spectrum, 0.9, 4, 0);

  //set the support and intensity
  poly.set_support(support,false);
  poly.set_intensity(data);

  //set the algorithm to hybrid input-output
  poly.set_algorithm(ER);

  //Initialise the current object ESW with a random numbers
  //"0" is the seed to the random number generator
  poly.initialise_estimate(7);

  //  poly.set_fftw_type(FFTW_ESTIMATE);

  //make a 2D object. This will be used to output the 
  //image of the current estimate.
  Double_2D result(nx,ny);
  object_estimate.get_2d(MAG,result);

  /******* for fun, let's get the autocorrelation *****/

  Double_2D autoc(nx,ny);
  //poly.get_intensity_autocorrelation(autoc);
  //write_image("test_autocorrelation.ppm", autoc, true); //"true" means log scale


  //  ProfilerStart("profile");


  /*** run the reconstruction ************/
  for(int a=0; a<cycles; a++){

    for(int i=0; i<er_iterations1; i++){

      cout << "iteration " << i << endl;

      //apply the set of poly CDI projections 
      poly.iterate(); 
      cout << "Current error is "<<poly.get_error()<<endl;

      //every "output_iterations" 
      //output the current estimate of the object
      if(i%output_iterations==0){

	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "poly_example_iteration_" << i +a*(hio_iterations+er_iterations1+er_iterations2+1)<< ".ppm";
	write_image(temp_str.str(), result);

	temp_str.clear();

	//uncomment to output the estimated diffraction pattern
	//poly.propagate_to_detector(object_estimate);
	//object_estimate.get_2d(MAG_SQ,result);
	//temp_str << "diffraction.ppm";
	//write_ppm(temp_str.str(), result, true);
	//poly.propagate_from_detector(object_estimate);
	//object_estimate.get_2d(MAG,result);

	//apply the shrinkwrap algorithm
	//1.5 is the gaussian width in pixels
	//0.1 is the threshold (10% of the maximum pixel).

      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	poly.apply_shrinkwrap(1.5,0.1);

    }

    //now change to the error reduction algorithm 
    poly.set_algorithm(HIO);

    for(int i=er_iterations1; i<(hio_iterations+er_iterations1+1); i++){

      cout << "iteration " << i << endl;

      poly.iterate(); 

      cout << "Current error is "<<poly.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "poly_example_iteration_" << i+a*(hio_iterations+er_iterations1+er_iterations2+1) << ".ppm";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//poly.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	poly.apply_shrinkwrap(1.5,0.1);

    }

    poly.set_algorithm(ER);

    for(int i=er_iterations1+hio_iterations; i<(hio_iterations+er_iterations1+er_iterations2+1); i++){

      cout << "iteration " << i << endl;

      poly.iterate();

      cout << "Current error is "<<poly.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "poly_example_iteration_" << i+a*(hio_iterations+er_iterations1+er_iterations2+1) << ".ppm";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//poly.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	poly.apply_shrinkwrap(1.5,0.1);

    }
  }


  //And we are done. "object_estimate" contained the final estimate of
  //the ESW.

  /** ignore the stuff below  
    Double_2D result2(nx,ny);

    double error=0;
    poly.get_best_result(0,error)->get_2d(MAG,result2);
    write_ppm("best_error.ppm", result2);
    cout << "Best error 0 is "<< error <<endl;

    poly.get_best_result(1,error)->get_2d(MAG,result2);
    write_ppm("best_error_1.ppm", result2);
    cout << "Best error 1 is "<< error <<endl;

    poly.get_best_result(2,error)->get_2d(MAG,result2);
    write_ppm("best_error_2.ppm", result2);
    cout << "Best error 2 is "<< error <<endl;

    poly.get_best_result(3,error)->get_2d(MAG,result2);
    write_ppm("best_error_3.ppm", result2);
    cout << "Best error 3 is "<< error <<endl; **/

  //  ProfilerStop();

  return 0;
}

