// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.

/**
 * @file PlanarCDI_example.c
 *
 * \a PlanarCDI_example.c This example reconstructs some planar
 * diffraction data (Lachie's data). The shrinkwrap algorithm is used
 * to improve the reconstruction. A combination of HIO and the
 * error-reduction algorithm are used.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include <io.h>
#include <Complex_2D.h>
#include <PlanarCDI.h>
#include <Double_2D.h>
#include <iomanip>

using namespace std;

int main(void){

  //Define some constants which will be used in the code.

  //the data file name
  string data_file_name = "image_files/planar_data.tif";

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = "image_files/planar_support.tiff";

  //number of cycles of ER and HIO to repeat
  const int cycles=2;

  //number of error reduction iterations to perform before the HIO.
  const int er_iterations1 = 50;

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 100;

  //number of error reduction iterations to perform after the HIO.
  const int er_iterations2 = 50;

  //output the current image every "output_iterations"
  int output_iterations = 10;

  //apply the shrinkwrap algorithm every "shrinkwrap iterations"
  int shrinkwrap_iterations = 50;

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
  read_image(support_file_name, support, nx, ny);

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

  //create the planar CDI object which will be used to
  //perform the reconstuction.
  PlanarCDI planar(object_estimate, 4);

  //set the support and intensity
  planar.set_support(support,false);
  planar.set_intensity(data);

  //Initialise the current object ESW with a random numbers
  //"0" is the seed to the random number generator
  planar.initialise_estimate(0);

  //planar.set_fftw_type(FFTW_ESTIMATE);
  //read_cplx("Planar_trans.cplx", object_estimate);

  //make a 2D object. This will be used to output the 
  //image of the current estimate.
  Double_2D result(nx,ny);
  object_estimate.get_2d(MAG,result);

  /******* for fun, let's get the autocorrelation *****/

  Double_2D autoc(nx,ny);
  planar.get_intensity_autocorrelation(autoc);
  write_image("test_autocorrelation.ppm", autoc, true); //"true" means log scale

  //set the algorithm to Error Reduction
  planar.set_algorithm(ER);

  //  ProfilerStart("profile");


  /*** run the reconstruction ************/
  for(int a=0; a<cycles; a++){

    for(int i=0; i<er_iterations1; i++){

      cout << "iteration " << i+a*(hio_iterations+er_iterations1+er_iterations2) << endl;

      //apply the set of planar CDI projections 
      planar.iterate(); 
      cout << "Current error is "<<planar.get_error()<<endl;

      //every "output_iterations" 
      //output the current estimate of the object
      if(i%output_iterations==0){

	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "planar_example_iteration_" << std::setw(3) << std::setfill('0') << i +a*(hio_iterations+er_iterations1+er_iterations2) << ".tiff";
	write_image(temp_str.str(), result);

	temp_str.clear();

	//uncomment to output the estimated diffraction pattern
	//planar.propagate_to_detector(object_estimate);
	//object_estimate.get_2d(MAG_SQ,result);
	//temp_str << "diffraction.ppm";
	//write_ppm(temp_str.str(), result, true);
	//planar.propagate_from_detector(object_estimate);
	//object_estimate.get_2d(MAG,result);

	//apply the shrinkwrap algorithm
	//1.5 is the gaussian width in pixels
	//0.1 is the threshold (10% of the maximum pixel).

      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	planar.apply_shrinkwrap(2.0,0.1);

    }

    //now change to the error reduction algorithm 
    planar.set_algorithm(HIO);

    for(int i=er_iterations1; i<(hio_iterations+er_iterations1); i++){

      cout << "iteration " << i+a*(hio_iterations+er_iterations1+er_iterations2) << endl;

      planar.iterate(); 

      cout << "Current error is "<<planar.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "planar_example_iteration_" << std::setw(3) << std::setfill('0') << i+a*(hio_iterations+er_iterations1+er_iterations2) << ".tiff";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//planar.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	planar.apply_shrinkwrap(2.0,0.1);

    }

    planar.set_algorithm(ER);

    for(int i=er_iterations1+hio_iterations; i<(hio_iterations+er_iterations1+er_iterations2); i++){

      cout << "iteration " << i+a*(hio_iterations+er_iterations1+er_iterations2) << endl;

      planar.iterate();

      cout << "Current error is "<<planar.get_error()<<endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << "planar_example_iteration_" << std::setw(3) << std::setfill('0')<< i+a*(hio_iterations+er_iterations1+er_iterations2) << ".tiff";
	write_image(temp_str.str(), result);
	temp_str.clear();

	//apply the shrinkwrap algorithm
	//planar.apply_shrinkwrap(1.5,0.1);
      }
      if(i%shrinkwrap_iterations==(shrinkwrap_iterations-1))
	planar.apply_shrinkwrap(2.0,0.1);

    }
  }


  //And we are done. "object_estimate" contained the final estimate of
  //the ESW.
  write_cplx("Planar_trans.cplx", object_estimate);


  return 0;
}

