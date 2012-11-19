// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

//author:  Nadia Davidson <nadiamd@unimelb.edu.au>
//date last modified: 17/01/2011

/**
 * @file FresnelCDI_example.c
 *
 *
 * \a FresnelCDI_example.c - This file provides an example of running
 * the Fresnel CDI reconstruction on real data (from Corey's
 * data). The white-field result from the FresnelCDI_wf_example.c is
 * required as input. The error-reduction algorithm is used. The phase
 * and magnitude of the transmission function are output to file.
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //read the data block in the file
  int nx = 1024;
  int ny = 1024;
  int status;
  Double_2D data;

  //read the diffraction data
  status = read_dbin("image_files/FCDI_data.dbin", nx, ny, data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }


  //read the reconstucted white field data
  Complex_2D wf(nx,ny);
  status = read_cplx("fwf_recovered.cplx", wf);
  if(!status){
    cout  << "Maybe you need to run ./FresnelCDI_WF_example.exe "
	 << "first... exiting"  << endl;
    //return(1);
  }

/*  Complex_2D illumination(nx, ny);
  status = read_cplx("fillumination_recovered.cplx", illumination);
  if(!status){
    cout  << "Maybe you need to run ./FresnelCDI_WF_example.exe "
      << "first... exiting"  << endl;
    return(1);
  }

*/

  /******* get the support from file and read it into an array *****/

  Double_2D support;
  status = read_tiff("image_files/FCDI_support.tiff", support);  
  if(!status){
    cout << "failed to get data from "
      <<".. exiting"  << endl;
    return(1);
  }
  if( support.get_size_x() != nx || support.get_size_y() != ny){
    cout << "dimensions of the support to not match ... exiting"  << endl;
    return(1);
  }

  /*******  set up the reconstuction *********************/

  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.

  Double_2D result(nx,ny);

  //create the projection object which will be used to
  //perform the reconstuction.
  Complex_2D object_estimate(nx,ny);
  FresnelCDI proj(object_estimate, //estimate
      wf, //white field 
      //illumination, //illumination
      4.892e-10, //wavelength
      0.909513 - 16.353e-3, //focal-to-detector
      18.513e-3 - 16.353e-3, //focal-to-sample
      13.5e-6, //pixel size
      .984729833); //normalisation between wf and image

  //set the support and intensity
  proj.set_support(support);

  proj.set_intensity(data);

  //set the algorithm to error reduction
  proj.set_algorithm(HIO);

  //Initialise the current object ESW
  proj.initialise_estimate(0);


  /*** run the reconstruction ************/

  for(int i=0; i<100; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
//    proj.iterate(); 
    cout << "Error: " << proj.get_error() << endl;

    //if(i%1==0){
      //output the current estimate of the object
      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,result);
      temp_str << "erfcdi_example_iter_" << i << ".ppm";
      write_ppm(temp_str.str(),result);
      temp_str.clear();

      //output the current estimate of the object
      Double_2D tmp(nx,ny);
      proj.get_illumination_at_sample().get_2d(MAG, tmp);
      ostringstream temp_stra ( ostringstream::out ) ;
      object_estimate.get_2d(PHASE,result);
      temp_stra << "fcdi_ill_example_iter_" << i << ".ppm";
      write_ppm(temp_stra.str(),data);//result);
      temp_str.clear();

      //apply the shrinkwrap algorithm
      proj.apply_shrinkwrap(2,0.1);
    //}
    proj.iterate(); 
  }

  //get the transmission function for the final iteration
  Complex_2D trans(nx,ny);
  proj.get_transmission_function(trans);

  //I can't see this without using log scale.
  trans.get_2d(MAG,result);
  write_ppm("fcdi_example_trans_mag.ppm",result);

  trans.get_2d(PHASE,result);
  write_ppm("fcdi_example_trans_phase.ppm",result);

  /* for(int i = 0; i< nx; i++){
     for(int j=0; j<ny; j++){
  //      printf("real: %lf \n",/*, imag: %lf, mag: %lf \n", result.get_real(i, j), result.get_imag(i, j), result.get(i, j));

  }}*/
  //trans.get_2d(PHASE,result);
  //write_ppm("fcdi_example_trans_phase.ppm",result);

  return 0;
}
