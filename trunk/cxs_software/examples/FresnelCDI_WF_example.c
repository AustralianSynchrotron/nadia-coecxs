//author:  Nadia Davidson <nadiamd@unimelb.edu.au>
//last modified: 17/01/2011

/**
 * @file FresnelCDI_WF_example.c
 *
 * \a FresnelCDI_WF_example.c - This example demonstrates how the
 * phase of a white field can be recovered and saved for use in
 * Fresnel CDI reconstruction. The data set comes from Corey.
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI_WF.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //read the detector data block in the file
  int nx = 1024;
  int ny = 1024;
  int status;
  Double_2D data;
  status = read_dbin("image_files/FCDI_wf_data.dbin", nx, ny, data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  /******* get the support from file and read it into an array *****/
  Double_2D support;
  status = read_tiff("image_files/FCDI_wf_support.tiff", support);  
  if(!status){
    cout << "failed to get data from "
	 <<".. exiting"  << endl;
    return(1);
  }
  //check that the dimensions are the same as the data
  if( support.get_size_x() != nx || support.get_size_y() != ny){
    cout << "dimensions of the support to not match ... exiting"  << endl;
    return(1);
  }


  /*******  set up the reconstuction *********************/

  //create the white field (wf) object which will store
  //the reconstucted illumination at the detector plane.
  Complex_2D wf(nx,ny); 
  
  //make the projection operator
  FresnelCDI_WF proj(wf, //inital estimate
		     4.892e-10, //wavelength
		     16.353e-3, //zone-to-focal length
		     0.909513 - 16.353e-3, //focal-to-detector
		     13.5e-6); //pixel size
 
  //set the support and intensity
  proj.set_support(support);
  
  proj.set_intensity(data);
  
  //no need to set the algorithm this time.
  //since there is only one for white field reconstruction.

  //initialise the current detector wavefield 
  //check FresnelCDI_WF.h to see what it is initialised to.
  proj.initialise_estimate(0);
  
  /*** run the reconstruction ************/

  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.
  Double_2D result(nx,ny);

  for(int i=0; i<15; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
    proj.iterate(); 
    cout << "Error: " << proj.get_error() << endl;

    if(i%5==0){
      //output the current estimate at the detector
      //(magnitude and phase)
      ostringstream temp_str ( ostringstream::out ) ;
      wf.get_2d(MAG,result);
      temp_str << "fcdi_wf_example_mag_iter_" << i << ".ppm";
      write_ppm(temp_str.str(),result);

      ostringstream temp_str2 ( ostringstream::out ) ;
      wf.get_2d(PHASE,result);
      temp_str2 << "fcdi_wf_example_phase_iter_" << i << ".ppm";
      write_ppm(temp_str2.str(),result);
    }
  }

  //save the result as a complex binary for later use
  status = write_cplx("wf_recovered.cplx", wf);
  if(status!=0)
    cout << "Successfully wrote out the reconstucted white field"
	 << " into wf_recovered.cplx"<<endl;


  return 0;
}
