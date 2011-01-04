#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "shrinkwrap.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //read the data block in the file
  int nx = 1024;
  int ny = 1024;
  Complex_2D wf(nx,ny);
  Complex_2D wf2(nx,ny);
  int status;
  Double_2D data;

  status = read_dbin("angelas_images/scan_029498.bin", nx, ny, data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  //read the data into an array
  status = read_cplx("angelas_images/white_field.bin", wf);

  Double_2D temp(nx,ny);
  wf.get_2d(MAG_SQ,temp);

  double data_sum = data.get_sum();
  double wf_sum = temp.get_sum();

  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }


  /******* get the support from file and read it into an array *****/

  Double_2D support;
  //  int nx_s, ny_s;

  //status = read_tiff("image_files/FCDI_support_B.tiff", support);  
  status = read_ppm("angelas_images/support.pgm", support);  
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
  FresnelCDI proj(object_estimate,
		  wf,
		  6.88795e-10,
		  0.466795+0.001573152, //check these numbers
		  0.001573152,
		  13.5e-6,
		  //0.72088);
		  data_sum/wf_sum); 
		  //0.984729833);
 
  //set the support and intensity
  proj.set_support(support);
  
  wf.get_2d(MAG,result);

  proj.set_intensity(data);//result);
  //set the algorithm to hybrid input-output
  proj.set_algorithm(ER);

  //Initialise the current object ESW with a random number
  proj.initialise_estimate(0);
  

  /*** run the reconstruction ************/

  for(int i=0; i<50; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
    proj.iterate(); 
    cout << "Error: " << proj.get_error() << endl;

    if(i%5==0){
      //output the current estimate of the object
      ostringstream temp_str ( ostringstream::out ) ;
      object_estimate.get_2d(MAG,result);
      temp_str << "fcdi_example_iter_" << i << ".ppm";
      write_ppm(temp_str.str(),result);
      temp_str.clear();
    
      //apply the shrinkwrap algorithm
      //proj.apply_shrinkwrap(2,0.1);
      //write_ppm("shrink.ppm", nx, ny, result);
    }
  }

  Complex_2D trans(nx,ny);
  proj.get_transmission_function(trans);

  trans.get_2d(MAG,result);
  write_ppm("trans_mag.ppm",result,true);
  trans.get_2d(PHASE,result);
  write_ppm("trans_phase.ppm",result);

  return 0;
}
