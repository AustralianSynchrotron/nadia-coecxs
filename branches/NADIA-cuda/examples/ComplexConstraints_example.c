/**
 * @file ComplexConstraints_example.c
 *
 * This example aims at demonstrating how complex costraints can be
 * inforced during a reconstruction. The example requires that
 * FresnelCDI_WF_example.exe and FresnelCDI_simulation_example.c
 * have been previously run. 
 *
 * Here, the complex constraint is apply to the transmission function
 * during Fresnel reconstruction, but constraints such as charge
 * flipping can be achieve in the same way on th exit-surface wave for
 * plane-wave reconstruction.
 *
 * @author Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * Last modified on 25/3/2011
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "TransmissionConstraint.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include <cstdlib>
#include <math.h>
#include <string>

using namespace std;

//Example of a custom constraint function.
//See example 4 below
void my_charge_flipping(Complex_2D & transmission){

  double max_thickness = 150e-9;
  double delta =  6.45e-4;
  double wavelength = 4.892e-10;
  double k = (2.0 * M_PI)/wavelength;

  double min_phase = -delta*k*max_thickness;
  double max_phase = 0;

  //we know the phase can't be out by more than
  for(int i=0; i<transmission.get_size_x(); i++){
    for(int j=0; j<transmission.get_size_y(); j++){

      if(transmission.get_phase(i,j)>max_phase)
	transmission.set_phase(i,j,max_phase);

      if(transmission.get_phase(i,j)<min_phase)
	transmission.set_phase(i,j,min_phase);
    }
  }

}


/**********************************/

int main(int argc, char * argv[]){

  /** All of this part is just reading files and setting constants **/

  //number of output iterations
  int output_iterations = 10;
  int total_iterations = 90;

  //load the diffraction image of the sample 
  //(made from the simulation example)
  Double_2D diffraction;
  int status = read_tiff("forward_projection.tiff", diffraction);
  if(!status){ //give an error if we couldn't open the file.
    cout  << "Maybe you need to run "
      << "./FresnelCDI_simulation_example.exe "
      << "first... exiting"  << endl;
    return(1);
  }

  //get the object dimensions
  int nx = diffraction.get_size_x();
  int ny = diffraction.get_size_y();

  //we will use the same white field as the Fresnel examples. For
  //this reason you will need to reconstruct it first using
  //FresnelCDI_WF_example.exe
  Complex_2D wf(nx,ny);
  status = read_cplx("wf_recovered.cplx", wf);
  if(!status){ //give an error if we couldn't open the file.
    cout  << "Maybe you need to run ./FresnelCDI_WF_example.exe "
      << "first... exiting"  << endl;
    return(1);
  }

  //load the support to use in the reconstruction
  Double_2D support;
  read_image("image_files/FCDI_support.tiff", support);
  //make sure the object and support are the same dimensions
  if(support.get_size_x() != nx || support.get_size_y() != ny){
    cout << "dimensions do not match....exiting" <<endl;
    return(1);
  }


  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current exit-surface-wave estimate.
  Double_2D result(nx,ny);

  //set the experimental parameters (all are in meters)
  double wavelength = 4.892e-10; //wavelength
  double fd = 0.8932; //focal to detector
  double fs = 2.5e-3; //focal to sample
  double ps = 13.5e-6; //pixel size

  //create the Complex_2D object which will be used to store
  //the reconstructed exit-surface-wave.
  Complex_2D object_estimate(nx,ny);


  /** 
   * set up the reconstuction 
   * - this part is the same as all the other examples
   */
  FresnelCDI proj(object_estimate, //estimate
      wf, //white field 
      wavelength, //wavelength
      fd, //focal-to-detector
      fs, //focal-to-sample
      ps //pixel size
      ); //normalisation of white-field to sample data

  //set the support and intensity
  proj.set_support(support);  
  proj.set_intensity(diffraction);
  proj.auto_set_norm();

  proj.initialise_estimate();

  //set the algorithm
  proj.set_algorithm(ER);


  /******************************************************************
   **          now set setup the complex contraint!                **
   ******************************************************************/

  //create a transmission constraint object. Here you
  //can set various constraints on the transmission function. 
  //For example charge-flipping, enforcing unity on the magnitude of
  //the transmission function, setting a constraint on the
  //refractive index of the material, custom constraints etc.

  //Example 1:
  //----------
  //Charge-flipping and unity enforced
  //This is set by default, so nothing extra is required
  TransmissionConstraint tc1;

  //Example 2:
  //----------
  //Charge-flipping and unity enforcment turned on
  //A complex constraint is set on the object.
  TransmissionConstraint tc2;

  //load an image which will be used to define the region
  //that the complex constraint is applied to.
  //in this case we just used the original object since the image
  //file is available.
  Double_2D area;
  read_image("image_files/FCDI_support.tiff", area);
  //image_files/FCDI_simulation_object.tiff",area);  
  //read_image("blah.tiff",area);  

  //create a ComplexConstraint using the constraint region we just
  //loaded from fial. All pixel values > 0 will be interpreted as part
  //of the constraint region. alpha1 = 1.0 and alpha2 = 0.0 
  //(see the J.Clark paper on complex contraints for information 
  //about these)
  ComplexConstraint c1(area,1.0,0.0);

  //add the complex constraint to the transmission constraint. You may
  //add a many complex constraint as you like (e.g. if the sample is 
  //composed of multiple materials).
  tc2.add_complex_constraint(c1);


  //Example 3:
  //----------
  //The same as Example 2, but this time we enforce a particular value
  //of c = beta/delta
  TransmissionConstraint tc3;

  ComplexConstraint c2(area,1.0,0.0);  
  double delta = 6.45e-4;
  double beta = 1.43e-4;

  tc3.add_complex_constraint(c2);    
  c2.set_fixed_c(beta/delta);

  //Example 4:
  //-----------
  //Setting a custom contraint on the transmission function.
  //The TransmissionConstraint class also allows a custom user
  //constraint to be invoked instead of or in addition to
  //the constraints provided. Do do this, you need to pass the
  //function name as a parameters to the object. The function
  //should be implemented in your code and need to be of the
  //form void function_name(Complex_2D & array);

  TransmissionConstraint tc4;
  tc4.set_enforce_unity(true);
  tc4.set_charge_flipping(false);
  tc4.set_custom_constraint(my_charge_flipping);


  //------------------------------------------------

  //Now set the contraint to use (here we use Example 2), but
  //you can change to another to see what happens.
  proj.set_complex_constraint(tc4);
  //just change tc2 in the line above.

  //------------------------------------------------


  /*** run the reconstruction ************/

  for(int i=0; i<  total_iterations+1; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
    proj.iterate(); 
    cout << "Error: " << proj.get_error() << endl;

    if(i%output_iterations==0){

      //output the current estimate of the object
      object_estimate.get_2d(MAG,result);
      ostringstream temp_str ( ostringstream::out ) ;
      temp_str << "fcdi_example_iter_" << i << ".tiff";
      write_tiff(temp_str.str(),result);

    }
  }

  //get the reconstructed transmission function for 
  //the final iteration
  Complex_2D trans(nx,ny);
  proj.get_transmission_function(trans);

  //and write the magnitude and phase of it to an image file
  trans.get_2d(MAG,result);
  write_tiff("trans_mag_recovered.tiff",result);

  trans.get_2d(PHASE,result);
  write_tiff("trans_phase_recovered.tiff",result);

  //print out the value of c = beta/delta recovered
  //(only useful if we ran example 2.
  cout << "The value of beta/delta recovered was " << c1.get_c_mean() << endl;
  cout << "The true value of beta/delta was " << c2.get_c_mean() << endl;


  return 0;

}

