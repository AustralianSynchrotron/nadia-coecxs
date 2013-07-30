// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.

/**
 * @file PhaseDiverse_example.c
 *
 * This example file shows how phase-diverse / ptychographic
 * reconstruction can be performed. The example uses data from Corey
 * Putkunz which can be found on osiris at
 * /data/cputkunz/phase_diverse_cdi/example_data.tar.gz
 * 
 *
 * @author T'Mir Julius <tdjulius@unimelb.edu.au> 
 *
 * Last modified on 26/07/2013
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <io.h>
#include <TransmissionConstraint.h>
#include <Complex_2D.h>
#include <Double_2D.h>
#include <FresnelCDI.h>
#include <PhaseDiverseCDI.h>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <utils.h>

using namespace std;

/**********************************/

int main(int argc, char * argv[]){

  ////////////////////////
  int nx = 1024;
  int ny = 1024;

  const int max_iterations = 15;

  //make the support
  Double_2D beam(nx,ny);
  double beam_fraction = 0.5;
  double i0 = (nx-1)/2.0;
  double j0 = (ny-1)/2.0;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(sqrt((i-i0)*(i-i0)+(j-j0)*(j-j0)) 
	 < beam_fraction*sqrt(j0*j0+i0*i0))
	beam.set(i,j,100);
      else{
	beam.set(i,j,0);
      }
    }
  }

  //Make a new PhaseDiverseCDI object
  PhaseDiverseCDI pd;

  //how many frames make up the reconstruction
  const int frames = 7;

  //initialise some arrays to hold data.
  FresnelCDI * proj[frames];
  Complex_2D * object_estimate[frames];
  Complex_2D * wf;
  Double_2D * diffraction;
  Double_2D * this_support;

  //set the experimental parameters (all are in meters)
  double wavelength = 4.892e-10; //wavelength
  double fz = 16.353e-3; //zone plate to focal distance
  double ps = 13.5e-6; //pixel size

  //zone plate to detector distances
  double fd[frames] = {0.909513, 0.909388, 0.909263, 0.909213,
		       0.909088, 0.909088, 0.909088};

  //sample to detector distances
  double fs[frames] = {18.513e-3, 18.388e-3, 18.263e-3, 18.213e-3,
		       18.088e-3, 18.088e-3, 18.088e-3};
  
  //white-field normalisation
  double norm[frames] = {0.984729833,0.97700431,0.986270638,0.967487825,
			 0.980945916, 0.97628279, 0.963066039 };

  
  //reconstructed white-field file names
  std::string wf_string[frames] = {"image_files/phasediverse_data/wf_A_1024.cplx",
				   "image_files/phasediverse_data/wf_B_1024.cplx",
				   "image_files/phasediverse_data/wf_C_1024.cplx",
				   "image_files/phasediverse_data/wf_D_1024.cplx",
				   "image_files/phasediverse_data/wf_E_1024.cplx",
				   "image_files/phasediverse_data/wf_F_1024.cplx",
				   "image_files/phasediverse_data/wf_G_1024.cplx"};

  //data file names
  std::string diff_string[frames] = {"image_files/phasediverse_data/A.dbin",
				     "image_files/phasediverse_data/B.dbin",
				     "image_files/phasediverse_data/C.dbin",
				     "image_files/phasediverse_data/D.dbin",
				     "image_files/phasediverse_data/E.dbin",
				     "image_files/phasediverse_data/F.dbin",
				     "image_files/phasediverse_data/G.dbin"};

  //these are the correct coordinates
  // double x_pos[frames] = {0,26,6,28,-14,34,38};
  // double y_pos[frames] = {-150,-149,-131,-127,-134,-95,211};

  //these are the coordinates smeared a bit
  double x_pos[frames] = {0,23,3,35,30,37,60};
  double y_pos[frames] = {-150,-145,-134,-123,-150,-93,210};

  TransmissionConstraint tc1;
  
  for(int n=0; n < frames; n++){

    //Set up the fresnel CDI in the same way you would
    //if you weren't doing phase-diversity
    object_estimate[n] = new Complex_2D(nx,ny);

    //read the white-field from file
    wf = new Complex_2D(nx,ny);
    read_cplx(wf_string[n], *wf);
    
    //read the data from file
    diffraction = new Double_2D(nx,ny);
    read_dbin(diff_string[n],nx,ny,*diffraction);
    
    //set-up the reconstruction for a single frame
    proj[n]= new FresnelCDI(*object_estimate[n],
			    *wf,
			    wavelength,
			    fd[n]-fz,
			    fs[n]-fz,
			    ps,
			    norm[n]);

    //set the support and intensity and initialise
    proj[n]->set_intensity(*diffraction);
    proj[n]->set_support(beam);
    proj[n]->initialise_estimate();

    //add a complex constraint
    proj[n]->set_complex_constraint(tc1);

    delete wf;
    delete diffraction;

    //New part.. Add the FresnelCDI to the PhaseDiverseCDI.
    pd.add_new_position(proj[n], x_pos[n], y_pos[n]);

  }

  //initialise the phase diverse transmission function
  pd.initialise_estimate();
  
  //now run the reconstruction
  for(int i=0; i< max_iterations; i++){

    //do some position alignment
    if(i==0)
      pd.adjust_positions(PhaseDiverseCDI::CROSS_CORRELATION);
    if(i==10)
      pd.adjust_positions(PhaseDiverseCDI::MINIMUM_ERROR);
    
    //iteration!
    pd.iterate();

  }


  Complex_2D * object = pd.get_transmission();

  Double_2D result(object->get_size_x(),object->get_size_y());

  //write the magnitude and phase of it to an image file
  object->get_2d(MAG,result);
  write_image("object_mag_recovered.tiff",result,false,0.4,1.0);

  object->get_2d(PHASE,result);
  write_image("object_phase_recovered.tiff",result,false, -1.2,0);
  
  //save the transmission function in case we want to use it later.
  write_cplx("trans.cplx",*object);


  return 0;
  
}
