/**
 * @file PhaseDiverse_example.c
 *
 *
 * @author Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * Last modified on 15/4/2011
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
#include "PhaseDiverseCDI.h"
#include <cstdlib>
#include <math.h>
#include <string>

#include <google/heap-profiler.h>

using namespace std;


void my_constraint(Complex_2D & object){

  //get the reconstructed transmission function for 
  //the final iteration
  
  int x = object.get_size_x();
  int y = object.get_size_y();

  for(int i=0; i<x; i++){
    for(int j=0; j<y; j++){

      /**      
      if(object.get_mag(i,j)>1)
	object.set_mag(i,j,1);
      
      if(object.get_mag(i,j)<0)
	object.set_mag(i,j,0);
      
      if(object.get_phase(i,j)>0)
	object.set_phase(i,j,-1*object.get_phase(i,j));
      
      if(object.get_phase(i,j)<-3.1415/2)
      object.set_phase(i,j,-3.1415/2);     
      **/

      double amp = object.get_mag(i,j);
      double phase = object.get_phase(i,j);

      if(amp>0){

	amp = log(amp);
	double c = fabs(amp/phase);

	if(c<0.05)
	  c=0.05;
	if(c>0.6)
	  c=0.6;

	double alpha2 = 0.5*c+0.5*c*c;
	
	object.set_mag(i,j,exp(c*phase));
	object.set_phase(i,j,(1-alpha2)*phase + alpha2*amp/c);
	
      }
    }
  }
}

/**********************************/

int main(int argc, char * argv[]){

  //  HeapProfilerStart("my_heap");

  //get the object dimensions
  int nx = 1024;
  int ny = 1024;
  
  double i0 = (nx-1)/2.0;
  double j0 = (ny-1)/2.0;

  Double_2D support(nx,ny);
  read_image("support_phase_diverse.tiff",support); //"image_files/FCDI_support.tiff",support);

  double beam_fraction = 0.5;
  Double_2D beam(nx,ny);
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


  Complex_2D object((nx+100), (ny+100));


  PhaseDiverseCDI pd(object,1,3,true);

  const int frames = 7;

  FresnelCDI * proj[frames];
  Complex_2D * object_estimate[frames];
  Complex_2D * wf;
  Double_2D * diffraction;
  Double_2D * this_support;

  //set the experimental parameters (all are in meters)
  double wavelength = 4.892e-10; //wavelength
  double fz = 16.353e-3;
  double ps = 13.5e-6; //pixel size
  double fd[frames] = {0.909513, 0.909388, 0.909263, 0.909213,
		       0.909088, 0.909088, 0.909088};

  double fs[frames] = {18.513e-3, 18.388e-3, 18.263e-3, 18.213e-3,
		       18.088e-3, 18.088e-3, 18.088e-3};
  //double fs[frames] = {18.088e-3, 18.388e-3, 18.263e-3, 18.213e-3};

  double norm[frames] = {0.984729833,0.97700431,0.986270638,0.967487825,
			 0.980945916, 0.97628279, 0.963066039 };
  //double norm[frames] = {0.981258,   0.977237,  0.978268,   0.970263};


  std::string wf_string[frames] = {"wf_A_1024.cplx",
				   "wf_B_1024.cplx",
				   "wf_C_1024.cplx",
				   "wf_D_1024.cplx",
				   "wf_E_1024.cplx",
				   "wf_F_1024.cplx",
				   "wf_G_1024.cplx"};

  std::string diff_string[frames] = {"A.dbin",
				     "B.dbin",
				     "C.dbin",
				     "D.dbin",
				     "E.dbin",
				     "F.dbin",
				     "G.dbin"};

  double x_pos[frames] = {0,26,6,28,-13,35,39};
  double y_pos[frames] = {0,-1,-19,-23,-16,-55,-362};

  TransmissionConstraint tc1;
  tc1.set_custom_constraint(my_constraint);

  for(int n=0; n < frames; n++){
    
    object_estimate[n] = new Complex_2D(nx,ny);
    wf = new Complex_2D(nx,ny);
    read_cplx(wf_string[n], *wf);
    
    diffraction = new Double_2D(nx,ny);
    read_dbin(diff_string[n],nx,ny,*diffraction);

    proj[n]= new FresnelCDI(*object_estimate[n],
			    *wf,
			    wavelength,
			    fd[n]-fz,
			    fs[n]-fz,
			    ps,
			    norm[n]);

    //set the support and intensity
    proj[n]->set_intensity(*diffraction);

    //make a temporary support which is shifted in x and y
    this_support = new Double_2D(nx,ny);

    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
	int i_ = i-x_pos[n];
	int j_ = j-y_pos[n];
	
	if(i_>0&&i_<nx&&j_>0&&j_<ny)
	  this_support->set(i,j,support.get(i_,j_));
      }
    }

    proj[n]->set_support(*this_support, true);  
    
    //    proj[n]->auto_set_norm();
    proj[n]->initialise_estimate();
    proj[n]->set_complex_constraint(tc1);

    //    if(n==0||n==1||n==2||n==6)
    pd.add_new_position(proj[n], x_pos[n],y_pos[n]); 
      
    delete wf;
    delete diffraction;
    delete this_support;

  }

  pd.initialise_estimate();

  //  HeapProfilerDump("mine1");

  for(int i=0; i< 100; i++){
    //proj[0]->iterate();
    pd.iterate();
    //    cout << proj[0]->get_error() << endl;
  }

  Double_2D result_2(1024,1024);
  object_estimate[0]->get_2d(MAG,result_2);
  write_image("temp.ppm",result_2);

  HeapProfilerDump("mine2");

  //  my_constraint(temp);
  Double_2D result(object.get_size_x(),object.get_size_y());

  //and write the magnitude and phase of it to an image file
  object.get_2d(MAG,result);
  write_image("object_mag_recovered.tiff",result,false,0.15,1.0);

  object.get_2d(PHASE,result);
  write_image("object_phase_recovered.tiff",result,false, -1.5,0);

  //  GetHeapProfile();
  // HeapProfilerStop();

  for(int n=0; n < 1; n++){
    delete proj[n];
    delete object_estimate[n];
  }


  return 0;
  
}
