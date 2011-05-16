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

//#include <google/heap-profiler.h>

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

  //  Double_2D support(nx,ny);
  //  read_image("support_phase_diverse.tiff",support); //"image_files/FCDI_support.tiff",support);




  //  Complex_2D object(nx+1000,ny+600);

  PhaseDiverseCDI pd(1,1,true);

  const int frames = 1; //72;
  const int total_frames = 441;

  FresnelCDI * proj[frames];
  Complex_2D * object_estimate[frames];
  Complex_2D * wf;
  Double_2D * diffraction;
  //Double_2D * this_support;

  //set the experimental parameters (all are in meters)
  double wavelength = 6.88795e-10; //wavelength
  double fz = 11.8e-3;
  double ps = 13.5e-6; //pixel size


  std::string wf_string = "white_field.bin";
  wf = new Complex_2D(nx,ny);
  read_cplx(wf_string, *wf);

  std::string diff_string[frames];
  int y_pos[frames];
  int x_pos[frames];

  char str[80];
  char diff_image_name[100];
  FILE * pFile = fopen ("angela_data_list_large.txt","r");


  int n=0;
  for(int i=0; i< total_frames; i++){
    int x_p = 0;
    int y_p = 0;

    fscanf (pFile, "%s %s %s %i %i", diff_image_name, str, str, &x_p, &y_p);

    if( y_p==3296 && x_p ==2952){ //y_p<2000 && y_p > -700 && x_p>500){ //500
      x_pos[n]=x_p;//+10*n;
      y_pos[n]=y_p;//-10*n;
      diff_string[n] = diff_image_name;
      cout << "file name "<<diff_string[n] <<" x pos="<<x_pos[n]<<" y pos= " << y_pos[n]<<endl;
      n++;
    } 
  }
  fclose (pFile);


  double beam_fraction = 0.5;
  Double_2D beam(nx,ny);
  
  for(int n=0; n < frames; n++){
    
    object_estimate[n] = new Complex_2D(nx,ny);
    
    diffraction = new Double_2D(nx,ny);
    read_dbin(diff_string[n],nx,ny,*diffraction);

    proj[n]= new FresnelCDI(*object_estimate[n],
			    *wf,
			    wavelength,
			    0.4918-fz,
			    0.0135-fz,
			    ps,
			    0.984729833);

    //set the support and intensity
    proj[n]->set_intensity(*diffraction);

    //make a temporary support which is shifted in x and y
    //    this_support = new Double_2D(nx,ny);

    /**    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
	int i_ = i-x_pos[n];
	int j_ = j-y_pos[n];
	
	if(i_>0&&i_<nx&&j_>0&&j_<ny)
	  this_support->set(i,j,support.get(i_,j_));
      }
      }**/

    proj[n]->get_illumination_at_sample().get_2d(MAG,beam);
    write_image("illumination_at_sample.tiff",beam);
    double max = beam.get_max();

    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
	if( beam.get(i,j) > beam_fraction*max  )
	  beam.set(i,j,1);
	else{
	  beam.set(i,j,0);
	}
      }
    }
    
    
    write_image("beam.tiff",beam);
    
    proj[n]->set_support(beam, true);  

    //    proj[n]->auto_set_norm();
    proj[n]->initialise_estimate();

    //   if(n==0)//||n==2||n==6)
      pd.add_new_position(proj[n], x_pos[n],y_pos[n]); 
      
    delete diffraction;
    //  delete this_support;

  }

  TransmissionConstraint tc1;
  tc1.set_charge_flipping(false);
  tc1.set_enforce_unity(true);
  
  ComplexConstraint cc(beam,1.0,0);
  //tc1.add_complex_constraint(cc);
  //tc1.set_custom_constraint(my_constraint);

  for(int i=0; i< frames; i++)
    proj[i]->set_complex_constraint(tc1);

  pd.initialise_estimate();
  
  // pd.get_transmission()->get_2d(MAG,temp_initial_ptyc);
  //write_image("temp_initial_ptyc.tiff",temp_initial_ptyc,false, 0,1.0);

  /**  Complex_2D guess(3776,3776);
  read_cplx("object_small.cplx",guess);
  pd.set_transmission(guess);**/

  //  HeapProfilerDump("mine1");
  pd.iterate();

  for(int i=0; i<15; i++){
    

    Double_2D temp_initial_ptyc(nx,ny);
    Complex_2D temp1(nx,ny);
    proj[0]->get_transmission_function(temp1);
    temp1.get_2d(MAG,temp_initial_ptyc);
    char buff[80];
    sprintf(buff,"temp_itr_%i_single.tiff",i);
    write_image(buff,temp_initial_ptyc,false,0,1.0);
    
    proj[0]->iterate();
    //    pd.iterate();
    cout << proj[0]->get_error() << endl;

    //    if(i==2||i==6){//||i==50){
    //   pd.adjust_positions(8);
    //   pd.adjust_positions(8,false);
    // }
  }

  
  Complex_2D trans(1024,1024);
  Double_2D result_2(1024,1024);
  object_estimate[0]->get_2d(MAG,result_2);
  write_image("temp0.ppm",result_2);

  proj[0]->get_transmission_function(trans);
  trans.get_2d(MAG,result_2);
  write_image("temp0_trans_mag.ppm",result_2,false,0,1.0);
  trans.get_2d(PHASE,result_2);
  write_image("temp0_trans_phase.ppm",result_2,false,-3,0);

  /**  object_estimate[1]->get_2d(MAG,result_2);
  write_image("temp1.ppm",result_2);

  object_estimate[2]->get_2d(MAG,result_2);
  write_image("temp2.ppm",result_2);**/


  //  HeapProfilerDump("mine2");

  //  my_constraint(temp);
  Complex_2D * object = pd.get_transmission();
  Double_2D result(object->get_size_x(),object->get_size_y());

  //and write the magnitude and phase of it to an image file
  object->get_2d(MAG,result);
  write_image("temp_mag.tiff",result,false,0.15,1.0); //"object_mag_5_pos_fixing.tiff",result,false,0.15,1.0);

  object->get_2d(PHASE,result);
  write_image("temp_phase.tiff",result,false,-1.0,0); //"object_phase_5_pos_fixing.tiff",result,false, -1.0,0);

  cout << "nx = "<<object->get_size_x()<<endl;
  cout << "ny = "<<object->get_size_y()<<endl;

  write_cplx("temp.cplx",*object); //"object_5_pos_fixing.cplx",*object);

  //  GetHeapProfile();
  // HeapProfilerStop();

  
  for(int f=0; f < frames; f++){
    cout << diff_string[f] <<" "
	 << "white_field.bin para.txt  "
	 << pd.get_final_x_position(f) << " "
	 << pd.get_final_y_position(f) << endl;
  };

  cout << "Difference with the orginal positions: " <<endl;
  for(int f=0; f < frames; f++){
    cout << x_pos[f] - pd.get_final_x_position(f) << " "
	 << y_pos[f] - pd.get_final_y_position(f) << endl;
  };
  

  for(int n=0; n < 1; n++){
    delete proj[n];
    delete object_estimate[n];
  }
  delete wf;


  return 0;
  
}
