#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "PlanarCDI.h"
#include "PhaseDiverseCDI.h"
#include "io.h" //
#include <sstream>
#include <typeinfo>

using namespace std;


void get_center(Double_2D & mag, int * x, int * y){
  
  double offset=0;
  double mean_x = 0;
  double mean_y = 0;
  double total = 0;

  for(int i=0; i<mag.get_size_x(); i++){
    for(int j=0; j<mag.get_size_y(); j++){

      mean_x += fabs(mag.get(i,j)-offset)*i;
      mean_y += fabs(mag.get(i,j)-offset)*j;
      
      total += fabs(mag.get(i,j)-offset);

    }
  }

  *x = mean_x/total;
  *y = mean_y/total;


}




PhaseDiverseCDI::PhaseDiverseCDI(Complex_2D & object):
  object(object){
  iterations_per_cycle = 1;
};


PhaseDiverseCDI::~PhaseDiverseCDI(){

  while(single_result.size()!=0){
    Complex_2D * temp = single_result.back();
    single_result.pop_back();
    delete temp;
  }

}

void PhaseDiverseCDI::add_new_position(PlanarCDI * local,
				       Double_2D & beam_shape,
				       double x, 
				       double y,
				       double probe_scaling){

  int nx = local->get_size_x() ;
  int ny = local->get_size_y() ;
    
  singleCDI.push_back(local);
  x_position.push_back(x);
  y_position.push_back(y);
  alpha.push_back(probe_scaling);
  single_result.push_back(new Complex_2D(nx,ny));
  this->beam_shape.push_back(&beam_shape);

  //do some check of the dimensions.

  //fix the weight..
  weight.push_back(1);

  get_result(local,*(single_result.back()));
  //update_to_object(single_result.size()-1);


  position_refine.push_back(4);
  position_refine_stable_itr.push_back(0);

};

void PhaseDiverseCDI::initialise_estimate(){
  
  for(int i=0; i<object.get_size_x(); i++){
    for(int j=0; j<object.get_size_y(); j++){
      object.set_real(i,j,1);
    }
  }

  for(int i=0; i<singleCDI.size(); i++){
    update_to_object(i);
  }
  
}


void PhaseDiverseCDI::iterate(){

  static int total_iterations = 0;

  bool parallel = false;  
    
  Double_2D result(1024,1024);

  for(int i=0; i<singleCDI.size(); i++){

    int x, y;

    /**
    char buf[50];
    sprintf(buf,"result_%i_abf.tiff",i);
    single_result.at(i)->get_2d(PHASE,result);
    write_image(buf,result);
    get_center(result, &x, &y);
    cout << "data "<<i<< ", pos before is " <<x<<","<<y<<endl;
    **/

    if(total_iterations>1 && total_iterations<20 && i!=0)
      update_from_object_with_shift(i);
    else
      update_from_object(i);
    set_result(singleCDI.at(i),*(single_result.at(i)));

    /**
    sprintf(buf,"result_%i_af.tiff",i);
    single_result.at(i)->get_2d(PHASE,result);
    write_image(buf,result);
    get_center(result, &x, &y);
    cout << "data "<<i<< ", pos after is " <<x<<","<<y<<endl; 
    **/

    for(int j=0; j<iterations_per_cycle; j++){
      singleCDI.at(i)->iterate();
      cout << "This Error="<<singleCDI.at(i)->get_error()<<endl;
    }

    get_result(singleCDI.at(i),*(single_result.at(i)));

    if(!parallel)
      add_to_object(i,weight.at(i),1-weight.at(i));
    //update_to_object(i);

  }
  
  if(parallel){
    add_to_object(0,0,0.2);
    for(int i=0; i<singleCDI.size(); i++){
      add_to_object(i,0.2,1);
    }
  }

  total_iterations++;

};

void PhaseDiverseCDI::get_result(PlanarCDI * local, Complex_2D & result){
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->get_transmission_function(result);
  }
  else
    local->get_exit_surface_wave(result);
};

void PhaseDiverseCDI::set_result(PlanarCDI * local, Complex_2D & result){
  
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->set_transmission_function(result);
  }
  else
    local->set_exit_surface_wave(result);
  
}



void PhaseDiverseCDI::update_from_object_with_shift(int n_probe){
 
  double step_size = position_refine.at(n_probe);
  double stable_itrs = position_refine_stable_itr.at(n_probe);

  const int stable_itrs_cut_off = 1;

  if(step_size<1)
    return;

  double x = x_position.at(n_probe);
  double y = y_position.at(n_probe);

  PlanarCDI * single = singleCDI.at(n_probe);
  
  double best_x=x;
  double best_y=y;
  double best_error=100;
   
  
  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){

      x_position.at(n_probe)=x+i*step_size;
      y_position.at(n_probe)=y+j*step_size;
      update_from_object(n_probe);
      set_result(single,*(single_result.at(n_probe)));
      single->iterate();

      if(single->get_error()<best_error){
	best_error = single->get_error();
	best_x = x+i*step_size;
	best_y = y+j*step_size;
      }
    }
  }
  
  cout << "moving probe "<< n_probe << " by "<<best_x-x<<" in x "
       << "and "<< best_y-y<<" in y." << endl;


  //setting to the best one :
  x_position.at(n_probe)=best_x;
  y_position.at(n_probe)=best_y;
  update_from_object(n_probe);
  
  if(best_x-x == 0 && best_y-y ==0){
    if(stable_itrs > stable_itrs_cut_off){
      position_refine.at(n_probe)=position_refine.at(n_probe)/2.0;
      position_refine_stable_itr.at(n_probe)=0;
    }
    else
      position_refine_stable_itr.at(n_probe)=stable_itrs+1;
  }

  //singleCDI.at(n_probe)->iterate();
  //cout << "This error="<<singleCDI.at(n_probe)->get_error()<<endl;  



}


void PhaseDiverseCDI::add_to_object(int n_probe, double weight, 
				    double old_weight){
  
  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D & small = *(single_result.at(n_probe));
  Double_2D & beam = *beam_shape.at(n_probe);
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small.get_size_x(); i_++){
    for(int j_=0; j_< small.get_size_y(); j_++){

      i = i_-x_offset;
      j = j_-y_offset;

      if(beam.get(i_,j_)>0 && i>=0 && j>=0 && 
	 i<object.get_size_x() && j<object.get_size_y()){

	double new_real = weight*small.get_real(i_,j_)
	  +old_weight*object.get_real(i,j);
	
	double new_imag = weight*small.get_imag(i_,j_)
	  +old_weight*object.get_imag(i,j);
	
	object.set_real(i,j,new_real);
	object.set_imag(i,j,new_imag);
      }
    }
  }

}

void PhaseDiverseCDI::update_from_object(int n_probe){

  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D * small = single_result.at(n_probe);
  Double_2D * beam = beam_shape.at(n_probe);
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small->get_size_x(); i_++){
    for(int j_=0; j_< small->get_size_y(); j_++){

      i = i_-x_offset;
      j = j_-y_offset;
      
      if(beam->get(i_,j_)>0 && i>=0 && j>=0 && 
	 i<object.get_size_x() && j<object.get_size_y()){

	double new_real = object.get_real(i,j);
	  //	  +(1-weight)*object.get_real(i,j);

	double new_imag = object.get_imag(i,j);
	  //  +(1-weight)*object.get_imag(i,j);
	
	small->set_real(i_,j_,new_real);
	small->set_imag(i_,j_,new_imag);
      }
    }
  }

};

void PhaseDiverseCDI::update_to_object(int n_probe){

  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D * small = single_result.at(n_probe);
  double w = weight.at(n_probe);
  Double_2D * beam = beam_shape.at(n_probe);
  

  /**  Double_2D out_s(small.get_size_x(),
		  small.get_size_y());

  Double_2D out_o(object.get_size_x(),
  object.get_size_y());**/

  //round off to a first approx. Will be fixed later.
  
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small->get_size_x(); i_++){
    for(int j_=0; j_< small->get_size_y(); j_++){

      //      cout << beam.get(i_,j_) <<endl;

      if(beam->get(i_,j_)>0){

	i = i_-x_offset;
	j = j_-y_offset;
	
	double new_real = w*small->get_real(i_,j_)
	  +(1-w)*object.get_real(i,j);
	double new_imag = w*small->get_imag(i_,j_)
	  +(1-w)*object.get_imag(i,j);
	
	object.set_real(i,j,new_real);
	object.set_imag(i,j,new_imag);
      }
    }
  }


};
