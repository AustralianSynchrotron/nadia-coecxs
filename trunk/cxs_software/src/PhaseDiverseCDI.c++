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




PhaseDiverseCDI::PhaseDiverseCDI(Complex_2D & object, int granularity):
  object(object), scale(granularity){
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

  //get_result(local,*(single_result.back()));
  //update_to_object(single_result.size()-1);

};

void PhaseDiverseCDI::initialise_estimate(){
  
  for(int i=0; i<object.get_size_x(); i++){
    for(int j=0; j<object.get_size_y(); j++){
      object.set_real(i,j,1);
    }
  }

  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i,0.5,0.5);
  }
  
}


void PhaseDiverseCDI::iterate(){

  static int total_iterations = 0;

  cout << "Iteration "<<total_iterations<<endl;

  bool parallel = true;  
    
  Double_2D result(1024,1024);

  for(int i=0; i<singleCDI.size(); i++){

    int x, y;

    if(i!=0 && total_iterations%7==6 && total_iterations<20)
      update_from_object_with_shift(i);
    update_from_object(i);

    /**    char buf[50];
    sprintf(buf,"result_%i_abf.tiff",i);
    single_result.at(i)->get_2d(PHASE,result);
    write_image(buf,result);**/
    
    for(int j=0; j<iterations_per_cycle; j++){
      singleCDI.at(i)->iterate();
      cout << "This Error="<<singleCDI.at(i)->get_error()<<endl;
    }

    /**    sprintf(buf,"result_%i_af.tiff",i);
    single_result.at(i)->get_2d(PHASE,result);    
    write_image(buf,result);**/

    //   get_center(result, &x, &y);
    // cout << "data "<<i<< ", pos before is " <<x<<","<<y<<endl;


    /**   single_result.at(i)->get_2d(PHASE,result);
    get_center(result, &x, &y);
    cout << "data "<<i<< ", pos after is " <<x<<","<<y<<endl; **/

    if(!parallel)
      add_to_object(i,weight.at(i),1-weight.at(i));

  }

  if(parallel){
    double fraction = 1.0/single_result.size();
    cout << fraction <<endl;
    add_to_object(0,fraction,0);
    for(int i=0; i<singleCDI.size(); i++){
      add_to_object(i,fraction,1);
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



int PhaseDiverseCDI::update_from_object_with_shift(int n_probe, double step_size, int tries){
  
  double x = x_position.at(n_probe);
  double y = y_position.at(n_probe);

  cout << "checking probe "<< n_probe << " position, " 
       << x << "," << y << " with step size " 
       << step_size << ". Try no.: "<<tries<< endl;

  //done, we found the best position.
  if(step_size<1.0/scale)
    return SUCCESS;

  //failed, we moved around a bit, but couldn't find a local minima
  //in the error metric.
  if(tries>10){
    cout << "giving up on probe "<< n_probe << ". Coule not find " 
	 << "it's position. Returning to the original. " 
	 << endl;
      return FAILURE;
  }
  
  PlanarCDI * single = singleCDI.at(n_probe);
  
  double best_x=x;
  double best_y=y;
  double best_error=100;
   
  double new_x;
  double new_y;
  double size_x = single_result.at(n_probe)->get_size_x();
  double size_y = single_result.at(n_probe)->get_size_y();
  
  //try the 9 positions around the current one
  //record the one with the lowest error metric.
  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){

      new_x=x+i*step_size;
      new_y=y+j*step_size;

      //checking whether the corrds are still within the 
      //boundary of the image
      //      if(new_x>0 && new_x<size_x && new_y>0 && new_y<size_y){

	x_position.at(n_probe)=new_x;
	y_position.at(n_probe)=new_y;
	
	update_from_object(n_probe);
	set_result(single,*(single_result.at(n_probe)));
	
	single->iterate();
	
	if(single->get_error()<best_error){
	  best_error = single->get_error();
	  best_x = x+i*step_size;
	  best_y = y+j*step_size;
	}
	  // }
    }
  }
  

  //set x and y to the best one :
  x_position.at(n_probe)=best_x;
  y_position.at(n_probe)=best_y;

  //recursively call this function with a smaller step size.
  if(best_x==x && best_y==y)
    step_size=step_size/2.0;
  
  int status = update_from_object_with_shift(n_probe, step_size, ++tries);

  if(status == FAILURE){
    //return to orginal coordinates
    x_position.at(n_probe) = x ;
    y_position.at(n_probe) = y;
    return FAILURE;
  }

  
  cout << "moving probe "<< n_probe << " by " 
       << x_position.at(n_probe)-x <<" in x "
       << "and "<< y_position.at(n_probe)-y<<" in y." 
       << endl;
  
  return SUCCESS;
  
}


void PhaseDiverseCDI::add_to_object(int n_probe, double weight, 
				    double old_weight){

  get_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));
  
  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D & small = *(single_result.at(n_probe));
  Double_2D & beam = *beam_shape.at(n_probe);
  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small.get_size_x(); i_++){
    for(int j_=0; j_< small.get_size_y(); j_++){

      i = (i_-x_offset)*scale;
      j = (j_-y_offset)*scale;

      //      cout << "i,j="<<i<<","<<j<<endl;

      if(beam.get(i_,j_)>0 && i>=0 && j>=0 && 
	 i<object.get_size_x() && j<object.get_size_y()){

	double new_real = weight*small.get_real(i_,j_)
	  +old_weight*object.get_real(i,j);
	
	double new_imag = weight*small.get_imag(i_,j_)
	  +old_weight*object.get_imag(i,j);
	
	update_to_object_sub_grid(i,j,new_real,new_imag); 

	//object.set_real(i,j,new_real);
	//object.set_imag(i,j,new_imag);
      }
    }
  }
}

void PhaseDiverseCDI::update_to_object_sub_grid(int i, int j, 
			       double real_value, 
			       double imag_value){
  
  for(int di=0; di < scale; di++){
    for(int dj=0; dj < scale; dj++){
      object.set_real(i+di,j+dj,real_value);
      object.set_imag(i+di,j+dj,imag_value);
    }
  }
}

void PhaseDiverseCDI::update_from_object_sub_grid(int i, int j, 
				 double & real_value, 
				 double & imag_value){
  
  real_value=0;
  imag_value=0;

  for(int di=0; di < scale; di++){
    for(int dj=0; dj < scale; dj++){
      real_value+=object.get_real(i+di,j+dj);
      imag_value+=object.get_imag(i+di,j+dj);
    }
  }

  real_value=real_value/(scale*scale);
  imag_value=imag_value/(scale*scale);

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

      i = (i_-x_offset)*scale;
      j = (j_-y_offset)*scale;
      
      if(beam->get(i_,j_)>0 && i>=0 && j>=0 && 
	 i<object.get_size_x() && j<object.get_size_y()){

	double new_real; //= object.get_real(i,j);
	double new_imag; //= object.get_imag(i,j);
	
	update_from_object_sub_grid(i, j, 
				    new_real, 
				    new_imag);

	small->set_real(i_,j_,new_real);
	small->set_imag(i_,j_,new_imag);
      }
    }
  }

  set_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));

};

