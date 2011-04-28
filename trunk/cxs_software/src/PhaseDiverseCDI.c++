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


/**void get_center(Double_2D & mag, int * x, int * y){
  
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

  }**/




PhaseDiverseCDI::PhaseDiverseCDI(
				 double beta, 
				 double gamma, 
				 bool parallel,
				 int nx,
				 int ny,
				 int min_x,
				 int min_y,
				 int granularity):scale(granularity),
						  beta(beta), 
						  gamma(gamma), 
						  parallel(parallel),
						  nx(nx),
						  ny(ny),
						  x_min(min_x),
						  y_min(min_y),
						  weights_set(false),
                                                  iterations_per_cycle(1),
                                                  object(0){
  
  if(nx!=0 && ny!=0)
     reallocate_object_memory(nx,ny);
  
};


PhaseDiverseCDI::~PhaseDiverseCDI(){

  while(single_result.size()!=0){
    Complex_2D * temp = single_result.back();
    single_result.pop_back();
    delete temp;

    Double_2D * temp_d = weights.back();
    weights.pop_back();
    delete temp_d;
  }

  if(object)
    delete object;
  
}

void PhaseDiverseCDI::reallocate_object_memory(int new_nx,int new_ny){

  if(object)
    delete object;

  object = new Complex_2D(new_nx,new_ny);
  
  nx = new_nx;
  ny = new_ny;

}

Complex_2D * PhaseDiverseCDI::get_transmission(){
  return object;
}

void PhaseDiverseCDI::set_transmission(Complex_2D & new_transmission){

  int new_nx = new_transmission.get_size_x();
  int new_ny = new_transmission.get_size_y();

  if(!object)
    reallocate_object_memory(new_nx,new_ny);

  if(new_nx!=nx||new_ny!=ny){
    cout << "Can not set the transmission function in "
	 << "PhaseDiverseCDI because the complex array "
	 << "given does not have the same dimensions as "
	 << "the current transmission function"<<endl;
    return;
  }

  object->copy(new_transmission);
  
}

void PhaseDiverseCDI::add_new_position(PlanarCDI * local,
				       double x, 
				       double y,
				       double alpha){
  
  weights_set = false;


  int lnx = local->get_size_x() ;
  int lny = local->get_size_y() ;
  
  singleCDI.push_back(local);
  x_position.push_back(x);
  y_position.push_back(y);
  this->alpha.push_back(alpha);
  single_result.push_back(new Complex_2D(lnx,lny));
  weights.push_back(new Double_2D(lnx,lny));

  cout << "Added position "<<singleCDI.size()-1<<endl;
  
  if(!object){
    x_min = -x;
    y_min = -y;
    reallocate_object_memory(lnx,lny);
  }

  //dynamically increase the object size to fit this frame in

  //what is the position (globally) of the first 
  //pixel if we don't increase the frame size  
  int extra_x=0;
  int extra_y=0;

  int global_x_min = -x-x_min;
  int global_y_min = -y-y_min;
  int global_x_max = lnx-x-x_min;
  int global_y_max = lny-y-y_min;

  if(global_x_min<0){
    x_min = -x;
    extra_x += -global_x_min;
  }
  if(global_y_min<0){
    y_min = -y;
    extra_y += -global_y_min;
  }
  if(global_x_max>nx)
    extra_x += global_x_max-nx;
  if(global_y_max>ny)
    extra_y += global_y_max-ny;

  if(extra_x || extra_y)
    reallocate_object_memory(nx+extra_x, ny+extra_y);

}




void PhaseDiverseCDI::initialise_estimate(){

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      object->set_real(i,j,1);
      object->set_imag(i,j,0);
    }
  }

  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i);
  }

}


void PhaseDiverseCDI::iterate(){

  static int total_iterations = 0;

  cout << "Iteration "<<total_iterations<<endl;

  bool parallel = true;  
    
  Double_2D result(1024,1024);

  for(int i=0; i<singleCDI.size(); i++){

    int x, y;

    //if(i!=0 && (total_iterations==2 || total_iterations==4))
    //  check_position(i);
    update_from_object(i);

    for(int j=0; j<iterations_per_cycle; j++){
      singleCDI.at(i)->iterate();
      cout << "This Error="<<singleCDI.at(i)->get_error()<<endl;
    }

    if(!parallel)
      add_to_object(i);

  }
  
  if(parallel){

    object->scale(1-beta);
       
    //  new_real += (1-beta)*object->get_real(i,j);
    //  new_imag += (1-beta)*object->get_imag(i,j);

    for(int i=0; i<singleCDI.size(); i++){
      add_to_object(i);
    }
  }

  total_iterations++;

};

void PhaseDiverseCDI::get_result(PlanarCDI * local, Complex_2D & result){
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->get_transmission_function(result);
  }
  else
    result.copy(local->get_exit_surface_wave());
};

void PhaseDiverseCDI::set_result(PlanarCDI * local, Complex_2D & result){
  
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->set_transmission_function(result);
  }
  else
    local->set_exit_surface_wave(result);
  
}



//frames need to be in order for this to work.
void PhaseDiverseCDI::adjust_positions(double step_size, bool forward){

  //make a copy of the tranmission function
  Complex_2D * pointer_to_object = object;
  int nx_c = nx;
  int ny_c = ny;
  double x_min_c = x_min; 
  double y_min_c = y_min; 
  bool parallel_c = parallel;  
  
  parallel = false;

  //leave the first frame, as everything will be aligned to it.
  int n=1;
  int limit=singleCDI.size();
  if(!forward){
    n=singleCDI.size()-2;
    limit=-1;
  }

  while(n!=limit){
    //for(int n=1; n<singleCDI.size(); n++){

    //store some before information for later use
    double before_x = x_position.at(n);
    double before_y = y_position.at(n);

    //replace the transmission function with
    //just the 1st frame information

    Complex_2D temp(1024,1024);
    object = &temp;
    nx = 1024;
    ny = 1024;
    x_min = -x_position.at(n);
    y_min = -y_position.at(n);
    
    for(int i=0; i<1024; i++){
      for(int j=0; j<1024; j++){
	temp.set_real(i,j,1);
	temp.set_imag(i,j,0);
      }
    }

    /**    char name[80];
    Double_2D pic(1024,1024);
    object->get_2d(MAG,pic);
    write_image("pic_0.tiff",pic,false,0,1);**/

    //if frame m overlaps with frame n
    //add it to the current array

    if(forward){
      for(int m=0; m < n; m++){
	add_to_object(m);
      }
    }
    else{
      for(int m=singleCDI.size()-1; m > n; m--){
	add_to_object(m);
      }
    }

    //    char name[80];
    //Double_2D pic(1024,1024);
    /**    object->get_2d(MAG,pic);
    sprintf(name,"pic_%i.tiff",n);
    write_image(name,pic,false, 0,1);**/

    check_position(n,step_size,0);

    cout << "moving prob " << n << " from ("
	 << before_x << "," << before_y << ") "
	 << "-> (" << x_position.at(n) 
	 << "," <<  y_position.at(n) << ")."
	 << endl;

    if(forward)
      n++;
    else
      n--;

  }
  
  object = pointer_to_object;
  nx = nx_c;
  ny = ny_c;
  x_min = x_min_c; 
  y_min = y_min_c; 
  parallel = parallel_c;

}

int PhaseDiverseCDI::check_position(int n_probe, double step_size, int tries){
  
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

  Complex_2D * object_copy = new Complex_2D(object->get_size_x(),object->get_size_y());
  object_copy->copy(*object);
  Complex_2D * single_copy = new Complex_2D(size_x,size_y);
  single_copy->copy(*single_result.at(n_probe));


  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){

      new_x=x+i*step_size;
      new_y=y+j*step_size;

      //checking whether the corrds are still within the 
      //boundary of the image
      //      if(new_x>0 && new_x<size_x && new_y>0 && new_y<size_y){

      x_position.at(n_probe)=new_x;
      y_position.at(n_probe)=new_y;
      
      add_to_object(n_probe);
      update_from_object(n_probe);

      /**     static int h = 0;
      char blah[80];
      Double_2D temp(1024,1024);
      single_result.at(n_probe)->get_2d(PHASE,temp);
      sprintf(blah,"blah_%i.tiff",h);
      write_image(blah,temp);
      h++;**/

      single->iterate();
	
      //      cout << "error at ("<<new_x<<","<<new_y
      //	   << ")"<<" is "<<single->get_error()<<endl;

      if(single->get_error()<best_error){
	best_error = single->get_error();
	best_x = x+i*step_size;
	best_y = y+j*step_size;
      }

      //reset the object estimate:
      object->copy(*object_copy);
      set_result(single,*single_copy);
      single_result.at(n_probe)->copy(*single_copy);

    }
  }
  delete object_copy;
  delete single_copy;

  //set x and y to the best one :
  x_position.at(n_probe)=best_x;
  y_position.at(n_probe)=best_y;

  //recursively call this function with a smaller step size.
  if(best_x==x && best_y==y)
    step_size=step_size/2.0;
  
  int status = check_position(n_probe, step_size, ++tries);

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



void PhaseDiverseCDI::set_up_weights(){
  
  if(weights_set)
    return;
  else
    weights_set=true;

  for(int n_probe=0; n_probe<singleCDI.size(); n_probe++){
    
    int nx = single_result.at(n_probe)->get_size_x();
    int ny = single_result.at(n_probe)->get_size_y();
    PlanarCDI * this_CDI = singleCDI.at(n_probe);
    double value;
    Double_2D illum_mag(nx,ny);
    double max;

    if(typeid(*this_CDI)==typeid(FresnelCDI)){
      const Complex_2D & illum = ((FresnelCDI*)(this_CDI))->get_illumination_at_sample();
      illum.get_2d(MAG,illum_mag);
      max = illum_mag.get_max();
    }
    
    for(int i=0; i< nx; i++){
      for(int j=0; j< ny; j++){
	
	if(typeid(*this_CDI)==typeid(FresnelCDI))
	  value = illum_mag.get(i,j)/max;
	else
	  value = 1;

	value=beta*alpha.at(n_probe)*pow(value,gamma);
	weights.at(n_probe)->set(i,j,value);
      }
    } 
  }
  
  if(parallel){
    
    Double_2D weight_norm(object->get_size_x(),object->get_size_y());
    
    for(int n=0; n<singleCDI.size(); n++){
      
      double x = x_position.at(n);
      double y = y_position.at(n);
      
      int this_size_x = single_result.at(n)->get_size_x();
      int this_size_y = single_result.at(n)->get_size_y();
      
      
      for(int i_=0; i_< this_size_x; i_++){
	for(int j_=0; j_< this_size_y; j_++){
	  
	  int i = (i_-x-x_min)*scale;
	  int j = (j_-y-y_min)*scale;
	  
	  if(i>=0&&j>=0&&i<nx&&j<ny){
	    
	    double new_weight = weights.at(n)->get(i_,j_);
	    double old_weight = weight_norm.get(i,j);
	    
	    if(n==0)
	      weight_norm.set(i,j,new_weight);
	    else
	      weight_norm.set(i,j,new_weight+old_weight);
	  }
	}
      }
    }
    
    for(int n=0; n<singleCDI.size(); n++){
      
      double x = x_position.at(n);
      double y = y_position.at(n);
      
      int this_size_x = single_result.at(n)->get_size_x();
      int this_size_y = single_result.at(n)->get_size_y();
      
      for(int i_=0; i_< this_size_x; i_++){
	for(int j_=0; j_< this_size_y; j_++){
	  
	  int i = (i_-x-x_min)*scale;
	  int j = (j_-y-y_min)*scale;
	  
	  if(i>=0&&j>=0&&i<nx&&j<ny){
	    
	    double old_weight = weights.at(n)->get(i_,j_);
	    double norm = weight_norm.get(i,j);

	    if(norm<=0)
	      weights.at(n)->set(i_,j_,0);
	    else
	      weights.at(n)->set(i_,j_,old_weight/norm);

	  }
	}
      }
    } 
    //    write_image("weight_norm.tiff",weight_norm);
 
  }

  //  write_image("weights_2.tiff",*(weights.at(2)));
 };


  


void PhaseDiverseCDI::add_to_object(int n_probe){
  
  get_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));
  
  set_up_weights();
  
  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D & small = *(single_result.at(n_probe));

  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small.get_size_x(); i_++){
    for(int j_=0; j_< small.get_size_y(); j_++){
      
      i = (i_-x_offset-x_min)*scale;
      j = (j_-y_offset-y_min)*scale;
      
      if(i>=0&&j>=0&&i<nx&&j<ny){
	
	double weight = weights.at(n_probe)->get(i_,j_);
	double new_real = weight*small.get_real(i_,j_);
	double new_imag = weight*small.get_imag(i_,j_);
	
	if(!parallel){
	  new_real += (1-weight)*object->get_real(i,j);
	  new_imag += (1-weight)*object->get_imag(i,j);
	}
	else{
	  new_real += object->get_real(i,j);
	  new_imag += object->get_imag(i,j);
	}
	
	update_to_object_sub_grid(i,j,new_real,new_imag); 
      }
      
    }
  } 
  
}

void PhaseDiverseCDI::update_to_object_sub_grid(int i, int j, 
			       double real_value, 
			       double imag_value){
  
  for(int di=0; di < scale; di++){
    for(int dj=0; dj < scale; dj++){
      object->set_real(i+di,j+dj,real_value);
      object->set_imag(i+di,j+dj,imag_value);
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
      real_value+=object->get_real(i+di,j+dj);
      imag_value+=object->get_imag(i+di,j+dj);
    }
  }

  real_value=real_value/(scale*scale);
  imag_value=imag_value/(scale*scale);

}


void PhaseDiverseCDI::update_from_object(int n_probe){

  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  Complex_2D * small = single_result.at(n_probe);

  
  //round off to a first approx. Will be fixed later.
  int i, j; //local coors (i,j) are global (object) coors

  //using the local coordinate system
  for(int i_=0; i_< small->get_size_x(); i_++){
    for(int j_=0; j_< small->get_size_y(); j_++){

      i = (i_-x_offset-x_min)*scale;
      j = (j_-y_offset-y_min)*scale;
      
      //      if(beam->get(i_,j_)>0 && i>=0 && j>=0 && 
      //	 i<object->get_size_x() && j<object->get_size_y()){

      if(i>=0 && j>=0 && i<nx && j<ny){

	double new_real; //= object->get_real(i,j);
	double new_imag; //= object->get_imag(i,j);
	
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

