/**
 * @file FresnelCDI.h
 * @class FresnelCDI
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief  The class which performs Fresnel CDI reconstruction.
 *
 * The class used for performing Fresnel CDI reconstruction (for
 * white-field reconstruction see FresnelCDI_WF). It inherits most
 * methods from BaseCDI, so please look at the documentation of this
 * class also. Although there are some differences in the underlying
 * code between this class and the planar case, the interface is
 * generally unchanged. Therefore users should refer to the
 * instructions for BaseCDI to understand how to use a FresnelCDI
 * object in their own code. Only the differences relevant to users
 * will be documented here.
 *
 */

#ifndef PHASED_H
#define PHASED_H

#include "BaseCDI.h"
#include <vector>

//forward declarations
class Complex_2D;
class PhaseDiverseCDI{

 protected:

  std::vector<BaseCDI*> singleCDI;
  std::vector<Complex_2D*> single_result;
  std::vector<Double_2D *> weights;
  std::vector<double> x_position;
  std::vector<double> y_position;
  
  //parameters controlling the feedback
  double beta;
  double gamma;
  std::vector<double> alpha;

  Complex_2D * object;

  int iterations_per_cycle;
  int scale;
  int nx,ny;
  
  int x_min;
  int y_min;
  bool parallel; 
  bool weights_set;

 public:

  PhaseDiverseCDI(double beta=1.0, 
		  double gamma=1.0,
		  bool parallel=false,
		  int granularity=1);
  
  ~PhaseDiverseCDI();

  void iterate();
  
  void add_new_position(BaseCDI * local, 
			double x=0, double y=0, 
			double alpha=1);
 
  void set_iterations_per_cycle(int iterations){
    iterations_per_cycle = iterations;
  }

  void set_feedback_parameter(double beta){
    this->beta = beta;    
    weights_set = false;

  };

  void set_amplification_factor(double gamma){
    this->gamma = gamma;
    weights_set = false;
 
  };

  void set_probe_scaling(int n_probe, double alpha){
    if(this->alpha.size() > n_probe)
      this->alpha.at(n_probe)=alpha;
    else{
      std::cout << "In PhaseDiverseCDI::set_probe_scaling, "
	   << "the probe number given is too large. "
	   << "Are you trying to set the value of alpha "
		<< "before calling add_new_position?" <<std::endl;
    }
    weights_set = false;

  };

  void initialise_estimate();
  
  Complex_2D * get_transmission();
  
  void set_transmission(Complex_2D & new_transmission);
  
  void adjust_positions(double step_size=4, bool forwards=true);
  
  double get_final_x_position(int n_probe){
    return x_position.at(n_probe);}
  
  double get_final_y_position(int n_probe){
    return y_position.at(n_probe);}



 private:

  void add_to_object(int n_probe);
  int check_position(int n_probe, double shift=4, int tries=0);
  void update_from_object(int n_probe);
  void get_result(BaseCDI * local, Complex_2D & result);
  void set_result(BaseCDI * local, Complex_2D & result);
  /**  void update_to_object_sub_grid(int i, int j, 
				 double real_value, 
				 double imag_value);
  void update_from_object_sub_grid(int i, int j, 
				 double & real_value, 
				 double & imag_value);**/

  void reallocate_object_memory(int new_nx,int new_ny);

  void scale_object(double factor);
  
  void get_object_sub_grid(Complex_2D & result,
			   double x_offset,
			   double y_offset);
  void set_up_weights();


  inline int get_global_x_pos(int x, double x_offset){
    return (x-x_offset-x_min); //*scale;
  }
 
  inline int get_global_y_pos(int y, double y_offset){
    return (y-y_offset-y_min);
  }

  inline int get_local_x_pos(int x, double x_offset){
    return x + x_offset + x_min;
  }

  inline int get_local_y_pos(int y, double y_offset){
    return y + y_offset + y_min; //y/scale + y_offset + y_min;
  }


};


#endif
