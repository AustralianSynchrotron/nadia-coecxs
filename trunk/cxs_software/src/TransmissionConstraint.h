#ifndef TRANS_H
#define TRANS_H

#include <vector>

class Double_2D;
class Complex_2D;

/**
 * @param phase_sign The sign of the phase to be fliped. ie. If you
 * pass the function 1 is will flip all positive phases to
 * negative. If you give it -1, it will flip all negative phases to
 * positive.
 *
 */

class ComplexConstraint{

 protected:
  double alpha1;
  double alpha2;
  Double_2D * region;
  double c_mean;
  bool fixed_c;

 public:
  ComplexConstraint(Double_2D & region,
		   double alpha1,
		   double alpha2){

    this->alpha1 = alpha1;
    this->alpha2 = alpha2;
    this->region = &region;
    fixed_c = false;
  };

  void set_fixed_c(double c_mean){
    fixed_c = true;
    this->c_mean = c_mean;
  };

  void set_c_mean(double c_mean){
    if(!fixed_c)
      this->c_mean = c_mean;
  };

  double get_c_mean(){
    return c_mean;
  }

  void set_alpha1(double alpha1){
    this->alpha1 = alpha1;
  };

  void set_alpha2(double alpha2){
    this->alpha2 = alpha2;
  };

  double get_new_mag(double old_mag, double old_phase){
    return exp((1-alpha1)*log(old_mag) + alpha1*c_mean*old_phase);
  };

  double get_new_phase(double old_mag, double old_phase){
    return (1-alpha2)*old_phase + alpha2*log(old_mag)/c_mean;
  };

  Double_2D * get_region(){
    return region;
  };

};


class TransmissionConstraint{

 protected:
  std::vector<ComplexConstraint*> complex_constraint_list;
  Double_2D * region_map;
  bool do_enforce_unity;
  bool do_charge_flip;
  int flip_sign;

  /** a function pointer to a cumstomized complex contraint */
  void (*custom_constraint)(Complex_2D&); 

 public:  
  
  TransmissionConstraint();

  ~TransmissionConstraint();

  virtual void add_complex_constraint(ComplexConstraint & new_constraint);

  void set_charge_flipping(bool enable, int flip_sign=1){
    do_charge_flip = enable;
    if(flip_sign>0)
      this->flip_sign=1;
    else
      this->flip_sign=-1;
  };

  void set_enforce_unity(bool enable){
    do_enforce_unity = enable;
  };

  void set_custom_constraint(void (*custom_constraint)(Complex_2D & tranmission)){
    this->custom_constraint = custom_constraint;
  };

  /** This is the fundamental method in this class */
  virtual void apply_constraint(Complex_2D & transmission);

};


#endif

