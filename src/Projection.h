#ifndef PROJECTION_H
#define PROJECTION_H

#define NALGORITHMS 9
#define NTERMS 5
#define NVECTORS 10


enum { PSF, PFS, PS, PF, PI };
enum { ER, BIO, BOO, HIO, DM, SF, ASR, HPR, RAAR};

class Complex_2D;
class FourierT;

class Projection{

 private:

  FourierT * fft;
  Complex_2D * complex;

  double ** support;
  double ** intensity_sqrt;
  double beta;
  double algorithm_structure[NTERMS];
  int algorithm;
  double current_error;
  // Complex_2D * x[NTERMS-1]; //one less since we don't need to 
  //contruct a temporary copy for the identity matrix

  Complex_2D * temp_complex_PFS;
  Complex_2D * temp_complex_PF;

 public:
  
  Projection(Complex_2D * complex);
  ~Projection();
  
  int iterate();
  
  Complex_2D * get_current_result(){
    return complex;
  };

  void set_support(double ** object_support);

  void set_intensity(double ** detector_intensity);

  void set_relaxation_parameter(double relaxation_parameter){
    beta = relaxation_parameter;
    set_algorithm(algorithm); //force algorithm constants to update
  };

  void initialise_estimate(double seed=0);
  
  double ** get_intensity_autocorrelation();


  void set_algorithm(int alg);
  void set_custom_algorithm(double m1, double m2, double m3, 
			    double m4, double m5, double m6, 
			    double m7, double m8,
			    double m9, double m10);
  
  void print_algorithm();
  double get_chi2_error(Complex_2D * c);
  double get_current_error();


 private:
  
  void project_support(Complex_2D * c);
  void project_intensity(Complex_2D * c);
  void scale_intensity(Complex_2D * c);
  bool is_support(int x, int y);
    
};

#endif
