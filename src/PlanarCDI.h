#ifndef PHASERETRIEVALBASE_H
#define PHASERETRIEVALBASE_H

#include <map>
#include <string>

/** The number of reconstruction algorithms */
#define NALGORITHMS 9 

/** The number of unique combination projectors (i.e. terms) */
#define NTERMS 5

/** The number of unique combination of terms (ie. vectors) */
#define NVECTORS 10

/** The terms in the iteration equation */
enum { PSF, PFS, PS, PF, PI }; 

/** The different algorithms to choose from */
enum { ER, BIO, BOO, HIO, DM, SF, ASR, HPR, RAAR};


//forward declarations
class Complex_2D;
class FourierT;

/**
 * @file PlanarCDI.h
 * @class PlanarCDI
 * @author  Nadia Davidson <> 
 *
 * @brief  The class which performs planar CDI reconstuction.
 *
 * The fundamental class used for performing planar CDI reconstuction.
 * A projection object should be created with Complex_2D object which
 * you wish the projection to act upon. The Complex_2D should already
 * be filled with an inital guess, or the
 * PlanarCDI::initialise_estimate() method can be used to set up an
 * inital guess with random numbers inside the support and zero
 * outside.
 *
 * The support and intensity in the detector plane must be set using
 * the approproate setter methods, and the algorithm and relaxation
 * parameters may also be set if something other than the default it
 * required. The PlanarCDI::set_custom_algorithm allows the
 * recostruction algorithm to be customised by the user. See page 22
 * of the H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
 * imaging using short wavelength light sources, Journal of Modern
 * Optics, 2010, DOI: 10.1080/09500340.2010.495459 for the formalism.
 *
 * The reconstruction can be the performed with the
 * PlanarCDI::iterate() method.  This is automatically update the
 * values in the Complex_2D estimate and the current error. The error
 * can be retrieved using PlanarCDI::get_current_error(). The
 * algorithm, relaxation parameter, support and intensity can be reset
 * at any point during reconstuction.
 */
class PlanarCDI{

 protected:

  /**a fourier transform objected used to perform 
     the forward and backward transformations */
  FourierT * fft; 

  /**a pointer to the Complex_2D object which is 
     altering during each iterations */
  Complex_2D * complex;

  /** a copy of the support */
  double ** support;

  /**a copy of the square root of the intensity at the detector plane. */
  double ** intensity_sqrt;

  /**the relaxation parameter */
  double beta;

  /**an array which holds the structure of the algorithm. i.e. the 
     coefficients for each term in the iteration equation.*/
  double algorithm_structure[NTERMS];

  /** the algorithm which is being used: ER, HIO etc. */
  int algorithm;

  /** the difference between the itensity in the detector plane 
      the current estimated intensity. */
  double current_error;

  /** temporary Complex_2Ds which are used in the computation of the
      PFS and PF terms for each iteration. */
  Complex_2D * temp_complex_PFS;
  Complex_2D * temp_complex_PF;

  static std::map<std::string,int> * algNameMap;

 public:
  
  PlanarCDI(Complex_2D * complex);
  ~PlanarCDI();
  
  virtual int iterate();
  virtual void initialise_estimate(int seed=0);
  
  /**  Complex_2D * get_current_result(){
    return complex;
    };**/

  void set_support(double ** object_support);

  void set_intensity(double ** detector_intensity);

  void set_relaxation_parameter(double relaxation_parameter){
    beta = relaxation_parameter;
    set_algorithm(algorithm); //force algorithm constants to update
  };


  
  double ** get_intensity_autocorrelation();

  void set_algorithm(int alg);

  void set_custom_algorithm(double m1, double m2, double m3, 
			    double m4, double m5, double m6, 
			    double m7, double m8,
			    double m9, double m10);
  
  void print_algorithm();

  double get_error();

  static int getAlgFromName(std::string algorithm_name){
    std::map<std::string,int>::iterator alg;
    alg = algNameMap->find(algorithm_name);
    if(alg == algNameMap->end() )
      return -1;
    return (alg->second);
  }

 protected:
  
  virtual void apply_support(Complex_2D * c);
  virtual void project_intensity(Complex_2D * c);
  virtual void scale_intensity(Complex_2D * c);
  static std::map<std::string,int> * set_up_algorithm_name_map();
  //bool is_support(int x, int y);
    
};

#endif
