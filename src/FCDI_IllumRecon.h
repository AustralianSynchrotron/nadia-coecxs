#ifndef FCDI_ILLUMRECON_H
#define FCDI_ILLUMRECON_H

#include <map>
#include <string>

//forward declarations
class Complex_2D;
class FourierT;

/**
 * @file FCDI_IllumRecon.h
 * @class FCDI_IllumRecon
 * @author  Nadia Davidson <> 
 *
 * @brief  The class which performs planar CDI reconstuction.
 *
 * The fundamental class used for performing planar CDI reconstuction.
 * A projection object should be created with Complex_2D object which
 * you wish the projection to act upon. The Complex_2D should already
 * be filled with an inital guess, or the
 * FCDI_IllumRecon::initialise_estimate() method can be used to set up an
 * inital guess with random numbers inside the support and zero
 * outside.
 *
 * The support and intensity in the detector plane must be set using
 * the approproate setter methods, and the algorithm and relaxation
 * parameters may also be set if something other than the default it
 * required. The FCDI_IllumRecon::set_custom_algorithm allows the
 * recostruction algorithm to be customised by the user. See page 22
 * of the H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
 * imaging using short wavelength light sources, Journal of Modern
 * Optics, 2010, DOI: 10.1080/09500340.2010.495459 for the formalism.
 *
 * The reconstruction can be the performed with the
 * FCDI_IllumRecon::iterate() method.  This is automatically update the
 * values in the Complex_2D estimate and the current error. The error
 * can be retrieved using FCDI_IllumRecon::get_current_error(). The
 * algorithm, relaxation parameter, support and intensity can be reset
 * at any point during reconstuction.
 */
class FCDI_IllumRecon{

 private:

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

  /** the difference between the itensity in the detector plane 
      the current estimated intensity. */
  double current_error;

  double zone_to_focal_length;

  double wavelength;

  double focal_to_detector_length;

  double pixel_length;

  Complex_2D * forward_coefficients;
  Complex_2D * backward_coefficients;

 public:
  
  FCDI_IllumRecon(Complex_2D * initial_guess,
		  double beam_wavelength,
		  double zone_focal_length,
		  double focal_detector_length,
		  double pixel_size);
  ~FCDI_IllumRecon();
  
  int iterate();
  
  /**  Complex_2D * get_current_result(){
    return complex;
    };**/

  void set_support(double ** object_support);

  void set_intensity(double ** detector_intensity);

  void initialise_estimate(int seed=0);
  
  double ** get_intensity_autocorrelation();

  double get_error();

 private:
  
  void project_support(Complex_2D * c);
  void project_intensity(Complex_2D * c);
  void scale_intensity(Complex_2D * c);
  
};

#endif
