/**
 * @file FresnelCDI.h
 * @class FresnelCDI
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief  The class which performs Fresnel CDI reconstruction.
 *
 * The class used for performing Fresnel CDI reconstruction (for
 * white-field reconstruction see FresnelCDI_WF). It inherits most
 * methods from PlanarCDI, so please look at the documentation of this
 * class also. Although there are some differences in the underlying
 * code between this class and the planar case, the interface is
 * generally unchanged. Therefore users should refer to the
 * instructions for PlanarCDI to understand how to use a FresnelCDI
 * object in their own code. Only the differences relevant to users
 * will be documented here.
 *
 */

#ifndef FCDI_H
#define FCDI_H

#include "PlanarCDI.h"

//forward declarations
class Complex_2D;
class FresnelCDI: public PlanarCDI{

 protected:

  /** the reconstructed white field */
  Complex_2D illumination;

  /** the white field in the sample plane */ 
  Complex_2D * illumination_at_sample;

  /** the transmission function */ 
  Complex_2D * transmission;

  /** the beam wavelength */
  double wavelength;

  /** the pixel size */
  double pixel_length;

  double focal_detector_length;
    
  double focal_sample_length;

  /** normalisation between the illumination with and without the
      sample placed in the beam */
  double norm;

  /** an array which holds a constants we use when propagating between
  difference planes */
  Complex_2D B_s; 

  /** an array which holds a constants we use when propagating between
  difference planes */
  //Complex_2D B_d;
  
 public:

  /** 
   * Create a FresnelCDI object. All of the dimensions which are given as
   * input should be in the same units. e.g metres.
   *
   * @param initial_guess The complex 2-D field which is modified by
   * this object. This represents the exit-surface-wave of the
   * sample. It may be pre-initialised to a best first guess (e.g
   * manually); it may be loaded from a file (e.g. the output from a
   * previous reconstruction job); or it may be unitialised, and the
   * default initialisation provided by FresnelCDI may be used.
   * @param white_field The reconstructed white-field (which has been
   * previously determined using FresnelCDI_WF or otherwise)
   * @param beam_wavelength Wavelength of the beam
   * @param focal_detector_length The distance from the focal point to
   * the detector
   * @param focal_sample_length The distance from the focal point to
   * the sample
   * @param pixel_size The diameter of one of the detector pixels.
   * @param normalisation The amount to scale the data before subtracting
   * the white field illumination. By default this is 1.0.
   *
   */
  FresnelCDI(Complex_2D & initial_guess,
	     Complex_2D & white_field,
	     double beam_wavelength,
	     double focal_detector_length,
	     double focal_sample_length,
	     double pixel_size,
	     double normalisation=1.0,
	     int n_best=0
	     );


  /**
   * The destructor for this class
   */
  virtual ~FresnelCDI();


  /**
   * Initialise the object esw (at the sample plane). The estimate is
   * constructed by subtracting the white-field magnitude from the
   * diffraction data. This forms the magnitude of the esw in the
   * detector plane. The phase is set randomly and will lie within the
   * range -0.1*2pi to 0.1*2pi. The esw is then propagated to the
   * sample plane and the support is applied.
   *
   * @param seed The seed of the random number generator 
   */
  virtual void initialise_estimate(int seed=0);

  virtual void auto_set_norm();

  /**
   * Get the transmission function from the current estimate of the
   * ESW. The ESW is divided by the illuminating wavefield in the
   * sample plane and unity is added.
   *
   * @param result The transmission function is copied into "result".
   * @param inforce_unit_mag Values above 1 (due to artifacts from the
   * reconstruction are reduced to 1). This helps when plotting the
   * magnitude. By default this is not done, but can be switched on by 
   * passing "true" for this parameter.
   */
  virtual void get_transmission_function(Complex_2D & result,
					 Complex_2D * esw = 0,
					 bool inforce_unit_mag=false);


  virtual void apply_support(Complex_2D & complex);

  virtual void set_transmission_function(Complex_2D & transmission,
					 Complex_2D * esw = 0);

  
  /**
   * This method overrides the one in PlanarCDI by adding/subtracting
   * the white-field before/after applying the intensity constraint.
   *
   * @param c The Complex_2D object to apply the intensity constraint
   * on.
   */
  virtual void scale_intensity(Complex_2D & c); 


  /**
   * Propagate to the sample plane using the paraxial free-space
   * equation.
   * @param c The Complex_2D field to propagate
   */
  virtual void propagate_from_detector(Complex_2D & c);

  /**
   * Propagate to the detector plane using the paraxial free-space
   * equation.
   * @param c The Complex_2D field to propagate
   */
  virtual void propagate_to_detector(Complex_2D & c);


  /**
   * Reset the normalisation factor which is used to scale the 
   * white-field data.
   * @param normalisation The new normalisation
   */
  void set_normalisation(double normalisation){
    norm = normalisation;
  };


   /**
    * Reset the experimental parameters. This is called when an object
    * of this type is constructed. I maybe called at any time during
    * reconstruction, for example to refine the experimental
    * parameters. This method resets the coefficient matricies which
    * are used in the propagation between the detector and sample
    * planes. All values should be given in the same length units.
    *
    * @param beam_wavelength The beam wavelength
    * @param focal_detector_length The distance between the focal plane
    *        and dector plane.
    * @param focal_sample_length The distance between the focal plane
    *        and sample plane.
    * @param pixel_size The size of one detector pixel.
    */ 
  void set_experimental_parameters(double beam_wavelength,
				   double focal_detector_length,
				   double focal_sample_length,
				   double pixel_size);
    


  void multiply_factors(Complex_2D & c, int direction);

};

#endif
