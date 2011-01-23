/**
 * @file PlanarCDI.h
 * @class PlanarCDI
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief  The class which performs planar CDI reconstruction.
 *
 * The fundamental class used for performing CDI
 * reconstruction; an object of type "PlanarCDI" allows planar CDI to
 * be performed, and it is also used by the Fresnel reconstruction
 * code (through inheritance).
 *
 * The reconstruction can be performed in the following way: 
 * <ol> 
 *
 * <li> The user should create an object to store the result (the exit
 * surface wave of type Complex_2D). For example using:
 * <br><kbd>Complex_2D esw_estimate(nx,ny)</kbd> <br>This can be
 * initialise to a best guess, or can be left for the PlanarCDI object
 * to initialise.</li>
 * 
 * <li>The PlanarCDI object should be created and the Complex_2D
 * previously created should be passed as input. For example:
 * <br><kbd>PlanarCDI my_planar_recon(esw_estimate)</kbd></li>
 *
 * <li>The support and intensity data should be set, and the exit
 * surface wave initialised to a first guess (if not already done):
 * <br><kbd>my_planar_recon.set_support(support) </kbd>
 * <br><kbd>my_planar_recon.set_intensity(data) </kbd>
 * <br><kbd>my_planar_recon.initialise_estimate() </kbd></li>
 *
 * <li>The reconstruction can be configured before, or at anytime
 * during the reconstruction using setter methods. For example the
 * algorithm choice and relaxation parameters:
 * <br><kbd>my_planar_recon.set_algorithm(HIO)</kbd>
 * <br><kbd>my_planar_recon.set_relaxation_parameter(0.95)</kbd>
 * <br>The algorithm is fully customisable through the
 * PlanarCDI::set_custom_algorithm method. See page 22 of the
 * H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive imaging
 * using short wavelength light sources, Journal of Modern Optics,
 * 2010, DOI: 10.1080/09500340.2010.495459 for the formalism.</li>
 *
 * 
 * <li>The reconstruction is then performed with multiple calls to the
 * PlanarCDI::iterate() method:
 * <br><kbd>my_planar_recon.iterate()</kbd> <br>This will
 * automatically update the values in the Complex_2D estimate and the
 * current error. The error can be retrieved using
 * PlanarCDI::get_current_error().
 *
 * </ol>
 */

#ifndef PHASERETRIEVALBASE_H
#define PHASERETRIEVALBASE_H

#include <map>
#include <string>
#include "FourierT.h"
#include "Double_2D.h"

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
class Double_2D;

class PlanarCDI{

 protected:

  /**a reference to the Complex_2D object which is altered during each
     iteration */
  Complex_2D & complex;

  /** The samplings in x */
  int nx;

  /** The samplings in y */
  int ny;

  /**a Fourier transform objected used to perform 
     the forward and backward transformations */
  FourierT fft; 

   /** temporary Complex_2Ds which are used in the computation of the
      PFS and PF terms for each iteration. */
  Complex_2D temp_complex_PFS;
  Complex_2D temp_complex_PF;
  Complex_2D temp_complex;
  Complex_2D temp_complex_PS;
  Complex_2D temp_complex_PSF;

  /**the relaxation parameter */
  double beta;

  /** a copy of the support */
  Double_2D support;

  /**a copy of the square root of the intensity at the detector plane. */
  Double_2D intensity_sqrt;

  /** the algorithm which is being used: ER, HIO etc. */
  int algorithm;

  /**an array which holds the structure of the algorithm. i.e. the 
     coefficients for each term in the iteration equation.*/
  double algorithm_structure[NTERMS];

  /** the difference between the intensity in the detector plane 
      the current estimated intensity. */
  double current_error;

  /** how many best estimates to store */
  int n_best;

  /** array of the best estimates */
  Complex_2D ** best_array;

  /** array of the error corresponding to each of the
      best estimates */
  double * best_error_array;

  /** a mapping between the algorithm name (string) and identification
      number */
  static std::map<std::string,int> * algNameMap;

 public:
  
  /**
   * Constructor. The default algorithms is set to HIO with a relaxation
   * parameter of 0.9.
   *
   * @param complex The complex 2-D field which is modified by this
   *   object. This represents the exit-surface-wave of the sample.
   *   It may be pre-initialised to a best first guess (e.g manually);
   *   it may be loaded from a file (e.g. the output from a previous
   *   reconstruction job); or it may be unitialised, and the default
   *   initialisation provided by PlanarCDI may be used.
   *
   * @param n_best The number of "best estimates" to keep. See
   *   PlanarCDI::get_best_result below for more detail. By default
   *   this option is switched off.
   */
  PlanarCDI(Complex_2D & complex, int n_best=0);

  /**
   * Destructor. This is not interesting for users.
   */
  virtual ~PlanarCDI();

 /**
  * The main method for running the reconstruction!  This method
  * transforms the Complex_2D given to it according to:
  * <br>\f$ x_{k+1} = g ( x_{k} ) \f$
  * <br> Where \f$g()\f$ depends on the algorithm used.
  * The current error is also updated when this method is called.
  */
  virtual int iterate();
  
  /**
   * Set the values in the Complex_2D to an initial estimate.  For
   * PlanarCDI, this is set to random numbers inside the support for
   * both the real and complex components, and zero outside the
   * support.
   *
   * @param seed The seed for the random number generator. By default
   * this is 0.
   */
  virtual void initialise_estimate(int seed=0);
  

  /**
   * While the reconstruction is running this object will store the
   * "n_best" number of best estimates according to the error metric.
   * So, for example, if HIO is used and the error enters and leaves
   * an error minima during reconstruction, the best (ie. lowest error)
   * Complex_2D results will be kept. They can then be retrieved using
   * this method.
   *
   * @param index get best result number "index". 0 is the best,
   *    1 is the 2nd best etc... 
   * @param error "error" is filled with the metric error of the 
   * selected best result
   */
  Complex_2D * get_best_result(int index, double & error);

   /**
    * Set the sample support. A "0" value in the Double_2D is
    * interpreted as being outside the support. All positive values
    * are considered inside the support. The support can be updated at
    * anytime during reconstruction by calling this method.
    *
    * @param object_support The sample support
    */ 
  void set_support(const Double_2D & object_support);
  
  /**
   * Set the detector diffraction image (i.e. the square of the
   * amplitude of the wavefield at the detector).
   *
   * @param detector_intensity The intensity at the detector 
   */ 
  void set_intensity(const Double_2D & detector_intensity);

  /**
   * Set the relaxation parameters of the reconstruction algorithm.
   * This is \f$ \beta \f$ used exactly as it is in Harry's review
   * paper. This maybe reset at anytime during the reconstruction by
   * calling this method. The default value is 0.9.
   *
   * @param relaxation_parameter The algorithm relaxation_parameter
   */ 
  void set_relaxation_parameter(double relaxation_parameter){
    beta = relaxation_parameter;
    set_algorithm(algorithm); //force algorithm constants to update
  };
  
  /**
   * Get the autocorrelation function of the intensity data.
   *
   * @param autoc The resulting autocorrelation.
   */ 
  void get_intensity_autocorrelation(Double_2D & autoc);

  /**
   * Set the algorithm. By default HIO is used.
   *
   * @param alg The algorithm. The options are:
   * <br> ER - error reduction
   * <br> BIO - basic input-output
   * <br> BOO - basic output-output
   * <br> HIO - hybrid input-output
   * <br> DM - difference map
   * <br> SF - solvent-flipping
   * <br> ASR - averaged successive reflections
   * <br> HPR - hybrid projection reflection
   * <br> RAAR - relaxed averaged alternating reflectors
   */ 
  void set_algorithm(int alg);

   /**
   * Iterative reconstruction algorithms can be expressed as a
   * combination of the following 5 operators:<br>\f$ \hat{P}_S
   * \hat{P}_F\f$, \f$\hat{P}_F \hat{P}_S\f$ , \f$ \hat{P}_S \f$ , \f$
   * \hat{P}_F \f$ and \f$ \hat{I} \f$ <br> Where \f$ \hat{P}_S \f$ is
   * the support projection and \f$ \hat{P}_F \f$ is the intensity
   * projection. These can be combined to form a basis of 10 vectors
   * (see Harray's review). The parameters to this method set the
   * coefficient for each of these vectors. Note: the vector set is
   * linearly dependant.
   *
   * Please see page 22 of the H.M. Quiney review: TUTORIAL REVIEW,
   * Coherent diffractive imaging using short wavelength light
   * sources, Journal of Modern Optics, 2010, DOI:
   * 10.1080/09500340.2010.495459
   *
   * @param m1 Coefficient to the \f$ [ \hat{P}_S \hat{P}_F -
   * \hat{P}_F \hat{P}_S ]x_k \f$ term
   * @param m2 Coefficient to the \f$ [ \hat{P}_S \hat{P}_F - 
   * \hat{P}_F ]x_k \f$ term
   * @param m3 Coefficient to the \f$ [ \hat{P}_S \hat{P}_F - 
   * \hat{P}_S ]x_k \f$ term
   * @param m4 Coefficient to the \f$ [ \hat{P}_S \hat{P}_F - 
   * \hat{I} ]x_k \f$ term
   * @param m5 Coefficient to the \f$ [ \hat{P}_F \hat{P}_S - 
   * \hat{P}_F ]x_k \f$ term
   * @param m6 Coefficient to the \f$ [ \hat{P}_F \hat{P}_S - 
   * \hat{P}_S ]x_k \f$ term
   * @param m7 Coefficient to the \f$ [ \hat{P}_F \hat{P}_S - 
   * \hat{I} ]x_k \f$ term
   * @param m8 Coefficient to the \f$ [ \hat{P}_F - \hat{P}_S ]x_k \f$ term
   * @param m9 Coefficient to the \f$ [ \hat{P}_F - \hat{I} ]x_k \f$ term
   * @param m10 Coefficient to the \f$ [ \hat{P}_S - \hat{I} ]x_k \f$ term
   */  
  void set_custom_algorithm(double m1, double m2, double m3, 
			    double m4, double m5, double m6, 
			    double m7, double m8,
			    double m9, double m10);
  
  
  /**
   * Print the current algorithm to standard output. The algorithm is
   * in the form: <br> \f$ x_{k+1} = x_k + (a \hat{P}_S \hat{P}_F + b
   * \hat{P}_F \hat{P}_S + c \hat{P}_S + d \hat{P}_F + e \hat{I} )x_k
   * \f$
   */
  void print_algorithm();

  /**
   * Get the difference between the estimated diffraction and the
   * actual diffraction pattern.  <br> \f$ \chi^2 = \frac {
   * \sum_{pixels} ( M_{pixel} - \sqrt{I_{pixel}} )^2 } { \sum_{pixels}
   * {I_{pixel}} } \f$.  <br> Where \f$ M \f$ is the magnitude of the
   * current estimate and \f$ I \f$ is the measured intensity.  
   *
   * NOTE: This method actually gives you the error for the previous
   * iteration. This is because it's much faster to calculate the
   * error during the intensity scaling step of the iteration rather
   * than performing a FFT just for the error calculation.
   *
   * @return The error metric
   */
  double get_error();

  /**
   * Given the name of an algorithm, e.g. "HIO", return the unique
   * integer identifier. This is used by the configuration file
   * reader, and probably won't be needed by users.
   *
   * @param algorithm_name The algorithm name 
   * @return The integer identifier
   */ 
  static int getAlgFromName(std::string algorithm_name){
    std::map<std::string,int>::iterator alg;
    alg = algNameMap->find(algorithm_name);
    if(alg == algNameMap->end() )
      return -1;
    return (alg->second);
  }
  
  /**
   * Apply the shrinkwrap algorithm. The magnitude of the current ESW
   * estimate will be convolved with a Gaussian, thresholded and set
   * as the new support.
   *
   * @param gauss_width The width ( \f$ \sigma \f$ ) of the Gaussian
   * in pixels. The default is 1.5.
   * @param threshold The threshold which is applied, 
   * as a fraction of the maximum pixel value. The default is 0.1 (10%).
   */  
  virtual void apply_shrinkwrap(double gauss_width=1.5, double threshold=0.1);  
  
  /**
   * Get the current support. This might be useful if shripwrap has
   * been applied, as it allows you to get the updated support.
   *
   * @param support The content of "support" is overridden with the
   * current support.
   */  
  void get_support(Double_2D & support);


  /**
   * Apply the support constraint.
   * 
   * @param c The complex field to apply the support constraint on
   */ 
  virtual void apply_support(Complex_2D & c);

  /**
   * Apply the intensity constraint. The ESW is projected to the
   * detector plane, the intensity is scaled to match the measured
   * data and the field is projected back into the sample plane.
   * 
   * @param c The complex field to apply the intensity constraint on
   */ 
  virtual void project_intensity(Complex_2D & c);


   /**
   * The intensity is scaled to match the measured data. The error is
   * updated in this method.
   * 
   * @param c The complex field to apply the scaling to
   */  
  virtual void scale_intensity(Complex_2D & c);


  /**
   * Propagate to the sample plane using a fast fourier transform.
   *
   * @param c The Complex_2D field to propagate
   */
  virtual void propagate_to_sample(Complex_2D & c);

  /**
   * Propagate to the detector plane using a fast fourier transform.
   * @param c The Complex_2D field to propagate
   */
  virtual void propagate_to_detector(Complex_2D & c);


 protected:

  /**
   * Sets up the mapping between algorithm name and number. This is
   * only called during the initialisation and should not be called
   * elsewhere.
   * 
   * @return The algorithm name-iteger mapping
   */  
  static std::map<std::string,int> * set_up_algorithm_name_map();


  /**
   * Apply a threshold to the values in "array".
   *
   * @param array The array to be thresholded
   * @param threshold The threshold which is applied, 
   * as a fraction of the maximum pixel value. The default is 0.1 (10%).
   */    
  void apply_threshold(Double_2D & array, double threshold);
  

  /**
   * Convolved "array" with a 2-D Gaussian function.  To speed up
   * computation we only convolve up to 4 pixels away from the
   * Gaussian peak.
   *
   * @param array The array to be convolved
   * @param gauss_width The width ( \f$ \sigma \f$ ) of the Gaussian
   * in pixels. The default is 1.5.
   * @param pixel_cut_off Convolve up to 4 pixels away from the
   * Gaussian peak.
   */  
  void convolve(Double_2D & array, double gauss_width, int pixel_cut_off=4);
  
    
};

#endif