/**
 * @file PhaseDiverseCDI.h
 * @class PhaseDiverseCDI
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief  The class which performs phase diverse and ptychographic 
 * reconstruction.
 *
 * This class can be used to perform phase diverse or ptychographic
 * reconstruction of either fresnel or plane-wave CDI data. Any number
 * of frames (also called local/single frames or probes in this
 * documentation) may be added to the reconstuction. In order to
 * perform a reconstuction, you will need to create either a new
 * FresnelCDI or PlanarCDI object for each of these 'local' datasets.
 * The FresnelCDI or PlanarCDI objects must then be passed to a
 * PhaseDiverseCDI object. Because the each sub-iteration (ie. each
 * iteration of a local frame) is performed using the code in
 * FresnelCDI/PlanarCDI, all the functionality available in these
 * classes is also available here. For example complex constraint can
 * be set, shrink-wrap can be used, the estimate can be initialised
 * using the results from a pervious reconstruction etc. An example
 * can be found in the /examples directory, demonstrating how to use
 * PhaseDiverseCDI.
 *
 * The displacement in position between 'local' frames may be either
 * transverse to the beam direction, or longitudinal, along the beam
 * direction. The longitudinal position (for FresnelCDI) is set when
 * the FresnelCDI objects are initially constructed. Transverse
 * positions are set when the FresnelCDI/PlanarCDI objects are add to
 * the PhaseDiverseCDI object. This class allows the transverse
 * positions to be automatically adjusted during reconstuction using
 * the "adjust_positions" function.
 *
 * The code allows for several options in the type of reconstruction
 * which is done:

 * - The reconstruction can be performed in series or parallel. In the
 *   case of series (see for example the paper....), a 'local' frame
 *   will undergo one or more iterations, the result will be updated
 *   to a 'global' estimate of the sample, and this estimate will form
 *   the starting point for the next frame. This process repeats. For
 *   a parallel reconstuction (see .....), each frame with
 *   independantly undergo one of more iteration, the result from all
 *   frames will be merged to form a new estimate of the sample, this
 *   estimate then becomes the starting point for the next iteration
 *   of all frames.
 *
 * - The number of local iterations to perform before updating the
 *   result to the 'global' function can be set.
 *
 * - The feedback parameter, beta, may be set. This quantity is used
 *   to set how much of the previous 'global' sample function will be
 *   left after the next 'global' iteration.
 *
 * - The amplification factor, gamma, and the probe scaling, alpha,
 *   may also be see. These parameters control the weighting of one
 *   frame (and pixels within a frame) with respect to each
 *   other. See the paper... for more detail.
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

  /* A vector holding a pointer to each of the FresnelCDI/PlanarCDI objects */
  std::vector<BaseCDI*> singleCDI;

  /* A vector holding a pointer to the transmission function/esw
     result for each FresnelCDI/PlanarCDI objects */
  std::vector<Complex_2D*> single_result;

  /** A vector holding the weighting function for each frame */  
  std::vector<Double_2D *> weights; 
  
  /* A vector of the transverse positions in x */
  std::vector<double> x_position;

  /* A vector of the transverse positions in y */
  std::vector<double> y_position;
  
  //parameters controlling the feedback
  double beta;
  double gamma;
  std::vector<double> alpha;

  /** The current estimate of the 'global' sample function */
  Complex_2D * object;

  /** The number of iterations to perform for each 'local' frame */
  int iterations_per_cycle;

  /** A factor which controls sub-pixel alignment. 
      1=regular, 2 = 2 'global' pixels for every 1 'local' pixel.
      i.e. 2 allows alignment to within half a pixel. */
  int scale;

  /** The size in x and y of the 'global' sample function. */
  int nx,ny;
  
  /** Minimum coordinates of the 'global' function based on the
   positions entered by the user. */
  int x_min;
  int y_min;

  /** a flag for running in either series or parallel mode */
  bool parallel; 

  /** a flag indicating whether the weights need to be recaculated */
  bool weights_set;

 public:

  /** 
   * Construct a PhaseDiverseCDI object for phase diverse or
   * ptychographic reconstruction. The data should be entered later
   * using the 'add_new_position' function.
   *
   * @param beta The feedback parameter. By default this is 1 (no feedback).
   * @param gamma The amplification factor. By default this is 1 (no amplification).
   * @param parallel true - run in parallel mode, false - run in series
   *        mode. By default series mode is set.
   * @param granularity  A factor which controls sub-pixel alignment. 
   *   1=regular, 2 = 2 'global' pixels for every 1 'local' pixel.
   *   i.e. 2 allows alignment to within half a pixel. 
   *   NOTE: This is not currently working properly!
   */
  PhaseDiverseCDI(double beta=1.0, 
		  double gamma=1.0,
		  bool parallel=false,
		  int granularity=1);
  
  /** 
   * Destructor for PhaseDiverseCDI
   */
  ~PhaseDiverseCDI();


  /**
   *
   *
   */
  void add_new_position(BaseCDI * local, 
			double x=0, double y=0, 
			double alpha=1);

  /**
   * Initialise the estimate of the 'global' object function. The
   * initialisation is performed using the current estimate for each
   * 'local' frame. For this reason you must first initialise each
   * FresnelCDI/PlanarCDI object individually before calling this
   * function.
   */
  void initialise_estimate();

  /** 
   * Perform one iteration. This with involve performing one or more
   * iterations for each sub-set of data. The exect implementation
   * depends on whether the reconstruction is being run in series or
   * parallel modes.
   */
  void iterate();
  
  /////////////////////////////////
  // Get and setter methods
  /////////////////////////////////

  /**
   * Set the number of 'local' frame iterations before updating the
   * result to the 'global' object function.
   *
   * @param iterations The number of 'local' iterations.
   */
  void set_iterations_per_cycle(int iterations){
    iterations_per_cycle = iterations;
  }

  /**
   * Set the feed-back parameter.
   *
   * @param beta The feedback parameter
   */
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
  
  /**
   * This function allows you to access the 'global' sample function.
   *
   * @return The current estimate of either the transmision (for
   * FresnelCDI) or exit-surface-wave (for PlanarCDI)
   */
  Complex_2D * get_transmission();
  
  /**
   * Set the 'global' sample function. This method could be used, for
   * example, instead of 'initialise_estimate' to initialise the
   * results to those from a previous reconstruction.
   *
   * @param new_transmission The transmission (for FresnelCDI) or
   * exit-surface-wave (for PlanarCDI) to copy.
   */
  void set_transmission(Complex_2D & new_transmission);
  
  ///////////////////////////////////////////
  // position adjustment
  //////////////////////////////////////////

  /**
   * Align the frames. The transverse positions will be automatically
   * adjusted. The new positions can be retrieved using the
   * "get_final_x(y)_position" functions.
   *
   * more to come....
   */
  void adjust_positions(double step_size=4, bool forwards=true);
  

  /**
   * Get the current position of the 'n_probe'th frame.
   *
   * @param n_probe The index of the frame you wish to retrieve the
   *                position for. 0=the first frame you added, 
   *                1=2nd frame added etc....
   * @return the transverse horizontal position
   */
  double get_final_x_position(int n_probe){
    return x_position.at(n_probe);}
  
  /**
   * Get the current position of the 'n_probe'th frame.
   *
   * @param n_probe The index of the frame you wish to retrieve the
   *                position for. 0=the first frame you added, 
   *                1=2nd frame added etc....
   * @return the transverse vertical position
   */
  double get_final_y_position(int n_probe){
    return y_position.at(n_probe);}


 private:
  
  /**
   * Update a 'local' frame result to the 'global' object.
   * The result is weighted before being added.
   *
   * @param n_probe The local frame number.
   */
  void add_to_object(int n_probe);

  /**
   * Update a 'local' frame result from the 'global' object.
   *
   * @param n_probe The local frame number.
   */
  void update_from_object(int n_probe);

  /**
   * This function is basically a wrapper to check the BaseCDI type so
   * that the correct type of sample function is retreive (either
   * transmission function or the exit-surface-wave.
   *
   * @param local A pointer to either a PlanarCDI or FresnelCDI object
   * @result result The estimate is copied from local to result
   */
  void get_result(BaseCDI * local, Complex_2D & result);

   /**
   * This function is basically a wrapper to check the BaseCDI type so
   * that the correct type of sample function is set (either
   * transmission function or the exit-surface-wave.
   *
   * @param local A pointer to either a PlanarCDI or FresnelCDI object
   * @result result The estimate is copied from result to local
   */ 
  void set_result(BaseCDI * local, Complex_2D & result);

  /**
   * This function is used to reallocate memory for the 'global'
   * sample object. It is only used when add_new_position is called,
   * in the case that the frame does not fit within the bounds of the
   * current object array.
   */
  void reallocate_object_memory(int new_nx, int new_ny);

  /**
   * This function is used during reconstruction in parallel mode. It
   * is similar to simply scaling all the elements of the 'global'
   * sample array, (e.g. with object.scale(factor). However, this
   * method preserves the value of the elements outside the boundary
   * of the sample. ie. it checks the weights of all the local frames
   * to work out which area to scale.
   * 
   * @param factor The scaled elements go from 'a' to 'factor*a'.
   *
   */
  void scale_object(double factor);
  
  /**  void get_object_sub_grid(Complex_2D & result,
			   double x_offset,
			   double y_offset); **/
 

  /**
   * Calculate the weights for each frame in the reconstruction.  See
   * a description above for how they are set.
   */
  void set_up_weights();

  /**
   * Needs some work......
   */
  int check_position(int n_probe, double shift=4, int tries=0);

  
  /**
   * get the global pixel "x" corrdinate using a local frame "x" and
   * the frame offset.
   *
   * @param x the local frame pixel position in x 
   * @param x_offset the offset of the local frame 
   *                 with respect to the other frames.
   */
  inline int get_global_x_pos(int x, double x_offset){
    return (x-x_offset-x_min); 
  }
 
  /**
   * See get_global_x_pos
   */
  inline int get_global_y_pos(int y, double y_offset){
    return (y-y_offset-y_min);
  }

  /**
   * See get_global_x_pos
   */
  inline int get_local_x_pos(int x, double x_offset){
    return x + x_offset + x_min;
  }
  
  /**
   * See get_global_x_pos
   */
  inline int get_local_y_pos(int y, double y_offset){
    return y + y_offset + y_min; 
  }


};


#endif
