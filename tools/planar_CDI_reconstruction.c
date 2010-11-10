/**
 * @file 
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au>
 *
 * @section DESCRIPTION
 *
 * @todo 
 * Allow input for the starting guess (phase and mag?)
 * 
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Projection.h"
#include "Config.h"
#include "FourierT.h"

using namespace std;

#define OUTPUT_MINIMAL 0
#define OUTPUT_ITER    1
#define OUTPUT_ERROR   2

int main(int argc, char * argv[]){

  /** work out which config file to use **/
  string config_file = "planar.config";

  //and set the seed of the inital guess
  int seed = 0;

  if(argc==1)
    cout << "No config file given, using "
	 << "the default: "<<config_file<<endl;

  if(argc==2){
    cout << "No seed value given... using the default (0)" <<endl;
    config_file=argv[1];
  }

  if(argc==3)
    seed = atoi(argv[2]);

  if(argc>3){
    cout << "Too many arguments given. Usage: "
	 << "planar_CDI_reconstuction <config filename> <seed>"<<endl;
      exit(0);
  }

  Config c(config_file);
  if(c.getStatus()==FAILURE){
    cout << "Could not open the file "<< config_file<< ". Exiting" <<endl;
    exit(0);
  }

  /** read the config file **/

  //the data file name
  string data_file_name = c.getString("data_file_name");

  string data_file_type = c.getString("data_file_type");

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = c.getString("support_file_name");

  string support_file_type = c.getString("support_file_type");

  //output filename prefix
  string output_file_name_prefix = c.getString("output_file_name_prefix");

  //output the current image every "output_iterations"
  int output_iterations = c.getDouble("output_iterations");

  //names of the algorithms to use in order
  list<string> * algorithms = c.getStringList("algorithm");
 
  //the number of iterations to perform for each algorithm
  list<int> * iterations = c.getIntList("iterations");
  
  //do some error checking. Were all the values we need present
  //in the config file?
  if(c.getStatus()==FAILURE){
    cout << "Problem reading the configuration file. Exiting" <<endl;
    exit(0);
  }
  if(algorithms->size()!=iterations->size()){
    cout << "The number of algorithms and the number "
	 <<"of iteractions do not match. Exiting"<< endl;
    exit(0);
  }

  /** get some optional configuration parameters **/
  int output_level = c.getInt("info_level");
  int output_diffraction_estimate = c.getInt("output_diffraction_estimate");
  int use_log_scale_for_diffraction = c.getInt("use_log_scale_for_diffraction");
  int use_log_scale_for_object = c.getInt("use_log_scale_for_object");

  /*** get the diffraction data from file and read into an array *****/
  int nx, ny;
  double ** data;
  int status;
  if(data_file_type=="tiff")
    status = read_tiff(data_file_name, &nx, &ny, &data);  
  else if(data_file_type=="ppm")
    status = read_ppm(data_file_name, &nx, &ny, &data);
  else{ //check that the file type is valid
    cout << "Can not process files of type \""<< data_file_type 
	 <<"\".. exiting"  << endl;
    return(1);
  }
  //and check that the file could be opened okay
  if(!status){
    cout << "failed to get data from "<< data_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }

  /******* get the support from file and read it into an array *****/
  double ** support;
  int nx_s=0, ny_s=0;
  if(support_file_type=="tiff")
    status = read_tiff(support_file_name, &nx_s, &ny_s, &support);  
  else if(support_file_type=="ppm")
    status = read_ppm(support_file_name, &nx_s, &ny_s, &support);
  else{ //check that the file type is valid
    cout << "Can not process files of type \""<< support_file_type 
	 <<"\".. exiting"  << endl;
    return(1);
  }
  if(!status){
    cout << "failed to get data from "<< support_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }
  if( nx_s != nx || ny_s != ny){
    cout << "Dimensions of the support to not match ... exiting"  << endl;
    return(1);
  }

  /*******  set up the reconstuction *********************/

  //create the projection object which will be used to
  //perform the reconstuction.
  Complex_2D object_estimate(nx,ny);
  Projection proj(&object_estimate);
 
  //set the support and intensity
  proj.set_support(support);
  proj.set_intensity(data);

  //Initialise the current object ESW with a random numbers
  proj.initialise_estimate(seed);

  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.
  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];

  //get up a temporary FFT in case we need to output
  //the guess of the diffraction pattern
  FourierT * fft;
  if(output_diffraction_estimate){
    fft = new FourierT(nx,ny);
  }

  /*** run the reconstruction ************/


  list<string>::iterator algorithms_itr = algorithms->begin();
  list<int>::iterator iterations_itr = iterations->begin();
  //loop over the algorithms
  int i=0;
  int cumulative_iterations = 0;
  while(algorithms_itr != algorithms->end()&&
	iterations_itr != iterations->end()){

    if(output_level!=OUTPUT_MINIMAL)
      cout << "Switching to the "<< (*algorithms_itr)
    	   <<" algorithm" << endl;
    
    //get the projection
    int alg = Projection::getAlgFromName(*algorithms_itr);
    if(alg == -1 ){
      std::cout << "Could not find reconstuction algorithm"
		<< " with the name "<< (*algorithms_itr)
		<< ". Exiting" << endl;
      exit(0);
    }

    proj.set_algorithm(alg);
    cumulative_iterations+=(*iterations_itr);

    for(; i < cumulative_iterations; i++){
      if(output_level!=OUTPUT_MINIMAL)
	cout << "Iteration " << i << endl;
     

      //apply the iterations  
      proj.iterate(); 
      if(output_level==OUTPUT_ERROR && i>0)
	cout << "Error for iteration "<<(i-1)<<" is " << proj.get_error() << endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,&result);
	temp_str << output_file_name_prefix << "_" << i << ".ppm" << flush;
	write_ppm(temp_str.str(), nx, ny, result, use_log_scale_for_object);
	//temp_str.clear();

	//output the estimation of the intensity in 
	//the detector plane if needed
	if(output_diffraction_estimate){
	  Complex_2D * temp = object_estimate.clone();
	  fft->perform_forward_fft(temp);
	  temp->get_2d(MAG_SQ,&result);
	  temp_str << output_file_name_prefix << "_diffraction_" << i << ".ppm" << flush;
	  write_ppm(temp_str.str(), nx, ny, result, use_log_scale_for_diffraction); 
	  delete temp;
	}
      }
    }
    iterations_itr++;
    algorithms_itr++;
  }

  //write out the final image
  object_estimate.get_2d(MAG,&result);
  ostringstream temp_str ( ostringstream::out ) ;
  temp_str << output_file_name_prefix << "_final.ppm" << flush;
  write_ppm(temp_str.str(), nx, ny, result, use_log_scale_for_object);
  
  if(output_diffraction_estimate){
    Complex_2D * temp = object_estimate.clone();
    fft->perform_forward_fft(temp);
    temp->get_2d(MAG_SQ,&result);
    temp_str << output_file_name_prefix << "_diffraction_final.ppm" << flush;
    write_ppm(temp_str.str(), nx, ny, result, use_log_scale_for_diffraction); 
    delete temp;
  }

  //clean up
  for(int i=0; i< nx; i++){
    delete[] result[i];
    delete[] support[i];
    delete[] data[i];
  }

  delete[] result;
  delete[] support;
  delete[] data;

  delete algorithms;
  delete iterations;
  
  if(output_diffraction_estimate)
    delete fft;
 
  return 0;
}

