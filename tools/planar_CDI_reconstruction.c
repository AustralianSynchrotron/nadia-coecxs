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

using namespace std;

int main(int argc, char * argv[]){

  /** work out which config file to use **/
  string config_file = "planar.config";
  
  if(argc==1)
    cout << "No config file given, using "
	 << "the default: "<<config_file<<endl;

  if(argc==2)
    config_file=argv[1];

  if(argc>2){
    cout << "Too many arguments given. Usage: "
	 << "planar_CDI_reconstuction <config filename>"<<endl;
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

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = c.getString("support_file_name");

  //output filename prefix
  string output_file_name_prefix = c.getString("output_file_name_prefix");

  //output the current image ever "output_iterations"
  int output_iterations = c.getDouble("output_iterations");

  //names of the algorithms to use in order
  list<string> * algorithms = c.getStringList("algorithm");
 
  //the number of iterations to perform for each algorithm
  list<int> * iterations = c.getIntList("iterations");
  
  //do some error checking. where all the values we need present
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
  

  /*** get the diffraction data from file and read into an array *****/
  int nx, ny;
  double ** data;
  int status = read_ppm(data_file_name, &nx, &ny, &data);  
  
  //check that the file could be opened okay
  if(!status){
    cout << "failed to get data from "<< data_file_name 
	 <<".. exiting"  << endl;
    return(1);
  }

  /******* get the support from file and read it into an array *****/
  double ** support;
  int nx_s, ny_s;
  status = read_tiff(support_file_name, &nx_s, &ny_s, &support);  
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
  proj.initialise_estimate(0);

  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.
  double ** result = new double*[nx];
  for(int i=0; i < nx; i++)
    result[i]= new double[ny];

 
  /*** run the reconstruction ************/

  list<string>::iterator algorithms_itr = algorithms->begin();
  list<int>::iterator iterations_itr = iterations->begin();
  //loop over the algorithms
  int i=0;
  int cumulative_iterations = 0;
  while(algorithms_itr != algorithms->end()&&
	iterations_itr != iterations->end()){

    //proj.set_algorithm((*algorithms_itr)); //<-- to be fixed    
    cout << "Switching to the "<<(*algorithms_itr)
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

      cout << "iteration " << i << endl;
      
      //apply the iterations  
      proj.iterate(); 
    
      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,&result);
	temp_str << output_file_name_prefix << "_" << i << ".ppm";
	write_ppm(temp_str.str(), nx, ny, result);
	temp_str.clear();

	//uncomment to output the estimated 
	/**Complex_2D * temp = object_estimate.clone();
	   fft.perform_forward_fft(temp);
	   temp->get_2d(MAG_SQ,&result);
	   temp_str << "diffraction.ppm";
	   write_ppm(temp_str.str(), nx, ny, result, true); 
	   delete temp;
	**/
      }
    }
    iterations_itr++;
    algorithms_itr++;
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

  return 0;
}

