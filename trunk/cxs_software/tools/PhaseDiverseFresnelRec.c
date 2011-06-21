/**
 * @file PhaseDiverseCDI_reconstruction
 *
 * A quick tool to run phase-diverse/ptychographic reconstruction
 * from the command line. The file format is similar to the one
 * used by Corey Putkunz's code 
 *
 * @author Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * Last modified on 21/6/2011
 *
 */
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "io.h"
#include "TransmissionConstraint.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "PhaseDiverseCDI.h"
#include "Config.h"
#include <string>
//#include "utils.h"
#include <vector>

using namespace std;

/**********************************/
void print_usage(){

  cout << "Usage: " << endl << endl
       << "PhaseDiverseFresnelRec.exe <data file filename> " 
       << "<iterations> <sub-iterations> " << endl
       << "optional: <beta> <gamma> <running-mode> <do alignment> <seed> " 
       << endl << endl;
}

int main(int argc, char * argv[]){

  //parameters (from corey's code):
  //data_list.txt ITERATIONS IT_SUB beta gamma

  if(argc<3){
    print_usage();
    exit(0);
  }

  string datalist_filename = argv[1];
  fstream filelist(datalist_filename.c_str(), fstream::in);
  if(!filelist.is_open()){
    cerr << "Error opening the data filename list.. exiting" 
	 << endl;
    exit(0);
  }  

  unsigned int iterations = atoi(argv[2]);
  if(iterations==0){
    cerr << "Invalid number of iterations: "<< iterations 
	 << " ... exiting" << endl;
    exit(0);
  }

  int subiterations = 1;
  double beta = 1.0;
  double gamma = 1.0;
  bool mode = false; //serial mode
  int seed = 0;
  bool alignment = false;

  if(argc>3)
    subiterations = atoi(argv[3]);

  if(argc>4)
    beta = atof(argv[4]); 

  if(argc>5)
    gamma = atof(argv[5]); 

  if(argc>6)
    mode = atoi(argv[6]); 

  if(argc>6)
    alignment = atoi(argv[7]); 
  
  if(argc>8)
    seed = atoi(argv[8]); 
  
  //initialise some arrays to hold data.
  vector<FresnelCDI *> proj;
  vector<Complex_2D *> object_estimate;

  //Make a complex contraint with the most
  //basic constraint
  TransmissionConstraint tc;
  tc.set_enforce_unity(true);
  tc.set_charge_flipping(false);

  //Make a new PhaseDiverseCDI object
  PhaseDiverseCDI pd(beta,gamma,mode);
  pd.set_iterations_per_cycle(subiterations);
  
  while(!filelist.eof()){
    
    int x_pos = 0;
    int y_pos = 0;

    string image_file_name = "";
    string white_field_file_name = "";
    string param_file_name = "";

    filelist >> image_file_name;
    filelist >> white_field_file_name;
    filelist >> param_file_name ;
    filelist >> x_pos ;
    filelist >> y_pos ;

    cout << image_file_name << " "<< white_field_file_name
	 << " "<< param_file_name << " " << x_pos << " " << y_pos <<endl;

    //read the parameter file 
    Config config_file(param_file_name);
    
    double wavelength = config_file.getDouble("LAMBDA");
    double z2 = config_file.getDouble("Z2");
    double z3 = config_file.getDouble("Z3");
    double zI = config_file.getDouble("ZI");

    double fs = zI - z2; //focal to sample distance
    double fd = z3 - z2; //focal to detector distance
    double ps = config_file.getDouble("IMAGE_WIDTH");
    double norm = config_file.getDouble("NORMALISATION");
    double nx = config_file.getInt("N");
    double ny = nx;

    if(config_file.getStatus()==FAILURE){
      cerr << "Could not read the parameter file "
	   << param_file_name << " correctly.."
	   << "exiting" << endl;
      exit(0);
    }
    
    Double_2D diffraction_image(nx,ny);
    Complex_2D white_field(nx,ny);
    read_image(image_file_name, diffraction_image, nx, ny);
    read_cplx(white_field_file_name, white_field);
    
    //make the support from a thresholded white-field
    Double_2D beam(nx,ny);
    white_field.get_2d(MAG,beam);
    double max = beam.get_max();
    double threshold = 0.5;

    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
	if( beam.get(i,j) > max*threshold )
	  beam.set(i,j,100);
	else
	  beam.set(i,j,0);
      }
    }
    
    //Set up the fresnel CDI in the same way you would
    //if you weren't doing phase-diversity
    object_estimate.push_back(new Complex_2D(nx,ny));
    
    //set-up the reconstruction for a single frame
    proj.push_back(new FresnelCDI(*object_estimate.back(),
				  white_field,
				  wavelength,
				  fd,
				  fs,
				  ps,
				  norm));
    
    //set the support and intensity and initialise
    proj.back()->set_intensity(diffraction_image);
    proj.back()->set_support(beam,true); //use fussy edges  
    proj.back()->initialise_estimate(seed);
    
    //add the most basic additional constraint
    proj.back()->set_complex_constraint(tc);

    //New part.. Add the FresnelCDI to the PhaseDiverseCDI.
    pd.add_new_position(proj.back(), x_pos, y_pos);
    
  }
  filelist.close();
  
  //initialise the phase diverse transmission function
  pd.initialise_estimate();
  
  //now run the reconstruction
  for(int i=0; i < iterations; i++){

    //do some position alignment
    if(alignment && i==0)
      pd.adjust_positions(PhaseDiverseCDI::CROSS_CORRELATION);
    if(alignment && i==10)
      pd.adjust_positions(PhaseDiverseCDI::MINIMUM_ERROR);
    
    //iteration!
    pd.iterate();
    
  }


  Complex_2D * object = pd.get_transmission();
  
  Double_2D result(object->get_size_x(),object->get_size_y());
  
  //write the magnitude and phase of it to an image file
  object->get_2d(MAG,result);
  write_image("trans_mag_recovered.tiff",result);

  object->get_2d(PHASE,result);
  write_image("trans_phase_recovered.tiff",result);
  
  //save the transmission function in case we want to use it later.
  write_cplx("trans.cplx",*object);


  return 0;
  
}
