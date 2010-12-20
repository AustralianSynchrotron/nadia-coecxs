#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>
#include "io.h"
#include "Double_2D.h"

using namespace std;

#define FAILURE 0
#define SUCCESS 1

/***************************************************************/

/***************************************************************/
int read_dbin(string file_name, int nx, int ny, Double_2D & data){
 
  //open the input file:
  FILE * file = fopen(file_name.c_str(), "r");

  //error check.
  if(!file){
    cout << "Could not open the file " << file_name << endl;
    return FAILURE;
  }

  double * buffer = new double [nx*ny];
  fread(buffer, sizeof(double), nx*ny, file);
  fclose(file);
  
  //Do a sanity check. Is the file size right 
  //for a nx by ny array with 16 bit pixel values?
  
  if(data.get_size_x()==0)
    data.allocate_memory(nx,ny);

  if(data.get_size_x()!=nx || data.get_size_y()!=ny ){
    cout << "The Double_2D object supplied has the wrong " 
	 << "dimensions" << endl;
    return FAILURE;
  }

  for(int i=0; i < nx; ++i){
    for(int j=0; j< ny; ++j){
      data.set(i,j,buffer[j*nx+i]);
    }
  }
  
  delete buffer;

  return SUCCESS; //success
    
}



/***************************************************************/

/***************************************************************/
int write_dbin(string file_name, const Double_2D & data){
  
  //open the input file:
  FILE * file = fopen(file_name.c_str(), "w");

  //error check.
  if(!file){
    cout << "Could not open the file " << file_name << endl;
    return FAILURE;
  }

  int nx = data.get_size_x();
  int ny = data.get_size_y();

  double * buffer = new double [nx*ny];
  for(int i=0; i < nx; ++i){
    for(int j=0; j< ny; ++j){
      buffer[j*nx+i] = data.get(i,j);
    }
  }

  fwrite(buffer, sizeof(double), nx*ny, file);
  fclose(file);
  
  delete buffer;

  return SUCCESS; //success
    
}
