#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>
#include "io.h"

using namespace std;

#define FAILURE 0
#define SUCCESS 1

/***************************************************************/

/***************************************************************/
int read_dbin(string file_name, int nx, int ny, double *** data){
 
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
  
  *data = new double*[nx];
  for(int i=0; i < nx; ++i){
    (*data)[i] = new double[ny];
    for(int j=0; j< ny; ++j){
      (*data)[i][j] = buffer[j*nx+i];
    }
  }
  
  delete buffer;

  return SUCCESS; //success
    
}

