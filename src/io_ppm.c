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
int write_ppm(string file_name, int nx, int ny, double ** data, bool log_scale){

  const int largest_pixel_value = pow(2.0,16.0)-1; //16 bit
  
   //find the maximum value of the image
   double array_maximum = 0;
   for(int i=0; i < nx ; ++i){
     for(int j=0; j < ny; ++j){
       if(array_maximum < data[i][j])
	 array_maximum = data[i][j];
     }
   }

   double scale_factor=1;
   int pixel_maximum=array_maximum;
   //   if(array_maximum > largest_pixel_value){
     scale_factor = ((double) largest_pixel_value)/array_maximum;
     pixel_maximum= largest_pixel_value;
     // }
   if(log_scale){
     scale_factor = ((double) largest_pixel_value)/log10(array_maximum*10);
     pixel_maximum= largest_pixel_value;
   }

   //cout << "array_maximum="<<array_maximum<<endl;
   //cout << "scale="<<scale_factor<<endl;

   //make the output file:
   ofstream new_file;
   new_file.open(file_name.c_str());
   
   //error check.
   if(!new_file.is_open()){
     cout << "Could not open the file " 
	  << file_name << "for writing" << endl;
     return FAILURE;
   }
   
   new_file << "P2" << endl;
   new_file << "#"<<file_name<<"\n";
   new_file << nx << " " << ny << endl;
   new_file << pixel_maximum << endl;

   for(int j=0; j < ny; ++j){
     for(int i=0; i < nx; ++i){
       if(log_scale){
	 if(data[i][j] >= 0.1)
	   new_file << (uint) (scale_factor*(log10((data[i][j])*10))) << " ";
	 else
	   new_file << "0 ";
       }
       else
	 new_file << (uint) (scale_factor*(data[i][j])) << " ";
     }
     new_file << endl;
   }
   new_file.close();
   
   return SUCCESS; //success
   
}


//if empty then ignore the line
//if it starts with a comment then ignore
//otherwise fill the data array with doubles.
void line_tokeniser(string line, vector<string> * data){
  
  string temp;
  std::istringstream iss(line);
  while ( getline(iss, temp, ' ') ){
    data->push_back(temp);
  }
  
  return;
  
}

/***************************************************************/

/***************************************************************/
int read_ppm(string file_name, int * nx, int * ny, double *** data){
  
  //open the input file:
  ifstream file(file_name.c_str());
  //  file.open(file_name);

  //error check.
  if(!file.is_open()){
    cout << "Could not open the file " << file_name << endl;
    return FAILURE;
  }

  string line;
  bool checked_type=false;
  bool checked_dimensions=false;
  bool checked_max=false;

  //make a temporary 1-d array;
  vector<string> string_data;
  
  while(!file.eof()){
    
    getline(file,line);
    line_tokeniser(line, &string_data);
    
    //only interested in lines which are non-empty and no
    if(string_data.size()>0 && (string_data.at(0)).at(0)!='#'){
      
      if(!checked_type){
	if((string_data.at(0)).compare("P2")!=0){
	  cout << "ppm wrong type. Only handling P2 at the moment" << endl;
	  return FAILURE; 
	}
	else{
	  checked_type=true;
	  string_data.clear();
	}
      }
      else if(!checked_dimensions){
	if(string_data.size()!=2){
	  cout << "Confused about the ppm dimensions.. exiting." << endl;
	  return FAILURE; 
	}
	else{
	  checked_dimensions=true;
	  std::istringstream(string_data.at(0)) >> *nx;
	  std::istringstream(string_data.at(1)) >> *ny;
	  string_data.clear();
	}
      }
      else if(!checked_max){
	if(string_data.size()!=1){
	  cout << "Confused about the ppm pixel maximum.. exiting." << endl;
	  return FAILURE; 
	}
	else{
	  checked_max=true;
	  string_data.clear();
	}
      }
    }
    else
      string_data.clear();
  }
  
  file.close();

  //do a sanity check
  if(string_data.size()!= (uint) (*nx)*(*ny)){
    cout << "Confused about ppm data. Dimensions"
	 << " don't match content... exiting." << endl;
    return FAILURE;
  }
  
  //fill the output values
  *data = new double*[*nx];
  for(int i=0; i < *nx; ++i)
    (*data)[i] = new double[*ny];
  

  for(int j=0,k=0; j< *ny; ++j)
    for(int i=0; i < *nx; ++i, ++k)
      std::istringstream(string_data.at(k)) >> (*data)[i][j];
  
  return SUCCESS; //success
 
}

