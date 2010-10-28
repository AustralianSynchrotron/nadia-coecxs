#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mfhdf.h"

using namespace std;

template <class T>
class ppmwriter {

  T * data_;
  int size_;
  int nx_;
  int ny_;

 public:

  ppmwriter(int nx, int ny){
    nx_ = nx;
    ny_ = ny;
    size_ = nx*ny;
    data_ = new T[size_];
  }
  
  
  ~ppmwriter(){
    delete[] data_;
  }

  T * data_pointer(){
    return data_;
  }

  void write_to_file(char * file_name){

    //find the maximum value of the image
    int max_value = 0;
    for(int i=0; i< size_ ; i++){
      if(max_value < data_[i])
	max_value = data_[i];
    }
    
    //make the output file:
    ofstream new_file;
    new_file.open(file_name);
    new_file << "P2" << endl;
    new_file << "#"<<file_name<<"\n";
    new_file << nx_ << " " << ny_ << endl;
    new_file << max_value << endl;
    for(int i=0; i< size_; i++){
      new_file << data_[i] << " ";
      if(i % nx_)
	new_file << endl;
    }
    new_file.close();
   
    cout << "Done" <<endl;
 
  }
};


/**************************************/
int main(int argc, char * argv[]){

  //check for 2 or 3 arguements
  if(argc<3 || argc>5 ){
    cout << "Wrong number of arguments" << endl;
    return 1;
  }

  //open the file
  int32 sd_id = SDstart(argv[1], DFACC_READ);
  if (sd_id == FAIL){
    cout << "Failed to open file:"<<argv[1]<< endl;
    return 1;
  }

  //read the data block in the file
  int32 sds_index;
  if(argc == 4) 
    sds_index = SDselect(sd_id,SDnametoindex(sd_id, argv[3]));
  else
    sds_index = SDselect(sd_id,SDnametoindex(sd_id, "data"));
  if (sd_id == FAIL){
    cout << "Failed to find data block in file:"<<argv[1]<< endl;
    return 1;
  }

  char sds_name[60];
  int32 rank, data_type, n_attrs;
  int32 dim_sizes[2];

  intn status = SDgetinfo(sds_index, sds_name, &rank, 
		     dim_sizes, &data_type, &n_attrs);

  if(status){
    cout << "Failed getting data info"<< endl;
    return 1;
  }
  if(rank!=2){
    cout << "Data block is not 2-dimensional"<< endl;
    return 1;
  }

  ppmwriter <uint16> mywriter(dim_sizes[0],dim_sizes[1]);
  

  //read the data into an array
  
  //  void * data2 = new uint16[dim_sizes[0]*dim_sizes[1]];
  
  

  int32 start[2] = {0};
  status = SDreaddata(sds_index, start, NULL, dim_sizes, (VOIDP) mywriter.data_pointer()); //(VOIDP)data );
  if (status == FAIL){
    cout << "Could not get data block"<< endl;
    return 1;
  }

  //close the SD readers
  status = SDendaccess(sds_index);
  status = SDend(sd_id);  

  mywriter.write_to_file(argv[2]);

  return 0;
};

