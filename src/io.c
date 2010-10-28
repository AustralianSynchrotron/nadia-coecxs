#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>

//some extra mess is required because the
//hdf and tiff libraries both want to define
//the same types, giving a conflict.
typedef char mfhdf_int8;
typedef long int mfhdf_int32;
typedef long unsigned int mfhdf_uint32;

typedef signed char tiffio_int8;
typedef int tiffio_int32;
typedef unsigned int tiffio_uint32;

#define int8 mfhdf_int8
#define int32 mfhdf_int32
#define uint32 mfhdf_uint32
#include "mfhdf.h"

#undef int8
#undef int32
#undef uint32

#define int8 tiffio_int8
#define int32 tiffio_int32
#define uint32 tiffio_uint32
#include "tiffio.h"

#undef int8
#undef int32
#undef uint32

#include "io.h"

using namespace std;

#define FAILURE 0
#define SUCCESS 1

/***************************************************************/
//bad code which allows the array type to be 
//decided at run-time.
/***************************************************************/
class anonomous_array{

  int type_;
  mfhdf_int8 * a_int8;
  uint8  * a_uint8;
  int16  * a_int16;
  uint16 * a_uint16;
  mfhdf_int32 * a_int32;
  mfhdf_uint32 * a_uint32;
  float32 * a_float32;
  float64 * a_float64;
  
 public:

  anonomous_array(int type, int size){
    type_ = type;

    a_int8=0;
    a_uint8=0;
    a_int16=0;
    a_uint16=0;
    a_int32=0;
    a_uint32=0;
    a_float32=0;
    a_float64=0;

    switch( type ){
    case DFNT_INT8:
      a_int8 = new mfhdf_int8[size];
      break;
    case DFNT_UINT8:
      a_uint8 = new uint8[size];
      break;
    case DFNT_INT16:
      a_int16 = new int16[size];
      break;
    case DFNT_UINT16:
      a_uint16 = new uint16[size];
      break;
    case DFNT_INT32:
      a_int32 = new mfhdf_int32[size];
      break;
    case DFNT_UINT32:
      a_uint32 = new mfhdf_uint32[size];
      break;
    case DFNT_FLOAT32:
      a_float32 = new float32[size];
      break;
    case DFNT_FLOAT64:
      a_float64 = new float64[size];
      break;
      
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
    }    
  }
  
  ~anonomous_array(){
    switch( type_ ){
    case DFNT_INT8:
      delete[] a_int8;
      break;
    case DFNT_UINT8:
      delete[] a_uint8;
      break;
    case DFNT_INT16:
      delete[] a_int16;
      break;
    case DFNT_UINT16:
      delete[] a_uint16;
      break;
    case DFNT_INT32:
      delete[] a_int32;
      break;
    case DFNT_UINT32:
      delete[] a_uint32;
      break;
    case DFNT_FLOAT32:
      delete[] a_float32;
      break;
    case DFNT_FLOAT64:
      delete[] a_float64;
      break;
      
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
    }    
  }

  void * return_array(){
    switch( type_ ){
    case DFNT_INT8:
      return a_int8;
    case DFNT_UINT8:
      return a_uint8;
    case DFNT_INT16:
      return a_int16;
    case DFNT_UINT16:
      return a_uint16;
    case DFNT_INT32:
      return a_int32;
    case DFNT_UINT32:
      return a_uint32;
    case DFNT_FLOAT32:
      return a_float32;
    case DFNT_FLOAT64:
      return a_float64;
     default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
      return FAILURE;
    }  
  }

  double return_array_value(int i){

    //coverting to double. Not ideal.
    switch( type_ ){
    case DFNT_INT8:
      return a_int8[i];
    case DFNT_UINT8:
      return a_uint8[i];
    case DFNT_INT16:
      return a_int16[i];
    case DFNT_UINT16:
      return a_uint16[i];
    case DFNT_INT32:
      return a_int32[i];
    case DFNT_UINT32:
      return a_uint32[i];
    case DFNT_FLOAT32:
      return a_float32[i];
    case DFNT_FLOAT64:
      return a_float64[i];
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
      return FAILURE;
    }          
  }


  
};


/***************************************************************/

/***************************************************************/
int write_ppm(string file_name, int nx, int ny, double ** data, bool log_scale){

  const int pixel_maximum = 1000;
    
   //find the maximum value of the image
   double array_maximum = 0;
   for(int i=0; i < nx ; ++i){
     for(int j=0; j < ny; ++j){
       if(array_maximum < data[i][j])
	 array_maximum = data[i][j];
     }
   }
   
   double scale_factor = pixel_maximum/array_maximum;
   if(log_scale)
     scale_factor = pixel_maximum/log10(array_maximum*10);

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
       
       //if((uint) pixel_maximum*(data[i][j]/array_maximum)!=0)
       //cout << (uint) (pixel_maximum*(data[i][j]/array_maximum)) << endl;

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
int read_tiff(string file_name, int * nx, int * ny, double *** data){
 
  //open the input file:
  TIFF* tif = TIFFOpen(file_name.c_str(), "r");
  
  if (tif) {
    int dircount = 0;
    for( ; TIFFReadDirectory(tif); dircount++);
    if(dircount > 1)
      cout << "Multiple directories in the file "<<file_name
	   << "... only using the first" <<endl;
  }
  else{
    cout << "Could not open the file "<<file_name<<endl;
    return FAILURE;
  }
  
  tiffio_uint32 w, h;
  uint16 bits_per_sample;
  uint16 samples_per_pixel;

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

  int strip_size = TIFFStripSize(tif);
  int bytes_per_pixel = bits_per_sample*samples_per_pixel/8;
  int pixels_per_strip = strip_size/bytes_per_pixel;
  
  double * grey_image = new double[w*h]; 
  tiffio_uint32 * colour_image = new tiffio_uint32[w*h];
  
  if(samples_per_pixel>1){ //see if the image is colour

    cout << "Processing colour image" << endl;
    TIFFReadRGBAImage(tif, w, h, colour_image, 0);
  }
  else{ //otherwise if the image is grey scale

    cout << "Processing grey scale image" << endl;

    anonomous_array * buffer;

    switch(bytes_per_pixel){ //otherwise set up the buffer for grey scale
    case(1):
      buffer = new anonomous_array(DFNT_UINT8,pixels_per_strip);
    break;
    case(2):
      buffer = new anonomous_array(DFNT_UINT16,pixels_per_strip);
    break;
    case(4):
      buffer = new anonomous_array(DFNT_UINT32,pixels_per_strip);
      break;
    case(8):
      buffer = new anonomous_array(DFNT_FLOAT64,pixels_per_strip);
      break;
    default:
      cout << "Confused about the tiff image.." <<endl;
      return FAILURE;
    }
  

    for ( tstrip_t strip = 0; strip < TIFFNumberOfStrips(tif); strip++){
      int bytes_read = TIFFReadEncodedStrip(tif, strip, buffer->return_array(), (tsize_t) - 1);
      for( int i=0; i< bytes_read/bytes_per_pixel; i++ ){
	grey_image[i+strip*pixels_per_strip] = buffer->return_array_value(i);
      }
    }
    
  }

  //copy to the image array
  *nx=w;
  *ny=h;
    
    *data = new double*[*nx];
    for(int i=0; i < *nx; ++i){
      (*data)[i] = new double[*ny];
      for(int j=0; j< *ny; ++j){
	if(samples_per_pixel>1){//if the image is colour we take the sum of colour value
	  int pixel = colour_image[(h-j)*w+i];
	  (*data)[i][j]=TIFFGetR(pixel)+TIFFGetG(pixel)+TIFFGetB(pixel);
	}
	else //grey scale:
	  (*data)[i][j]= grey_image[j*w+i];
      }
    }

  delete[] grey_image;
  delete[] colour_image;

  TIFFClose(tif);
  
  return SUCCESS; //success
    
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


/***************************************************************/
int read_hdf4(string file_name, int * nx, int *ny, 
	      double *** data, char * data_name){

  //open the file
  mfhdf_int32 sd_id = SDstart(file_name.c_str(), DFACC_READ);
  if (sd_id == FAIL){
    cout << "Failed to open file:"<<file_name<< endl;
    return FAILURE;
  }
  
   //read the data block in the file
   mfhdf_int32 sds_index;
   sds_index = SDselect(sd_id,SDnametoindex(sd_id, data_name));
   if (sd_id == FAIL){
     cout << "Failed to find data block in file:"<<file_name<< endl;
     return FAILURE;
   }
   
   char sds_name[60];
   mfhdf_int32 rank, data_type, n_attrs;
   mfhdf_int32 dim_sizes[2];

   int status = SDgetinfo(sds_index, sds_name, &rank, 
			   dim_sizes, &data_type, &n_attrs);

   if(status){
     cout << "Failed getting data info"<< endl;
     return FAILURE;
   }
   if(rank!=2){
     cout << "Data block is not 2-dimensional"<< endl;
     return FAILURE;
   }

   anonomous_array array(data_type,dim_sizes[0]*dim_sizes[1] );

   mfhdf_int32 start[2] = {0};
   status = SDreaddata(sds_index, start, NULL, 
		       dim_sizes, array.return_array()  );
   if (status == FAIL){
     cout << "Could not get data block"<< endl;
     return FAILURE;
   }

   //close the SD readers
   status = SDendaccess(sds_index);
   status = SDend(sd_id);  


   //Fill return parameters
   *data = new double*[dim_sizes[0]];
   for(int i=0; i< dim_sizes[0]; ++i){
     (*data)[i] = new double[dim_sizes[1]];
     for(int j=0; j< dim_sizes[1]; ++j){
       (*data)[i][j] = array.return_array_value(dim_sizes[1]*j+i);
     }
   }
   
   *nx = dim_sizes[0];
   *ny = dim_sizes[1];

   return SUCCESS;

 }
 
