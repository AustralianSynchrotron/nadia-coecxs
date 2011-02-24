#include <string>
#include <cmath>

class Complex_2D;
class Double_2D;

/**
 * Read a ppm file. Returns a 2D of the data. 
 *
 * @param file_name The name of the file to read from
 * @param data The array to be filled with data
 */
int read_ppm(std::string file_name, Double_2D & data);

/**
 * Read a tiff file. Returns a 2D of the data.
 * 
 * @param file_name The name of the file to read from
 * @param data The array to be filled with data
 */
int read_tiff(std::string file_name, Double_2D & data);

/**
 * Read a HDF4 file. Returns a 2D of the data. 
 *
 * @param file_name The name of the file to read from
 * @param data The array to be filled with data
 * @param data_name The name of the block in the HDF4 where the data is.
 * By default it looks for the "data" block.
 */
int read_hdf4(std::string file_name, Double_2D & data, 
	      char * data_name="data");



int read_dbin(std::string file_name, int nx, int ny, Double_2D & data);

int read_cplx(std::string file_name, Complex_2D & complex);



/** 
 * Write a 2D array to a ppm file. The data will be saved as a 16 bit
 * grey-scale image. Please note that the data will be scaled to fit
 * within the range 0 - 2^16. Hence this method is only useful for
 * viewing the final output and should not be used if you plan to
 * reopen the file at a later date and pass the data back to the
 * reconstruction algorithm.
 *
 * @param file_name The name of the file to write to
 * @param data The array to be written to file
 * @param log_scale Output on log scale? true/false. Default is false.
 */ 
int write_ppm(std::string file_name, const Double_2D & data, 
	      bool log_scale=false);

/** 
 * Write a 2D array into a binary file. The data will be saved in 64
 * bit format. Because of their high precision and because they are
 * not scaled, files produces in this way can used to save data from
 * the reconstruction which you plan to reload into the algorithm at a
 * later date.
 *
 * @param file_name The name of the file to write to
 * @param data The array to be written to file
 */ 
int write_dbin(std::string file_name, const Double_2D & data);

/** 
 * Write a 2D complex array into a binary file (the array of type
 * fftw3 is written out). This method is useful for saving data from
 * the reconstruction which you plan to reload into the algorithm at a
 * later date. For example saving the current estimate of the object
 * exist surface wave for later use, saving the white field
 * reconstructed from 3-plane propogation for use in Fresnel
 * reconsturction, etc.
 *
 * @param file_name The name of the file to write to 
 * @param data The complex array to be written to file
 */ 
int write_cplx(std::string file_name, const Complex_2D & complex);

/** 
 * Write a 2D array to a tiff file. The data will be saved as a 16 bit
 * grey-scale image. Please note that the data will be scaled to fit
 * within the range 0 - 2^16. Hence this method is only useful for
 * viewing the final output and should not be used if you plan to
 * reopen the file at a later date and pass the data back to the
 * reconstruction algorithm.
 *
 * @param file_name The name of the file to write to
 * @param data The array to be written to file
 * @param log_scale Output on log scale? true/false. Default is false.
 */ 
int write_tiff(std::string file_name, const Double_2D & data);


//used to transform an array of doubles from -x_min .... x_max
//to fall between 0... pixel_max (which is usually 2^16).
inline unsigned int io_scale_value(double min, double max, 
			  int pixel_max, 
			  double value, bool log_scale){
  
  double grad;

  if(log_scale){

    //first scale to be between 1 and some large number (no. of pixels?):
    grad = (double)(pixel_max-1)/(double)(max-min);
    value = grad*(value-min)+1;
    
    //adjust to log scale
    min = 0;
    max = log10(pixel_max);
    value = log10(value);
  }
  
  grad = (double)pixel_max/(double)(max-min);
  return grad*(value-min);
  
}

//int write();
//int read();
