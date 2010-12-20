#include <string>

class Complex_2D;
class Double_2D;

/** Write a 2D array to a ppm file */
//int write_ppm(std::string file_name, Double_2D data, bool log_scale=false);
int write_ppm(std::string file_name, const Double_2D & data, bool log_scale=false);

/** Read a ppm file. Returns a 2D of the data. */
int read_ppm(std::string file_name, Double_2D & data);

/** Read a tiff file. Returns a 2D of the data. */
int read_tiff(std::string file_name, Double_2D & data);

/** Read a HDF4 file. Returns a 2D of the data. */
int read_hdf4(std::string file_name, Double_2D & data, 
	      char * data_name="data");


int read_dbin(std::string file_name, int nx, int ny, Double_2D & data);

int write_dbin(std::string file_name, const Double_2D & data);


int read_cplx(std::string file_name, Complex_2D & complex);

int write_cplx(std::string file_name, const Complex_2D & complex);

