


void apply_shrinkwrap(int nx, int ny, double *** recon,
		      double gauss_width, double threshold);

void apply_threshold(int nx, int ny, double *** array, 
	       double threshold);

void convolve(int nx, int ny, double *** array, double gauss_width);

double gauss_2d(double x, double y, 
		double sigma_x, double sigma_y, 
		double amp);
