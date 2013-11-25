#ifndef _TYPES_H_
#define _TYPES_H_
#include <fftw3.h>
#include <Complex_2D.h>
/////////////////////////////////////////////////////
#ifndef DOUBLE_PRECISION
#define FFTW_PLAN fftwf_plan
#define FFTW_COMPLEX fftwf_complex
#define FFTW_EXECUTE fftwf_execute
#define FFTW_PLAN_DFT_2D fftwf_plan_dft_2d
#define FFTW_PLAN_DFT_R2C_2D fftwf_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D fftwf_plan_dft_c2r_2d
#define FFTW_MPI_INIT fftwf_mpi_init
#define FFTW_MPI_LOCAL_SIZE_2D fftwf_mpi_local_size_2d
#define FFTW_MPI_PLAN_DFT_2D fftwf_mpi_plan_dft_2d
#define FFTW_MPI_EXECUTE_DFT fftwf_mpi_execute_dft
#define FFTW_DESTROY_PLAN fftwf_destroy_plan
#define FFTW_MALLOC fftwf_malloc
#define FFTW_FREE fftwf_free
#define MPI_MYREAL MPI_FLOAT
#else //DOUBLE
#define FFTW_PLAN fftw_plan
#define FFTW_COMPLEX fftw_complex
#define FFTW_EXECUTE fftw_execute
#define FFTW_PLAN_DFT_2D fftw_plan_dft_2d
#define FFTW_PLAN_DFT_R2C_2D fftw_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D fftw_plan_dft_c2r_2d
#define FFTW_MPI_INIT fftw_mpi_init
#define FFTW_MPI_LOCAL_SIZE_2D fftw_mpi_local_size_2d
#define FFTW_MPI_PLAN_DFT_2D fftw_mpi_plan_dft_2d
#define FFTW_MPI_EXECUTE_DFT fftw_mpi_execute_dft
#define FFTW_DESTROY_PLAN fftw_destroy_plan
#define FFTW_MALLOC fftw_malloc
#define FFTW_FREE fftw_free
#define MPI_MYREAL MPI_DOUBLE
#endif //DOUBLE_PRECISION
//////////////////////////////////////////////////

/** the function failed */
#define FAILURE 0

/** the function finished successfully */
#define SUCCESS 1

/** real component */
#define REAL 0

/** imaginary component */
#define IMAG 1

/** magnitude */
#define MAG 2

/** phase */
#define PHASE 3

/** magnitudes squared */
#define MAG_SQ 4
#endif
