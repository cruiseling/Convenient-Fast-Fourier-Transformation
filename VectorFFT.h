/// @file VectorFFT.h

/////////////////////////////////
// Fast Fourier Transformation 
/////////////////////////////////

#ifndef VectorFFT_h
#define VectorFFT_h 1

// usual header
#include <complex>
#include <fftw3.h>
#include <iostream>

/// Forward FFT from real to complex.
class VectorFFT {
 private:
  fftw_plan planfor;  ///< plan for forward FFT.
  fftw_plan planback;  ///< plan for backward FFT.

  double* x_;
  fftw_complex* y_;
  int nr;
 public:
  int n_;         ///< Size of input array.
  
  /// Constructor.
  VectorFFT(int);
  /// Destructor.
  ~VectorFFT();
  /// Perform the FFT.
  void fft(std::complex<double>*, double*);
  /// Perform the inverse FFT.
  void ifft(double*, std::complex<double>*);
  /// Return the pointer of x_
  double* get_ptr_x_();
  /// Return the pointer of y_
  std::complex<double>* get_ptr_y_();

};

/// @param[in] n Size of the input array.
inline VectorFFT::VectorFFT(int n) {
  n_ = n;
  nr = ceil((double)(n + 1) / 2);
  x_ = new double[n];
  std::fill(x_, x_ + n, 0);
  y_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
  planfor = fftw_plan_dft_r2c_1d(n, x_, y_, FFTW_ESTIMATE);
  planback = fftw_plan_dft_c2r_1d(n, y_, x_, FFTW_ESTIMATE);
  return;
}

inline VectorFFT::~VectorFFT() {
  delete[] x_;
  fftw_free(y_);
  fftw_destroy_plan(planfor);
  fftw_destroy_plan(planback);
}

/// @param[in] x Input double array
/// @param[out] y Output complex array. y = FFT(x)
inline void VectorFFT::fft(std::complex<double>* y, double* x) {
	std::copy(x, x + n_, x_);
	fftw_execute(planfor);
	for (int ii = 0; ii < nr; ++ii) {
		y[ii] = std::complex<double>(y_[ii][0], y_[ii][1]);
	}
	return;
}

/// @param[out] y Output double array. y = IFFT(x)
/// @param[in] x Input complex array
inline void VectorFFT::ifft(double* y, std::complex<double>* x) {
	for (int ii = 0; ii < nr; ++ii) {
		y_[ii][0] = real(x[ii]);		
		y_[ii][1] = imag(x[ii]);
	}

	fftw_execute(planback);
	for (int ii = 0; ii < n_; ++ii) {
		y[ii] = x_[ii] / n_;
	}
  return;
}

inline double* get_ptr_x_() {
 return x_;
}

inline std::complex<double>* get_ptr_y_() {
 return y_;
}

/// Product between two complex arrays.
///
/// @param[out] z Output complex array. z = x * y.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of input array.
inline void VectorFFT_mult(std::complex<double>* z, std::complex<double>* x, 
                           std::complex<double>* y, int n) {
  for (int ii = 0; ii < n; ++ii) {
    z[ii] = x[ii] * y[ii];
  }
  return;
}

/// In-place subtract the product of two complex arrays.
///
/// @param[out] z Output complex array. z = z - x * y.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of input array.
inline void VectorFFT_diff(std::complex<double>* z, std::complex<double>* x, 
                           std::complex<double>* y, int n) {
  for (int ii = 0; ii < n; ++ii) {
	  z[ii] -= x[ii] * y[ii];
  }
  return;
}

/// In-place plus the product of two complex arrays.
///
/// @param[out] z Output complex array. z = z + x * y.
/// @param[in] x First input complex array.
/// @param[in] y Second input complex array.
/// @param[in] n Length of input array.
inline void VectorFFT_sum(std::complex<double>* z, std::complex<double>* x,
                          std::complex<double>* y, int n) {
  for (int ii = 0; ii < n; ++ii) {
	  z[ii] += x[ii] * y[ii];
  }
  return;
}
