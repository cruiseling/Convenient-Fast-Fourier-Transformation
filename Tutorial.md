# Fast Fourier Transformation (FFT) and Inverse Transformation (IFFT)

For a length $n$ real vector 

$$\bm x = [x_1,\ldots, x_n],$$

to compute its FFT 

$$y = \text{FFT}(x),$$

where $y$ is a length $n$ complex vector, we can run the following code

``` cpp
int n;                                           /// vector length
double* x;                                       /// real array of x
complex<double>* y;                              /// complex array that stores FFT(x)

VectorFFT vf = new VectorFFT(n);                 
vf->fft(y, x);                                   /// get the result in y
```



For a length $n$ complex vector

$$\bm x = [x_1,\ldots, x_n],$$

to compute its IFFT

$$ y = \text{IFFT}(x),$$

where $y$ is a length $n$ real vector, we can run the following code

``` cpp
int n;                                           /// vector length
complex<double>* x;                              /// complex array of x
double* y;                                       /// real array that stores IFFT(x)

VectorFFT vf = new VectorFFT(n);                
vf->ifft(y, x);                                  /// get the result in y
```

# Fast Polynomial Convolution

For two polynomials 

$$ a(x) = a_0 + a_1 x + \ldots + a_p x^p, \qquad b(x) = b_0 + b_1 x + \ldots + b_q x^q ,$$

their product results in a new polynomial

$$ c(x) = a(x) \times b(x) = c_0 + c_1 x + \ldots + c_n x^n, \quad n = p+q. $$

The computation of coefficients $\bm c = [c_0, c_1,\ldots, c_n]$, which is the convolution between vectors $\bm a = [a_0, a_1,\ldots, a_p]$ and $\bm b = [b_0, b_1,\ldots, b_q]$, can be obtained in two approaches: the usual method which takes $O(n^2)$ steps, and the FFT approach which takes only $O(n \log n)$ steps. The fast FFT approach can be easily implemented with our `VectorFFT` class.

``` cpp
int p, q;                                                   /// be aware that length of vector a is p+1, length of vector b is q+1, length of vector c is p+q+1
double *a, *b, *c;                                          /// real array of a, b and array for storing array c

double* a0 = new double[p+q+1];
copy(a, a+p+1, a0);
fill(a0+p+1, a0+p+q+1, 0);
complex<double>* a_fft = new complex<double>[p+q+1];
double* b0 = new double[p+q];
copy(b, b+q+1, a0);
fill(b0+q+1, b0+p+q+1, 0);
complex<double>* b_fft = new complex<double>[p+q+1];
complex<double>* c_fft = new complex<double>[p+q+1];
VectorFFT vf = new VectorFFT(n);
vf->fft(a_fft, a0);
vf->fft(b_fft, b0);
VectorFFT_mult(c_fft, a_fft, b_fft, p+q+1);
vf->ifft(c, c_fft);                                          /// get the result in c
```
