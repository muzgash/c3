#ifndef C3_H
#define C3_H

#include<math.h>
#include<fftw3.h>
#include<stdlib.h>



void wavnum(int N,double deltat,double *w);
double gaussian(double x,double a,double c);
double ricker(double x,double to, double vfrq);
double xricker(double x, double xo, double t, double to, double f, double c);
void fftc2c(fftw_complex *in,fftw_complex *out, int Nx);
void ifftc2c(fftw_complex *in,fftw_complex *out, int Nx);
void fftc2c2d(fftw_complex *in,fftw_complex *out, int Nz, int Nx);
void ifftc2c2d(fftw_complex *in,fftw_complex *out, int Nz, int Nx);
void double2fftw_complex(double *real, fftw_complex *out, int size);
void transpose(float *P, float *Pt,int Nx, int Nt);
void float2fftw_complex(float *real, fftw_complex *out, int size);
void rickerkw(int Nx, int Nt, double xo, double dx, double dt,double c,double f, fftw_complex *out);

#endif
