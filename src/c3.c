#include<math.h>
#include<c3.h>
#include<stdlib.h>
#include<stdio.h>



//=======================================================================
/*__________________________   Frecuencias-Numeros de onda   ___________________________*/
void wavnum(int N,double deltat,double *w)
{
  int i;
  for(i=0; i<N/2+1; i++)   /*positive frequencies*/ 
  {
    w[i] = i*2*M_PI/(N*deltat);
  }
  for(i=N/2+1; i<N; i++)   /*negative frequencies*/ 
  {
    w[i] = (i-N)*2*M_PI/(N*deltat);
  }
}

double gaussian(double x,double a,double c)
{
	double result;
	result=a*exp(-(x*x)/(2*c*c));     
	return result;
}

void float2fftw_complex(float *real, fftw_complex *out, int size)
{
    int i;
	for(i=0; i<size;i++) {
	    out[i][0]=real[i];
		out[i][1]=0.0;
	}
}

void double2fftw_complex(double *real, fftw_complex *out, int size)
{
    int i;
	for(i=0; i<size;i++) {
	    out[i][0]=real[i];
		out[i][1]=0.0;
	}
}

void transpose(float *P,float *Pt, int N1, int N2)
{
    int i,j;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
		    Pt[j*N1+i]=P[i*N2+j];
		    Pt[j*N1+i]=P[i*N2+j];
		}
	}
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
		    P[i*N2+j]=P[i*N2+j];
		    P[i*N2+j]=P[i*N2+j];
		}
	}
}

void fftc2c(fftw_complex *in,fftw_complex *out, int Nx)
{
	fftw_plan plan_forward;
	plan_forward = fftw_plan_dft_1d (Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute ( plan_forward );
	fftw_destroy_plan ( plan_forward );
}

void ifftc2c(fftw_complex *in,fftw_complex *out, int Nx)
{
	fftw_plan plan_backward;
	plan_backward = fftw_plan_dft_1d (Nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute ( plan_backward );
	fftw_destroy_plan ( plan_backward );
}

void fftc2c2d(fftw_complex *in,fftw_complex *out, int Nz, int Nx)
{
	fftw_plan plan_forward;
	plan_forward = fftw_plan_dft_2d (Nz, Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute ( plan_forward );
	fftw_destroy_plan ( plan_forward );
}

void ifftc2c2d(fftw_complex *in,fftw_complex *out, int Nz, int Nx)
{
	fftw_plan plan_backward;
	plan_backward = fftw_plan_dft_2d (Nz, Nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute ( plan_backward );
	fftw_destroy_plan ( plan_backward );
}
 
//=======================================================================
/*__________________________   pulso de Ricker   ___________________________*/


double ricker(double t,double to, double vfrq)
{
	double result;
	double alpha;
	alpha=M_PI*M_PI*vfrq*vfrq*(t-to)*(t-to);
	result=(1-2*alpha)*exp(-alpha);     
	return result;
}

double xricker(double x, double xo, double t, double to, double f, double c)
{
    return ricker(t,to,f)*exp(-0.5*pow((x-xo)/c,2));
}

void rickerkw(int Nx, int Nt, double xo, double dx, double dt,double c,double f, fftw_complex *out)
{
    int i,j;
    fftw_complex *in;
    in = fftw_malloc(Nx*Nt*sizeof(fftw_complex));
	//fprintf(stderr,"%f\n",f);
	for(j=0;j<Nt;j++){
	    for(i=0;i<Nx;i++) {
		    in[j*Nx+i][0]=xricker(i*dx,xo,j*dt,0.1,f,c);
			in[j*Nx+i][1]=0.0;
			//printf("%d %d %f\n",i,j,in[j*Nx+i][0]);
		}
		//printf("\n");
	}
	fftc2c2d(in,out,Nt,Nx);
	/*for(j=0;j<Nt;j++){
	    for(i=0;i<Nx;i++) {
			printf("%d %d %f\n",i,j,in[j*Nx+i][0]);
		}
		printf("\n");
	}*/
	fftw_free(in);
}

