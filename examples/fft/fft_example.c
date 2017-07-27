#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<c3.h>

int main()
{
    int i,j;
	fftw_complex *in=fftw_malloc(100*sizeof(fftw_complex));
	fftw_complex *out=fftw_malloc(100*sizeof(fftw_complex));

    FILE *fft=fopen("fft.txt","w");
	for(i=0;i<100;i++){
	    in[i][0]=sin(3*M_PI*i/700.0)*cos(i);
	    in[i][1]=0.0;
	    out[i][0]=0.0;
	    out[i][1]=0.0;
	}
	//Realizamos la transformada de Fourier de complejos a complejos
	fftc2c(in+25,out+50,50);

	//Se escribe a disco
	for(i=0;i<100;i++){
	    fprintf(fft,"%d %lf\n",i,out[i][0]);
	}
	fclose(fft);
    system("gnuplot -e \"p 'fft.txt' w l; pause -1\"");
    return 0;
}
