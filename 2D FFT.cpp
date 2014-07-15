// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <ctime>
#include <iostream>
using namespace std;

#define PI    M_PI    /* pi to machine precision, defined in math.h */
#define TWOPI    (2.0*PI)

/*
 FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)

 Inputs:
    data[] : array of complex* data points of size 2*NFFT+1.
        data[0] is unused,
        * the n'th complex number x(n), for 0 <= n <= length(x)-1, is stored as:
            data[2*n+1] = real(x(n))
            data[2*n+2] = imag(x(n))
        if length(Nx) < NFFT, the remainder of the array must be padded with zeros

    nn : FFT order NFFT. This MUST be a power of 2 and >= length(x).
    isign:  if set to 1, 
                computes the forward FFT
            if set to -1, 
                computes Inverse FFT - in this case the output values have
                to be manually normalized by multiplying with 1/NFFT.
 Outputs:
    data[] : The FFT or IFFT results are stored in data, overwriting the input.
*/

void four1(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
    if (j > i) {
        tempr = data[j];     data[j] = data[i];     data[i] = tempr;
        tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
    }
    j += m;
    }
    mmax = 2;
    while (n > mmax) {
    istep = 2*mmax;
    theta = TWOPI/(isign*mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
        for (i = m; i <= n; i += istep) {
        j =i + mmax;
        tempr = wr*data[j]   - wi*data[j+1];
        tempi = wr*data[j+1] + wi*data[j];
        data[j]   = data[i]   - tempr;
        data[j+1] = data[i+1] - tempi;
        data[i] += tempr;
        data[i+1] += tempi;
        }
        wr = (wtemp = wr)*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
    }
}

/********************************************************
* The following is a test routine that generates a ramp *
* with 10 elements, finds their FFT, and then finds the *
* original sequence using inverse FFT                   *
********************************************************/

int main(int argc, char * argv[])
{
    int i;
    int j;
    int h;
    int Nx;
    int NFFT;
    double input[16384]={0};
    double output[2*16384]={0};
    double final[2*16384]={0};
    double *y;
    double *Y;
    double *x;
    double *X;

    /* generate a ramp with 10 numbers */
    Nx = 128;

        for(i=0; i<16384;i++)
        {
            input[i]=i;
        }

        float t=clock();

    x=&input[0];

    /* calculate NFFT as the next higher power of 2 >= Nx */
    NFFT = (int)pow(2.0, ceil(log((double)Nx)/log(2.0)));
    //printf("NFFT = %d\n", NFFT);

    /* allocate memory for NFFT complex numbers (note the +1) */
    X = (double *) malloc((2*NFFT+1) * sizeof(double));

    /* Storing x(n) in a complex array to make it work with four1. 
    This is needed even though x(n) is purely real in this case. */
    for(j=0; j<Nx; j++){
        for(i=0; i<Nx; i++)
        {
            X[2*i+1] = *(x+i+NFFT*j);
            X[2*i+2] = 0.0;
        }
    /* pad the remainder of the array with zeros (0 + 0 j)
    for(i=Nx; i<NFFT; i++)
    {
        X[2*i+1] = 0.0;
        X[2*i+2] = 0.0;
    }*/





    /* calculate FFT */
        four1(X, NFFT, 1);

    //printf("\nFFT:\n");
    /*for(i=0; i<NFFT; i++)
    {
        printf("X[%d] = (%.2f + j %.2f)\n", (i+Nx*j), X[2*i+1], X[2*i+2]);
    }*/
        
        for(h=0;h<(2*NFFT);h++){
            output[h+(2*NFFT)*j]=*(X+h+1);
        }

    }

    /*    for(i=0; i<(4*NFFT); i++)
    {
        printf("Output[%d] = (%.2f + j %.2f)\n", i, output[2*i], output[2*i+1]);
    }*/

    
        y=&output[0];

        Y = (double *) malloc((2*NFFT+1) * sizeof(double));


    for(j=0;j<NFFT;j++)
    {
        for(i=0;i<NFFT;i++)
        {
            
        Y[2*i+1] = *(y+(2*NFFT)*i+2*j);
        Y[2*i+2] = *(y+(2*NFFT)*i+2*j+1);
        }
        four1(Y, NFFT, 1);
        /*for(i=0; i<NFFT; i++)
        {
            printf("Y[%d] = (%.2f + j %.2f)\n", i, Y[2*i+1], Y[2*i+2]);
        }*/

        for(h=0;h<NFFT;h++)
        {
            final[(2*NFFT)*h+2*j]=Y[2*h+1];
            final[(2*NFFT)*h+2*j+1]=Y[2*h+2];
        }

    }

        /*for(i=0; i<128; i++)
        {
            printf("Final[%d] = (%.2f + j %.2f)\n", i, final[2*i], final[2*i+1]);
        }*/
    


    t=clock()-t;

    cout<<"Time(mil secs) :"<<t<<endl;
    
    char key;
    cin>>key;
    return 0;
}