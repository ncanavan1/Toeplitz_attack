#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795  // Value of Pi 

void DFT(double complex *x, double complex *X, int N){
  
  for(int n = 0; n < N; n++){
    for(int k = 0; k < N; k++){
      double complex e_val = cexp(-2*I*PI*n*k/N);
      X[k] = X[k] + x[n]*e_val;
    }
  }
}



void slice_arr_C(double complex *arr1, double complex *arr2, int p1, int N){
  int j = 0;
  for(int i = p1; i < 2*N; i +=2){
    arr2[j] = arr1[i];
    j++;
  }
}

void FFT_recursive(double complex *x, double complex *X, int N){
  //should pad before FFT initiates to N = power of 2

  if (N <= 16){
    DFT(x, X, N);
  }

  else{
    double complex x_even[N/2];
    double complex x_odd[N/2];
    slice_arr_C(x, x_even, 0, N/2);
    slice_arr_C(x, x_odd, 1, N/2);


    double complex X_even[N/2];
    double complex X_odd[N/2];


    //FFT wont work without thisbut it seems unneeded
    for(int i=0; i <N/2; i++){
      X_even[i] = 0 + 0*I;
      X_odd[i] = 0 + 0*I;
    }


    FFT_recursive(x_even, X_even, N/2);
    FFT_recursive(x_odd, X_odd, N/2);

    
    
    double complex wn = cexp(-2*I * PI / N);
    double complex w = 1;
    for (int k = 0; k < N / 2; k++) {
        X[k] = X_even[k] + w*X_odd[k];    // First half
        X[k + N / 2] = X_even[k] - w*X_odd[k]; // Second half
        w = w*wn;
    }
  }
}


void conj_array(double complex *arr, int N){
  for(int i = 0; i < N; i++){
    arr[i] = conj(arr[i]);
  }
}


void embed_toeplitz_on_circulant(double complex *circ, int* row1, int* col1, int rowLen, int colLen, int circulant_size){
  int m = 0;
  for(int i = 0; i < colLen; i++){
    circ[m] = col1[i] + 0*I;
    m++;
  }
  for(int i = rowLen-1; i >= 1; i--){
    circ[m] = row1[i] + 0*I;
    m++;
  }  
}

void pad_key(complex double *padded_key, int *key, int padLen, int keyLen){
  for(int i = 0; i < keyLen; i++){
    padded_key[i] = key[i] + 0*I;
  }
  for(int i = keyLen; i < padLen; i++){
    padded_key[i] = 0 + 0*I;
  }
}




void Toeplitz_hash_fft(int* input_key, int* output_key, int* row1, int *col1, int rowLen, int colLen){

    int circulant_size = rowLen + colLen - 1;
    double complex circ[circulant_size];
    double complex padded_key[circulant_size];
    double complex v[circulant_size];
    double complex y[circulant_size];
    double complex u[circulant_size];
    double complex cx[circulant_size];

    for(int i=0; i < circulant_size; i++){
      circ[i] = 0 + 0*I;
      padded_key[i] = 0 + 0*I;
      v[i] = 0 + 0*I;
      y[i] = 0 + 0*I;
      u[i] = 0 + 0*I;
      cx[i] = 0 + 0*I;

    }

    embed_toeplitz_on_circulant(circ, row1, col1, rowLen, colLen, circulant_size); //embed to next power of 2
    pad_key(padded_key, input_key, circulant_size, rowLen);

    FFT_recursive(circ, v, circulant_size);
    FFT_recursive(padded_key, y, circulant_size);
    for(int i = 0; i < circulant_size; i++){
      u[i] = v[i] * y[i];
    }

    conj_array(u,circulant_size);
    FFT_recursive(u, cx, circulant_size);
    conj_array(cx,circulant_size);

    for(int i = 0; i < colLen; i++){
      cx[i] = cx[i]/circulant_size ;
      int bit_rounded = round(creal(cx[i]));
      output_key[i] = bit_rounded%2; 
    }
}

int main(){
  int row1[] = {1,1,0,0,1};
  int col1[] = {1,0,0,1};
  int rowLen = 5;
  int colLen = 4;

  int input_key[] = {0,0,1,1,0};
  int output_key[colLen];

  Toeplitz_hash_fft(input_key, output_key, row1, col1, rowLen, colLen);
  for(int i = 0; i < colLen; i++){
    printf("%d,",output_key[i]);
  }
}