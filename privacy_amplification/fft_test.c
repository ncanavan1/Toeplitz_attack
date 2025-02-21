#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795  // Value of Pi 

void DFT(double complex *x, double complex *X, int N){
  
  for(int n = 0; n < N; n++){
    for(int k = 0; k < N; k++){
      double complex e_val = cexp(-2*I*PI*k*n/N);
      X[k] = X[k] + x[n]*e_val;
    }
  }
}


void IDFT(double complex *x, double complex *X, int N){
  
  for(int n = 0; n < N; n++){
    for(int k = 0; k < N; k++){
      double complex e_val = cexp(-2*I*PI*k*n/N)/N;
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

    
    
    double complex factor;
    for (int k = 0; k < N / 2; k++) {
        factor = cexp(-2.0 * I * PI * k / N);  // Twiddle factor e^(-2πik/N)
        X[k] = X_even[k] + factor * X_odd[k];    // First half
        X[k + N / 2] = X_even[k] - factor * X_odd[k]; // Second half
    }
  }
}


void IFFT_recursive(double complex *x, double complex *X, int N){
  //should pad before FFT initiates to N = power of 2

  if (N <= 16){
    IDFT(x, X, N);
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

    
    
    double complex factor;
    for (int k = 0; k < N / 2; k++) {
        factor = cexp(-2.0 * I * PI * k / N);  // Twiddle factor e^(-2πik/N)
        X[k] = X_even[k] + factor * X_odd[k];    // First half
        X[k + N / 2] = X_even[k] - factor * X_odd[k]; // Second half
    }
  }
}





void embed_toeplitz_on_circulant(double complex *circ, int* row1, int* col1, int rowLen, int colLen, int circulant_size){
  int m = 0;
  for(int i = 0; i < rowLen; i++){
    circ[m] = row1[i] + 0*I;
    m++;
  }
  for(int i = colLen-1; i >= 1; i--){
    circ[m] = col1[i] + 0*I;
    m++;
  }
  while (m < circulant_size)
  {
    circ[m] = 0 + 0*I;
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

int calc_circulant_size(int rowLen, int colLen){
  int circulant_size = rowLen + colLen - 1;
  if(circulant_size % 2 > 0){
    int p = 0;
    while(pow(2,p) < circulant_size){
      p++;
    }
    circulant_size = pow(2,p);
  }
}

// Function to perform FFT (Recursive)
void fft(complex double *X, int n) {
  if (n <= 1) return;

  // Divide
  complex double *X_even = malloc(n / 2 * sizeof(complex double));
  complex double *X_odd = malloc(n / 2 * sizeof(complex double));

  for (int i = 0; i < n / 2; i++) {
      X_even[i] = X[i * 2];
      X_odd[i] = X[i * 2 + 1];
  }

  // Recursion
  fft(X_even, n / 2);
  fft(X_odd, n / 2);

  // Combine
  for (int k = 0; k < n / 2; k++) {
      complex double t = cexp(-2.0 * I * PI * k / n) * X_odd[k];
      X[k] = X_even[k] + t;
      X[k + n / 2] = X_even[k] - t;
  }

  free(X_even);
  free(X_odd);
}

// Function to perform Inverse FFT (Recursive)
void ifft(complex double *X, int n) {
  // Conjugate the complex numbers
  for (int i = 0; i < n; i++) {
      X[i] = conj(X[i]);
  }

  // Forward FFT
  fft(X, n);

  // Conjugate again and scale
  for (int i = 0; i < n; i++) {
      X[i] = conj(X[i]) / n;
  }
}



void Toeplitz_hash_fft(int* input_key, int* output_key, int* row1, int *col1, int rowLen, int colLen){

    int circulant_size = calc_circulant_size(rowLen, colLen);
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


    //FFT_recursive(circ, v, circulant_size);
    //FFT_recursive(padded_key, y, circulant_size);
    fft(circ,circulant_size);
    fft(padded_key,circulant_size);
    for(int i = 0; i < circulant_size; i++){
      u[i] = circ[i] * padded_key[i];
    }
    
    //IFFT_recursive(u, cx, circulant_size);
    ifft(u,circulant_size);


    for(int i = 0; i < colLen; i++){
      int bit_rounded = round(creal(cx[i]));
      output_key[i] = bit_rounded%2; 
    }
}

int main(){
  int row1[] = {0,1,1,0,1};
  int col1[] = {0,1,0};
  int rowLen = 5;
  int colLen = 3;

  int input_key[] = {1,1,0,1,1};
  int output_key[3];

  Toeplitz_hash_fft(input_key, output_key, row1, col1, rowLen, colLen);
  printf("Done");

  complex double test[] = {1,1,1,0,0,1,1,1};
  fft(test, 8);
  ifft(test, 8);
  printf("Done");
}










int main1(void){
  
  int N = 32;
  int x_int[] = {0,1,1,0,1,1,0,1,1,1,0,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,0,0,1,0,0,0};
  double complex x[32];
  double complex X[N];
  double complex X_op[N];

  for(int i=0; i < N; i++){
    x[i] = x_int[i] + 0*I;
    X[i] = 0 + 0*I;
    X_op[i] = 0 + 0*I;
  }


  FFT_recursive(x, X, N);
  IFFT_recursive(X, X_op, N);
  for(int i = 0; i < N; i ++){
    X_op[i] = X_op[i]/N;
    int rounded_mod_x = round(creal(X_op[i]));
    printf("%d,", rounded_mod_x%2);
  }



  printf("Done");

}