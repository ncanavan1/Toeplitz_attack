#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.1415926535897932384626433832795  // Value of Pi 


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


void conj_array(double complex *arr, int N){
  for(int i = 0; i < N; i++){
    arr[i] = conj(arr[i]);
  }
}


int reverse_bit(int num, int s) {
    int res = 0;
    for (int i = 0; i < s; i++) {
        if (num & (1 << i))
            res |= 1 << (s - 1 - i);
    }
    return res;
}

void reverse_array(double complex *X, int N, int s){
    for(int i = 0; i < N; i ++){
        int rb = reverse_bit(i,s);
        if(i < rb){
            double complex tmp = X[i];
            X[i] = X[rb];
            X[rb] = tmp;
        }
    }
}


//https://github.com/Swati-Verma671/Computation-of-DFT-using-Radix-2-DIT-FFT-algorithm/blob/main/code.m
void DIT_FFT(double complex *X, int N){
    int s = round(log2(N));
    reverse_array(X,N,s);
    for(int stage=1; stage <= s; stage++){
        int p = 0;
        int q = 0 + pow(2,stage-1);
        int n = 0;
        while (n<=pow(2,(stage-1))-1 && q <=N)
        {
            double complex w = cexp(-2*I*PI*n/(pow(2,stage)));
            double complex y = X[p] + w*X[q];
            double complex z = X[p] - w*X[q];
            X[p] = y;
            X[q] = z;
            p++;
            q++;
            n++;
            if(q % (int)pow(2,stage) == 0){
                p = p + pow(2,stage-1);
                q = q + pow(2,stage-1);
                n=0;
            }
        }
    }
}

void DIT_IFFT(double complex *X, int N){
    conj_array(X, N);
    DIT_FFT(X,N);
    conj_array(X,N);

    for(int i = 0; i < N; i++){
      X[i] = X[i]/N;
    }
}


void Toeplitz_hash_fft_butterfly(int* input_key, int* output_key, int* row1, int *col1, int rowLen, int colLen){

    int circulant_size = rowLen + colLen - 1;
    double complex circ[circulant_size];
    double complex padded_key[circulant_size];
    double complex u[circulant_size];
    double complex cx[circulant_size];

    for(int i=0; i < circulant_size; i++){
      circ[i] = 0 + 0*I;
      padded_key[i] = 0 + 0*I;
      u[i] = 0 + 0*I;
      cx[i] = 0 + 0*I;

    }

    embed_toeplitz_on_circulant(circ, row1, col1, rowLen, colLen, circulant_size); //embed to next power of 2
    pad_key(padded_key, input_key, circulant_size, rowLen);

    
    DIT_FFT(circ,circulant_size);
    DIT_FFT(padded_key,circulant_size);
    for(int i = 0; i < circulant_size; i++){
      u[i] = circ[i] * padded_key[i];
    }

    DIT_IFFT(u, circulant_size);


    for(int i = 0; i < colLen; i++){
      int bit_rounded = round(creal(u[i]));
      output_key[i] = bit_rounded%2; 
    }
    int k=7;

}

void naive_mult_result(int *row1, int*col1, int *input_key, int * output_key, int rowLen, int colLen){

  printf("\n\n");

  for(int i = 0; i < colLen; i++){
    output_key[i] = 0;
    for(int j = 0; j < rowLen; j++){
      if (i >= j){
        printf("%d,", col1[i-j]);
        output_key[i] = output_key[i] + (col1[i-j]*input_key[j]);
      }
      else{
        printf("%d,", row1[j-i]);
        output_key[i] = output_key[i] + (row1[j-i]*input_key[j]);
      }
    }
    printf("\n");
  }

    printf("Final Result\n");
  for(int i = 0; i < colLen; i++){
    output_key[i] = (int)output_key[i]%2;
    printf("%d,", output_key[i]);
  }
}


void generate_binary_list(int arr[], int length) {
    for (int i = 0; i < length; i++) {
        arr[i] = rand() % 2;  // Generates either 0 or 1
    }
}

int main(){

    //double complex test[] = {0,1,2,3,4,5,6,7};
    //reverse_array(tes1,0,1,0,0,1,0,1,0,1t,8,3);

    int rowLen = 80;
    int colLen = 49;
    int row1[rowLen];
    int col1[colLen];
    int input_key[rowLen];
    generate_binary_list(row1, rowLen);
    generate_binary_list(col1, colLen);
    generate_binary_list(input_key, rowLen);

    //int row1[] = {1,0,1,1,1,1,0,0,1,0};
    //int col1[] = {1,0,1,1,0,0,0};
    //int input_key[] = {0,1,1,1,0,1,0,1,1,1};
    int output_key_fft[colLen];
    int output_key_naive[colLen];



    Toeplitz_hash_fft_butterfly(input_key, output_key_fft, row1, col1, rowLen, colLen);
    naive_mult_result(row1, col1, input_key, output_key_naive, rowLen, colLen);

    int err = 0;
    for(int i = 0; i < colLen; i++){
      if(output_key_fft[i] != output_key_naive[i]){
        printf("Error on Bit %d \n",i);
        err++;
      }
    }
    if(err == 0){
      printf("\nFFT matches Naive Implimentation\n");
    }
}


