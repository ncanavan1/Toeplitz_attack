#include "hal.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define IDLE 0
#define KEY 1
#define PLAIN 2
#define PI 3.1415926535897932384626433832795  // Value of Pi 
#define BUFLEN 64

uint8_t memory[BUFLEN];
uint8_t tmp[BUFLEN];
char asciibuf[BUFLEN];
uint8_t pt[16];

static void delay_2_ms(void);


void my_puts(char *c)
{
  do {
    putch(*c);

  } while (*++c);
}

static void delay_2_ms()
{
  for (volatile unsigned int i=0; i < 0xfff; i++ ){
    ;
  }
}

void my_read(char *buf, int len)
{
  for(int i = 0; i < len; i++) {
    while (buf[i] = getch(), buf[i] == '\0');

    if (buf[i] == '\n') {
      buf[i] = '\0';
      return;
    }
  }
  buf[len - 1] = '\0';
}


int char_to_int(char *c_arr, int len){
  int num = 0;
  for(int i = 0; i < len; i++){
    num = num + (c_arr[i] - '0')*pow(10,i);
  }
  return num;
}

void int_to_char_arr(int *i_arr, char *c_arr, int len){
  for(int i = 0; i < len; i++){
    c_arr[i] = i_arr[i] + '0';
  }
}

void long_delay(){
  for(int i = 0; i < 1000; i++){
    asm volatile(
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    "nop"       "\n\t"
    ::
    );
  }
}


void print_int_array_serial(const int *arr, int len) {
    char buffer[16]; // Temporary buffer to hold string representations of integers
    for (int i = 0; i < len; i++) {
        // Convert integer to string
        sprintf(buffer, "%d ", arr[i]);
        
        // Send the string over serial
        for (char *ptr = buffer; *ptr != '\0'; ptr++) {
            putch(*ptr); // Send each character
        }
    }
    putch('\n'); // Add a newline for better readability
}


////////////////////////
/// FFT STUFF //////////

void embed_toeplitz_on_circulant(float complex *circ, int* row1, int* col1, int rowLen, int colLen, int circulant_size){
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

void pad_key(complex float *padded_key, int *key, int padLen, int keyLen){
  for(int i = 0; i < keyLen; i++){
    padded_key[i] = key[i] + 0*I;
  }
  for(int i = keyLen; i < padLen; i++){
    padded_key[i] = 0 + 0*I;
  }
}


void conj_array(float complex *arr, int N){
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

void reverse_array(float complex *X, int N, int s){
    for(int i = 0; i < N; i ++){
        int rb = reverse_bit(i,s);
        if(i < rb){
          float complex tmp = X[i];
            X[i] = X[rb];
            X[rb] = tmp;
        }
    }
}


//https://github.com/Swati-Verma671/Computation-of-DFT-using-Radix-2-DIT-FFT-algorithm/blob/main/code.m
void DIT_FFT(float complex *X, int N){

    //int collected = 0;

    int s = round(log2(N));
    reverse_array(X,N,s);
    trigger_high();
    for(int stage=1; stage <= s; stage++){
        int p = 0;
        int q = 0 + pow(2,stage-1);
        int n = 0;
        while (n<=pow(2,(stage-1))-1 && q <=N)
        {

            if (stage == 2){
              trigger_low();
            }

            float complex w = 1;//cexp(-2*I*PI*n/(pow(2,stage)));

            float complex y = X[p];
            float complex z = X[q];


            // for(int i = 0; i < 20; i++){
            //     __asm__ volatile ("nop");
            // }

            z *= w;
            X[p] = y+z;
            X[q] = y-z;


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

void DIT_IFFT(float complex *X, int N){
    conj_array(X, N);
    DIT_FFT(X,N);
    conj_array(X,N);

    for(int i = 0; i < N; i++){
      X[i] = X[i]/N;
    }
}


void Toeplitz_hash_fft_butterfly(int* input_key, int* output_key, int* row1, int *col1, int rowLen, int colLen){

    int circulant_size = rowLen + colLen - 1;
    float complex *circ = (float complex*)malloc(circulant_size * sizeof(float complex));
    float complex *padded_key = (float complex*)malloc(circulant_size * sizeof(float complex));
   // double complex u[circulant_size];
   // double complex cx[circulant_size];

    for(int i=0; i < circulant_size; i++){
      circ[i] = 0 + 0*I;
      padded_key[i] = 0 + 0*I;
    //  u[i] = 0 + 0*I;
    //  cx[i] = 0 + 0*I;

    }

    embed_toeplitz_on_circulant(circ, row1, col1, rowLen, colLen, circulant_size); //embed to next power of 2
    pad_key(padded_key, input_key, circulant_size, rowLen);

    DIT_FFT(circ,circulant_size);
    DIT_FFT(padded_key,circulant_size);
    for(int i = 0; i < circulant_size; i++){
      circ[i] = circ[i] * padded_key[i];
    }

    DIT_IFFT(circ, circulant_size);


    for(int i = 0; i < colLen; i++){
      int bit_rounded = round(creal(circ[i]));
      output_key[i] = bit_rounded%2; 
    }
    free(circ);
    free(padded_key);
}

void decimalToBinary(int num, int *binaryNum, int len) {   
  if (num == 0) {
      return;
  }
 int i=len-1;
 
 for ( ;num > 0; ){
    binaryNum[i--] = num % 2;
    num /= 2;
 }
}

void binary_char_to_int(char *char_arr, int *int_arr, int len){
  for(int i = 0; i < len; i++){
    int_arr[i] = char_arr[i] - '0';
  }
}

/*
void generate_binary_list(int arr[], int length) {
  for (int i = 0; i < length; i++) {
      arr[i] = rand() % 2;  // Generates either 0 or 1
  }
}
*/
//////////////////////////
/// END FFT STUFF ////////

//Should not need num
void take_base_10_input(int *arr, int arrLen, int num){

    char arr_char[32];
    my_read(arr_char, 32);
    my_puts("\n");
    //int num = atoi(arr_char);
    //num = char_to_int(arr_char,5);
    decimalToBinary(num, arr, arrLen);
    print_int_array_serial(arr, arrLen);
}


void take_base_2_input(int *arr_int, int arrLen){
  int maxlen = 1024;
  char arr_char[maxlen];
  my_read(arr_char, maxlen);
  my_puts("\n");
  binary_char_to_int(arr_char, arr_int, arrLen);
  print_int_array_serial(arr_int, arrLen);

}

void Toeplitz_setup(int *row1, int *col1, int rowLen, int colLen){
    delay_2_ms();
    my_puts("Toeplitz Privacy Amplifification\n");
    my_puts("\n\n");

    //Give them one last warning
    my_puts("Define Toeplitz Matrix\n");

   // trigger_low();

    //Get password
    int row_int = 89435;
    my_puts("\nPlease enter row spec (Max 32 Bits): \n");
    take_base_2_input(row1, rowLen);


   // trigger_low();

    int col_int = 5547;
    my_puts("\nPlease enter col spec (Max 32 Hex Bits): \n");
    take_base_2_input(col1, colLen);
}



int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    //trigger_low();


    //Memory allocation Uncomment if want full PA working

    /*
    int rowLen = 5;
    int colLen = 4;
    int *row1 = (int*)malloc(rowLen * sizeof(int));
    int *col1 = (int*)malloc(colLen * sizeof(int));

    for(int i = 0; i < rowLen; i++){
      row1[i] = 0;
    }
    for(int i = 0; i < colLen; i++){
      col1[i] = 0;
    }

    Toeplitz_setup(row1, col1, rowLen, colLen);
   */

    int rowLen = 8;
    int colLen = 1;
    while(1){



      int *input_key = (int*)malloc(rowLen * sizeof(int));
      //int *output_key_fft = (int*)malloc(colLen * sizeof(int));
      for(int i = 0; i < rowLen; i++){
        input_key[i] = 0;
      }
      //for(int i = 0; i < colLen; i++){
      //  output_key_fft[i] = 0;
      //}

      my_puts("\nPlease enter input key (Max 32 Hex Bits): \n");
      take_base_2_input(input_key, rowLen);

      //Toeplitz_hash_fft_butterfly(input_key, output_key_fft, row1, col1, rowLen, colLen);


      ///// JUST FFT ON SECRET KEY ///////////////



      int circulant_size = rowLen + colLen - 1;
      float complex *padded_key = (float complex*)malloc(circulant_size * sizeof(float complex));

      for(int i=0; i < circulant_size; i++){
        padded_key[i] = 0 + 0*I;
        padded_key[i] += input_key[i]; 
      }
      //pad_key(padded_key, input_key, circulant_size, rowLen);
      
      //trigger_high();
      DIT_FFT(padded_key,circulant_size);


      /////////////////////////////////////////////


      //my_puts("\nCalculated P.A. key: \n");
      //print_int_array_serial(output_key_fft, colLen);  

      free(input_key);
      //free(output_key_fft);
      //trigger_low();

      }
      return 1;
  }


