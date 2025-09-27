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


void take_base_2_input(int *arr_int, int arrLen){
  int maxlen = 1024;
  char arr_char[maxlen];
  my_read(arr_char, maxlen);
  my_puts("\n");
  binary_char_to_int(arr_char, arr_int, arrLen);
  print_int_array_serial(arr_int, arrLen);

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

void run_xor_test(){

  int repeats = 100;
  int xor_res = 0;
  int a = 0;
  int b = 0;

  for(int i = 0; i < repeats; i++){
    if(i%4 == 0){
      a = 0;
      b = 0;
    }
    else if (i%4 == 1){
      a = 0;
      b = 1;
    }
    else if (i%4 == 2){
      a = 1;
      b = 0;
    }
    else if (i%4 == 3){
      a = 1;
      b = 1;
    }

    trigger_high();
    xor_res = a ^ b;
    trigger_low();
    my_puts("\nDone");

  }

}


int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    int count = 0;
    int a = 0;
    int b = 0;

    my_puts("Start?");

    while(1){
      my_puts("\nPress any key to start\n");
      int maxlen = 1024;
      char arr_char[maxlen];
      my_read(arr_char, maxlen);  
      run_xor_test();  

    }
    return 1;
  }


