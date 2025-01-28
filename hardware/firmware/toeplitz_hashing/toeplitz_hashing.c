#include "hal.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#define IDLE 0
#define KEY 1
#define PLAIN 2

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


void char_to_int_arr(char *c_arr, int *i_arr, int len){
  for(int i = 0; i < len; i++){
    i_arr[i] = c_arr[i] - '0';
  }
}

void int_to_char_arr(int *i_arr, char *c_arr, int len){
  for(int i = 0; i < len; i++){
    c_arr[i] = i_arr[i] + '0';
  }
}


void toeplitz_hash(int *row_spec, int *col_spec, int *ec_key, int *pa_key, int pa_len, int l_key){
  for(int j = 0; j < pa_len; j++){
    for(int i = 0; i < l_key; i++){
      int intermediate_val = 0;
      if(i >= j){
        intermediate_val = row_spec[i-j]*ec_key[i];
        pa_key[j] = (pa_key[j] + intermediate_val)%2;
      }
      else{
        intermediate_val = col_spec[j-i]*ec_key[i];
        pa_key[j] = (pa_key[j] + intermediate_val)%2;
      }
    }
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

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    int pa_len = 8; //place holder sizes
    int l_key = 10;

    char row_spec[l_key+1];
    char col_spec[pa_len+1];
    char ec_key[l_key+1];
    char pa_key[pa_len+1];

    int row_spec_i[l_key];
    int col_spec_i[pa_len];
    int ec_key_i[l_key];
    int pa_key_i[pa_len];

    for(int i = 0; i < pa_len; i ++){
      pa_key[i] = '0';
      pa_key_i[i] = 0;
    }

    while(1){

        my_puts("\n\n****Provide Toeplitz Matrix Specification*****\n");
        delay_2_ms();
        my_puts("\n\n");

        trigger_low();

        //Get Row spec
        my_puts("Please enter row spec (size l_key): ");
        my_read(row_spec, l_key+1);

        trigger_high();

        my_puts(row_spec);
        char_to_int_arr(row_spec, row_spec_i, l_key);
        print_int_array_serial(row_spec_i,l_key);
        my_puts("\n");


        trigger_low();

        //Get coloumn spec
        my_puts("Please enter col spec (size pa_len): ");
        my_read(col_spec, pa_len+1);

        trigger_high();

        my_puts(col_spec);
        my_puts("\n");

        trigger_low();

        //Get EC key
        my_puts("Please enter Error Corrected Key (size l_key): ");
        my_read(ec_key, l_key+1);

        trigger_high();

        my_puts(ec_key);
        my_puts("\n");


        char_to_int_arr(col_spec, col_spec_i, pa_len);
        char_to_int_arr(ec_key, ec_key_i, l_key);

        toeplitz_hash(row_spec_i, col_spec_i, ec_key_i, pa_key_i, pa_len, l_key);

        int_to_char_arr(pa_key_i, pa_key, pa_len);

        //my_puts("PA KEY: ");
        //my_puts(pa_key);

        int_to_char_arr(row_spec_i, row_spec, l_key);
        my_puts(row_spec);

    }

    return 1;
}

