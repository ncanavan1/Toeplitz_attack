#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>



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

int main(void)
{

    int pa_len = 8; //place holder sizes
    int l_key = 10;

    char row_spec[l_key];
    char col_spec[pa_len];
    char ec_key[l_key];
    char pa_key[pa_len];

    int row_spec_i[l_key];
    int col_spec_i[pa_len];
    int ec_key_i[l_key];
    int pa_key_i[pa_len];

    for(int i = 0; i < pa_len; i ++){
      pa_key[i] = '0';
      pa_key_i[i] = 0;
    }

    while(1){

        printf("\n\n****Provide Toeplitz Matrix Specification*****\n");
        printf("\n\n");


        //Get Row spec
        printf("Please enter row spec (size l_key): ");
        scanf("%s",row_spec);


        printf(row_spec);
        char_to_int_arr(row_spec, row_spec_i, l_key);
        //print_int_array_serial(row_spec_i,l_key);
        printf("\n");



        //Get coloumn spec
        printf("Please enter col spec (size pa_len): ");
        scanf("%s", col_spec);


        printf(col_spec);
        printf("\n");


        //Get EC key
        printf("Please enter Error Corrected Key (size l_key): ");
        scanf("%s", ec_key);


        printf(ec_key);
        printf("\n");


        char_to_int_arr(col_spec, col_spec_i, pa_len);
        char_to_int_arr(ec_key, ec_key_i, l_key);

        toeplitz_hash(row_spec_i, col_spec_i, ec_key_i, pa_key_i, pa_len, l_key);

        int_to_char_arr(pa_key_i, pa_key, pa_len);

        printf("PA KEY: ");
        printf(pa_key);


    }

    return 1;
}

