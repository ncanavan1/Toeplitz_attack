#include "hal.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>

void my_puts(char *c)
{
  do {
    putch(*c);

  } while (*++c);
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


static volatile float volatile_sink;

float take_base_10_input(){

  char arr_char[32];
  my_read(arr_char, 32);
  return atof(arr_char);
}

void print_float_hex(float f) {
  union {
      float f;
      uint32_t u;
  } u;
  u.f = f;

  char buf[16];
  sprintf(buf, "0x%08lX", (unsigned long)u.u);
  my_puts(buf);
}

float randomFloat()
{
    return (float)(rand()) / (float)(RAND_MAX);
}

float xor_floats(float a, float b) {
  uint32_t ua, ub, uc;

  // reinterpret float as 32-bit integer
  memcpy(&ua, &a, sizeof(float));
  memcpy(&ub, &b, sizeof(float));

  // XOR the integer representations
  uc = ua ^ ub;

  // reinterpret result back to float
  float c;
  memcpy(&c, &uc, sizeof(float));
  return c;
}

int main(void) {
  platform_init();
  init_uart();
  trigger_setup();

  /* Simple prompt (optional) */
  my_puts("FP_MUL_READY\n");

  char line[128];
  while (1) {

      my_puts("\nGive two floats: \n");

      float a = take_base_10_input();
      float b = take_base_10_input();

      my_puts("Input a: ");
      print_float_hex(a);
      my_puts("\nInput b: ");
      print_float_hex(b);

      float mask_a = randomFloat();
      float mask_b = randomFloat();

      float a_mask = xor_floats(a, mask_a);
      float b_mask = xor_floats(b, mask_b);


      /* trigger -> multiply -> untrigger */

      trigger_high();
      // for(int i = 0; i < 10; i++){
      //     __asm__ volatile ("nop");
      // }

      /* The multiply itself */

      /////////////////////////
      /////////////////////////
      //UNMASKED VERSION
      
      //volatile float c = a * b;

      /////////////////////////
      /////////////////////////
      //MASKED VERSION
      volatile float c_mask = a_mask * b_mask;
      volatile float mask_c = mask_a * mask_b;

      volatile float c = xor_floats(c_mask, mask_c);

      /*
        * Small software fence to ensure result is materialised and not optimized away.
        * (volatile_sink forces the compiler to keep the result)
        */
      volatile_sink = c_mask;
      volatile_sink = mask_c;
      volatile_sink = c;

      // for(int i = 0; i < 10; i++){
      //   __asm__ volatile ("nop");
      // }

      trigger_low();

      /* Optionally perform a tiny delay if you want measurable time window:
        * for (volatile int i=0;i<100;i++);  // uncomment if needed
        */

      /* Print back the result */
      my_puts("\nOutout c: ");
      print_float_hex(c);

  }

    return 0;
}


