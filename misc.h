#include <stdio.h> 
#include <sys/resource.h> 
#include <sys/time.h> 
#include <unistd.h> 

/////////////////////
//    misc.c 
/////////////////
/*
  input: nonnegative integer x
  output: x!! = x(x-2)...1            */
long double doublefactorial(int x);
long double log_doublefactorial(int x);
void print_cpu_time(); 
