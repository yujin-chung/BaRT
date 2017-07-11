#include <stdio.h> 
#include <sys/resource.h> 
#include <sys/time.h> 
#include <unistd.h> 

/*
  input: nonnegative integer x
  output: x!! = x(x-2)...1            */
long double doublefactorial(int x)
{
  int i;
  long double fac=1.0;
  
  if(x==0)
    return 0;
  else
    {
      for(i=x;i>=1;i=i-2)
	{
	  fac*=(long double) i;
	}
      //printf("doublefactorial: res=%Le\n",fac);
      return fac;
    }
}

/*
  input: nonnegative integer x
  output: log(x!!) = x+(x-2)+...+1            */
long double log_doublefactorial(int x)
{
  int i;
  long double fac=0.0;
  
  if(x==0)
    return 0;
  else
    {
      for(i=x;i>=1;i=i-2)
	{
	  fac +=(long double) i;
	}
      //printf("doublefactorial: res=%Le\n",fac);
      return fac;
    }
}


void print_cpu_time() 
{

  struct rusage usage, cusage; 
  int current_utime_tv_sec, current_utime_tv_usec, current_stime_tv_sec, current_stime_tv_usec;
  getrusage (RUSAGE_SELF, &usage); 

  current_utime_tv_sec = usage.ru_utime.tv_sec;
  current_utime_tv_usec = usage.ru_utime.tv_usec;
  current_stime_tv_sec = usage.ru_stime.tv_sec;
  current_stime_tv_usec = usage.ru_stime.tv_usec;
  printf ("cumulative CPU time: %ld.%06ld sec user, %ld.%06ld sec system\n", 
          usage.ru_utime.tv_sec, usage.ru_utime.tv_usec, 
          usage.ru_stime.tv_sec, usage.ru_stime.tv_usec); 
  printf ("CPU time: %ld.%06ld sec user, %ld.%06ld sec system\n", 
          usage.ru_utime.tv_sec-current_utime_tv_sec, usage.ru_utime.tv_usec-current_utime_tv_usec, 
          usage.ru_stime.tv_sec-current_stime_tv_sec, usage.ru_stime.tv_usec-current_stime_tv_usec); 
  printf("max resident set size (physical memory allocated in kilobytes): %ld\n",usage.ru_maxrss);
  printf("memory size for code: %ld\n",usage.ru_ixrss);
  printf("memory size for statistics, globals, and new/malloc: %ld\n",usage.ru_isrss);

  /*

  getrusage (RUSAGE_CHILDREN, &cusage); 

  printf ("CPU time: %ld.%06ld sec user, %ld.%06ld sec system\n", 
         cusage.ru_utime.tv_sec, cusage.ru_utime.tv_usec, 
          cusage.ru_stime.tv_sec, cusage.ru_stime.tv_usec); 
  */
} 

 
