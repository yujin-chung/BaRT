#include "run_sampler.h"

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  
  time0 = clock ();
  
  if (argc != 2)
    biomc2_error ( " USAGE: %s <file>", basename (argv[0]));

 
  // Yujin Chung 5/19/2014
  // Added the below "printf" commend and removed the following commend.
  printf ("Recombination detection program based on segment-wise topology distance\n");
  // printf ("|\t %s\n| Recombination detection program based on segment-wise topology distance\n", ProgramVersion);

  /*---------------------------------------
    By Yujin Chung on 2/12/2013
    -- changed print-out comments
  ----------------------------------------- */
  printf ("|\t Yujin Chung and Cecile Ane\n|\n");
  printf ("|\t This program is modefied from biomc2 by Leonardo de Oliveira Martins and Hirohisa Kishino\n|\n");
  printf ("| This program is protected by the GPL license\n\n");
  
  /*---------------------------------------
     Added by Yujin Chung 8/5/2014
    ----------------------------------------- */
printf("Warning:: the software requires the joint distribution of a ranked tree and a distance given the other ranked tree.\n\n");

  random_number = new_biomc2_rng ();

  /*---------------------------------------
    By Yujin Chung on 2/12/2013
    start running the new program
  ----------------------------------------- */
  #ifdef DEBUG
  // REMOVE - YC 8/5/2014
 printf("REMOVE: Calling 'run_sampler_gibbs()'.\n\n");
 // END of REMOVE
#endif //DEBUG
 run_sampler_gibbs(argv[1]); // BaRT
  // run_sampler (argv[1]); // boimc2

  #ifdef DEBUG
  // REMOVE - YC 8/5/2014
 printf("REMOVE: Done with 'run_sampler_gibbs()'.\n\n");
 // END of REMOVE
#endif // DEBUG
  
  time1 = clock ();
  fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  
  del_biomc2_rng (random_number);
  return EXIT_SUCCESS;
}



/*----------------------------------------------
  By Yujin Chung on 2/12/2013
  The following function is from biomc2-1.9
------------------------------------------------ */
/*
int
main (int argc, char **argv)
{
  clock_t time0, time1;
  
  time0 = clock ();
  
  if (argc != 2)
    biomc2_error ( " USAGE: %s <file>", basename (argv[0]));
  
  printf ("|\t %s\n| Recombination detection program based on segment-wise topology distance\n", ProgramVersion);
  printf ("|\t Leonardo de Oliveira Martins and Hirohisa Kishino\n|\n");
  printf ("| This program is protected by the GPL license\n\n");
  
  random_number = new_biomc2_rng ();
  
  run_sampler (argv[1]);
  
  time1 = clock ();
  fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  
  del_biomc2_rng (random_number);
  return EXIT_SUCCESS;
}
*/
