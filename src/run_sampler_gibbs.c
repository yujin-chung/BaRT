/*!
 * Copyright (C) 2016 Yujin Chung
 *
 */

#include "run_sampler.h"

/* local function prototypes */

/*---------------------------------
  By Yujin Chung on 2/6/2013:
  New functions were created because of 
  new data type "chain_data_gibbs".
----------------------------------- */
void print_initial_info_gibbs(chain_data_gibbs chain, int interval);
void print_header (void);
void set_frequency_to_zero (acceptance_frequency *afreq);
void estimate_time (int seconds, double ntimes);
void print_acceptance_rates_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq);
void mcmc_setup_initial_values_gibbs (chain_data_gibbs chain);
void mcmc_set_temperature_gibbs (chain_data_gibbs chain1, chain_data_gibbs chain2, int *i, int *j);
void mcmc_forward_backward_gibbs (chain_data_gibbs chain, int *iter, acceptance_frequency *afreq, convergence C);
void mcmc_forward_gibbs  (chain_data_gibbs chain, acceptance_frequency *afreq);
void mcmc_backward_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq);

void clean_linked_topology_gibbs (chain_data_gibbs chain, int *iter);

void sample_to_output_file_gibbs (chain_data_gibbs chain_file, chain_data_gibbs chain, int *j);
void open_output_files_gibbs (chain_data_gibbs chain);
void close_output_files_gibbs (chain_data_gibbs chain);


/*-------------------------------
  By Yujin Chung on 2/6/2013:
  This function is to trigger an analysis containing
  -- call functions to read input files
  -- call functions to initiate MCMC parameters
  -- call functions to print out results

  It replaces "run_sampler".
---------------------------------- */
void run_sampler_gibbs (char *filename)
{
  #ifdef DEBUG
  // REMOVE - YC 8/5/2014
  printf("Calling 'run_sampler_gibbs()'.\n\n");
  // END of REMOVE
#endif // DEBUG

  int i, j, time_estim = 0, interval;
  //  int n_mini_bkp, n_cycles_bkp;
  //  double kT_mini_bkp, kT_cycles_bkp;
  clock_t time0, time1;
  acceptance_frequency cold_freq, heat_freq;
  chain_data_gibbs chain1, chain2;
  convergence C_cold, C_heat;
  
  /*-----------------------------------------
    By Yujin Chung on 2/28/2013
    read input files and create objects of data type "chain_data_gibbs" 
    containing all the parameters in the model and used in MCMC simulations 
  ----------------------------------------------*/
  #ifdef DEBUG
  printf("REMOVE: Reading input files - calling 'new_chain_data_from_file_gibbs()'.\n\n");
#endif// DEBUG
  chain1 = new_chain_data_from_file_gibbs (filename);
  chain2 = new_chain_data_from_file_gibbs (filename);

  #ifdef DEBUG
  printf("REMOVE: Done with new_chain_data_from_file_gibbs()\n\n");
  // printf("REMOVE: Set up settings.\n\n");
#endif// DEBUG
  setup_acceptance_names_gibbs ();
/*------------------------------------
  By Yujin Chung on 3/4/2013
  the following functions calculate the partition functions
  % fixit: double-check
---------------------------------------- */
  #ifdef DEBUG
  printf("REMOVE: calling mcmc_setup_initial_values_gibbs(chain1) .\n\n");
#endif //DEBUG
  mcmc_setup_initial_values_gibbs (chain1);
  mcmc_setup_initial_values_gibbs (chain2);
   #ifdef DEBUG
  printf("REMOVE: Done with mcmc_setup_initial_values_gibbs(chain1) .\n\n");
#endif //DEBUG
  // REMOVE - YC 8/5/2014
  //  printf("\t REMOVE: Done with setup - initial values.\n\n");
  // END of REMOVE
  time1 = clock ();
  
  interval = chain1->ngen/chain1->nsamples;
  chain2->nsamples = chain1->nsamples;
  
  print_initial_info_gibbs (chain1, interval);
  
  set_frequency_to_zero ( &cold_freq);
  set_frequency_to_zero ( &heat_freq);
  
  /* Posterior topologies */
  open_output_files_gibbs (chain1);
  
  /* stdout header */
  printf ("WARM-UP STAGE (simulated annealing)\n");
  
  // fixit: what does "i" stand for? -- YC 8/5/2014
  for (i = 1; i <= 5; i++)
  {
	  /*	 pre- burn-in to unlink topologies on heated chain */
	  for (j = 1; j <= chain1->warmup; j++)  // chain->warmup = 100; -- YC 8/5/2014
	  {
	      #ifdef DEBUG
      printf("cycle i = %d  iter j=%d\n",i,j);
#endif // DEBUG
	    if (!(j%chain1->ratio_spr)) // ratio_spr: iterations at which we update topologies by SPR instead of NNI -- YC 9/9/2016
		  {
			  branch_swap = &(spr_apply_random_move_spr);
			  chain1->logD = &(chain1->logDspr);
			  chain2->logD = &(chain2->logDspr);

			  // YC 9/9/2016
			  // logDspr: SPR-neighborhood size \f$\log(2(n-3)(2n-7))\f$
			  //          used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$.
			  // logD: pointer to chain_data::logDspr or chain_data::logDnni,
			  //       according to the operation
		  }
		  else
		  {
			  branch_swap = &(spr_apply_random_move_nni);
			  chain1->logD = &(chain1->logDnni);
			  chain2->logD = &(chain2->logDnni);

			  // YC 9/9/2016
			  // logDnni: NNI-neighborhood size \f$\log(2(n-3))\f$, 
			  //          used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$.
			  
		  }
      
		  mcmc_set_temperature_gibbs (chain1, chain2, &i, &j);
		  mcmc_forward_backward_gibbs (chain1, &j, &cold_freq, NULL);
		  mcmc_forward_backward_gibbs (chain2, &j, &heat_freq, NULL);
	  }
    
    time0 = time1; time1 = clock (); 
    
    printf ("\n[warm-up] iteration %5d, %8.4f secs, ", i*(j-1), (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    printf ("cooling scheme (1/temperature): %f - %f \n", chain1->beta_zero, chain1->beta_n);
    print_header ();
    printf ("ch 1: "); print_acceptance_rates_gibbs (chain1, &cold_freq);
    printf ("ch 2: "); print_acceptance_rates_gibbs (chain2, &heat_freq);
    set_frequency_to_zero ( &cold_freq);
    set_frequency_to_zero ( &heat_freq);
    
    time_estim += (int)(time1 - time0);
  } // End of Warm-up
  
  /*  
      chain1->n_mini    = chain2->n_mini    = n_mini_bkp;
      chain1->kT_mini   = chain2->kT_mini   = kT_mini_bkp;
      chain1->n_cycles  = chain2->n_cycles  = n_cycles_bkp;
      chain1->kT_cycles = chain2->kT_cycles = kT_cycles_bkp;
      chain1->acceptance_probability = &(acceptance_probability_bayesian);
      chain2->acceptance_probability = &(acceptance_probability_bayesian);
  */
  
  printf ("\nEstimated time to completion : "); 
  estimate_time ( time_estim, (double)(chain1->ngen + (5 * chain1->burnin))/(double)(5 * chain1->warmup) );
  printf ("\n\n");
  time_estim = 0;
  chain1->kT = chain2->kT = 1.;
  /* Setup SRQ initial states based on independent cold chains */
  printf ("BURN-IN STAGE (no swap between cold and heated chains): ");
  printf ("Sampling initial state for convergence statistics\n");
  for (i = 1; i <= 5; i++) {
    for (j = 1; j <= chain1->burnin; j++) {
      #ifdef DEBUG
      printf("cycle i = %d  iter j=%d\n",i,j);
#endif // DEBUG
      if (!(j%chain1->ratio_spr)) {
	branch_swap = &(spr_apply_random_move_spr);
	chain1->logD = &(chain1->logDspr);
	chain2->logD = &(chain2->logDspr);
      }
      else {
	branch_swap = &(spr_apply_random_move_nni);
	chain1->logD = &(chain1->logDnni);
	chain2->logD = &(chain2->logDnni);
      }
      
      mcmc_forward_backward_gibbs (chain1, &j, &cold_freq, NULL);
      mcmc_forward_backward_gibbs (chain2, &j, &heat_freq, NULL);
    }
    time0 = time1; time1 = clock (); 
    printf ("\n[burn-in] iteration %5d, %8.4f secs\n", i*(j-1), (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    print_header ();
    printf ("ch 1: "); print_acceptance_rates_gibbs (chain1, &cold_freq);
    printf ("ch 2: "); print_acceptance_rates_gibbs (chain2, &heat_freq);
    set_frequency_to_zero ( &cold_freq);
    set_frequency_to_zero ( &heat_freq);
    time_estim += (int)(time1 - time0);
  }
  C_cold = prepare_convergence_gibbs (chain1);
  C_heat = prepare_convergence_gibbs (chain2);
  
  printf ("\nEstimated time to completion : "); 
  //	estimate_time ( time_estim, (double)(chain1->ngen)/(double)(10 * chain1->warmup) );
  estimate_time ( time_estim, (double)(chain1->ngen)/(double)(5 * chain1->burnin) );
  printf ("\n\n");
  
  /* main sampler iteration */
  printf ("\nMAIN SAMPLER (saving prior/posterior parameters/topologies to file)\n");
  chain1->kT = 1.;              /* cold chain */
  chain2->kT = chain2->kT_heat; /* heated chain */
  
  for (j = 1; j <= chain1->ngen; j++) {
    
    #ifdef DEBUG
    printf("ngen =%d iter j= %d\n",chain1->ngen, j);
    for(int ii=0; ii<chain1->n_segments; ii++)
      printf("seg %d: dist %d  ",ii,chain1->segment[ii]->dist);
    printf("\n");
    int count_block = 1;
    for (topology tt=chain1->topol->first; (tt->next); tt = tt->next)
      { 
	printf ("block=%d first seg=%d last seg=%d tree=%s dist=%d\n",count_block, tt->first_segment, tt->last_segment,topology_to_string(tt), chain1->segment[tt->last_segment]->dist);
	count_block++;
      }
    topology tt=chain1->topol->last;
    printf ("block=%d first seg=%d last seg=%d tree=%s dist=%d\n",count_block, tt->first_segment, tt->last_segment,topology_to_string(tt), chain1->segment[tt->last_segment]->dist);
#endif //DEBUG
    
    if (!(j%chain1->ratio_spr)) {
      branch_swap = &(spr_apply_random_move_spr);
      chain1->logD = &(chain1->logDspr);
      chain2->logD = &(chain2->logDspr);
    }
    else {
      branch_swap = &(spr_apply_random_move_nni);
			chain1->logD = &(chain1->logDnni);
			chain2->logD = &(chain2->logDnni);
    }
    
    if (chain1->kT == 1.) {
      mcmc_forward_backward_gibbs (chain1, &j, &cold_freq, C_cold);
      mcmc_forward_backward_gibbs (chain2, &j, &heat_freq, C_heat);
    }
    else {
      mcmc_forward_backward_gibbs (chain2, &j, &cold_freq, C_cold);
      mcmc_forward_backward_gibbs (chain1, &j, &heat_freq, C_heat);
    }
    
    if (!(j%chain1->swap_interval)) propose_swap_chains_gibbs (chain1, chain2, &cold_freq, &heat_freq);
    
    if(!(j%interval)) {
      /* stream IO info is stored only in chain1 */
      if (chain1->kT == 1.) 
	sample_to_output_file_gibbs (chain1, chain1, &j);
      else 
	sample_to_output_file_gibbs (chain1, chain2, &j);
    }
    
    /* show accept_prob 20 times */
    if (!(j%(chain1->ngen/20))) {
      time0 = time1;
      time1 = clock ();
      
      printf ("\n[main sampler] iteration %d, remaining time : ", j);
      estimate_time ( (int)(time1 - time0), (double)((chain1->ngen - j) * 20)/(double)(chain1->ngen) );
      printf ("\n");
      print_header ();
      
      if (chain1->kT == 1.) {
	printf ("cold: "); print_acceptance_rates_gibbs (chain1, &cold_freq);
	printf ("heat: "); print_acceptance_rates_gibbs (chain2, &heat_freq);
      }
      else {
	printf ("cold: "); print_acceptance_rates_gibbs (chain2, &cold_freq);
	printf ("heat: "); print_acceptance_rates_gibbs (chain1, &heat_freq);
      }
      set_frequency_to_zero ( &cold_freq);
      set_frequency_to_zero ( &heat_freq);
    }
  }
  printf ("\n [cold chain] "); report_convergence (C_cold);
  printf ("\n [heat chain] "); report_convergence (C_heat);
  
  close_output_files_gibbs (chain1);
  del_chain_data_gibbs (chain1);
  del_chain_data_gibbs (chain2);
}


/*-----------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "print_initial_info"
-------------------------------------- */
void
print_initial_info_gibbs(chain_data_gibbs chain, int interval)
{
  printf ("  number of segments = %d\n", chain->n_segments);
  printf ("  ngen = %d burnin = %d(x10) sample_interval (after burnin) = %d warm-up = %d(x10) \n", chain->ngen, chain->burnin, interval, chain->warmup);
  printf ("  period, in number of iterations, between SPR moves (otherwise NNI): %d \n",  chain->ratio_spr);
  printf ("- in Al-Awadhi updates (changing the break-point structure):\n");
  printf ("            temperature = %f x normal, mini-sampler size = %d \n", chain->kT_cycles, chain->n_cycles);
  printf ("- in topology updates (without changing the break-point structure):\n");
  printf ("            temperature = %f x normal, mini-sampler size = %d \n", chain->kT_mini, chain->n_mini);
  printf ("            how many times, per iteration, to update topologies: %d \n", chain->ratio_topol_update);  
  printf ("  heated chain: heat factor 1/kT = %f, swap interval = %d\n\n", chain->kT_heat, chain->swap_interval);
  printf (" the acceptance rates below are (number of accepted)/(number of tries);\n");
  printf (" prop means (number of tries for this move)/(total number of tries for moves)\n\n");
}	


/*-----------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "print_acceptance_rates"
-------------------------------- */
void
print_acceptance_rates_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)
{
  int i;
  double total = 0., fraction = 0.;

  for (i=0; i<3;i++) 
    total += (double)(afreq->propose[i]);
  if (total < 1.) 
    total=1.;

  for (i=0; i < 3; i++) {	
    if (afreq->propose[i]) 
      fraction = (double)(afreq->accept[i])/(double)(afreq->propose[i]);
    else 
      fraction = 0.;
    printf ("%.5f %.2f|", fraction, (double)(afreq->propose[i])/total);
  } 
  printf (" ");
  for (i=3; i < AFREQSIZE; i++) {	
    if (afreq->propose[i]) 
      fraction = (double)(afreq->accept[i])/(double)(afreq->propose[i]);
    else 
      fraction = 0.;
    printf ("%.5f ", fraction);
  }
  printf (" [%d %d]\n", chain->topol->nSPR, chain->topol->nCOP);
}

/*---------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "mcmc_setup_initial_values"
------------------------------ */
void mcmc_setup_initial_values_gibbs (chain_data_gibbs chain)
{
  // int i;
 
  // for (i=0; i < chain->n_segments - 1; i++) {
    /*----------------------------------
      By Yujin Chung on 2/6/2013:
      penalty(weight) w 
    ------------------------------------- */
    // w = chain->segment[i]->penalty + WEIGHT;

    /*-------------------------------
      By Yujin Chung on 2/6/2013:
      The partition function is replaced.
By Yujin Chung on 3/4/2013
  In BaRT, lambda is the same across the segments and it is a member of chain_data_gibbs
    ------------------------------------- */
    partitionFunction_gibbs (&(chain->lnN), &(chain->lambda), chain);   
    // lnTruncate_gibbs (&(chain->lnN[i]), &(chain->segment[i]->lambda), chain);
    // *** YC 5/24/2017 ***
    // lnTruncate() computes the partition (or normalizing) function of truncated Poisson prior distribution on gene trees
    // }  
} 


/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "mcmc_set_temperature".
---------------------------------------- */
void
mcmc_set_temperature_gibbs (chain_data_gibbs chain1, chain_data_gibbs chain2, int *i, int *j)
{ 
/*
	chain1->beta_n = chain1->beta_zero + (log ( (double)((((*i)-1)*chain1->warmup) + (*j)) ))/4.;
	chain2->beta_n = exp (chain1->beta_n);
	chain1->beta_n = chain2->beta_n;
 */
  chain1->beta_n = chain1->beta_zero * log ( (double)( (*j)) + EXP_1 );
  chain1->kT = chain2->kT = chain1->beta_n;
}


/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "mcmc_forward_backward".
---------------------------------------- */
void mcmc_forward_backward_gibbs (chain_data_gibbs chain, int *iter, acceptance_frequency *afreq, convergence C)
{
  #ifdef DEBUG
  printf("\n***** Calling mcmc_foward_backward_gibbs() *****\n");
  printf("This mcmc_forward_backward_gibbs() is the key function for mcmc updates for all parameters.\n");

  printf("iter = %d\n",(*iter));
  for(int ii=0; ii<chain->n_segments; ii++)
    printf("seg %d: dist %d  ",ii,chain->segment[ii]->dist);
  printf("\n");
  int count_block = 1;
  for (topology tt=chain->topol->first; (tt->next); tt = tt->next)
    { 
      printf ("block=%d first seg=%d last seg=%d tree=%s dist=%d\n",count_block, tt->first_segment, tt->last_segment,topology_to_string(tt), chain->segment[tt->last_segment]->dist);
      count_block++;
    }
  topology tt=chain->topol->last;
  printf ("block=%d first seg=%d last seg=%d tree=%s dist=%d\n\n",count_block, tt->first_segment, tt->last_segment,topology_to_string(tt), chain->segment[tt->last_segment]->dist);
  /*
  printf("n_segments = %d\n", chain->n_segments);
  printf("'breakpoint' is ");
  for(int y=0; y<chain->n_segments; y++)
    printf("%d ",chain->breakpoint[y]);
  printf("\n");
  */
#endif // DEBUG
  
  int i;
  topology t;

  #ifdef DEBUG
  printf("[5/31/2016 YC] This function mcmc_forward_gibbs() is for updating topologies or segments.\n");
#endif // DEBUG
  mcmc_forward_gibbs (chain, afreq);
  clean_linked_topology_gibbs (chain, iter);
  
 #ifdef DEBUG
  printf("Another update: mcmc_topol_update_outside_gibbs().\n");
#endif // DEBUG
  for (i = 0; i < chain->ratio_topol_update; i++)
    {
      // YC 9/9/2016
      // ratio_topol_update: the number of topology updates per iteration,
      //                     to reduce correlation between samples.
      
      for (t=chain->topol->first; (t); t = t->next)
	mcmc_topol_update_outside_gibbs (chain, t, afreq);
      clean_linked_topology_gibbs (chain, iter);
    }
  if (C) 
    update_convergence_gibbs (chain, C);
  
 #ifdef DEBUG
  printf("Another update: mcmc_backward_gibbs().\n");
#endif // DEBUG
  mcmc_backward_gibbs (chain, afreq);
  clean_linked_topology_gibbs (chain, iter);
  
 #ifdef DEBUG
  printf("Again: mcmc_topol_update_outside_gibbs().\n");
#endif // DEBUG
  for (i = 0; i < chain->ratio_topol_update; i++)
    {
    for (t=chain->topol->last; (t); t = t->prev)
      mcmc_topol_update_outside_gibbs (chain, t, afreq);
    clean_linked_topology_gibbs (chain, iter);
  }
  if (C) 
    update_convergence_gibbs (chain, C);
  
 
  if (!((*iter)%5)) {
      #ifdef DEBUG
  printf("substitution parameter updates.\n");
#endif // DEBUG
    for (t=chain->topol->first; (t); t = t->next) 
      for (i=t->first_segment; i <= t->last_segment; i++)
	{
	  proposal_update_substitution_rate_gibbs  (chain, i, t, afreq);
	  proposal_update_substitution_kappa_gibbs (chain, i, t, afreq);
	}
  }
  
  /*------------------------------------
By Yujin Chung on 2/6/2013:
      Weight w and related parameters were removed.
  By Yujin Chung on 3/6/2013
  lambda is the same across the segments in BaRT.
  ---------------------------------------- */
#ifdef DEBUG
  printf("Calling proposal_update_lambda_gibbs().\n");
#endif // DEBUG
  proposal_update_lambda_gibbs  (chain, afreq);
  /*
  for (i=0; i < chain->n_segments - 1; i++) {
    proposal_update_lambda_gibbs  (chain, i, afreq);
    proposal_update_penalty (chain, i, afreq);
  }
  */

  if (!((*iter)%2))
    {
      #ifdef DEBUG
  printf("Updating hyperparameters.\n");
#endif // DEBUG
    proposal_update_prior_rate_gibbs  (chain, afreq);
    proposal_update_prior_kappa_gibbs (chain, afreq);
  }

  
  #ifdef DEBUG
  printf("\n***** Exiting mcmc_foward_backward_gibbs() *****\n\n");
#endif // DEBUG
  return;
}


/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "mcmc_forward".
---------------------------------------- */
void mcmc_forward_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)
{
  #ifdef DEBUG
  printf("Calling mcmc_forward_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)\n");
#endif //DEBUG
  topology t = chain->topol->first;
  double u, log2 = chain->log2;

  if (t->last_segment > t->first_segment) {
    if (!(t->next)) log2 = 0.;   
    proposal_update_topol_birth_forward_gibbs (t, chain, afreq, &log2);  
  }
  t = t->next;
  
  while (t) {
    log2 = 0.;
    
    if (t->last_segment > t->first_segment) {
      u = biomc2_rng_uniform (random_number);
      if (u < 0.4) {
	//				if (!(t->next)) log2 = chain->log2;
	proposal_update_topol_birth_forward_gibbs (t, chain, afreq, &log2);
      }
      else if (u < 0.8) {
	if (!(t->next)) log2 = -chain->log2;
	proposal_update_topol_death_forward_gibbs (t, chain, afreq, &log2);
      }
      else {
	proposal_update_topol_shift_forward_gibbs (t, chain, afreq);
      }
    }
    
    else {
      if (!(t->next)) log2 = -chain->log2;
      proposal_update_topol_death_forward_gibbs (t, chain, afreq, &log2);
    }
    
    if ( (t->next) && (!(chain->segment[t->last_segment]->dist)) ) { t = t->next; }
    t = t->next;
  }
}


/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "mcmc_backward".
---------------------------------------- */
void
mcmc_backward_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)
{
  topology t = chain->topol->last;
  double u, log2 = chain->log2;
  
  if (t->last_segment > t->first_segment) {
    if (!(t->prev)) log2 = 0.;
    proposal_update_topol_birth_backward_gibbs (t, chain, afreq, &log2);
  }
  t = t->prev;
  
  while (t) {
    log2 = 0.;
    
    if (t->last_segment > t->first_segment) {
      u = biomc2_rng_uniform (random_number);
      if (u < 0.4) {
	//				if (!(t->prev)) log2 = chain->log2;
	proposal_update_topol_birth_backward_gibbs (t, chain, afreq, &log2);
      }
      else if (u < 0.8) {
	if (!(t->prev)) log2 = -chain->log2;
	proposal_update_topol_death_backward_gibbs (t, chain, afreq, &log2);
      }
      else {
	proposal_update_topol_shift_backward_gibbs (t, chain, afreq);
      }
    }
    
    else {
      if (!(t->prev)) log2 = -chain->log2;
      proposal_update_topol_death_backward_gibbs (t, chain, afreq, &log2);
    }
    
    if ( (t->prev) && (!(chain->segment[t->prev->last_segment]->dist)) ) { t = t->prev; }
    t = t->prev;
  }
}

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "clean_linked_topology".
---------------------------------------- */
void
clean_linked_topology_gibbs (chain_data_gibbs chain, int *iter)
{
  int r_lt_l = 0;
  topology t, redundant, keep;
  
  chain->topol->nSPR = chain->topol->nCOP = 0;
  for (t = chain->topol->first; (t->next); ) {
    if (chain->segment[t->last_segment]->dist) {
      chain->topol->nCOP++; chain->topol->nSPR += chain->segment[t->last_segment]->dist;
      t = t->next;
    }
    else { /* the two segments will be merged (since dist = 0) */
      /* r_lt_l = right topol larger (more segments) than left topol */
      r_lt_l = t->next->last_segment - t->last_segment + t->first_segment - t->next->first_segment;
      if (r_lt_l > 0) 
	{ redundant = t; keep = t->next; }
      else 
	{ redundant = t->next; keep = t; }
      
      /* swap likelihood nodes */
      copy_topology_from_topology_mapping (redundant, keep, chain->split);
      swap_likelihood_gibbs (redundant, chain->segment, keep);
/*
  redundant->likelihood_proposal = redundant->likelihood_accepted;
  if (debug_topol (redundant, chain)) { 
  printf ("redundant\n");
  }
  redundant->likelihood_accepted = redundant->likelihood_proposal;
*/
      keep->likelihood_accepted += redundant->likelihood_accepted;
      if (r_lt_l > 0) 
	keep->first_segment = redundant->first_segment;
      else 
	keep->last_segment = redundant->last_segment;
			
      remove_topology (redundant, chain->topol);
      t = keep;
/*		
		keep->likelihood_proposal = keep->likelihood_accepted;
		if (debug_topol (keep, chain)) {
		printf ("keep\n");
		}
		keep->likelihood_accepted = keep->likelihood_proposal;
 */
    }
  }
}
  

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "sample_to_output_file".
---------------------------------------- */
void
sample_to_output_file_gibbs (chain_data_gibbs chain_file, chain_data_gibbs chain, int *j)
{
  int i=0;
  char *s;
  double lnLik = 0.;
  topology t;
  
  for (t=chain->topol->first; (t); t = t->next) {
    create_topology_bipartition_l (t);
    order_topology_by_bipartition (t->root->left);
    s = topology_to_string (t);
    fprintf (chain_file->treefile, "tree rep%ds%d = %s\n",(*j), i++, s);
    fflush (chain_file->treefile);
    free (s);
  }
  
  fprintf (chain_file->distfile, "%-8d ", (*j));
  
  full_likelihood_gibbs (&lnLik, chain);
  fprintf (chain_file->distfile, "%12f ", lnLik);
  fprintf (chain_file->distfile, "%4d ", chain->topol->nSPR);
  fprintf (chain_file->distfile, "%4d ", chain->topol->nCOP);
  fprintf (chain_file->distfile, "%12f ", chain->prior_rate);
  fprintf (chain_file->distfile, "%12f ", chain->prior_kappa);

  for (t=chain->topol->first; (t->next); t = t->next) { 
    fprintf (chain_file->distfile, "| %4d %2d ", t->last_segment, chain->segment[t->last_segment]->dist);
  }
  fprintf (chain_file->distfile, "\n");
  fflush (chain_file->distfile);
  
}
     

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "open_output_files".
---------------------------------------- */
void
open_output_files_gibbs (chain_data_gibbs chain)
{
  int i;
  
  /* posterior parameters */
  if (chain->prior) {
    chain->distfile = biomc2_fopen ("prior.dist", "w");
    chain->treefile = biomc2_fopen ("prior.tree", "w");
  } else {
    chain->distfile = biomc2_fopen ("post.dist", "w");
    chain->treefile = biomc2_fopen ("post.tree", "w");
  }
  
  fprintf (chain->treefile, "#NEXUS\n\nBegin trees;\n Translate\n");
  fprintf (chain->treefile, "\t1  %s", chain->taxlabel[0]);
  for (i=1; i < chain->ntax; i++) 
    fprintf (chain->treefile, ",\n\t%d  %s", i+1, chain->taxlabel[i]);
  fprintf (chain->treefile, "\n;\n");
  
  
  /* file Header */ 
  fprintf (chain->distfile, "Nseg %d Nsamples %d \n",chain->n_segments, chain->nsamples);
  for (i=0; i<chain->n_segments-1; i++) fprintf (chain->distfile, "%d ", chain->breakpoint[i]);
  fprintf (chain->distfile, "%d\n", chain->nsites);
  // YC 5/31/2017
  fprintf(chain->distfile,"iter\tloglik\tnSPR\tnCOP\tprior_rate\tprior_kappa |last_segment  dist | ... \n");
}

  

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "close_output_files".
---------------------------------------- */
void
close_output_files_gibbs (chain_data_gibbs chain)
{
  fprintf (chain->treefile, "\nEnd;\n");
  fclose (chain->treefile); 
  fclose (chain->distfile);
}


