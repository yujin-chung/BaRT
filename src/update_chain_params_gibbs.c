#include "update_chain.h"
#include "shape.h"

/* local function prototypes */

void update_variable_logscale (double *result, double *current, double *step_size);
void update_variable (double *result, double *current, double *step_size);

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "full_likelihood"
---------------------------------------- */
void
full_likelihood_gibbs (double *r, chain_data_gibbs chain)
{
  int i;
  double d1 = 0.;
  double sumD = 0.; // the sum of RF distances by Yujin Chung on 3/6/2013
  topology t;
  
  /*------------------------------------
    By Yujin Chung on 2/21/2013
    related to substitution rate
    ---------------------------------------- */
  (*r) = ( (chain->prior_rate/chain->hp_rate) * log (chain->hp_rate) );

  for (t=chain->topol->first; (t); t = t->next) 
    (*r) += t->likelihood_accepted;
  for (i=0; i < chain->n_segments; i++) 
    d1 += chain->segment[i]->rate;

  (*r) += ( (d1/chain->prior_rate) * log (chain->prior_rate) );

  /*-----------------------------------
    By Yujin Chung on 2/6/2013:
    Through this loop, prior distribution is calculated.
    Penalty w and related parameters were removed
      By Yujin Chung on 3/4/2013
      In biomc2, lambda could be different across the segments
      In BaRT, lambda is the same across the segments
  --------------------------------------- */
  // By Yujin Chung on 3/6/2013
  // computing the prior distribution on gene trees
  for (i=0; i < chain->n_segments - 1; i++) //{
    sumD += chain->segment[i]->dist;
  
  (*r) -= chain->lambda * sumD;
  (*r)-= chain->lnN;
  
  // By Yujin Chung on 3/6/2013
  // computing the hyper prior distribution on lambda (beta~(1/recombination rate) in Chung et al.)
  (*r) += (chain->hp_lambda_alpha - 1.) * log(chain->lambda);
  (*r) -= chain->lambda * chain->hp_lambda_beta;

    /*--------- 
       d1 = chain->segment[i]->penalty + WEIGHT;
    (*r) -= ( chain->lnN[i] + d1 * chain->factorial[ chain->segment[i]->dist ] ) ;
    (*r) += ( log (chain->segment[i]->lambda) * (chain->hp_lambda_alpha - 1. + (chain->segment[i]->dist * d1)) );
    (*r) += ( (chain->hp_penalty_alpha - 1.) * log (chain->segment[i]->penalty) );
    (*r) -= ( chain->segment[i]->lambda * (d1 + chain->hp_lambda_beta) + 
	      chain->hp_penalty_beta * chain->segment[i]->penalty );
    ------------------*/
    // }
}


/*------------------------------------
  By Yujin Chung on 2/22/2013
  This function calculates the partition(normalizing) function 
  of the Gibbs distribution for gene trees.
  It is on log scale.
---------------------------------------- */
void partitionFunction_gibbs (double *r, double *lambda, chain_data_gibbs chain)
{
  #ifdef DEBUG
	// REMOVE - YC 8/5/2014
  //	printf("\t In partitionFunction_gibbs.\n\n");
	// END of REMOVE
#endif //DEBUG
  

  int i,j,k; // index used in "for" loop

  int ntaxa = chain->ntax;
   #ifdef DEBUG
	// REMOVE - YC 8/5/2014
  // printf("ntaxa = %d\n",ntaxa);
	// END of REMOVE
#endif

  /*------- computing the number of tree shapes------------*/
  int nshape = no_urshapes(ntaxa);  
   #ifdef DEBUG
	// REMOVE - YC 8/5/2014
  // printf("nshape = %d\n",nshape);
	// END of REMOVE
#endif
  
 /*-------- calculating the full joint distribution ------------*/    
  long double***  propDistance_2trees_full;
  // if(ntaxa>15)
    {
      propDistance_2trees_full =  malloc(nshape * sizeof(long double**));
      for(i=0;i<nshape;i++)
	{
	  propDistance_2trees_full[i] = malloc(nshape * sizeof(long double*));
	  for(j=0;j<nshape;j++)
	    {	     
	      propDistance_2trees_full[i][j] = malloc((ntaxa-2) * sizeof(long double));  
	    }
	}

      #ifdef DEBUG
        printf("Computing the full joint distribution of RF distance and a tree: calling print_fullJointDistribution(ntaxa,nshape,propDistance_2trees_full);\n");
#endif // DEBUG
      print_fullJointDistribution(ntaxa,nshape,propDistance_2trees_full);
      #ifdef DEBUG
       printf("Done with print_fullJointDistribution(ntaxa,nshape,propDistance_2trees_full);.\n");
#endif //DEBUG
    }

  /*---------- Read joint distribution from a text file --------------*/
  // printf("Reading jointDistribution....\n");
  FILE *fopen(), *read_jointDist;
  char filename[10000];
  // sprintf(filename,"JointDistributions/JointDistribution_%dtaxa.txt",ntaxa);
  sprintf(filename,"JointDistribution_%dtaxa.txt",ntaxa);
  read_jointDist = fopen(filename,"r"); 
  long double***  propDistance_2trees;
  propDistance_2trees =  malloc(nshape * sizeof(long double**));
  for(i=0;i<nshape;i++)
    {
      //printf("i=%d\t",i);
      propDistance_2trees[i] = malloc(nshape * sizeof(long double*));
      for(j=0;j<nshape;j++)
	{	    
	  propDistance_2trees[i][j] = malloc((ntaxa-2) * sizeof(long double));  
	  for(k=ntaxa-3;k>=0;k--)
	    {
	      fscanf(read_jointDist, "%Le", &propDistance_2trees[i][j][k]);
	    }
	}
      //printf("i=%d j=%d k=%d\n",i,j,k);
    }
  fclose(read_jointDist);
  // printf("Done.\n");
  

  /*---------- calculating the independence approximation 
    to the normalizing constant (log scale) --------------*/
  (*r) = 0.;
  (*r) = log_normalizingConstant_indep(propDistance_2trees,ntaxa, chain->n_segments, nshape, (*lambda));
}


/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "lnTruncate".
  This function calculates the partition function 
  of truncated Poisson prior distribution on gene trees
  (originally used in biomc2)

  List of changes:
  -- It has three arguments
  -- penalty w and related parameters were removed.
---------------------------------------- */
/*
void 
lnTruncate_gibbs (double *r, double *l, chain_data_gibbs chain)
{
  /*------------------------------------
    By Yujin Chung on 2/22/2013
    l: lambda, the parameter of the truncated Poisson distribution
  ---------------------------------------- */
/*
  int j;
  double lg, lnl = log (*l);
  
  (*r) = 0.;
  for (j=0; j <= chain->max_distance; j++) {
    /*------------------------------------
      By Yujin Chung on 2/22/2013
      chain->factorial[k] has the value of log(k!)      
      ---------------------------------------- */
/*
    lg = (double)(j) * lnl - (*l) -  chain->factorial[j];
    /*------------------------------------
      By Yujin Chung on 2/6/2013:
      penalty w and related parameters were removed.
    ---------------------------------------- */
/*
    (*r) += exp (lg);
    // (*r) += exp ((*w) * lg);
  }
  (*r) = log ((*r));
}
*/









/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "calculate_moments".

  List of changes:
  -- weight w and related parameters were removed

  By Yujin Chung on 3/6/2013:
  This function is not used in BaRT and biomc2.
  So, we stop modifying it for BaRT.
---------------------------------------- */
void calculate_moments_gibbs (double *r1, double *r2, int i, chain_data_gibbs chain)
{
  // By Yujin Chung on 3/6/2013:
  // Warning message is added in case this function is called.
  printf("Warning!! calculate_moments_gibbs() is called and it may need to be modified.\n");

  //  int j;
  // double lg;
/*------------------------------------
  By Yujin Chung on 3/4/2013
  In BaRT, lambda is the same across the segments
  ---------------------------------------- */
  // double lnl = log (chain->lambda); 
  //  double lnl = log (chain->segment[i]->lambda); 

  /*------------------------------------
  By Yujin Chung on 2/7/2013
  weight w and related parameters were removed
  ---------------------------------------- */
  // double w = chain->segment[i]->penalty + WEIGHT;
  
  //(*r1) = (*r2) = 0.;
  /* first moment about zero (E[X]) */
  //for (j=1; j <= chain->max_distance; j++) {
    /*------------------------------------
  By Yujin Chung on 3/4/2013
  In BaRT, lambda is the same across the segments
  ---------------------------------------- */
  //lg = (double)(j) * lnl - (chain->lambda) - chain->factorial[j];
    // lg = (double)(j) * lnl - (chain->segment[i]->lambda) - chain->factorial[j];

    /*------------------------------------
      By Yujin Chung on 2/7/2013
      weight w and related parameters were removed
      ---------------------------------------- */
    //lg = exp (lg);
    // lg = exp (w * lg);
    //(*r1) += ( (double)(j)     * lg ); 
  //}
  
  /*------------------------------------
  By Yujin Chung on 3/4/2013
  lnN contains the partition function of the full Gibbs prior distribution on gene trees
  ---------------------------------------- */
  //(*r1) /= exp (chain->lnN);
  // (*r1) /= exp (chain->lnN[i]);
  /* second moment about E[X], also called second central moment (Var[X]) */
  //for (j=1; j <= chain->max_distance; j++) {
/*------------------------------------
  By Yujin Chung on 3/4/2013
  In BaRT, lambda is the same across the segments
  ---------------------------------------- */
  // lg = (double)(j) * lnl - (chain->lambda) - chain->factorial[j];
    // lg = (double)(j) * lnl - (chain->segment[i]->lambda) - chain->factorial[j];

    /*------------------------------------
      By Yujin Chung on 2/7/2013
      weight w and related parameters were removed
      ---------------------------------------- */
  // lg = exp (lg);
    // lg = exp (w * lg);
  //  (*r2) += ( ((double)(j)-(*r1)) * ((double)(j)-(*r1)) * lg ); 
  // }
/*------------------------------------
  By Yujin Chung on 3/4/2013
  lnN contains the partition function of the full Gibbs prior distribution on gene trees
  ---------------------------------------- */
  // (*r2) /= exp (chain->lnN);
  // (*r2) /= exp (chain->lnN[i]);
}










/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "propose_swap_chains".
---------------------------------------- */
void
propose_swap_chains_gibbs (chain_data_gibbs chain1, chain_data_gibbs chain2, acceptance_frequency *cold_freq, acceptance_frequency *heat_freq)
{
  double lnl_1, lnl_2;

  cold_freq->propose[FREQswapchain]++;
  heat_freq->propose[FREQswapchain]++;
  
  full_likelihood_gibbs (&lnl_1, chain1);
  full_likelihood_gibbs (&lnl_2, chain2);
  
  if (log (biomc2_rng_uniform_pos (random_number)) < ((chain2->kT - chain1->kT) * (lnl_1 - lnl_2)) ) {
    lnl_1      = chain1->kT;
    chain1->kT = chain2->kT;
    chain2->kT = lnl_1;
    
    cold_freq->accept[FREQswapchain]++;
    heat_freq->accept[FREQswapchain]++;
  }
}

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "proposal_update_lambda".
  By Yujin Chung on 3/6/2013:
  In biomc2, this function update one lambda of a vector of lambda's at a time.
  In BaRT, lambda is scalar and we have the same lambda across the segment.
  Therefore, biomc2 has one more argument "int n" to indicate which lambda to update, which is not necessary in BaRT.
---------------------------------------- */
// By Yujin Chung on 3/6/2013:
// To propose a new value of lambda, calculate the acceptance probability, and replace lambda if required
void proposal_update_lambda_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)
{
  
  double nlambda; // proposed lambda by Yujin Chung 3/6/2013
  double ln_r, transition;
  double sumD=0.; // the sum of RF distances by Yujin Chung 3/6/2013
  int i; // index for a loop
  
  /*------------------------------------
    By Yujin Chung on 2/7/2013
    weight w and related parameters were removed
  ---------------------------------------- */
  // double w = chain->segment[n]->penalty + WEIGHT;
	
  afreq->propose[FREQlambda]++;

/*------------------------------------
  By Yujin Chung on 3/4/2013
  In BaRT, lambda is the same across the segments
---------------------------------------- */
  // By Yujin Chung on 3/6/2013:
  // This function proposes a new value of lambda
  update_variable (&nlambda, &(chain->lambda), &(chain->update_lambda));
  //  update_variable (&nlambda, &(chain->segment[n]->lambda), &(chain->update_lambda));

  /*-----------------------------------------------------------
    By Yujin Chung on 2/28/2013
   The following calculates the partition function of the Gibbs prior distribution with a new lambda  
   Note that "lambda" is used for beta which scales with 1/recombination rate (see Chung et al.)
    By Yujin Chung on 3/4/2013
    lnN_proposal contains the partition function of the Gibbs prior distribution of all gene trees.
  -------------------------------------------------------*/
 // By Yujin Chung on 3/6/2013:
  // to calculate the partition function with the new lambda
  partitionFunction_gibbs (&(chain->lnN_proposal), &nlambda, chain);
  // lnTruncate_gibbs (&(chain->lnN_proposal[n]), &nlambda, chain);

  /*------------------------------------
    By Yujin Chung on 2/7/2013
    weight w and related parameters were removed
    By Yujin Chung on 3/4/2013
    lambda is the same across the segments
    lnN_proposal contains the partition function of the Gibbs prior distribution of all gene trees.
  ---------------------------------------- */
  for (i=0; i < chain->n_segments - 1; i++)
    sumD += chain->segment[i]->dist;
  ln_r  = (chain->lambda - nlambda) * (sumD+ chain->hp_lambda_beta);
  ln_r += ( ( chain->hp_lambda_alpha - 1.) * (log (nlambda) - log (chain->lambda)) );
  ln_r += ( chain->lnN - chain->lnN_proposal ); 
  /*---
   ln_r  = (chain->segment[n]->lambda - nlambda) * (w + chain->hp_lambda_beta);
  ln_r += ( (w * (double)(chain->segment[n]->dist) + chain->hp_lambda_alpha - 1.) * (log (nlambda) - log (chain->segment[n]->lambda)) );
ln_r += ( chain->lnN - chain->lnN_proposal );
  ------*/

  ln_r = (chain->kT * ln_r);
 /*------------------------------------
    By Yujin Chung on 3/4/2013
    lambda is the same across the segments
  ---------------------------------------- */
  transition = log (nlambda) - log (chain->lambda);
  // transition = log (nlambda) - log (chain->segment[n]->lambda);
  
  if (chain->acceptance_probability_gibbs (chain, &ln_r, &transition)) {
    /*------------------------------------
    By Yujin Chung on 3/4/2013
    lambda is the same across the segments
    lnN_proposal contains the partition function of the Gibbs prior distribution of all gene trees.
    ---------------------------------------- */
    chain->lnN = chain->lnN_proposal;
    chain->lambda = nlambda;
    /*
    chain->lnN[n] = chain->lnN_proposal[n];
    chain->segment[n]->lambda = nlambda;
    */
    afreq->accept[FREQlambda]++;
  }
}

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "proposal_update_substitution_rate"
---------------------------------------- */
void 
proposal_update_substitution_rate_gibbs (chain_data_gibbs chain, int n, topology t, acceptance_frequency *afreq)
{
  double new_rate, a2, diff_likelihood, transition;
  
  afreq->propose[FREQrate]++;
  
  update_variable (&new_rate, &(chain->segment[n]->rate), &(chain->update_rate));
  update_Q_matrix_gibbs (chain->segment[n], new_rate);

  chain->ln_likelihood_proposal_gibbs (chain->segment[n], t);

  diff_likelihood = chain->segment[n]->likelihood_proposal - chain->segment[n]->likelihood_accepted;

  a2 = chain->kT * (diff_likelihood + (chain->segment[n]->rate - new_rate)/chain->prior_rate);
  
  transition = log (new_rate) - log (chain->segment[n]->rate);
  
  if (chain->acceptance_probability_gibbs (chain, &a2, &transition)) {
    accept_likelihood_proposal_gibbs (chain->segment[n], t);
    chain->segment[n]->rate = new_rate;
    t->likelihood_accepted += diff_likelihood;
    afreq->accept[FREQrate]++;
  }
  else
    update_Q_matrix_gibbs (chain->segment[n], chain->segment[n]->rate);
}

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "proposal_update_substitution_kappa"
---------------------------------------- */
void 
proposal_update_substitution_kappa_gibbs (chain_data_gibbs chain, int n, topology t, acceptance_frequency *afreq)
{
  double new_kappa, a2, diff_likelihood, transition;
  
  afreq->propose[FREQkappa]++;
  
  update_variable (&new_kappa, &(chain->segment[n]->kappa), &(chain->update_kappa));
  update_Q_eigenvalues_gibbs (chain->segment[n], new_kappa);
  update_Q_matrix_gibbs (chain->segment[n], chain->segment[n]->rate);

  chain->ln_likelihood_proposal_gibbs (chain->segment[n], t);
  diff_likelihood = chain->segment[n]->likelihood_proposal - chain->segment[n]->likelihood_accepted;
  
  a2 = chain->kT * (diff_likelihood + (chain->segment[n]->kappa - new_kappa )/chain->prior_kappa);
  
  transition = log (new_kappa) - log (chain->segment[n]->kappa);

  if (chain->acceptance_probability_gibbs (chain, &a2, &transition)) {
    accept_likelihood_proposal_gibbs (chain->segment[n], t);
    chain->segment[n]->kappa = new_kappa;
    t->likelihood_accepted += diff_likelihood;
    afreq->accept[FREQkappa]++;
  }
  else {
    update_Q_eigenvalues_gibbs (chain->segment[n], chain->segment[n]->kappa);
    update_Q_matrix_gibbs (chain->segment[n], chain->segment[n]->rate);
  }
}

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "proposal_update_prior_rate".
---------------------------------------- */
void
proposal_update_prior_rate_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)
{
  double new_prior, diff, sum = 0., transition;
  int i;
  
  afreq->propose[FREQrateprior]++;
  
  update_variable (&new_prior, &(chain->prior_rate), &(chain->update_rate));

  for (i=0; i < chain->n_segments; i++) 
    sum += chain->segment[i]->rate;

  diff  = sum * (1./chain->prior_rate - 1./new_prior);
  diff += ( 1./chain->hp_rate * (chain->prior_rate - new_prior) );
  diff += ( (double)(chain->n_segments) * (log (chain->prior_rate) - log (new_prior)) );
  
  diff = (chain->kT * diff);
  transition = log (new_prior) - log (chain->prior_rate);

  if (chain->acceptance_probability_gibbs (chain, &diff, &transition)) {
    chain->prior_rate = new_prior;
    afreq->accept[FREQrateprior]++;
  }
}


/*------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "proposal_update_prior_kappa".
---------------------------------------- */
void
proposal_update_prior_kappa_gibbs (chain_data_gibbs chain, acceptance_frequency *afreq)
{
  double new_prior, diff, sum = 0., transition; 
  int i;
  
  afreq->propose[FREQkappaprior]++;
  
  update_variable (&new_prior, &(chain->prior_kappa), &(chain->update_kappa));
  
  for (i=0; i < chain->n_segments; i++) 
    sum += chain->segment[i]->kappa;
  
  diff  = sum * (1./chain->prior_kappa - 1./new_prior);
  diff += ( 1./chain->hp_kappa * (chain->prior_kappa - new_prior) );
  diff += ( (double)(chain->n_segments) * (log (chain->prior_kappa) - log (new_prior)) );
  
  diff = (chain->kT * diff);
  transition = log (new_prior) - log (chain->prior_kappa);
  
  if (chain->acceptance_probability_gibbs (chain, &diff, &transition)) {
    chain->prior_kappa = new_prior;
    afreq->accept[FREQkappaprior]++;
  }
}


