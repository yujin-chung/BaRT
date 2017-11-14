/* Copyright (C) 2006  Leonardo de Oliveira Martins
 * 
 * leo at lbm ab a u-tokyo ac jp 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the 
 * Free Software Foundation, Inc., 51 Franklin Street, 
 * Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*! \file chain_data.c 
 *  \brief 
 */

#include "chain_data.h"

/* local function prototypes */

/*-----------------------------------------------
  By Yujin Chung on 2/6/2013:
  The following functions are created 
  because of new structure chain_data_gibbs and phylogeny_gibbs
------------------------------------------------- */
void read_biomc_file_gibbs (char *filename, chain_data_gibbs chain);
void initialize_segment_chain_data_gibbs (chain_data_gibbs chain, nexus_alignment a, nexus_treespace t);
void initialize_segment_equilibrium_frequencies_gibbs (chain_data_gibbs chain);
void initialize_breakpoint_chain_data_gibbs (chain_data_gibbs chain, nexus_alignment a);
void initialize_poisson_chain_data_gibbs (chain_data_gibbs chain);
void store_chain_data_site_pattern_gibbs (chain_data_gibbs chain, nexus_alignment a);
void initialize_phylogeny_idx_fixed_size_gibbs (chain_data_gibbs chain, int n, int nchar, double *siteweight, int *idx, int *nidx);

void ln_likelihood_moved_branches_nodata_gibbs (phylogeny_gibbs *p, topology t);
void ln_likelihood_proposal_nodata_gibbs (phylogeny_gibbs p, topology t);

void initialize_eigenvector_gibbs (chain_data_gibbs chain);



/*-------------------------------------------------------------
  By Yujin Chung on 2/6/2013:
  Structure "new_chain_data_from_file_gibbs"
  -- contains all parameters we need to keep track through MCMC simulations.
  -- set default values of parameters
  -- replaces "new_chain_data_from_file" from biomc2.

  List of changes:
  -- penalty w and related objects was removed
  -- "chain_data_gibbs" is used instead of "chain_data", "chain_data_struct"
----------------------------------------------------------------- */
chain_data_gibbs 
new_chain_data_from_file_gibbs (char *filename)
{
  int i;
  nexus_alignment a = NULL;
  nexus_treespace t = NULL;
  chain_data_gibbs chain;
  
  chain = 
    (chain_data_gibbs) biomc2_malloc (sizeof (struct chain_data_gibbs));
 
  /* default values */

  chain->update_rate       = 0.5; // Perturbation value \f$\xi_\mu\f$ for substitution rate. -- Yujin Chung 2/6/2013 
  chain->update_kappa      = 0.5; // Perturbation value \f$\xi_\kappa\f$ for ti/tv ratio.  -- Yujin Chung 2/6/2013 
  /*--------------------------------------------
    2/6/2013
    Perturbation value \f$\xi_\lambda\f$ for prior topology distance \lambda
    Warning(notation)!! we used \beta instead of \lambda in our papers and \beta scales with the inverse of recombination rate.
  ----------------------------------------------*/
  chain->update_lambda     = 0.5; 
  
  /*--------------------------------
    By Yujin Chung on 2/6/2013:
    chain_data_gibbs "chain" does not have members of the followings.
    Penalty w and related parameters were removed.
  ------------------------------------ */
  // chain->update_penalty    = 0.5; 
  // chain->hp_penalty_alpha  = 1.;
  // chain->hp_penalty_beta   = 1.;

  /* ---------------------------------
     By Yujin Chung on 2/20/2013:
     parameters of hyperprior distribution 
     for lambda (beta in Gibbs distribution, see our papers)
  -------------------------------- */
  chain->hp_lambda_alpha   = 1.;
  chain->hp_lambda_beta    = 1.;
  
  /*------------------------------------
    By Yujin Chung on 2/20/2013
    --fixit: don't need "max_distance"?
    ---------------------------------------- */
  chain->max_distance      = -1;
  chain->ratio_spr         = 5; // Iterations at which we update topologies by SPR instead of NNI -- Yujin Chung 2/20/2013 

  chain->hp_rate           = 1.; // Hyperprior \f$\mathcal{M}\f$ for substitution rate (exponential). -- Yujin Chung 2/6/2013 
  chain->hp_kappa          = 1.; // Hyperprior \f$\mathcal{K}\f$ for ti/tv ratio (exponential).  -- Yujin Chung 2/6/2013 

  chain->ratio_topol_update = 1;
  chain->beta_final        = 1.;
  chain->kT_heat           = 0.9;
  chain->swap_interval     = 2;
  chain->n_cycles          = 2;
  chain->kT_cycles         = 0.8;
  chain->n_mini            = 2.;
  chain->kT_mini           = 0.8;
  chain->prior    = false;
  chain->ngen     = 10000; 
  chain->nsamples = 100; 
  chain->warmup   = 100;
  chain->burnin   = 1000;
  chain->kT        = 1.;
  chain->tree_filename[0] = '\0';
  chain->seq_filename[0]  = '\0';

  chain->n_segments = 0;
  chain->treefile = NULL;

  #ifdef DEBUG
  printf("REMOVE: Calling read_biomc_file_gibbs (filename, chain);\n");
#endif //DEBUG
  read_biomc_file_gibbs (filename, chain);
  #ifdef DEBUG
  printf("REMOVE: Done with read_biomc_file_gibbs (filename, chain);\n\n");
#endif //DEBUG
  if (strlen (chain->tree_filename) ) t = read_nexus_treespace_file (chain->tree_filename);
  else biomc2_error ( "no TREEFILE command in control file\n");

  if (strlen (chain->seq_filename) ) a = read_nexus_alignment_file (chain->seq_filename);
  else biomc2_error ( "no SEQFILE command in control file\n");

  if ((chain->n_segments < 2) || (chain->n_segments > a->nchar)) {
    if (a->n_charset > 1) chain->n_segments = a->n_charset;
    else { chain->n_segments = a->nchar/2; a->n_charset = 0; }
  }
  else a->n_charset = 0;

  if (a->ntax != t->nleaves) {
    del_nexus_alignment (a);
    del_nexus_treespace (t);
    biomc2_error ( 
		  "ntax disagree between sequence and tree files\n");
  }
  if (a->ntax > MaxLeaves) {
    del_nexus_alignment (a);
    del_nexus_treespace (t);
    biomc2_error ( 
		  "number of leaves > %d, cannot calculate tree distance\n", MaxLeaves);
  }

  chain->nsites = a->nchar;
  chain->ntax = a->ntax;
  chain->taxlabel = a->taxlabel; 
  /*-------------------------
    By Yujin Chung on 2/20/2013
    Compare two topologies based on their splits, here represented as
    bit-strings (e.g. the pair 00110001, 11001110). 
    
    split_space_struc is defined in "topology.h"     
    --------------------------*/
  chain->split = new_split_space (t->nleaves);

  if (chain->ratio_spr < 1)    
    chain->ratio_spr    = 2 * chain->ngen;
  if (chain->ratio_topol_update < 1) 
    chain->ratio_topol_update = 0;
  if (chain->ratio_topol_update > chain->ngen/100) 
    chain->ratio_topol_update = chain->ngen/100;
  
  if (chain->kT_heat < 1.e-6) chain->kT_heat = 1.e-6;
  if (chain->kT_heat > (1. - 1.e-8))    chain->kT_heat = 1. - 1.e-8;
  if (chain->swap_interval > chain->ngen/100) chain->swap_interval = chain->ngen/100;
  /* If ngen < 100 the above line will make swap_interval=0 */
  if (chain->swap_interval < 1) chain->swap_interval = 1;
  
  if ((chain->ngen) < 20) chain->ngen = 20;
  
  if (chain->warmup < 1) chain->warmup = 1;
  if (chain->burnin < 1) chain->burnin = 1;

  if (chain->beta_final < 1.e-4) chain->beta_final = 1.e-4;
  if (chain->beta_final > 1.e8)  chain->beta_final = 1.e8;
  chain->beta_zero = log( (double)(chain->warmup) + EXP_1 );
  chain->beta_zero = chain->beta_final/chain->beta_zero;
  chain->acceptance_probability_gibbs = &(acceptance_probability_bayesian_gibbs);

  /* Absurd limits */
  if (chain->kT_cycles < 1.e-6) chain->kT_cycles = 1.e-6;
	if (chain->kT_cycles > 1.e6)  chain->kT_cycles = 1.e6;
	if (chain->n_cycles < 1) { chain->n_cycles = 0; chain->kT_cycles = 1.; }
  if (chain->n_cycles > 1e4) chain->n_cycles = 1e4;
	if (chain->kT_mini < 1.e-6) chain->kT_mini = 1.e-6;
	if (chain->kT_mini > 1.e6)  chain->kT_mini = 1.e6;
	if (chain->n_mini < 1)   chain->n_mini = 1;
	if (chain->n_mini > 1e4) chain->n_mini = 1e4;

	if(chain->max_distance > chain->split->inodes) chain->max_distance = chain->split->inodes;
	
	if (chain->max_distance < 1) { chain->max_distance   = 1; chain->max_dist_split = 0; }
	else chain->max_dist_split = chain->max_distance; 

	chain->log2 = log (2.);
	/* SPR neighbourhood of a topology */
	chain->logDspr = log (2.*(double)(a->ntax - 3)) + log ((double)(2*a->ntax - 7));
	/* NNI neighbourhood of a topology */
	chain->logDnni = log (2.*(double)(a->ntax - 3));
	chain->logD = &(chain->logDnni);

	// printf("\nREMOVE: Initialize the breakpoints \n\n");
	initialize_breakpoint_chain_data_gibbs (chain, a);
	// Yujin Chung Aug 19, 2016
	// The function below initializes the trees of segments
	// printf("\nREMOVE: Initializ the segment and chain data \n\n");
	initialize_segment_chain_data_gibbs (chain, a, t);
  // printf("\nREMOVE: Initializ the segments' equilibrium and frequencies  \n\n");
  initialize_segment_equilibrium_frequencies_gibbs (chain);  
  // printf("\nREMOVE: Initializ the poisson and chain data \n\n");
  initialize_poisson_chain_data_gibbs (chain);
  
 
  if (chain->prior) { // simulation of prior -- YC
    #ifdef DEBUG
    printf("REMOVE: chain->prior = 1 (simulation of prior)\n");
#endif //DEBUG
    chain->ln_likelihood_moved_branches_gibbs = &(ln_likelihood_moved_branches_nodata_gibbs);
    chain->ln_likelihood_proposal_gibbs = &(ln_likelihood_proposal_nodata_gibbs);
  }
  else {	// simulation of posterior -- YC 
    #ifdef DEBUG
    printf("REMOVE: chain->prior = 0 (simulation of posterior)\n");
#endif //DEBUG
    chain->ln_likelihood_moved_branches_gibbs = &(ln_likelihood_moved_branches_gibbs);
    chain->ln_likelihood_proposal_gibbs = &(ln_likelihood_proposal_gibbs);
  }

  /* assuming there is only one topology */
  chain->topol->first->likelihood_accepted = 0.;
  if (!chain->prior) // simulation of posterior -- YC
    {
      #ifdef DEDUG
      printf("REMOVE: chain->prior = 0 (simulation of posterior)\n");
      //  printf("** REMOVE: YC Aug 19, 2016 **\nLook at the function accept_likelihood_proposal_gibbs() where the prior is computed.. probably.\n ");
      // printf("** REMOVE: YC Sep 9, 2016 **\nThe function seems to just compute the likelihood. I guess the tree prior isn't be computed yet.\n ");
#endif // DEBUG
      for (i=chain->topol->first->first_segment; i <= chain->topol->first->last_segment; i++) {
	ln_likelihood_proposal_gibbs (chain->segment[i], chain->topol->first);
	accept_likelihood_proposal_gibbs (chain->segment[i], chain->topol->first);
	//			printf ("segment %d \t lnLik: %f\n", i , chain->segment[i]->likelihood_accepted);
	chain->topol->first->likelihood_accepted += chain->segment[i]->likelihood_accepted;
      }
    }

  printf ("initial topology likelihood : %f\n", chain->topol->first->likelihood_accepted);
  
  del_nexus_alignment (a);
  del_nexus_treespace (t);
 	
  return chain;
}


/*-----------------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "initialize_breakpoint_chain_data"
--------------------------------------------------- */
void
initialize_breakpoint_chain_data_gibbs (chain_data_gibbs chain, nexus_alignment a)
{
  int i;

  chain->breakpoint = (int*) biomc2_malloc ((chain->n_segments - 1) * sizeof (int));
  if (a->n_charset) 
    for (i=0; i < chain->n_segments-1; i++) 
      chain->breakpoint[i] = (int)((double)(a->charset_end[i] + a->charset_start[i+1])/2.);
  else 
    for (i=0; i < chain->n_segments-1; i++) 
      chain->breakpoint[i] = (int)((double)((i+1) * chain->nsites)/(double)(chain->n_segments)); 
}

/*-----------------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "initialize_segment_chain_data"
--------------------------------------------------- */
void
initialize_segment_chain_data_gibbs (chain_data_gibbs chain, nexus_alignment a, nexus_treespace t)
{ 
  int i, j, lik_vector_size;
  double rate, site_weight[a->nchar]; 
  int  idx[a->nchar], n_idx, label_order[t->nleaves];

  if (chain->n_mini > chain->n_cycles) 
    lik_vector_size = chain->n_mini;
  else 
    lik_vector_size = chain->n_cycles;

  /* vector label_order will receive positions based on a->taxlabel */
  map_hashtable_order (label_order, a->taxlabel_hash, t->taxlabel, t->nleaves);
  /* trees in nexus_treespace will have leaves ordered according to alignment
   * and rooted at a->taxlabel[0] */
  order_nexus_treespace_id (t, label_order);

  chain->site_pattern = (int*) biomc2_malloc (chain->nsites * sizeof (int));
  for (i=0; i < chain->nsites; i++) 
    chain->site_pattern[i] = i;
  /* this function will modify nexus_alignment (shrink it to distinct patterns) */
  store_chain_data_site_pattern_gibbs (chain, a);


  chain->segment = (phylogeny_gibbs*) biomc2_malloc (chain->n_segments * sizeof (phylogeny));

  for (i=0; i < 6; i++) 
    chain->pi[i] = 0.;

  for (i=0; i < chain->n_segments; i++) {
    initialize_phylogeny_idx_fixed_size_gibbs (chain, i, a->nchar, site_weight, idx, &n_idx);
    chain->segment[i] = new_phylogeny_from_nexus_gibbs (a, site_weight, idx, n_idx, lik_vector_size);
    chain->segment[i]->pi = chain->pi;
    for (j=0; j < 4; j++) {
      chain->segment[i]->z1[j] = chain->z1[j];
      chain->segment[i]->z2[j] = chain->z2[j];
    }
  }
  chain->topol = new_linked_topol (chain->n_segments, chain->ntax);
  j = biomc2_rng_uniform_int (random_number, t->ntrees);
  for (i=0; i < chain->topol->size; i++)
    copy_topology_from_nexus_tree (chain->topol->tree[i], t->T[j]); 

  if (t->T[j]->has_branches) {
    rate = average_rate_nexus_tree (t->T[j]);
    for (i=0; i < chain->n_segments; i++) {
      chain->segment[i]->rate = rate;
    }
  }
  else {
    for (i=0; i < chain->n_segments; i++) {
      chain->segment[i]->rate = average_rate_distance_method_gibbs (chain->segment[i]);
      /* if segment has no segregating sites our initial guess comes from hyperprior  */
      if (chain->segment[i]->rate <= 0)
	chain->segment[i]->rate = chain->hp_rate;
    }
  }
  
  /* random topology generator */
  //	j = biomc2_rng_uniform_int (random_number, t->ntrees);
  //	copy_topology_from_nexus_tree (chain->topol->tree[chain->topol->size - 2], t->T[j]);
}


/* -------------------------------------
   By Yujin Chung on 2/6/2013
   It replaces "initialize_segment_equilibrium_frequencies".
----------------------------------------- */
void
initialize_segment_equilibrium_frequencies_gibbs (chain_data_gibbs chain)
{
  int i;
  double average = 0., total_freq = 0.;

  for (i=0; i < chain->n_segments; i++) {
    update_equilibrium_frequencies_gibbs (chain->pi, chain->segment[i]);
  }
  for (i=0; i < 4; i++) {
    if (chain->pi[i] < 1.) chain->pi[i] = 0.5; /* if data does not have i-th base (throw error?) */
    total_freq += chain->pi[i];
  }
  // chain->pi contains the empirical frequencies of ACGT -YC
  printf("Empirical frequencies of ACGT: ");
  for (i=0; i < 4; i++) { 
    chain->pi[i] /= total_freq;
    printf ("[ %f ] ", chain->pi[i]); 
  }
  printf("\n");
  
  chain->pi[4] = chain->pi[1] + chain->pi[3]; // PI[4] = PI_Y = PI_C + PI_T
  chain->pi[5] = chain->pi[0] + chain->pi[2]; // PI[5] = PI_R = PI_A + PI_G
  initialize_eigenvector_gibbs (chain);

  for (i=0; i < chain->n_segments; i++) {
    chain->segment[i]->kappa = chain->hp_kappa;
    update_Q_eigenvalues_gibbs (chain->segment[i], chain->segment[i]->kappa);
    update_Q_matrix_gibbs (chain->segment[i], chain->segment[i]->rate);
    average += chain->segment[i]->rate;
  }
	
  chain->prior_rate  = average/(double)(i);
  chain->prior_kappa = chain->hp_kappa;
}

/*-----------------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "initialize_poisson_chain_data"

  Log of chaines
  -- penalty w and related parameters were defunctionalized


  By Yujin Chung on 2/6/2013
  --fixit: related to the prior distribution on gene trees
------------------------------------------------- */
void
initialize_poisson_chain_data_gibbs (chain_data_gibbs chain)
{
  int i;
  
  /*------------------------------------
    By Yujin Chung on 3/4/2013
    In biomc2, lnN and lnN_proposal are array, containing partition functions for two segments.
    In BaRT, lnN and lnN_proposal contains the partition function of the Gibbs prior distribution on all gene trees.
    ---------------------------------------- */
  // chain->lnN = (double*) biomc2_malloc ((chain->n_segments - 1) * sizeof (double));
  // chain->lnN_proposal = (double*) biomc2_malloc ((chain->n_segments - 1) * sizeof (double));

  for (i=0; i < chain->n_segments - 1; i++) {
    /*---------------------------------------------------
      By Yujin Chung on 2/6/2013:
      It is related to the prior on trees
      This should be revised
    ------------------------------------------------------ */
    chain->lambda  = chain->hp_lambda_alpha / chain->hp_lambda_beta;
    // chain->segment[i]->lambda  = chain->hp_lambda_alpha / chain->hp_lambda_beta;

    /*-----------------------------------
      By Yujin Chung on 2/6/2013
      Penalty w and related parameters are removed
     -----------------------------------------*/
    // chain->segment[i]->penalty = chain->hp_penalty_alpha / chain->hp_penalty_beta;
  }

  chain->factorial = (double*) biomc2_malloc ((chain->ntax) * sizeof (double));

  /*------------------------------------
    By Yujin Chung on 2/22/2013
    chain->factorial[k] stores log(k!).
    fixit: we don't need factorial as a member of chain
  ---------------------------------------- */
  chain->factorial[0] = chain->factorial[1] = 0.;
  for (i=2; i < chain->ntax; i++) 
    chain->factorial[i] = chain->factorial[i-1] + log ((double)(i));
}

/*----------------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "store_chain_data_site_pattern"
------------------------------------------------- */
void
store_chain_data_site_pattern_gibbs (chain_data_gibbs chain, nexus_alignment a)
{
  bool equal;
  int s1, s2, seq, lastsite = a->nchar;
  int chunked_gain, index[a->nchar];

  for (s1=0; s1 < lastsite; s1++) 
    index[s1] = s1;
  
  for (s1 = 0; s1 < lastsite - 1; s1++) {
    for (s2 = s1 + 1; s2 < lastsite; s2++) {
      
      for (equal = true, seq = 0; (equal == true) && (seq < a->ntax); seq++)
        if (a->character[seq][s1] != a->character[seq][s2]) 
	  equal = false;
      
      /* s1 and s2 have identical patterns */
      if (equal == true) {
	chain->site_pattern[index[s2]] = s1;
	lastsite--;
	if (s2 < lastsite) {
	  for (seq = 0; seq < a->ntax; seq++)
	    a->character[seq][s2] = a->character[seq][lastsite];
	  index[s2] = index[lastsite];
	}
	s2--;/* compare again, since character[][s2] now contains last site */
      }// if (equalpattern)
    }// for (s2)
    chain->site_pattern[index[s1]] = s1;
  }// for (s1)
  chain->site_pattern[index[s1]] = s1;
  
  chunked_gain = a->nchar - lastsite;
	
  if (chunked_gain > 0) {
    a->nchar -= chunked_gain;
    for (seq = 0; seq < a->ntax; seq++) {
      a->character[seq] =
	(char *) biomc2_realloc ((char*) a->character[seq], (a->nchar + 1) * sizeof (char));
      a->character[seq][a->nchar] = '\0';
    }
  }
}

/*----------------------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "initialize_phylogeny_idx_fixed_size".
------------------------------------------------------ */
void 
initialize_phylogeny_idx_fixed_size_gibbs (chain_data_gibbs chain, int n, int nchar, double *siteweight, int *idx, int *nidx)
{
  int i, start, end;
  (*nidx) = 0;
  
  for (i = 0; i < nchar; i++) siteweight[i] = 0.;
	
  if (n==0) start=0;
  else start = chain->breakpoint[n-1];
  
  if (n==chain->n_segments-1) end = chain->nsites;  
  else end = chain->breakpoint[n]; 
  
  /* frequency of pattern in segment */
  for (i = start; i < end; i++)
    siteweight[chain->site_pattern[i]] += 1.;
  
  /* vector pointing to patterns */
  for (i = 0; i < nchar; i++)
    if (siteweight[i] > 0.) idx[(*nidx)++] = i;
  /* DEBUG */
//	printf ("charset segment%d \t = %d \t- %d;\n", n+1, start+1, end);
}


/*-------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "ln_likelihood_moved_branches_nodata".
--------------------------------------- */
void
ln_likelihood_moved_branches_nodata_gibbs (phylogeny_gibbs *p, topology t) { return; }

/*-------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "ln_likelihood_proposal_nodata".
--------------------------------------- */
void 
ln_likelihood_proposal_nodata_gibbs (phylogeny_gibbs p, topology t) { return; }


/*-------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "acceptance_probability_bayesian".
--------------------------------------- */
bool
acceptance_probability_bayesian_gibbs (chain_data_gibbs chain, double *acceptance, double *transition)
{
  if ( log (biomc2_rng_uniform (random_number)) < ((*acceptance) + (*transition))) 
    return true;
  else 
    return false;
}

/*-------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "acceptance_probability_bayesian".
--------------------------------------- */
bool
acceptance_probability_annealing_gibbs (chain_data_gibbs chain, double *acceptance, double *transition)
{
  if ( biomc2_rng_uniform (random_number) < (1./(1. - (chain->beta_n * ((*acceptance) + (*transition))))) ) 
    return true;
  else 
    return false;
}


/*-------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "initialize_eigenvector".
--------------------------------------- */
void
initialize_eigenvector_gibbs (chain_data_gibbs chain)
{

	chain->z1[0][0] = chain->pi[0];
	chain->z1[0][1] = chain->pi[1];
	chain->z1[0][2] = chain->pi[2];
	chain->z1[0][3] = chain->pi[3];
	chain->z1[1][0] = -chain->pi[0] * chain->pi[4];
	chain->z1[1][1] =  chain->pi[1] * chain->pi[5];
	chain->z1[1][2] = -chain->pi[2] * chain->pi[4];
	chain->z1[1][3] =  chain->pi[3] * chain->pi[5];
	chain->z1[3][0] = chain->z1[3][2] = chain->z1[2][1] = chain->z1[2][3] = 0.;
	chain->z1[3][3] = chain->z1[2][0] = 1.;
	chain->z1[3][1] = chain->z1[2][2] = -1.;

	chain->z2[0][0] =	 chain->z2[0][1] =	chain->z2[0][2] = chain->z2[0][3] = 1.; 
	chain->z2[1][0] = -1./chain->pi[5];
	chain->z2[1][1] =  1./chain->pi[4];
	chain->z2[1][2] = -1./chain->pi[5];
	chain->z2[1][3] =  1./chain->pi[4];
	chain->z2[3][0] =  chain->z2[3][2] = chain->z2[2][1] = chain->z2[2][3] = 0.;
	chain->z2[3][1] = -chain->pi[3]/chain->pi[4];
	chain->z2[3][3] =  chain->pi[1]/chain->pi[4];
	chain->z2[2][0] =  chain->pi[2]/chain->pi[5];
	chain->z2[2][2] = -chain->pi[0]/chain->pi[5];

}


/*-------------------------------------
  By Yujin Chung on 2/6/2013
  It replaces "read_biomc_file".
--------------------------------------- */
void
read_biomc_file_gibbs (char *filename, chain_data_gibbs chain) 
{
  FILE *file;
  char *line=NULL, *needle_tip;
  size_t linelength = 0;

  file = biomc2_fopen (filename, "r");

  /* the variable *line should point always to the same value (no line++ or
   * alike) */
  while (biomc2_getline (&line, &linelength, file) != -1) {
    line = remove_nexus_comments (&line, &linelength, file);
    /* do not parse empty lines ( = "\n\0") */	
    if (strlen (line) > 1) { 
      /* read hyperprior information */
      
      if (strcasestr (line, "RATE") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++; // Exp(rate) => E[x] = rate
	if ((sscanf (needle_tip, " %lf ", &chain->hp_rate)) != 1)	
          biomc2_error ( "could not read RATE information (it should be a single float or integer number)\n");
      }
      
      if (strcasestr (line, "KAPPA") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++; // Exp(kappa) => E[x] = kappa
	if ((sscanf (needle_tip, " %lf ", &chain->hp_kappa)) != 1)	
          biomc2_error ( "could not read KAPPA information (it should be a single float or integer number)\n");
      }
      
      else if (strcasestr (line, "DISTANCE") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++; // gamma (a,b) -> E[x] = a/b 
        if ((sscanf (needle_tip, " %lf %lf ", &(chain->hp_lambda_alpha), &(chain->hp_lambda_beta))) != 2)	
          biomc2_error ( "could not read DISTANCE information (it should be _two_ float or integer numbers, one for alpha and one for beta)\n");
      }
      
      else if (strcasestr (line, "TREEFILE") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
        if ((sscanf (needle_tip, " %s ", chain->tree_filename)) != 1)	
	  biomc2_error ( "could not read TREEFILE information (it should be a single file name without spaces or non-standard characters)\n");
      }
      
      else if (strcasestr (line, "SEQFILE") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	if ((sscanf (needle_tip, " %s ", chain->seq_filename)) != 1)	
	  biomc2_error ( "could not read SEQFILE information (it should be a single file name without spaces or non-standard characters)\n");
      }
      
      else if (strcasestr (line, "NGEN") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	if ((sscanf (needle_tip, " %d ", &chain->ngen)) != 1)	
          biomc2_error ( "could not read NGEN information (it should be a single, integer value)\n");
      }
      
      else if (strcasestr (line, "WARMUP") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	if ((sscanf (needle_tip, " %d ", &chain->warmup)) != 1)	
	  biomc2_error ( "could not read WARMUP information (it should be a single, integer value)\n");
      }
      
      else if (strcasestr (line, "NSAMPLES") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	if ((sscanf (needle_tip, " %d ", &(chain->nsamples))) != 1)	
	  biomc2_error ( "could not read NSAMPLES information (it should be a single, integer value)\n");
      }
      
      else if (strcasestr (line, "BURNIN") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
        if ((sscanf (needle_tip, " %d ", &(chain->burnin))) != 1)	
          biomc2_error ( "could not read BURNIN information (it should be a single, integer value)\n");
      }
      
      else if (strcasestr (line, "NSEGMENTS") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	if ((sscanf (needle_tip, " %d ", &(chain->n_segments))) != 1)	
          biomc2_error ( "could not read NSEGMENTS information (it should be a single, integer number)\n");
      }
      
      else if (strcasestr (line, "PRIOR") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	if ((sscanf (needle_tip, " %d ", (int*)&(chain->prior))) != 1)	
	  biomc2_error ( "could not read PRIOR information (it should be a single, integer number)\n");
      }
      
      else if (strcasestr (line, "UPDATE") && (needle_tip = strcasestr (line, "="))) {
	needle_tip++;
	/*--------------------------------------
	  By Yujin Chung on 2/6/2013:
	  penalty w and related parameters were removed
	--------------------------------------------- */
	if ((sscanf (needle_tip, " %lf %lf %lf ", &(chain->update_rate), &(chain->update_kappa), &(chain->update_lambda))) != 3)	
          biomc2_error ( "could not read UPDATE information (it should be _three_ real numbers); Try removing the line with \"update\")\n");
	// if ((sscanf (needle_tip, " %lf %lf %lf %lf ", &(chain->update_rate), &(chain->update_kappa), &(chain->update_lambda), &(chain->update_penalty))) != 4)	
        //  biomc2_error ( "could not read UPDATE information (it should be _four_ real numbers); Try removing the line with \"update\")\n");
      }
      
      /*---------------------------------------------
	By Yujin Chung on 2/6/2013
	Penalty w and related parameters were defunctionalized
      -------------------------------------------------------- */
      /*
      else if (strcasestr (line, "PENALTY") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
        if ((sscanf (needle_tip, " %lf %lf ", &(chain->hp_penalty_alpha), &(chain->hp_penalty_beta))) != 2)	
          biomc2_error ( "could not read PENALTY information (it should be _two_ float or integer numbers, one for alpha and one for beta)\n");
      }
      */
      
      /*-----------------------------------------------
	By Yujin Chung on 2/6/2013:
	It is related to prior on trees.
      -----------------------------------------------*/
      else if (strcasestr (line, "MAXDIST") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
        if ((sscanf (needle_tip, " %d ", &(chain->max_distance)) ) != 1)	
          biomc2_error ( "could not read MAXDIST information\n");
			}
 
      else if (strcasestr (line, "PROPOSAL") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
				if ((sscanf (needle_tip, " %d %d ", &(chain->ratio_spr), &(chain->ratio_topol_update) ) ) != 2)
          biomc2_error ( "could not read PROPOSAL information\n");
			}

      else if (strcasestr (line, "ANNEAL") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
	if ((sscanf (needle_tip, " %lf ", &(chain->beta_final)) ) != 1)	
          biomc2_error ( "could not read ANNEAL information\n");
      }

      else if (strcasestr (line, "MC3") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
	if ((sscanf (needle_tip, " %lf %d ", &(chain->kT_heat), &(chain->swap_interval)) ) != 2)	
          biomc2_error ( "could not read MC3 information\n");
      }
      
      else if (strcasestr (line, "AWADHI") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
	if ((sscanf (needle_tip, " %lf %d ", &(chain->kT_cycles), &(chain->n_cycles)) ) != 2)	
          biomc2_error ( "could not read AWADHI information\n");
      }
      
      else if (strcasestr (line, "MINISAMPLER") && (needle_tip = strcasestr (line, "="))) {
        needle_tip++;
	if ((sscanf (needle_tip, " %lf %d ", &(chain->kT_mini), &(chain->n_mini)) ) != 2)	
          biomc2_error ( "could not read MINISAMPLER information\n");
      }
      
    } // if (line)
  } //while (biomc2_getline)
  fclose (file);
  if (line) free (line);
}

/*-------------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "del_chain_data".
--------------------------------------------- */
void
del_chain_data_gibbs (chain_data_gibbs p)
{
  int i;
  
  if (p) {
    if (p->factorial)    free (p->factorial);
    /*------------------------------------
      By Yujin Chung on 3/4/2013
      In biomc2, lnN and lnN_proposal are pointers of array.
      In BaRT, they are "double" (scalar).
      ---------------------------------------- */
    // if (p->lnN)          free (p->lnN);
    // if (p->lnN_proposal) free (p->lnN_proposal);
    if (p->breakpoint)   free (p->breakpoint);
    if (p->site_pattern) free (p->site_pattern);
    
    if (p->segment) {
      for (i=p->n_segments-1; i >=0; i--) del_phylogeny_gibbs (p->segment[i]);
      free (p->segment);
    }
    
    if (p->taxlabel) {
      for (i=p->ntax-1; i>=0; i--)
        if (p->taxlabel[i]) free (p->taxlabel[i]);
      free (p->taxlabel);
    }
    
    del_split_space (p->split);
    del_linked_topol (p->topol);
    
    
    free (p);
  }
}

/*---------------------------------------------
  By Yujin Chung on 2/6/2013:
  It replaces "setup_acceptance_names".
  
------------------------------------------------*/
void setup_acceptance_names_gibbs (void)
{
	acceptance_names[FREQadd_gibbs]        = "breakpoint_addition";
	acceptance_names[FREQremove_gibbs]     = "breakpoint_removal";
	acceptance_names[FREQshift_gibbs]      = "breakpoint_shift_pos";
	acceptance_names[FREQcycle_gibbs]      = "al_awadhi_cycle_sampler";
	acceptance_names[FREQheat_gibbs]       = "heated minisampler";
	acceptance_names[FREQtopol_gibbs]      = "nonrecomb_topology";
	/*----------------------------------------------------
	  By Yujin Chung on 2/28/2013
	  Although \lambda is used in this software for \beta in the Gibbs prior distribution, "beta" is printed out on monitor and in result files.
	  % fixit: double-check that  "beta" is printed out on monitor and in result files.
	-------------------------------------------------------- */
	acceptance_names[FREQlambda_gibbs]     = "beta";
	// acceptance_names[FREQlambda_gibbs]     = "lambda";

	/*------------------------------------
	  By Yujin Chungon 2/6/2013:
	  penalty w and related parame were removed.
	-------------------------------------- */
	// acceptance_names[FREQpenalty]    = "penalty";

	acceptance_names[FREQrate_gibbs]       = "rate";
	acceptance_names[FREQrateprior_gibbs]  = "rate_prior";
	acceptance_names[FREQkappa_gibbs]      = "kappa";
	acceptance_names[FREQkappaprior_gibbs] = "kappa_prior";
	acceptance_names[FREQswapchain_gibbs]  = "swap_chains";

	a_names[FREQadd_gibbs]        = " Tadd  ";
	a_names[FREQremove_gibbs]     = "Tremove";
	a_names[FREQshift_gibbs]      = "Tshift ";
	a_names[FREQcycle_gibbs]      = "Tcycle ";
	a_names[FREQheat_gibbs]       = " Theat ";
	a_names[FREQtopol_gibbs]      = "Tno.bkp";
	/*----------------------------------------------------
	  By Yujin Chung on 2/28/2013
	  Although \lambda is used in this software for \beta in the Gibbs prior distribution, "beta" is printed out on monitor and in result files.
	  % fixit: double-check that  "beta" is printed out on monitor and in result files.
	  -------------------------------------------------------- */
	a_names[FREQlambda_gibbs]     = "beta   ";
	// a_names[FREQlambda_gibbs]     = "lambda ";

	/*------------------------------------
	  By Yujin Chungon 2/6/2013:
	  penalty w and related parame were removed.
	-------------------------------------- */
	// a_names[FREQpenalty]    = "pnalty ";

	a_names[FREQrate_gibbs]       = " rate  ";
	a_names[FREQrateprior_gibbs]  = "rate.p ";
	a_names[FREQkappa_gibbs]      = " kappa ";
	a_names[FREQkappaprior_gibbs] = "kappa.p";
	a_names[FREQswapchain_gibbs]  = " swap  ";
}




