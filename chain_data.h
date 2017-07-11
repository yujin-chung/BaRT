/*! 
 * Copyright (C) 2006	Leonardo de Oliveira Martins
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

/*! \file chain_data.h 
 *
 *  \brief Header file for chain_data.c
 */

#ifndef _biomc2_chain_data_h
#define _biomc2_chain_data_h

#include "linked_topol.h"

#define AFREQSIZE 13

/*-----------
  By Yujin Chung on 2/6/2013:
  A new enum name "Acceptance_gibbs" has been created.
  It does not have "FREQpenalty" (see "Acceptance")
  -----------------------*/
/*! \brief Vector positions of acceptance rates (counts). */
typedef enum 
{
	FREQadd_gibbs = 0,    /*!< \brief Addition (birth) of break-point. */
	FREQremove_gibbs,     /*!< \brief Removal (death) of break-point. */
	FREQshift_gibbs,      /*!< \brief Break-point location shift. */
	FREQcycle_gibbs,      /*!< \brief Al-Awadhi update within birth/death move. */
	FREQheat_gibbs,       /*!< \brief Al-Awadhi-like update (heat/cold swap) inside FREQtopol move (does not change break-point structure). */
	FREQtopol_gibbs,      /*!< \brief Topology update without changing break-points location. */
	/*------------------------------------
	  By Yujin Chung on 2/20/2013
	  related with hyperparameter of the Poisson-like prior distribution on gene trees
	  
	  -- fixit: is it better to change it to "FREQbeta_gibbs"?
	            maybe not. beta is used for other things.
	---------------------------------------- */
	FREQlambda_gibbs,     /*!< \brief Lambda (\f$\lambda\f$) update. */
	/*------------------------
	  By Yujin Chung on 2/6/2013:
	  "FREQpenalty" is removed
	-----------------------*/
	// FREQpenalty,    /*!< \brief Penalty (\f$w\f$) update. */
	FREQrate_gibbs,       /*!< \brief Average rate per segment (\f$\mu_i\f$). */
	FREQrateprior_gibbs,  /*!< \brief Prior for rates (\f$\mu_0\f$). */
	FREQkappa_gibbs,      /*!< \brief Transition:transversion ratio (\f$\kappa_i\f$) per segment. */
	FREQkappaprior_gibbs, /*!< \brief Prior for ti/tv ratios (\f$\kappa_0\f$).  */
	FREQswapchain_gibbs   /*!< \brief Swap between cold/heated chains. */
} Acceptance_gibbs;


/*-----------
  By Yujin Chung on 2/6/2013:
  The original new enum name "Acceptance" in biomc2,
  which is replaced by "Acceptance_gibbs" 
  for the software with gibbs prior distribution.
  -----------------------*/
/*! \brief Vector positions of acceptance rates (counts). */
typedef enum 
{
	FREQadd = 0,    /*!< \brief Addition (birth) of break-point. */
	FREQremove,     /*!< \brief Removal (death) of break-point. */
	FREQshift,      /*!< \brief Break-point location shift. */
	FREQcycle,      /*!< \brief Al-Awadhi update within birth/death move. */
	FREQheat,       /*!< \brief Al-Awadhi-like update (heat/cold swap) inside FREQtopol move (does not change break-point structure). */
	FREQtopol,      /*!< \brief Topology update without changing break-points location. */
	FREQlambda,     /*!< \brief Lambda (\f$\lambda\f$) update. */
	FREQpenalty,    /*!< \brief Penalty (\f$w\f$) update. */
	FREQrate,       /*!< \brief Average rate per segment (\f$\mu_i\f$). */
	FREQrateprior,  /*!< \brief Prior for rates (\f$\mu_0\f$). */
	FREQkappa,      /*!< \brief Transition:transversion ratio (\f$\kappa_i\f$) per segment. */
	FREQkappaprior, /*!< \brief Prior for ti/tv ratios (\f$\kappa_0\f$).  */
	FREQswapchain   /*!< \brief Swap between cold/heated chains. */
} Acceptance;

/*! \brief Vectors with acceptance rates (accepted moves over tried moves) per chain.  */
typedef struct
{
  int propose[AFREQSIZE], /*!< \brief Proposal counter. */
    accept[AFREQSIZE];  /*!< \brief Acceptance counter. */
} acceptance_frequency;

/*! \brief Strings with names of acceptance moves (long version). */
char *acceptance_names[AFREQSIZE]; 
/*! \brief Strings with names of acceptance moves (short version, for headers). */
char *a_names[AFREQSIZE]; 

/*-------------------
   By Yujin Chung on 2/6/2013:
   A new structure name "chain_data_gibbs" has been created.
----------------------*/
typedef struct chain_data_gibbs* chain_data_gibbs;

/* -------------------
   By Yujin Chung on 2/6/2013:
   The following structure name "chain_data" was used in biomc2, 
   but is not used any more and replaced by "chain_data_gibbs" 
   for the new software with gibbs distributions 
----------------------*/
typedef struct chain_data_struct* chain_data;


/*-------------------
   By Yujin Chung on 2/4/2013:
   A new structure "chain_data_gibbs" has been created.
   This structure basically contains the same members in "chain_data_struct" 
   except for "hp_penalty_alpha", "hp_penalty_beta".
----------------------*/
/*! \brief MCMC chain data (topologies, likelihoods, sampler parameters, etc.). */
struct chain_data_gibbs
{
  /*! \brief Vector with gene segments. */ 
  phylogeny_gibbs *segment;
  /*! \brief Linked list of topologies. */
  linked_topol topol;
  /*! \brief Number of segments. */
  int n_segments;
  /*! \brief Pattern vector for each site, where each element describes the
   * pattern to which the site belongs. This is the only place where the actual
   * gene segments (and not the "chunked data" pattern) are stored. */
  int *site_pattern;
  /*! \brief Size of chain_data_struct::site_pattern vector (equal to NCHAR in nexus file). */
  int nsites;
  /*! \brief Number of taxa. */
  int ntax;
  /*! \brief Vector of break-point locations for every segment. */
  int *breakpoint;
  /*! \brief Sequence names, as in nexus_alignment_struct. */
  char **taxlabel;
  /*-------------------------
    By Yujin Chung on 2/20/2013
    Compare two topologies based on their splits, here represented as
    bit-strings (e.g. the pair 00110001, 11001110). 
    
    -- fixit: split_space_struc is defined in "topology.h" 
    and it seems to be related with SPR distances.
  --------------------------*/
  /*! \brief split_partition_struct information. */
  split_space split;
  
  /*! \brief Base (A,G,T,C) equilibrium frequencies. */
  double pi[6];
  /*! \brief Eigenvectors for HKY model. */
  double z1[4][4], z2[4][4];
  
  /*! \brief Prior \f$\mu_0\f$ for substitution rate (exponential). */
  double prior_rate;
  /*! \brief Prior \f$\kappa_0\f$ for ti/tv ratio (exponential). */
	double prior_kappa;
 /*------------------------------------
  By Yujin Chung on 2/20/2013
  parameter of the prior distribution on gene trees
  In our papers, this should be "beta" from Gibbs prior distributions
  "beta" scales with the inverse of the average recombination rate per segment.

  By Yujin Chung on 3/4/2013
  In biomcs2, chain->segment[i] has different lambda as member.
  In BaRT, chain contains only one lambda as member.
  ---------------------------------------- */
  /*! \brief Prior \f$\lambda_i\f$ for topology distances (modified truncated Poisson). */
  double lambda;
  /*------------------------------------
    By Yujin Chung on 2/20/2013
    lnN contains the partition function 
    of the current truncated Poisson distribution on gene trees
    lnN_proposal contains the partition function 
    of the proposed truncated Poisson distribution on gene trees

    By Yujin Chung on 3/4/2013
    In biomc2, chain has a vector of partition functions (product of truncated Poisson dist)
    In BaRT, we compute the partition function of the full prior Gibbs distribution on gene trees
---------------------------------------- */
  /*! \brief Poisson truncation term. */
  double lnN;
  // double *lnN;
  /*! \brief Poisson truncation term in proposal step. */
  double lnN_proposal;
  // double *lnN_proposal;

  /*! \brief Temperature of the chain (namelly \f$1/kT\f$) */
	double kT;

  /* MCMC parameters (from control file) */

  /*! \brief File name with tree topologies. */
  char tree_filename[512];
  /*! \brief File name with sequence data. */
  char seq_filename[512];

   /*------------------------------------
    By Yujin Chung on 2/20/2013
    --fixit: with the new prior distribution on gene trees,
    we don't need "max_distance" and "max_dist_split"
---------------------------------------- */
  /*! \brief Maximum distance between topologies in the modified Poisson distribution of Bayesian model. */
  int max_distance;
  /*! \brief Maximum distance to be calculated by the split method, less or
   *  equal than chain_data_struct::max_distance. */
  int max_dist_split;
	/*! \brief Number of iterations of the main MCMC sampler. */
  int ngen;
  /*! \brief Number of samples from the posterior distribution. */
  int nsamples;
  /*! \brief Number of iterations of the burn-in stage. */
  int burnin;
  /*! \brief Number of iterations of the warm-up (to sample initial values). */
	int warmup;
	/*! \brief Simulation of prior distribution ( 1 => prior, 0 => posterior). */
	bool prior;
  /*! \brief Hyperprior \f$\mathcal{M}\f$ for substitution rate (exponential). */
  double hp_rate;
  /*! \brief Hyperprior \f$\mathcal{K}\f$ for ti/tv ratio (exponential). */
  double hp_kappa;
  /*------------------------------------
    By Yujin Chung on 2/20/2013
    hyperparameters of prior distribution on gene trees
    --fixit: need to change their names?
    ---------------------------------------- */
  /*! \brief Hyperprior \f$\alpha_\lambda\f$ for topology distances (assuming gamma). */
  double hp_lambda_alpha;
  /*! \brief Hyperprior \f$\beta_\lambda\f$ for topology distances (assuming gamma). */
  double hp_lambda_beta;

  /*-----------------------
    By Yujin Chung on 2/4/2013
     "hp_penalty_alpha" and "hp_penalty_beta" are removed.
  --------------------------- */
  /*! \brief Hyperprior \f$\alpha_w\f$ for penalty (assuming gamma). */
  // double hp_penalty_alpha;
  /*! \brief Hyperprior \f$\beta_w\f$ for penalty (assuming gamma) */
  // double hp_penalty_beta;

  /*! \brief Perturbation value \f$\xi_\mu\f$ for substitution rate. */	
	double update_rate;
  /*! \brief Perturbation value \f$\xi_\kappa\f$ for ti/tv ratio. */
	double update_kappa;
 /*------------------------------------
    By Yujin Chung on 2/20/2013
    related to the prior distribution on gene trees
    --fixit: need to change its name?
    ---------------------------------------- */
  /*! \brief Perturbation value \f$\xi_\lambda\f$ for prior topology distance \f$\lambda\f$.  */
  double update_lambda;

  /*-------------------------------------
    By Yujin Chung on 2/4/2013
    "update_penalty_alpha" is removed.
  ----------------------------------------- */
  /*! \brief Perturbation value \f$\xi_w\f$ for prior topology distance \f$w\f$. */
  //double update_penalty;
  
  double beta_zero,  /*! \brief initial temperature (\f$\beta_0\f$) in \f$\beta_n=\beta_0\log(n+e)\f$ 
											 for the warm-up annealing stage. */
         beta_final,  /*! \brief final temperature in \f$\beta_n=\beta_0\log(n+e)\f$ for the warm-up annealing stage. */
         beta_n;     /*! \brief current temperature (\f$\beta_n\f$) in \f$\beta_n=\beta_0\log(n+e)\f$ 
                       for the warm-up annealing stage. */
  /*! \brief Interval between chain swap updates (MC-MCMC). */
  int swap_interval;
  /*! \brief Temperature of heated chain (MC-MCMC). */
  double kT_heat;

  /*! \brief Number of cycles in Al-Awadhi update (when changing number of breakpoints). */
	int n_cycles;
  /*! \brief Temperature \f$1/kT\f$ of Al-Awadhi update (when changing number of breakpoints). */
  double kT_cycles;
  /*! \brief Number of cycles in Al-Awadhi-like update (heated minisampler). */
  int n_mini;
  /*! \brief Temperature \f$1/kT\f$ of Al-Awadhi-like update (heated minisampler). */ 
  double kT_mini;
  /*! \brief Factorial (\f$\log(n!)\f$) from zero to chain_data_struct::ntax */
  double *factorial;
  FILE *treefile, /*! \brief Output topology file for all segments (whenever topology changes) and samples. */
       *distfile; /*! \brief Output parameters file, one line per sample, inspired by dualbrothers. */
 /*-------------------------------------
    By Yujin Chung on 2/20/2013
    Related to the prior distribution on gene trees
    --fixit: need to change its name?
  ----------------------------------------- */
  int ratio_spr,          /*! \brief Iterations at which we update topologies by SPR instead of NNI. */	
      ratio_topol_update; /*! \brief Number of topology updates per iteration, to reduce correlation between samples. */	
  double log2, /*! \brief \f$\log(2)\f$, used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$. */
/*-------------------------------------
    By Yujin Chung on 2/20/2013
    Related to proposing and updating trees
    -- fixit: double-check it if it isn't related to prior distribution on gene trees
  ----------------------------------------- */
         logDspr, /*! \brief SPR-neighborhood size \f$\log(2(n-3)(2n-7))\f$, 
										used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$. */
         logDnni, /*! \brief NNI-neighborhood size \f$\log(2(n-3))\f$, 
										used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$. */
				 *logD;   /*! \brief pointer to chain_data::logDspr or chain_data::logDnni, according to the operation */
  /*! \brief Pointer to likelihood_moved_branches() (posterior sampling) or
   * likelihood_moved_branches_nodata() (prior sampling) */
  void (*ln_likelihood_moved_branches_gibbs) (phylogeny_gibbs *p, topology t);
  /*! \brief Pointer to likelihood_proposal() (posterior sampling) or
   * likelihood_proposal_nodata() (prior sampling) */
  void (*ln_likelihood_proposal_gibbs) (phylogeny_gibbs p, topology t);
	/*! \brief Pointer to acceptance probability calculation: standard Bayesian
	 * procedure or new simulated annealing */
  bool (*acceptance_probability_gibbs) (chain_data_gibbs chain, double *acceptance, double *transition);
};



/*------------------------------------
  By Yujin Chung on 2/20/2013
  The following structure is from biomc2-1.9
---------------------------------------- */
/*! \brief MCMC chain data (topologies, likelihoods, sampler parameters, etc.). */
struct chain_data_struct
{
  /*! \brief Vector with gene segments. */ 
  phylogeny *segment;
  /*! \brief Linked list of topologies. */
  linked_topol topol;
  /*! \brief Number of segments. */
  int n_segments;
  /*! \brief Pattern vector for each site, where each element describes the
   * pattern to which the site belongs. This is the only place where the actual
   * gene segments (and not the "chunked data" pattern) are stored. */
  int *site_pattern;
  /*! \brief Size of chain_data_struct::site_pattern vector (equal to NCHAR in nexus file). */
  int nsites;
  /*! \brief Number of taxa. */
  int ntax;
  /*! \brief Vector of break-point locations for every segment. */
  int *breakpoint;
  /*! \brief Sequence names, as in nexus_alignment_struct. */
  char **taxlabel;
	/*! \brief split_partition_struct information. */
	split_space split;

	/*! \brief Base (A,G,T,C) equilibrium frequencies. */
	double pi[6];
	/*! \brief Eigenvectors for HKY model. */
	double z1[4][4], z2[4][4];

	/*! \brief Prior \f$\mu_0\f$ for substitution rate (exponential). */
	double prior_rate;
  /*! \brief Prior \f$\kappa_0\f$ for ti/tv ratio (exponential). */
	double prior_kappa;
 
 
  /*! \brief Poisson truncation term. */
   double *lnN;
  /*! \brief Poisson truncation term in proposal step. */
   double *lnN_proposal;

  /*! \brief Temperature of the chain (namelly \f$1/kT\f$) */
	double kT;

  /* MCMC parameters (from control file) */

  /*! \brief File name with tree topologies. */
  char tree_filename[512];
  /*! \brief File name with sequence data. */
  char seq_filename[512];

  /*! \brief Maximum distance between topologies in the modified Poisson distribution of Bayesian model. */
  int max_distance;
  /*! \brief Maximum distance to be calculated by the split method, less or
   *  equal than chain_data_struct::max_distance. */
  int max_dist_split;
	/*! \brief Number of iterations of the main MCMC sampler. */
  int ngen;
  /*! \brief Number of samples from the posterior distribution. */
  int nsamples;
  /*! \brief Number of iterations of the burn-in stage. */
  int burnin;
  /*! \brief Number of iterations of the warm-up (to sample initial values). */
	int warmup;
	/*! \brief Simulation of prior distribution ( 1 => prior, 0 => posterior). */
	bool prior;
  /*! \brief Hyperprior \f$\mathcal{M}\f$ for substitution rate (exponential). */
  double hp_rate;
  /*! \brief Hyperprior \f$\mathcal{K}\f$ for ti/tv ratio (exponential). */
  double hp_kappa;
  /*! \brief Hyperprior \f$\alpha_\lambda\f$ for topology distances (assuming gamma). */
  double hp_lambda_alpha;
  /*! \brief Hyperprior \f$\beta_\lambda\f$ for topology distances (assuming gamma). */
  double hp_lambda_beta;
  /*! \brief Hyperprior \f$\alpha_w\f$ for penalty (assuming gamma). */
  double hp_penalty_alpha;
  /*! \brief Hyperprior \f$\beta_w\f$ for penalty (assuming gamma) */
  double hp_penalty_beta;
  /*! \brief Perturbation value \f$\xi_\mu\f$ for substitution rate. */	
	double update_rate;
  /*! \brief Perturbation value \f$\xi_\kappa\f$ for ti/tv ratio. */
	double update_kappa;
  /*! \brief Perturbation value \f$\xi_\lambda\f$ for prior topology distance \f$\lambda\f$.  */
  double update_lambda;
  /*! \brief Perturbation value \f$\xi_w\f$ for prior topology distance \f$w\f$. */
  double update_penalty;
  double beta_zero,  /*! \brief initial temperature (\f$\beta_0\f$) in \f$\beta_n=\beta_0\log(n+e)\f$ 
											 for the warm-up annealing stage. */
         beta_final,  /*! \brief final temperature in \f$\beta_n=\beta_0\log(n+e)\f$ for the warm-up annealing stage. */
         beta_n;     /*! \brief current temperature (\f$\beta_n\f$) in \f$\beta_n=\beta_0\log(n+e)\f$ 
                       for the warm-up annealing stage. */
  /*! \brief Interval between chain swap updates (MC-MCMC). */
  int swap_interval;
  /*! \brief Temperature of heated chain (MC-MCMC). */
  double kT_heat;

  /*! \brief Number of cycles in Al-Awadhi update (when changing number of breakpoints). */
	int n_cycles;
  /*! \brief Temperature \f$1/kT\f$ of Al-Awadhi update (when changing number of breakpoints). */
  double kT_cycles;
  /*! \brief Number of cycles in Al-Awadhi-like update (heated minisampler). */
  int n_mini;
  /*! \brief Temperature \f$1/kT\f$ of Al-Awadhi-like update (heated minisampler). */ 
  double kT_mini;
  /*! \brief Factorial (\f$\log(n!)\f$) from zero to chain_data_struct::ntax */
  double *factorial;
  FILE *treefile, /*! \brief Output topology file for all segments (whenever topology changes) and samples. */
       *distfile; /*! \brief Output parameters file, one line per sample, inspired by dualbrothers. */
  int ratio_spr,          /*! \brief Iterations at which we update topologies by SPR instead of NNI. */	
      ratio_topol_update; /*! \brief Number of topology updates per iteration, to reduce correlation between samples. */	
  double log2, /*! \brief \f$\log(2)\f$, used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$. */
         logDspr, /*! \brief SPR-neighborhood size \f$\log(2(n-3)(2n-7))\f$, 
										used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$. */
         logDnni, /*! \brief NNI-neighborhood size \f$\log(2(n-3))\f$, 
										used in topology transition term \f$q(x'\mid x)/q(x\mid x')\f$. */
				 *logD;   /*! \brief pointer to chain_data::logDspr or chain_data::logDnni, according to the operation */
  /*! \brief Pointer to likelihood_moved_branches() (posterior sampling) or
   * likelihood_moved_branches_nodata() (prior sampling) */
  void (*ln_likelihood_moved_branches) (phylogeny *p, topology t);
  /*! \brief Pointer to likelihood_proposal() (posterior sampling) or
   * likelihood_proposal_nodata() (prior sampling) */
	void (*ln_likelihood_proposal) (phylogeny p, topology t);
	/*! \brief Pointer to acceptance probability calculation: standard Bayesian
	 * procedure or new simulated annealing */
  bool (*acceptance_probability) (chain_data chain, double *acceptance, double *transition);
};

/* Global function prototypes */

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  The following functions were created 
  because of new structure "chain_data_gibbs".
------------------------------------ */
chain_data_gibbs new_chain_data_from_file_gibbs(char *filename);
void del_chain_data_gibbs (chain_data_gibbs p);
/*! \brief Bayesian acceptance probability. */
bool acceptance_probability_bayesian_gibbs (chain_data_gibbs chain, double *acceptance, double *transition);
/*! \brief Simulated annealing acceptance probability. */
bool acceptance_probability_annealing_gibbs (chain_data_gibbs chain, double *acceptance, double *transition);

/*------------------------------------
  By Yujin Chung on 2/6/2013:
  The following functions are from biomc2
------------------------------------ */
chain_data new_chain_data_from_file (char *filename);
void del_chain_data (chain_data p);
/*! \brief Bayesian acceptance probability. */
bool acceptance_probability_bayesian (chain_data chain, double *acceptance, double *transition);
/*! \brief Simulated annealing acceptance probability. */
bool acceptance_probability_annealing (chain_data chain, double *acceptance, double *transition);
void setup_acceptance_names (void);
void setup_acceptance_names_gibbs (void); // YC 5/24/2017
#endif
