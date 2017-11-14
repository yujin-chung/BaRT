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

/*! \file likelihood.h 
 *
 *  \brief likelihood calculation functions 
 */

#ifndef _biomc2_likelihood_h
#define _biomc2_likelihood_h

#include "phylogeny.h"

/* Global function prototypes */

/*---------------------------------------
  By Yujin Chung on 2/6/2013:
  New functions were created because of
  new structure type phylogeny gibbs
  ----------------------------------------*/
void update_Q_eigenvalues_gibbs (phylogeny_gibbs p, double kappa);
void update_Q_matrix_gibbs (phylogeny_gibbs p, double lambda);

void     ln_likelihood_gibbs (phylogeny_gibbs *p, topology t);
void accept_likelihood_gibbs (phylogeny_gibbs *p, topology t);
void     ln_likelihood_moved_branches_gibbs (phylogeny_gibbs *p, topology t);
void accept_likelihood_moved_branches_gibbs (phylogeny_gibbs *p, topology t);
void     ln_likelihood_proposal_gibbs (phylogeny_gibbs p, topology t);
void accept_likelihood_proposal_gibbs (phylogeny_gibbs p, topology t);

/*--------------------------------------------
  By Yujin Chung on 2/6/2013:
  The following functions are from biomc2-1.9
------------------------------------------- */
void update_Q_eigenvalues (phylogeny p, double kappa);
void update_Q_matrix (phylogeny p, double lambda);

void     ln_likelihood (phylogeny *p, topology t);
void accept_likelihood (phylogeny *p, topology t);
void     ln_likelihood_moved_branches (phylogeny *p, topology t);
void accept_likelihood_moved_branches (phylogeny *p, topology t);
void     ln_likelihood_proposal (phylogeny p, topology t);
void accept_likelihood_proposal (phylogeny p, topology t);

#endif
