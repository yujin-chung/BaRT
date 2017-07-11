/* --------------------------------------
 Modified on June 15th, 2012
 -----------------------------------*/

#include <math.h>
#include <stdio.h>
#include "misc.h"
#include "shape.h"
#include <stdlib.h>

// #include <stdlib.h>  // need it for malloc()


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491
#endif

/*
  Return the number of all possible trees
  input: the number of taxa, x
  output: (2*x-5)!!                       */
long double treeSpaceSize(int ntaxa)
{
  long double res = doublefactorial(2*ntaxa-5);
  //printf("treeSpaceSize:: res = %Le\n",res);
  return res;
}

long double log_treeSpaceSize(int ntaxa)
{
  long double res = log_doublefactorial(2*ntaxa-5);
  //printf("treeSpaceSize:: res = %Le\n",res);
  return res;
}

/*
  E[d(T1,T2)]
 */
long double expected_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu=0.0, cond_E;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_E=0.0; 
      ptrees_givenShape=0.0;
      for(count_dist=0;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_E += count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      mu += cond_E * ptrees_givenShape;
    }
  //printf("mu = %Le\n",mu);
  return mu;
}

/*
  E[d(T1,T2)|d(T1,T2)>=2]
 */
long double expected_dist2(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu=0.0, cond_E;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_E=0.0; 
      ptrees_givenShape=0.0;
      for(count_dist=2;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_E += count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      mu += cond_E * ptrees_givenShape;
    }
  //printf("mu = %Le\n",mu);
  return mu;
}

long double var_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu = expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double cond_Esquared;
  long double E_squareDist = 0.0, var;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_Esquared = 0.0;
      ptrees_givenShape = 0.0;
      for(count_dist=0;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_Esquared +=  count_dist*count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      E_squareDist += cond_Esquared * ptrees_givenShape;
    }
  var = E_squareDist - mu*mu;

  //printf("var_d = %Le\n",var);
  return var;
}

void condi_mu(int ntaxa, int nshape,long double*** propDistance_twoTrees, long double conditional_mean[])
{
  int maxdist_twoTrees = ntaxa-3;
  int count_shape1, count_shape2, count_dist;
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      conditional_mean[count_shape1] = 0;
      for(count_shape2=0;count_shape2<nshape;count_shape2++)
	for(count_dist=0;count_dist<=maxdist_twoTrees;count_dist++)
	  conditional_mean[count_shape1] += count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
    }  
}

long double cov_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu = expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double conditional_mean[nshape];
  condi_mu(ntaxa,nshape,propDistance_twoTrees,conditional_mean);
  long double cond_E;
  long double E_squaredE = 0.0,cov;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_E=0.0;
      ptrees_givenShape=0.0;
      for(count_dist=0;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_E += count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist] * conditional_mean[count_shape2];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      E_squaredE += cond_E * ptrees_givenShape;
    }
  cov = E_squaredE - mu*mu;
  //printf("E_squareE = %Le\t mu*mu=%Le\t cov=%Le\n",E_squaredE,mu*mu,cov);
  return cov;
}

long double corr_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  long double cov = cov_dist(ntaxa,nshape,propDistance_twoTrees);
  long double var = var_dist(ntaxa,nshape,propDistance_twoTrees);
  long double corr = (long double) cov/var;
  return corr;
}

long double cov_dist_original(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu = expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double cond_E;
  long double E_squaredE = 0.0,cov;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_E=0.0;
      ptrees_givenShape=0.0;
      for(count_dist=0;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_E += count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      E_squaredE += cond_E*cond_E * ptrees_givenShape;
    }
  cov = E_squaredE - mu*mu;
  //printf("E_squareE = %Le\t mu*mu=%Le\t cov=%Le\n",E_squaredE,mu*mu,cov);
  return cov;
}

long double var_dist2(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu=expected_dist2(ntaxa, nshape,propDistance_twoTrees);
  long double cond_Esquared;
  long double E_squareDist=0.0, var;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_Esquared = 0.0;
      ptrees_givenShape=0.0;
      for(count_dist=2;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_Esquared +=  count_dist*count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      E_squareDist += cond_Esquared * ptrees_givenShape;
    }
  var = E_squareDist - mu*mu;

  //printf("var_d = %Le\n",var);
  return var;
}

long double cov_dist2(int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double ptrees_givenShape; // the number of trees with the given tree shape
  int count_shape1, count_shape2, count_dist;
  long double mu=expected_dist2(ntaxa, nshape,propDistance_twoTrees);
  long double cond_E;
  long double E_squaredE=0.0,cov;
    // mu = 2.0/3; var = 2.0/9; cov=0.0; // ntaxa = 4
    // mu = 8.0/5, var = 28.0/75, cov=0.0; // ntaxa = 5
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    {
      cond_E=0.0;
      ptrees_givenShape=0.0;
      for(count_dist=2;count_dist<=maxdist_twoTrees;count_dist++)
	{
	  ptrees_givenShape += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
	  //printf("propDistance_twoTrees[count_shape1][count_shape1][count_dist]=%Le\t  ptrees_givenShape=%Le\n",propDistance_twoTrees[count_shape1][count_shape1][count_dist],ptrees_givenShape);
	  for(count_shape2=0;count_shape2<nshape;count_shape2++)
	    {
	      cond_E += count_dist * propDistance_twoTrees[count_shape1][count_shape2][count_dist];
	    }
	}
      //printf("count_shape1 = %d\t cond_E=%Le\t ptrees_givenShape=%Le\n",count_shape1,cond_E,ptrees_givenShape);
      E_squaredE += cond_E*cond_E * ptrees_givenShape;
    }
  cov = E_squaredE - mu*mu;
  //printf("E_squareE = %Le\t mu*mu=%Le\t cov=%Le\n",E_squaredE,mu*mu,cov);
  return cov;
}

/*
  mu_L = E[D_1+...+D_{L-1}] = (L-1)E[D]
 */
long double mu_L(int nsites, int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  long double mu = expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double mu_L = (long double) (nsites-1)*mu;
  //printf("mu_L2 = %Le\n",mu_L);
  return mu_L;
}

/*
  var_L = Var[D_1+...+D_{L-1}] = (L-1)Var[D] + L(L-1)Cov[D_1,D_2]
 */
long double var_L(int nsites, int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  long double var_L = (long double) (nsites-1)*var_dist(ntaxa, nshape,propDistance_twoTrees)+nsites*(nsites-1)*cov_dist(ntaxa, nshape,propDistance_twoTrees);
  //printf("var_L2 = %Le\n",var_L);
  return var_L;
}

long double mu_L2(int nsites, int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  long double Z1 = treeSpaceSize(ntaxa);
  long double mu = expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double mu_L = (long double) (nsites-1)*(mu - 2*Z1*(ntaxa-3)/pow(Z1,nsites));
  //printf("mu_L2 = %Le\n",mu_L);
  return mu_L;
}

long double var_L2(int nsites, int ntaxa, int nshape,long double*** propDistance_twoTrees)
{
  long double Z1 = treeSpaceSize(ntaxa);
  long double mu = expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double sigma2 = var_dist(ntaxa, nshape,propDistance_twoTrees)+cov_dist(ntaxa, nshape,propDistance_twoTrees);
  long double p1 =  (long double) 2*Z1*(ntaxa-3)*(nsites-1)/pow(Z1,nsites);
  long double var_L = (long double) (nsites-1)*sigma2 + p1*(2*mu-1-p1);
  //printf("var_L2 = %Le\n",var_L);
  return var_L;
}

////////////////////////////////////
//
//  Approximations
//
/////////////////////////////////////

long double log_normalizingConstant_Normalapprox(int nsites,int ntaxa, int nshape, long double b, long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double maxdist = (long double)(nsites-1)*maxdist_twoTrees; 
  long double Z1 = treeSpaceSize(ntaxa), logZ1 = (long double) log(Z1);

  int count_shape1, count_shape2, count_dist;

  long double ptrees_givenShape[nshape]; // the number of trees with the given tree shape

   
  /*
  long double mu=expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double var = var_dist(ntaxa, nshape,propDistance_twoTrees);
  long double cov = cov_dist(ntaxa, nshape,propDistance_twoTrees);
  */
  long double mu=expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double var = var_dist(ntaxa, nshape,propDistance_twoTrees);
  long double cov = cov_dist(ntaxa, nshape,propDistance_twoTrees);
  long double mu_dist = (nsites-1)*mu, var_dist = (nsites-1)*(var+cov);
  long double mu_L = mu_L2(nsites,ntaxa, nshape,propDistance_twoTrees);
  long double var_L = var_L2(nsites,ntaxa, nshape,propDistance_twoTrees);
  long double x1 = (2-0.5-mu_dist+b*var_dist)/(long double)sqrt(2*var_dist), x2= (maxdist+0.5-mu_dist+b*var_dist)/(long double)sqrt(2*var_dist);
  //printf("x1=%Le\t x2=%Le\n",x1,x2);
  //printf("maxdist=%d\t mu_dist=%f\t var_dist=%f\t z1=%f\t z2=%f\n",maxdist,mu_dist,var_dist,z1,z2);
  //double prob1 = (1-(z1/sqrt(pow(z1,2))*sqrt(1-exp(-(pow(z1,2)/2)*(4/M_PI+a*pow(z1,2)/2)/(1+a*pow(z1,2)/2)))))/2;
  //double prob2 = (1-(z2/sqrt(pow(z2,2))*sqrt(1-exp(-(pow(z2,2)/2)*(4/M_PI+a*pow(z2,2)/2)/(1+a*pow(z2,2)/2)))))/2;
  //double prob1 = (1-sqrt(2)*x1/(sqrt(M_PI)*exp(-x1*x1)+sqrt(M_PI*exp(-x1*x1)+8*x1*x1)))/2;
  //double prob2 = (1-sqrt(2)*x2/(sqrt(M_PI)*exp(-x2*x2)+sqrt(M_PI*exp(-x2*x2)+8*x2*x2)))/2;
  long double prob1 =erfl(x1), prob2 = erfl(x2); // erfcl() in math.h
  //printf("prob1=%Le\t prob2 =%Le\n",prob1,prob2);

  /* when nshape =1
  double z_apprx= Z1+Z1*propDistance_twoTrees[1]*Z1*(nsites-1)*exp(-b);
  if(prob1!=0.0)
      z_apprx+=(pow(Z1,nsites)-Z1-Z1*propDistance_twoTrees[1]*Z1*(nsites-1))*prob1*exp(-b*mu_dist+var_dist*(b*b)/2); 
  //z_apprx+=(pow(Z1,nsites)-Z1-Z1*propDistance_twoTrees[1]*Z1*(nsites-1))*(prob1-prob2)*exp(-b*mu_dist+var_dist*(b*b)/2); 
  */
  for(count_shape1=0;count_shape1<nshape;count_shape1++){
    ptrees_givenShape[count_shape1]=0.0;
    for(count_dist = 0; count_dist <= maxdist_twoTrees; count_dist++){
      ptrees_givenShape[count_shape1] += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
      //printf("count_shape=%d\tcount_dist=%d\tptrees_givenShape[count_shape]=%f\n",count_shape,count_dist,ptrees_givenShape[count_shape]);
      }
  }   
  long double logz_apprx, term1,term2=0.0, term3,term3_0,term3_1,term3_2,g1,g2;
  long double multiplier = -nsites*logZ1;
  /*
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    for(count_shape2=0;count_shape2<nshape;count_shape2++)
      term2 += propDistance_twoTrees[count_shape1][count_shape2][1]*ptrees_givenShape[count_shape1]*Z1;
  term2 *= (nsites-1);
  */
  term1 = logZ1 + multiplier;
  term2 = (long double) (log(2) + logZ1 + log(ntaxa-3) + log(nsites-1) - b + multiplier);
  //term3 = (long double) exp(log(pow(Z1,nsites) - Z1-term2)+log(prob2-prob1)-b*mu_dist+b*b*var_dist/2 -log(2));
  term3_0 = (long double) log(1 - exp(logZ1-nsites*logZ1) - exp(log(2)+logZ1+log(ntaxa-3)+log(nsites-1)-nsites*logZ1));
  term3_1 = (long double) (log(prob2-prob1) - log(2));
  term3_2 = -b*mu_dist+b*b*var_dist/2;
  term3 = term3_0 + term3_1 + term3_2;
  g1 = (long double) term1+log(1+exp(term2-term1));
  g2 = term3;
  if(g2-g1<0)
    logz_apprx = (long double) g1+log(1+exp(g2-g1)) - multiplier;
  else
    logz_apprx = (long double) g2+log(1+exp(g1-g2)) - multiplier;
  //printf("term1=%Le term2=%Le term3=%Le (term3_0=%Le term3_1=%Le term3_2=%Le)\n",term1,term2,term3,term3_0,term3_1,term3_2);
  //printf("g1=%Le g2=%Le g2-g1=%Le logz_approx=%Le\n",g1,g2,g2-g1,logz_apprx);
  //printf("in term3:: ((long double) pow(Z1,nsites) - Z1-term2)=%Le\t (prob2-prob1)=%Le\n",((long double) pow(Z1,nsites) - Z1-term2),(prob2-prob1));
   /*
   if(prob1!=0.0){
     term3=exp( nsites*log(Z1) - b*mu_dist+var_dist*(b*b)/2 + log(prob1-prob2) ) - (Z1+term2)*exp(-b*mu_dist+var_dist*(b*b)/2)*(prob1-prob2); 
     z_apprx += term3;
     //z_apprx+=(pow(Z1,nsites)-Z1-term2)*prob1*exp(-b*mu_dist+var_dist*(b*b)/2); 
   }
   */

  return logz_apprx;
}



/* Return the log density of normal distribution at x
 */
long double log_pdfNormalDist(long double x, long double mean, long double var){
  return (long double) -log(2*M_PI*var)/2-pow(x-mean,2)/(2*var);
}

long double log_normalizingConstant_Normalapprox2(int nsites,int ntaxa, int nshape, long double b, long double*** propDistance_twoTrees)
{
  int maxdist_twoTrees = ntaxa-3;  // the maximum (RF-distance/2) between two trees, (RF distance: 0,1, 2, ..., (ntaxa-3))
  long double maxdist = (long double)(nsites-1)*maxdist_twoTrees; 
  long double Z1 = treeSpaceSize(ntaxa), logZ1 = (long double) log(Z1);

  int count_shape1, count_shape2, count_dist;

  long double ptrees_givenShape[nshape]; // the number of trees with the given tree shape

   
  /*
  long double mu=expected_dist(ntaxa, nshape,propDistance_twoTrees);
  long double var = var_dist(ntaxa, nshape,propDistance_twoTrees);
  long double cov = cov_dist(ntaxa, nshape,propDistance_twoTrees);
  */
  long double mu=expected_dist2(ntaxa, nshape,propDistance_twoTrees);
  long double var = var_dist2(ntaxa, nshape,propDistance_twoTrees);
  long double cov = cov_dist2(ntaxa, nshape,propDistance_twoTrees);
  long double mu_dist = (nsites-1)*mu, var_dist = (nsites-1)*(var+cov);
  long double mu_L = mu_L2(nsites,ntaxa, nshape,propDistance_twoTrees);
  long double var_L = var_L2(nsites,ntaxa, nshape,propDistance_twoTrees);
  long double x1 = (2-0.5-mu_dist+b*var_dist)/(long double)sqrt(2*var_dist), x2= (maxdist+0.5-mu_dist+b*var_dist)/(long double)sqrt(2*var_dist);
 
  long double prob1 =erfl(x1), prob2 = erfl(x2); // erfcl() in math.h
  //printf("prob1=%Le\t prob2 =%Le\n",prob1,prob2);


  for(count_shape1=0;count_shape1<nshape;count_shape1++){
    ptrees_givenShape[count_shape1]=0.0;
    for(count_dist = 0; count_dist <= maxdist_twoTrees; count_dist++){
      ptrees_givenShape[count_shape1] += propDistance_twoTrees[count_shape1][count_shape1][count_dist];
      //printf("count_shape=%d\tcount_dist=%d\tptrees_givenShape[count_shape]=%f\n",count_shape,count_dist,ptrees_givenShape[count_shape]);
      }
  }   
  long double logz_apprx, term1,term2=0.0, term3,term3_0,term3_1,term3_2,g1,g2,g3;
  long double multiplier = -nsites*logZ1;
  /*
  for(count_shape1=0;count_shape1<nshape;count_shape1++)
    for(count_shape2=0;count_shape2<nshape;count_shape2++)
      term2 += propDistance_twoTrees[count_shape1][count_shape2][1]*ptrees_givenShape[count_shape1]*Z1;
  term2 *= (nsites-1);
  */
  term1 = logZ1 + multiplier;
  term2 = (long double) (log(2) + logZ1 + log(maxdist_twoTrees) + log(nsites-1) - b + multiplier);
  //term3 = (long double) exp(log(pow(Z1,nsites) - Z1-term2)+log(prob2-prob1)-b*mu_dist+b*b*var_dist/2 -log(2));
  term3_0 = (long double) log(1 - exp(logZ1+multiplier) - exp(log(2)+logZ1+log(ntaxa-3)+log(nsites-1)+multiplier));
  
  long double x;

  if(maxdist>2)
    {
      term3_1 = (long double) log(0);
      for(x=2;x<=maxdist;x++)
	{
	  g1 = term3_1;
	  g2 = (long double) -b*x+ log_pdfNormalDist(b*x, b*mu_dist, b*b*var_dist);//+log(var_dist)/2;
	  
	  /*
	  g2 = (long double) log_pdfNormalDist(x+1/exp(b*b*x), mu_dist, var_dist);
	  g3 = (long double) log_pdfNormalDist(x-1/exp(b*b*x), mu_dist, var_dist);
	  printf("g2=%Le g3=%Le g2-g3=%Le\n",g2,g3,g2-g3);
	  if(g2-g3>0)
	    g2 =(long double) -b*x + g2+log(1- exp(g3-g2)) -b*b*x;
	  else
	    g2 =(long double) -b*x + g3+log(1- exp(g2-g3)) - b*b*x;
	  printf("g1=%Le g2=%Le\n",g1,g2);
	  */
	  
	  //g2 = (long double) -b*x+ (log(erfl(x+0.00000001-mu_dist/sqrt(var_dist))- erfl(x-0.00000001-mu_dist/sqrt(var_dist)))-log(2));
	    //log_pdfNormalDist(x, mu_dist, var_dist);
	  if(exp(g1)+exp(g2)!=0)
	    {
	      if(g2-g1<0)
		term3_1 = (long double) g1+log(1+exp(g2-g1));
	      else
		term3_1 = (long double) g2+log(1+exp(g1-g2));
	    }	  
	  //printf("x=%Le -b*x=%Le g1=%Le g2=%Le term3_1=%Le\n",x,-b*x,g1,g2,term3_1);
	}
      term3 = term3_0 + term3_1;
    }
  else    
    term3 = (long double) log(0);
  g1 = (long double) term1+log(1+exp(term2-term1));
  g2 = term3;
  if(g2-g1<0)
    logz_apprx = (long double) g1+log(1+exp(g2-g1)) - multiplier;
  else
    logz_apprx = (long double) g2+log(1+exp(g1-g2)) - multiplier;
  // printf("term1=%Le term2=%Le term3=%Le (term3_0=%Le term3_1=%Le )\n",term1,term2,term3,term3_0,term3_1);
  //printf("g1=%Le g2=%Le g2-g1=%Le multiplier = %Le logz_approx=%Le\n",g1,g2,g2-g1,multiplier,logz_apprx);
  //printf("in term3:: ((long double) pow(Z1,nsites) - Z1-term2)=%Le\t (prob2-prob1)=%Le\n",((long double) pow(Z1,nsites) - Z1-term2),(prob2-prob1));
   /*
   if(prob1!=0.0){
     term3=exp( nsites*log(Z1) - b*mu_dist+var_dist*(b*b)/2 + log(prob1-prob2) ) - (Z1+term2)*exp(-b*mu_dist+var_dist*(b*b)/2)*(prob1-prob2); 
     z_apprx += term3;
     //z_apprx+=(pow(Z1,nsites)-Z1-term2)*prob1*exp(-b*mu_dist+var_dist*(b*b)/2); 
   }
   */

  return logz_apprx;
}

/*------------------------------
  Independence approximation
-------------------------------- */
// %fixit: need to revise this function, so that it uses the marginal distribution of RF distant given a tree shape, instead of the joint distribution.
// 
// argument: nsites represents the number of segments if the length of a segment is larger than 1.
long double log_normalizingConstant_indep(long double*** propDistance_twoTrees, int ntaxa, unsigned int nsites, int nshape, long double b)
{
  long double z=0.0, logz;
  long double Z1 = treeSpaceSize(ntaxa);
  int arrSize = (ntaxa-3)+1;

  int d, count_shape1, count_shape2;
  /* when nshape=1
  for(d=0;d<arrSize;d++){
    z+=propDistance_twoTrees[d]*exp(-b*d);
  }
  */
  double maxdist_twoTrees = ntaxa-3;
  long double ptrees_givenShape[nshape];
  unsigned int count_shape, count_dist;
  for(count_shape=0;count_shape<nshape;count_shape++){
    ptrees_givenShape[count_shape]=0.0;
    for(count_dist = 0; count_dist <= maxdist_twoTrees; count_dist++){
      ptrees_givenShape[count_shape] += propDistance_twoTrees[count_shape][count_shape][count_dist];
      //printf("count_shape=%d\tcount_dist=%d\tptrees_givenShape[count_shape]=%f\n",count_shape,count_dist,ptrees_givenShape[count_shape]);
    }
  }

  for(d=0;d<arrSize;d++)
    for(count_shape1=0;count_shape1<nshape;count_shape1++)
      for(count_shape2=0;count_shape2<nshape;count_shape2++)
	z+= propDistance_twoTrees[count_shape1][count_shape2][d] * exp(-b*d) * ptrees_givenShape[count_shape1];   
  
  logz=(nsites-1)*(long double)log(z)+nsites*(long double)log(Z1);
  return logz;
}


////////////////////////////////////
//
//  Normalizing Constants
//
/////////////////////////////////////


void log_Z_L_T(long double logZ[MAX_ARR], int ntaxa, int nshapes, int nsites, int rankT, long double b, int maxD, long double*** propDistance_twoTrees)
{

  //printf("logZ_L_T:: start\n");
  //printf("nshapes=%d nsites=%d maxD=%d\n",nshapes,nsites,maxD+1);
  long double Z1 = treeSpaceSize(ntaxa),logZ1 = (long double) log(Z1);
  long double Z=0, sublogZ=0, sumMaxD,temp;
  int s,d,x,l, s2;
  long double log_multiplier = (long double) b*(ntaxa-3);
  //long double tmp[4][699][139];
  //printf("logZ_L_T:: end\n");
  long double matlogZ[nshapes][maxD+1][nsites+1];
  //printf("logZ_L_T:: end\n");

  //printf("rankT = %d maxD = %d nsites = %d\n",rankT,maxD, nsites);
  for(l=0;l<=nsites;l++)
    for(d=0;d<=maxD;d++)
      for(s=0;s<nshapes;s++)
	{
	  //if(d>(l-1)*(ntaxa-3))
	  // {
	      //printf("here: d = %d l=%d (l-1)*(ntaxa-3)=%d\n",d,l, (l-1)*(ntaxa-3));
	  //  matZ[s][d][l] = 0;
	  // }
	  //else{
	    if(l<2)
	      matlogZ[s][d][l] = (long double) log(0);
	    else if(l==2)
	      {		
		if(d<=ntaxa-3)		     
		  {
		    matlogZ[s][d][l] = 0;
		    for(x=0;x<=d;x++)
		      {
			/*
			  Compute log(f1+f2) based on gi=log(fi)
			  log(f1+f2) = log(exp(g1)*(1+exp(g2-g1)) = g1 + log(1+exp(g2-g1))
			  Therefore,
			  double logsum = g[0];
			  for(i=1; i<=k; i++)
			     logsum += log(1+exp(g[i]-logsum));
			 */
			for(s2=0;s2<nshapes;s2++)
			  {
			    if(s2==0 || exp(sublogZ) == 0)
			      {
				sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]);				
				//printf("logZ1=%Le\t log(propDistance_twoTrees[s][s2][x])=%f\n",logZ1,log(propDistance_twoTrees[s][s2][x]));
			      }
			    else
			      {
				temp = sublogZ;
				sublogZ += (long double) log(1+ exp(logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) - temp));				
				//printf("log(propDistance_twoTrees[s][s2][x])=%f\n",log(propDistance_twoTrees[s][s2][x]));
			      }
			    //printf("l=%d\t d=%d\t s=%d\t x=%d\t s2=%d\t sublogZ=%Le\n",l,d,s,x,s2,sublogZ);
			  }
			if(x==0)
			  matlogZ[s][d][l] = sublogZ -b*x + log_multiplier;
			else
			  {
			    temp = matlogZ[s][d][l];
			    matlogZ[s][d][l] += (long double) log(1+exp(sublogZ - b*x + log_multiplier - temp));
			  }
			//printf("when nsites=2: s=%d, d=%d l=%d b=%Le| x = %d subZ=%Le exp(-b*x) = %Le exp(b*(ntaxa-3-x))=%Le\n",s,d,l,b,x,subZ,(long double) exp(-b*(double) x),(long double) exp(b*(long double) (ntaxa-3-x)));
			//printf("multiplier!\n");
		      }
		    sublogZ = matlogZ[s][d][l] - log_multiplier;
		    matlogZ[s][d][l] = sublogZ;
		  }
		else
		  matlogZ[s][d][l] = matlogZ[s][ntaxa-3][l];
		
		
		//printf("matZ[s][d][l] =%Le\n",matZ[s][d][l]);
	      }
	    else // if l>2
	      {
		if(d<=ntaxa-3)
		  sumMaxD = d;
		else
		  sumMaxD = ntaxa-3;
		for(x=0;x<=sumMaxD;x++)
		  {
		    for(s2=0;s2<nshapes;s2++)
		      {
			//printf("s = %d x = %d\n",s,x);
			if(s2==0 || exp(sublogZ)==0)
			  {
			    if(d-2*x<0)
			      sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) + matlogZ[s2][0][l-1];
			    else
			      sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) + matlogZ[s2][d-x][l-1];
			  }
			else
			  {
			    temp = sublogZ;
			    if(d-2*x<0)
			      sublogZ += (long double) log(1+ exp(logZ1 + log(propDistance_twoTrees[s][s2][x]) + matlogZ[s2][0][l-1] - temp));
			    else
			      sublogZ += (long double) log(1+ exp(logZ1+log(propDistance_twoTrees[s][s2][x])+matlogZ[s2][d-x][l-1] - temp));
			  }
			//printf("when nsites>2: Z1 = %Le propDistance_twoTrees[s][s2][x]=%Le matZ[s2][d-2*x][l-1]=%Le\n",Z1,propDistance_twoTrees[s][s2][x],matZ[s2][d-x][l-1]);
		      }
		    if(x==0)
		      matlogZ[s][d][l] = sublogZ - b*x + log_multiplier;
		    else	
		      {		   
			temp = matlogZ[s][d][l];
			matlogZ[s][d][l] += (long double) log(1+ exp(sublogZ - b*x + log_multiplier - temp)) ; 
		      }
		    
		    //printf("when nsites>2: s=%d, d=%d l=%d| x = %d subZ=%Le exp(-b*x)=%Le exp(b*(ntaxa-3-x))= %Le matZ[s][d][l]=%Le \n",s,d,l,x,subZ,(long double) exp(- b*(double) x),(long double) exp(b*(double) (ntaxa-3-x)),matZ[s][d][l]);
		    //printf("ntaxa-3-x=%d ntaxa-3-x=%Le\n",ntaxa-3-x,(long double)(ntaxa-3-x));
		    //printf("multiplier!\n");
		  }
		sublogZ = matlogZ[s][d][l] - log_multiplier;
		matlogZ[s][d][l] = sublogZ;
		//matZ[s][d][l] /= multiplier;
		//printf("s=%d\t d=%d\t l=%d\t matZ[s][d][l] = %Le\n",s,d,l,matlogZ[s][d][l]);
	      }
	  }
  
  for(s=0;s<nshapes;s++)
    logZ[s] = matlogZ[s][maxD][nsites];
  //	}
  /*  
  for(l=0;l<=nsites;l++)
    for(d=0;d<=maxD;d++)
      for(s=0;s<nshapes;s++)
	printf("s=%d d=%d l=%d matZ[s][d][l]=%Le\n",s,d,l,matlogZ[s][d][l]);
 //*/
  //printf("ntaxa = %d, rankT = %d maxD = %d nsites = %d\n",ntaxa,rankT,maxD, nsites);
  //printf("final:  matZ[0][D2][2]=%Le\n", matZ[0][ntaxa-3][2]);
  //printf("final:  matZ[rankT][maxD][nsites]=%Le\n", matZ[rankT][maxD][nsites]);
      //*/
  //return Z = matlogZ[rankT][maxD][nsites];

}




void log_Z_L_T_lessMemory(long double logZ[MAX_ARR], int ntaxa, int nshapes, int nsites, int rankT, long double b, int maxD, long double*** propDistance_twoTrees)
{
  
  //printf("logZ_L_T:: start\n");
  //printf("nshapes=%d nsites=%d maxD=%d\n",nshapes,nsites,maxD+1);
  /*int i,j,k;
  for(i=0;i<nshapes;i++)
    {
      for(j=0;j<nshapes;j++)
	{	    
	  for(k=0;k<=ntaxa-3;k++)
	    {
	      printf("%d %d %d %Le\t",i,j,k, propDistance_twoTrees[i][j][k]);
	    }
	  printf("\n");
	}
    }
  */
  long double Z1 = treeSpaceSize(ntaxa),logZ1 = (long double) log(Z1);
  long double Z=0, sublogZ=0, sumMaxD,temp;
  int s,d,x,l, s2;
  long double log_multiplier = (long double) b*(ntaxa-3);
  //long double tmp[4][699][139];
  //printf("maxD=%d\n",maxD);
  //printf("logZ_L_T:: end\n");
  //long double matlogZ_present[nshapes][maxD+1], matlogZ_past[nshapes][maxD+1];
  long double** matlogZ_present;
  matlogZ_present = malloc(nshapes * sizeof(long double*));
  for(s=0;s<nshapes;s++)
    {
      matlogZ_present[s] = malloc((maxD+1) * sizeof(long double));
    }

  long double matlogZ_past[nshapes][ntaxa-3+1], matlogZ_past_atZero[nshapes];
  //printf("logZ_L_T:: end\n");

  //printf("rankT = %d maxD = %d nsites = %d\n",rankT,maxD, nsites);
  for(l=2;l<=nsites;l++)
    {
      for(d=0;d<=maxD;d++)
	{
	  if(d<=ntaxa-3)
	    sumMaxD = d;
	  else
	    sumMaxD = ntaxa-3;
	  if(l>2)
	    {
	      for(s=0;s<nshapes;s++)
		{
		  matlogZ_past_atZero[s] = matlogZ_present[s][0];
		  for(x=0;x<sumMaxD;x++)
		    matlogZ_past[s][x] = matlogZ_past[s][x+1];
		  matlogZ_past[s][(int) sumMaxD] = matlogZ_present[s][d];
		}
	    }
	  for(s=0;s<nshapes;s++)
	    {
	      if(l==2)
		{		
		  if(d<=ntaxa-3)		     
		    {
		      matlogZ_present[s][d] = 0;
		      for(x=0;x<=d;x++)
			{
			  /*
			    Compute log(f1+f2) based on gi=log(fi)
			    log(f1+f2) = log(exp(g1)*(1+exp(g2-g1)) = g1 + log(1+exp(g2-g1))
			    Therefore,
			    double logsum = g[0];
			    for(i=1; i<=k; i++)
			    logsum += log(1+exp(g[i]-logsum));
			  */
			  for(s2=0;s2<nshapes;s2++)
			    {
			      if(s2==0 || exp(sublogZ) == 0)
				{
				  sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]);				
				  //printf("logZ1=%Le\t log(propDistance_twoTrees[s][s2][x])=%f\n",logZ1,log(propDistance_twoTrees[s][s2][x]));
				}
			      else
				{
				  temp = sublogZ;
				  sublogZ += (long double) log(1+ exp(logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) - temp));				
				  //printf("log(propDistance_twoTrees[s][s2][x])=%f\n",log(propDistance_twoTrees[s][s2][x]));
				}
			      //printf("l=%d\t d=%d\t s=%d\t x=%d\t s2=%d\t sublogZ=%Le\n",l,d,s,x,s2,sublogZ);
			    }
			  if(x==0)
			    matlogZ_present[s][d] = sublogZ -b*x + log_multiplier;
			  else
			    {
			      temp = matlogZ_present[s][d];
			      matlogZ_present[s][d] += (long double) log(1+exp(sublogZ - b*x + log_multiplier - temp));
			    }
			  //printf("when nsites=2: s=%d, d=%d l=%d b=%Le| x = %d subZ=%Le exp(-b*x) = %Le exp(b*(ntaxa-3-x))=%Le\n",s,d,l,b,x,subZ,(long double) exp(-b*(double) x),(long double) exp(b*(long double) (ntaxa-3-x)));
			  //printf("multiplier!\n");
			}
		      sublogZ = matlogZ_present[s][d] - log_multiplier;
		      matlogZ_present[s][d] = matlogZ_past[s][d] = sublogZ;
		    }
		  else
		    matlogZ_present[s][d] = matlogZ_present[s][ntaxa-3];
		  
		  
		  //printf("matZ[s][d][l] =%Le\n",matZ[s][d][l]);
		}
	      else // if l>2
		{
		  if(d<=ntaxa-3)
		    sumMaxD = d;
		  else
		    sumMaxD = ntaxa-3;
		  for(x=0;x<=sumMaxD;x++)
		    {
		      for(s2=0;s2<nshapes;s2++)
			{
			  //printf("s = %d x = %d\n",s,x);
			  if(s2==0 || exp(sublogZ)==0)
			    {
			      if(d-2*x<0)
				sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) + matlogZ_past_atZero[s2];
			      else if(x==0)
				sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) + matlogZ_present[s2][d];
			      else
				sublogZ = logZ1 + (long double) log(propDistance_twoTrees[s][s2][x]) + matlogZ_past[s2][ntaxa-3-x];
			    }
			  else
			    {
			      temp = sublogZ;
			      if(d-2*x<0)
				sublogZ += (long double) log(1+ exp(logZ1 + log(propDistance_twoTrees[s][s2][x]) + matlogZ_past_atZero[s2] - temp));
			      else
				sublogZ += (long double) log(1+ exp(logZ1+log(propDistance_twoTrees[s][s2][x])+matlogZ_past[s2][ntaxa-3-x] - temp));
			    }
			  //printf("when nsites>2: Z1 = %Le propDistance_twoTrees[s][s2][x]=%Le matZ[s2][d-2*x][l-1]=%Le\n",Z1,propDistance_twoTrees[s][s2][x],matZ[s2][d-x][l-1]);
			}
		      if(x==0)
			matlogZ_present[s][d] = sublogZ - b*x + log_multiplier;
		      else	
			{		   
			  temp = matlogZ_present[s][d];
			  matlogZ_present[s][d] += (long double) log(1+ exp(sublogZ - b*x + log_multiplier - temp)) ; 
			}
		      
		      //printf("when nsites>2: s=%d, d=%d l=%d| x = %d subZ=%Le exp(-b*x)=%Le exp(b*(ntaxa-3-x))= %Le matZ[s][d][l]=%Le \n",s,d,l,x,subZ,(long double) exp(- b*(double) x),(long double) exp(b*(double) (ntaxa-3-x)),matZ[s][d][l]);
		      //printf("ntaxa-3-x=%d ntaxa-3-x=%Le\n",ntaxa-3-x,(long double)(ntaxa-3-x));
		      //printf("multiplier!\n");
		    }
		  sublogZ = matlogZ_present[s][d] - log_multiplier;
		  matlogZ_present[s][d] = sublogZ;
		  //matZ[s][d][l] /= multiplier;
		  //printf("s=%d\t d=%d\t l=%d\t matZ[s][d][l] = %Le\n",s,d,l,matlogZ[s][d][l]);
		}
	    }   
	}
      /*
      for(d=0;d<=maxD;d++)
	for(s=0;s<nshapes;s++)
	  matlogZ_past[s][d] = matlogZ_present[s][d];
      */
    }
  for(s=0;s<nshapes;s++)
    logZ[s] = matlogZ_present[s][maxD];
  //	}
  /*  
      for(l=0;l<=nsites;l++)
      for(d=0;d<=maxD;d++)
      for(s=0;s<nshapes;s++)
      printf("s=%d d=%d l=%d matZ[s][d][l]=%Le\n",s,d,l,matlogZ[s][d][l]);
  //*/
  //printf("ntaxa = %d, rankT = %d maxD = %d nsites = %d\n",ntaxa,rankT,maxD, nsites);
  //printf("final:  matZ[0][D2][2]=%Le\n", matZ[0][ntaxa-3][2]);
  //printf("final:  matZ[rankT][maxD][nsites]=%Le\n", matZ[rankT][maxD][nsites]);
  //*/
  //return Z = matlogZ[rankT][maxD][nsites];
  
}





long double logZ_L(int ntaxa, int nshapes, int nsites, long double b, long double*** propDistance_twoTrees)
{
  long double Z1 = treeSpaceSize(ntaxa), logZ1 = (long double) log(Z1);
  long double logZ, Z=0.0, subZ = 0.0, recurlogZ[MAX_ARR], temp
;
  int s,d;
  long double multiplier = exp(b);

  //printf("nshapes = %d, ntaxa = %d, nsites = %d Z1 = %Le log(Z1)=%Le logZ1 = %Le\n",nshapes,ntaxa,nsites,Z1,(long double)log(Z1),logZ1);
  
  //printf("logZ_L:: start\n");
  log_Z_L_T_lessMemory(recurlogZ,ntaxa,nshapes,nsites,s,b,(nsites-1)*(ntaxa-3),propDistance_twoTrees);
  //printf("logZ_L:: end\n");
  for(s=0;s<nshapes;s++)
    {
      subZ = 0;
      // computing the number of tree topologies with the given tree shape
      for(d=0;d<=ntaxa-3;d++)
	{
	  //printf("s=%d d=%d propDistance_twoTrees[0][s][d]=%Le\n",s,d,propDistance_twoTrees[0][s][d]);
	  subZ += propDistance_twoTrees[0][s][d];
	}
      //recurZ = Z_L_T(ntaxa,nshapes,nsites,s,b,(nsites-1)*(ntaxa-3),propDistance_twoTrees);
      //printf("s = %d b=%Le maxD =%d subZ=%Le recurZ =%Le\n",s,b,(nsites-1)*(ntaxa-3),subZ,recurZ);
      if(s==0)	
	Z = (long double) log(subZ) + recurlogZ[s];
      else
	{
	  temp = Z;
	  Z += (long double) log(1+exp(log(subZ)+recurlogZ[s] - temp));
	}
    }
  logZ = Z + logZ1;
  //printf("final log(Z) = %Le Z1 = %Le log(Z1)=%Le logZ=%Le \n",(long double)log(Z),Z1,logZ1,logZ);
  return logZ;  
}
   
