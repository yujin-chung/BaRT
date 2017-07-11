///////////////////////////////////////////
//
//  Bryant and Steel (2009)
//   The calculation of RF metric distribution given a tree shape
//
////////////////////////////////////////////////////

// Modified for Cecile Ane on July 14, 2011 by Yujin Chung 
// updated on March 14, 2012 by Yujin Chung

#include <math.h>
#include <stdio.h>
#include <stdlib.h> // malloc
#include "shape.h"


// Return # internal children nodes
int no_intchildrenNodes(struct node *root,struct node *nd)
{
  int i; 
  int n=0;
  
  for(i=0;i<3;i++)
    if(nd->listnodes[i] != NULL)
      if(nd->listnodes[i]->id != root->id)
	if(nd->listnodes[i]->terminal == 0)
	  n++;

  return n;
}


// Return # internal children nodes of the root
int no_childrenNodes_forRooted(struct node *nd)
{
  int i; 
  int n=0;
  
  for(i=0;i<3;i++)
    if(nd->listnodes[i] != NULL)
      n++;
  
  return n;
}


// Return # internal nodes
int no_intnod(struct node *nd, struct node *childnd)
{
  int i;
  int n=0;
  
  for(i=0;i<3;i++)
    if(childnd->listnodes[i] != NULL)
      if(childnd->listnodes[i]->terminal == 0)
	if(childnd->listnodes[i]->id != nd->id)
	  {
	    n++; //printf("id=%d\t terminal=%d\t n=%d\n",nd->listnodes[i]->id,nd->listnodes[i]->terminal,n);
	    n += no_intnod(childnd,childnd->listnodes[i]);
	  }
  
  return n;
}

// Return # internal nodes of the root
int no_intnod_root(struct node *root)
{
  int i;
  int n=0;
  
  for(i=0;i<3;i++)
    if(root->listnodes[i] != NULL)
      if(root->listnodes[i]->terminal == 0)
	{
	  n++;
	  n += no_intnod(root,root->listnodes[i]);
	}
  return n;
}

// Return # tree topologies with a given "nint" internal edges
//  The input argument -- the number of internal edges
double no_trees(int nint)
{
  int i;
  double n=1;
  if(nint !=0)    
    for(i=1;i<=(2*nint+1);i+=2)
      n *= i;
  return n;
}

double R(struct node *root, struct node *nd, int s, int k)
{
  // Bryant and Steel (2009)
  int i, findChildnodes=0, sub_k1, sub_k2, sub_s;
  double r=0, sub_r1, sub_r2; 
  double sub_r3;
  double temp;
  int n_e = no_intnod(root,nd); // the number of internal edges
  int n_childrenNod = no_intchildrenNodes(root,nd); // the number of children nodes
  struct node *childNode1, *childNode2;
  
  
  
  //printf("start R: n_e=%d\t nd->id=%d s=%d  no_intnod(nd)=%d\t k=%d\t n_trees=%d n_children=%d\n",n_e,nd->id,s,no_intnod(root,nd),k,no_trees(k),n_childrenNod);

  if(s==0) // Lemma 1
    {
      if(k==n_e)
	{
	  r = no_trees(k); 
	}
    }
  else
    {
      if(k>n_e) // Lemma 2-1
	{
	  r = 0;
	}
      else
	{
	  if(n_childrenNod == 0) // Lemma 2-2
	    {
	      r = 0;
	    }
	  else if(n_childrenNod == 1) // Lemma 2-3
	    {	      
	      findChildnodes=0;
	      i=0;
	      while(findChildnodes==0)
		{
		  if(nd->listnodes[i] != NULL)
		    {
		      //printf("in R: i=%d\t nd->id=%d\t nd->listnodes[i]->id=%d\t nd->listnodes[i]->terminal=%d\n",i,nd->id,nd->listnodes[i]->id,nd->listnodes[i]->terminal);
		      //if(nd->listnodes[i]->id != nd->id)
		      if(nd->listnodes[i]->terminal == 0 && nd->listnodes[i]->id != root->id)
			{
			  childNode1 = nd->listnodes[i];
			  findChildnodes=1;
			}
		    }
		  i++;
		}
	      
	      //printf("in R: no_children_nodes=%d, childrenNode1=%d\t \n",n_childrenNod,childNode1->id);
	      
	      
	      if(k==0)		
		{
		  r=0;
		  //printf("sub_k1's max =%d\n",(no_intnod(nd,childNode1)-s+1));
		  for(sub_k1=0;sub_k1<=(no_intnod(nd,childNode1)-s+1);sub_k1++)
		    {
		      r += R(nd,childNode1,s-1,sub_k1);
		      //printf("in R: sub_k1=%d\t max_k1=%d\t R=%d, no_intnod(nd,childNode1)=%d\n",sub_k1,(no_intnod(nd,childNode1)-s+1),R(nd,childNode1,s-1,sub_k1),no_intnod(nd,childNode1));
		    }
		}
	      else
		r = R(nd,childNode1,s,k-1) * (2*k+1);
	      
	    }
	  else if(n_childrenNod == 2)
	    {	      
	      findChildnodes=0; i=0;
	      while(findChildnodes<=1)
		{		      
		  if(nd->listnodes[i] != NULL)
		    if(nd->listnodes[i]->id != root->id)
		      if(nd->listnodes[i]->terminal == 0)
			{
			  if(findChildnodes==0)
			    childNode1 = nd->listnodes[i];
			  else if(findChildnodes==1)
			    childNode2 = nd->listnodes[i];
			  findChildnodes++;
			}
		  i++;
		}
	      
	      //printf("in R: nodeID=%d\t no_children_nodes=%d, childrenNode1=%d\t childrenNode2=%d\n",nd->id,n_childrenNod,childNode1->id,childNode2->id);

	      if(k==0) // Lemma 2-4
		if(s <=1)
		  r = 0;
		else
		  {
		    r=0;
		    for(sub_s=0;sub_s<=s-2;sub_s++)
		      {
			if(sub_s <= no_intnod(nd,childNode1))
			  {
			    sub_r1=sub_r2=0;
			    for(sub_k1=0;sub_k1<=no_intnod(nd,childNode1)-sub_s;sub_k1++)
			      sub_r1 += R(nd,childNode1,sub_s,sub_k1);
			    if(no_intnod(nd,childNode2)-(s-2-sub_s)>=0)
			      {
				for(sub_k2=0;sub_k2<=no_intnod(nd,childNode2)-(s-2-sub_s);sub_k2++)
				  {
				    sub_r2 += R(nd,childNode2,s-2-sub_s,sub_k2);
				    //printf("sub_k2=%d\t R=%d\t sub_r2=%d\n",sub_k2,R(nd,childNode2,s-2-sub_s,sub_k2),sub_r2);
				  }
			      }
			    temp = log(sub_r1) + log(sub_r2);
			    r += exp(temp);
			  }
			
			//printf("in R: s=%d\t k=%d\t sub_s=%d\t r1=%d\t r2=%d\t R=%d\n",s,k,sub_s,sub_r1,sub_r2,r);
		      }	
		  }
	      else // if k >=1, Lemma 2-5
		{
		  r=0;
		  for(sub_s=0;sub_s<=s;sub_s++)
		    {
		      
		      if(sub_s==s && k>=2)
			{
			  // the third term
			  sub_r3=0;
			  for(sub_k1=0;sub_k1<=k-2;sub_k1++)
			    {			      		   
			      //printf("in R: 3rd term:: (childNode1->id=%d,sub_s=%d,sub_k1=%d),(childNode2->id=%d,s-sub_s=%d,k-2-sub_k1=%d)\n",childNode1->id,sub_s,sub_k1,childNode2->id,s-sub_s,k-2-sub_k1);			      			     
			      temp =  log(R(nd,childNode1,sub_s,sub_k1))+log(R(nd,childNode2,s-sub_s,k-2-sub_k1)) + log(no_trees(k)) - log(no_trees(sub_k1))-log(no_trees(k-2-sub_k1));
			      //printf("in R: 3rd term(1):: temp =%f\n",temp);
			      sub_r3 += exp(temp);
			    }
			  r += sub_r3; 
			  //printf("in R: s=%d\t k=%d\t sub_s=%d\t r1=%d\t r2=%d\t r3=%d\n",s,k,sub_s,sub_r1,sub_r2,sub_r3);
			}
		      else if(sub_s !=s)
			{
			  if(no_intnod(nd,childNode1)-sub_s >= 0)
			    {	
			      // the first term
			      sub_r1=0;
			      for(sub_k1=0;sub_k1<=no_intnod(nd,childNode1)-sub_s;sub_k1++)
				{				  	   
				  //printf("in R: 1st term:: (childNode1->id=%d,sub_s=%d,sub_k1=%d),(childNode2->id=%d,s-1-sub_s=%d,k-1=%d)\n",childNode1->id,sub_s,sub_k1,childNode2->id,s-1-sub_s,k-1);
				   sub_r1 += R(nd,childNode1,sub_s,sub_k1);
				}
			      temp = log(sub_r1) +log( R(nd,childNode2,s-1-sub_s,k-1))-log(no_trees(k-1))+log( no_trees(k));	
			      //printf("in R: 1st term:: exp(temp) =%f\n",exp(temp)); 
			      r += exp(temp);
			      //printf("Lemma 2-5, case1:\n\t id=%d\t s-1-sub_s=%d\t k-1=%d\n",childNode2->id,s-1-sub_s,k-1);
			    }

			  if(no_intnod(nd,childNode2)-sub_s >= 0)
			    {
			      // the second term
			      sub_r2 = 0;
			      for(sub_k1=0;sub_k1<=no_intnod(nd,childNode2)-sub_s;sub_k1++)
				{		   
				  //printf("in R: 2nd term:: (childNode2->id=%d,sub_s=%d,sub_k1=%d),(childNode1->id=%d,s-1-sub_s=%d,k-1=%d)\n",childNode2->id,sub_s,sub_k1,childNode1->id,s-1-sub_s,k-1);
				  sub_r2 += R(nd,childNode2,sub_s,sub_k1);
				}
			      temp = log(sub_r2)+log( R(nd,childNode1,s-1-sub_s,k-1))+log(no_trees(k))-log(no_trees(k-1));
			      //printf("in R: 2nd term:: temp =%d\n",temp);
			      /*
			      if(sub_r2<0)
				printf("in R: s=%d\t k=%d\t sub_r2=%d\n",s,k,sub_r2);
			      */
			      r += exp(temp);
			      //printf("Lemma 2-5, case2\n");//:\n\t id=%d\t s-1-sub_s=%d\t k-1=%d\n",childNode2->id,s-1-sub_s,k-1);
			    }

			  if(k>=2)
			    {
			      // the third term
			      sub_r3=0;
			      for(sub_k1=0;sub_k1<=k-2;sub_k1++)
				{			      		   
				  //printf("in R: 3rd term(2):: (childNode1->id=%d,sub_s=%d,sub_k1=%d),(childNode2->id=%d,s-sub_s=%d,k-2-sub_k1=%d)\n",childNode1->id,sub_s,sub_k1,childNode2->id,s-sub_s,k-2-sub_k1);
				  //printf("in R:sub_k1=%d s=%d\t k=%d\t sub_r3=%f\n",sub_k1,s,k,sub_r3);
				  //temp = no_trees(k);
				  //printf("part of 3rd term(2) :no_trees(k)=%f\n",temp);
				  //temp = log(R(nd,childNode1,sub_s,sub_k1))+log(R(nd,childNode2,s-sub_s,k-2-sub_k1))+log(no_trees(k))-log(no_trees(sub_k1)*no_trees(k-2-sub_k1));
				  temp = log(R(nd,childNode1,sub_s,sub_k1))+log(R(nd,childNode2,s-sub_s,k-2-sub_k1))-log(no_trees(sub_k1))-log(no_trees(k-2-sub_k1))+log(no_trees(k));
				  sub_r3 += exp(temp);
				   //printf("in R: 3rd term(2):: temp =%f\n",temp);
				  /*
				  if(sub_r3<0)
				    {
				      printf("Warning!::in R:sub_k1=%d s=%d\t k=%d\t sub_r3=%f\n",sub_k1,s,k,sub_r3);
				      printf("\t R(nd,childNode1,sub_s,sub_k1) = %d\n ", R(nd,childNode1,sub_s,sub_k1));
				      printf("\t R(nd,childNode2,s-sub_s,k-2-sub_k1) = %d\n ",R(nd,childNode2,s-sub_s,k-2-sub_k1));
				      printf("\t no_trees(k)/(no_trees(sub_k1)*no_trees(k-2-sub_k1))=%d\n",no_trees(k)/(no_trees(sub_k1)*no_trees(k-2-sub_k1)));
				    }
				  */
				}
			      r += sub_r3;
			    }
			  //if(s==1)
			  //printf("in R: s=%d\t k=%d\t sub_s=%d\t r1=%d\t r2=%d\t r3=%d\n",s,k,sub_s,sub_r1,sub_r2,sub_r3);
			}		      
		    }
		  
		}
	    }
	}
    }
  
  //printf("End R: nd->id=%d\t s=%d\t k=%d\t r=%f\n",nd->id,s,k,r); 
  //  /*
  //    if(k==6 && s==1)
  // if(r>0)
  //printf("in R: nd->id=%d\t s=%d\t k=%d\t r=%f\n",nd->id,s,k,r); 
  //*/
  return r;
}

double r_s(struct node *nd, int s)
{
  // Eq. (8) in Bryant and Steel (2009)
  int k;
  double r=0, temp;
  for(k=0;k<=no_intnod_root(nd)-s;k++)
    {
      temp =R(nd->listnodes[0],nd,s,k);
      r += temp;
      //      if(s==1|| s==0)
      //printf("in r_s:: s=%d\t k=%d\t nd->listnodes[0]->id=%d\t R(nd,s,k)=%f\n",s,k,nd->listnodes[0]->id, temp);
    }
  return r;
}

int factorial(int x)
{
  int i,fac=1;
  if(x==0)
    return 1;
  else
    {
      for(i=1;i<=x;i++)
	fac *= i;
      return fac;
    }
}

double log_factorial(int x)
{
  int i;
  double logfac=0;
  if(x==0)
    return 0;
  else
    {
      for(i=1;i<=x;i++)
	logfac += log(i);
      return logfac;
    }
  
}

void gen_r(double gen_r_s[], int ntaxa, struct node *nd)
{
  int s;
  int n_intnod = ntaxa-3;
  for(s=0;s<=n_intnod;s++)
    {
      gen_r_s[s] = r_s(nd,s);
      //      printf("in RFdistribution: s=%d\t gen_r_s=%f\n",s, gen_r_s[s]);
    }
}


double RFdistribution(struct node *nd, int ntaxa, int d)
{
  int n_intnod = ntaxa-3;
  double gen_r_s[n_intnod+1];
  int s;
  double q_s,temp;
  int iseven;
  
  
  gen_r(gen_r_s,ntaxa,nd);

  q_s = gen_r_s[d];
  //  printf("in RFdistribution: d=%d\t gen_r_s=%f q_s=%f\n",d,gen_r_s[d],q_s);
  if(d<n_intnod)
    {
      for(s=d+1;s<=n_intnod;s++)
	{
	  iseven = (s-d)/2;
	  if(2*iseven == s-d)
	    {
	      temp = log(gen_r_s[s])+log_factorial(s)-log_factorial(d)-log_factorial(s-d);
	      q_s += exp(temp);
	      //printf("even\t q_s=%d\n",factorial(s-d));
	    }
	  else
	    {
	      temp = log(gen_r_s[s])+log_factorial(s)-log_factorial(d)-log_factorial(s-d);
	      //printf("odds:: temp = %f exp(temp)=%f\n",temp,exp(temp));
	      q_s -= exp(temp);
	      //printf("odd\n");
	    }
	  //  printf("in RFdistribution: d=%d\t s=%d\t gen_r_s=%f q_s=%f\n",d,s,gen_r_s[s],q_s);
	}
    }
  return q_s;  
}
