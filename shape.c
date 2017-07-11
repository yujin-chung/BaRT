/*
  Created Feb. 7, 2011 by Yujin Chung

  Updated Feb. 14, 2012 by Yujin Chung
*/



/* 
   file name: shape.c

   brief:
          list tree shapes   
 */



#include <math.h>
#include <stdio.h>
#include <stdlib.h> // malloc
#include "shape.h"


#define MAX_TAXA 15


//int no_rshapes(int ntaxa);
int listrtshape(int ntaxa);
struct node *talloc(void);
int F(int i, int j, int k);
int no_trptr(int ntaxa);
//int no_urshapes(int ntaxa);
int no_dbltr(int ntaxa);
int delta2(int m);
int delta3(int m);
int U(int alpha, int beta, int ntaxa);
int invW1(int z, int w);
int invW2(int z, int w);
double derv_i1(int i, int z, int alpha, int beta);
double derv_i2(int i, int z, int alpha, int beta);
void rRanks(int ntaxa, int rank, int subrank[2]);
void printNode(struct node *node, struct node *childnode);

/* talloc: make a node */
struct node *talloc(void)
{
  return (struct node *) malloc(sizeof(struct node));
}

//////////////////////////////////////////////////
//
//  No. tree topologies on a specific shape
//
//////////////////////////////////////////////////

int size_subtree(struct node *tr)
{
  int i, n=0;

  if(tr->terminal == 1)
    n=1;
  else
    {
      for(i=0;i<3;i++)
	if(tr->listnodes[i] != NULL)
	  {
	    if(tr->terminal == 0)
	      n += size_subtree(tr->listnodes[i]);
	    else
	      n++;
	  }
    }

  return n;
}

int identicalShape(struct node *tr1, struct node *tr2)
{
  int size1 = size_subtree(tr1);
  int size2 = size_subtree(tr2);
  //printf("size1 = %d, size2 = %d\n",size1,size2);
  if(tr1->rank==tr2->rank && size1==size2)
    return 1;
  else
    return 0;
}

double log_noSymmetries(struct node *tr)
{
  int identity1, identity2;
  int nchild = no_childrenNodes_forRooted(tr); //printf("nchild = %d\n",nchild);
  double res = 0.0;
  //printf("res= %f\n",res);
  if(tr->terminal==0)
    {
      if(nchild==3)
	{
	  identity1 = identicalShape(tr->listnodes[0],tr->listnodes[1]);
	  identity2 = identicalShape(tr->listnodes[1],tr->listnodes[2]);
	  //printf("identity1 = %d, identity2 = %d\n",identity1,identity2);
	  if(identity1==1 && identity2==1)
	    {
	      res += log(3)+log(2)+3*log_noSymmetries(tr->listnodes[0]);
	      //printf("res= %f\n",res);
	    }
	  else if(identity1==1)
	    res += log(2) + 2*log_noSymmetries(tr->listnodes[0]) + log_noSymmetries(tr->listnodes[2]);
	  else if(identity2==1)
	    res += log(2) + log_noSymmetries(tr->listnodes[0]) + 2*log_noSymmetries(tr->listnodes[2]);
	  else
	    res += log_noSymmetries(tr->listnodes[0]) + log_noSymmetries(tr->listnodes[1]) + log_noSymmetries(tr->listnodes[2]);
	}
      else if(nchild==2)
	{
	  //printf("nchild is 2\n");
	  identity1 = identicalShape(tr->listnodes[0],tr->listnodes[1]);
	  if(identity1==1)
	    res += log(2) + 2*log_noSymmetries(tr->listnodes[0]);
	  else
	    res += log_noSymmetries(tr->listnodes[0]) + log_noSymmetries(tr->listnodes[1]);
	}
      
    }
  //printf("final res= %f\n",res);
  return res;
}

double log_noTrees_onShape(int ntaxa, struct node *tr)
{
  double log_numerator = 0;
  double log_denominator = log_noSymmetries(tr); //printf("log_denominator = %f\n",log_denominator);

  int i=1;
  for(i=1;i<=ntaxa;i++)
    {
      log_numerator += log(i);
    }

  return log_numerator - log_denominator;
}

/* 
   transform tree from LLC ordering form to unrooted tree
 */

struct node *LLC2urtree(struct node *root)
{
  struct node *ur;

  if(root->listnodes[2]==NULL) // double truncated tree
    {
      root->listnodes[0]->listnodes[2] = root->listnodes[1];
      root->listnodes[1]->listnodes[2] = root->listnodes[0];

      /*
      int i; 
      for(i=0;i<3;i++)
	if(root->listnodes[0]->listnodes[i] !=NULL)
	  printf("i=%d\t id=%d\n",i,root->listnodes[0]->listnodes[i]->id);
      */
      ur = root->listnodes[0];
    }
  else
    ur = root;
  
  return ur;
}

/*
  rearrange node (fill out NULL node)
 */

struct node *noNULL(struct node *nd, struct node *childnd)
{
  int i;
  
  if(childnd->terminal==1)
    childnd->listnodes[0] = nd;
  else
    {
      for(i=0;i<3;i++)
	{
	  if(childnd->listnodes[i]==NULL)
	    childnd->listnodes[i] = nd;
	  else if(childnd->listnodes[i]->id != nd->id)
	    childnd->listnodes[i] = noNULL(childnd,childnd->listnodes[i]);
	}
    }
  
  return childnd;
}


struct node *noNULL_node(struct node *root)
{
  int i;
  for(i=0;i<3;i++)
    root->listnodes[i] = noNULL(root,root->listnodes[i]);

  return root;
}


/* 
   reroot to a node nearest to a tip 
*/
struct node *reroot(struct node *root)
{
  struct node *nd;
  int i, isthenode;

  // rerooting
  if(root->terminal==1)
    {
      i=0; isthenode=0;
      while(isthenode==0)
	{
	  if(root->listnodes[i] != NULL)
	    {
	      nd = root->listnodes[i];
	      isthenode=1;
	    }
	  i++;
	}
    }
  else
    {
      if(root->listnodes[0]->terminal!=1)
	{
	  nd = reroot(root->listnodes[0]);
	}
      else
	nd = root;
    }

  return nd; 
}

/*
struct node *pseudo_rooted(struct node *nd, struct node *childnd)
{
  printf("%d\n",childnd->terminal);
  if(childnd->terminal==0)
    {
      printf("yes\n");
      childnd->listnodes[0]=NULL;
      printf("yes\n");
    }
  else
    {
      int i;
      for(i=0;i<3;i++)
	{
	  if(childnd->listnodes[i]->id == nd->id)
	    {
	      childnd->listnodes[i] = NULL;
	    }
	}
    }
  return childnd;
}
  
struct node *pseudo_rooted_Root(struct node *root)
{
  // root =  unrooted tree
  int i; 
  for(i=0;i<3;i++)
    root->listnodes[i] = pseudo_rooted(root,root->listnodes[i]);

  return root;
}
*/

void printAllShapes(int ntaxa)
{
  int st_sizeRank[no_urshapes(ntaxa)][6];
  int nshapes = no_urshapes(ntaxa);  
  int r;
  int intnod = ntaxa+1,ternod=0;
  struct node *tree;

  subtree_sizeRank(ntaxa, st_sizeRank); 

   for(r=1;r<=nshapes;r++)
    {   
      tree = invcrank_original(ntaxa,intnod,ternod,r,st_sizeRank[r-1][0],st_sizeRank[r-1][1],st_sizeRank[r-1][2],st_sizeRank[r-1][3],st_sizeRank[r-1][4],st_sizeRank[r-1][5]);      
      printRoot_tips(tree);
    }
}

void printRoot_tips(struct node *node)
{ 
  /*
    to print out a tree in newick format
    printing leaf label only
  */
  
  printf("(");
  if(node->terminal==1) //if terminal node
    {
      printf("%d",node->id);
    }  
  else
    {  
      printNode_tips(node,node->listnodes[0]);
      printf(",");
      //printf("yes, node->listnodes[1]->id=%d\n",node->listnodes[1]->id);
      printNode_tips(node,node->listnodes[1]);
      if(node->listnodes[2] != NULL)// && node->listnodes[2]->id != node->id)
	{
	  printf(",");
	  printNode_tips(node,node->listnodes[2]);
	}
    }
  printf(")");
  printf("\n");
}

void printNode_tips(struct node *node, struct node *childnode)
{
  if(childnode->terminal==1) //if terminal node
    printf("%d",childnode->id);
  else
    {
      printf("("); 
      //printf("nd->id=%d\t childnd->id=%d\t childnode->listnodes[0]->id=%d\n",node->id, childnode->id,childnode->listnodes[0]->id);
      if(childnode->listnodes[0] != NULL && childnode->listnodes[0]->id != node->id)	
	if(childnode->listnodes[0]->id != node->id)
	  {
	    printNode_tips(childnode,childnode->listnodes[0]);
	    printf(",");
	  }
      //printf("nd->id=%d\t childnd->id=%d\t childnode->listnodes[1]->id=%d\n",node->id, childnode->id,childnode->listnodes[1]->id);
      printNode_tips(childnode, childnode->listnodes[1]);
      if(childnode->listnodes[2] != NULL && childnode->listnodes[2]->id != node->id)
	{	  
	  if(childnode->listnodes[2]->id != node->id)
	    {
	      //printf(" childnode->id=%d\t node->id=%d\n", childnode->id,node->id);
	      printf(",");
	      printNode_tips(childnode, childnode->listnodes[2]);
	    }
	}
      printf(")");
    }
}


void printRoot(struct node *node)
{ 
  /*
    to print out a tree in newick format
  */

  printf("(");
  if(node->terminal==1) //if terminal node
    {
      printf("%d",node->id);
      int i;
      for(i=0;i<3;i++)
	if(node->listnodes[i] !=NULL)
	  {
	    printNode(node,node->listnodes[i]);
	  }
    }  
  else
    {  
      printNode(node,node->listnodes[0]);
      printf(",");
      //printf("yes, node->listnodes[1]->id=%d\n",node->listnodes[1]->id);
      printNode(node,node->listnodes[1]);
      if(node->listnodes[2] != NULL)// && node->listnodes[2]->id != node->id)
	{
	  printf(",");
	  printNode(node,node->listnodes[2]);
	}
    }
  printf(")");
  if(node->terminal!=1)
    printf("%d:%d",node->id,node->rank);
  printf("\n");
}

void printNode(struct node *node, struct node *childnode)
{
  if(childnode->terminal==1) //if terminal node
    printf("%d",childnode->id);
  else
    {
      printf("("); 
      //printf("nd->id=%d\t childnd->id=%d\t childnode->listnodes[0]->id=%d\n",node->id, childnode->id,childnode->listnodes[0]->id);
      if(childnode->listnodes[0] != NULL && childnode->listnodes[0]->id != node->id)	
	if(childnode->listnodes[0]->id != node->id)
	  {
	    printNode(childnode,childnode->listnodes[0]);
	    printf(",");
	  }
      //printf("nd->id=%d\t childnd->id=%d\t childnode->listnodes[1]->id=%d\n",node->id, childnode->id,childnode->listnodes[1]->id);
      printNode(childnode, childnode->listnodes[1]);
      if(childnode->listnodes[2] != NULL && childnode->listnodes[2]->id != node->id)
	{	  
	  if(childnode->listnodes[2]->id != node->id)
	    {
	      //printf(" childnode->id=%d\t node->id=%d\n", childnode->id,node->id);
	      printf(",");
	      printNode(childnode, childnode->listnodes[2]);
	    }
	}
      printf(")");
      printf("%d:%d",childnode->id,childnode->rank);
    }
}


/*

  Used the algorithms in

       Furnas (1984), The Generation of Random, Binary Unordered Trees (especially page 210 - 223 and 230-231)

  to list all shapes.  

*/

void subtree_sizeRank(int ntaxa, int sizeRank[][6])
{
  // return an array containing sizes and ranks of subtrees
  // sizeRank[][0..2] : sizes from left subtree to right subtree
  // sizeRank[][3..5] : ranks from left subtree to right subtree
  
  int r,rL, rM, rR, str_rM, str_rR;
  int alpha, beta, gamma, p;
  int i;
 
  r=0;
  for(alpha=1;alpha<=ntaxa/2;alpha++)
    for(beta=alpha;beta<=ntaxa/2;beta++)
      {
	gamma = ntaxa - alpha - beta;
	if(gamma<=ntaxa/2 && gamma>=beta && gamma < alpha+beta) // triple truncated tree
	  {
	    // decide corresponding rank
	    for(rL=1;rL<=no_rshapes(alpha);rL++)
	      {
		if(alpha<beta)	       
		  str_rM = 1;
		else
		  str_rM = rL;
		for(rM=str_rM;rM<=no_rshapes(beta);rM++)
		  {
		    if(beta<gamma)	       
		      str_rR = 1;
		    else
		      str_rR = rM;		  
		    for(rR=str_rR;rR<=no_rshapes(ntaxa-alpha-beta);rR++)
		      {
			sizeRank[r][0]=alpha;
			sizeRank[r][1]=beta;
			sizeRank[r][2]=gamma;
			sizeRank[r][3] = rL;
			sizeRank[r][4] = rM;
			sizeRank[r][5] = rR;
			r++;		    
		      }
		  }
	      }
	  }
	else if(gamma==0) // double truncated tree
	  {
	    // decide corresponding rank
	    for(rL=1;rL<=no_rshapes(alpha);rL++)
	      for(rM=rL;rM<=no_rshapes(beta);rM++)
		{
		  sizeRank[r][0]=sizeRank[r][1]=ntaxa/2;
		  sizeRank[r][2]=0;
		  sizeRank[r][3] = rL;
		  sizeRank[r][4] = rM;
		  sizeRank[r][5] = 0;
		  /*
		  for(i=0;i<6;i++)
		    printf("sizeRank=%d\t",r,sizeRank[r][i]);
		  */
		  r++;
		}
	  }
      }
}

void rRanks(int ntaxa, int rank, int subrank[4])
{
  int alpha, p, V,x,Rrank,Lrank;
  
  // Decide sizes (no. terminal nodes) of left and right subtrees
  // sizes: alpha, ntaxa-alpha
  alpha=0; p=0;
  while(p < rank)
    {
      alpha++;
      //printf("Node: id=%d\t terminal=%d\t, int_nod=%d\t ter_nod=%d\n",nd->id,nd->terminal,int_nod,ter_nod);
      p += no_rshapes(alpha)*no_rshapes(ntaxa-alpha);
      //printf("Node: id=%d\t terminal=%d\t, int_nod=%d\t ter_nod=%d\n",nd->id,nd->terminal,int_nod,ter_nod);
    }
  
  // Decide ranks of subtrees based on their sizes
  p -= no_rshapes(alpha)*no_rshapes(ntaxa-alpha);
  V = rank-p-1;
  if(ntaxa == 2*alpha)
    {
      x = no_rshapes(alpha);
      Lrank = 1 + fmax( (2*x-1-sqrt(pow(2*x-1,2)-8*(V-x+1)))/2 , (2*x+1-sqrt(pow(2*x+1,2)-8*V))/2 );
      Rrank = 1 + V - x*(x+1)/2 + (x-(Lrank-1))*(x-(Lrank-1)+1)/2 + (Lrank-1);
    }
  else
    {
      Lrank = 1 + V/no_rshapes(ntaxa-alpha);
      Rrank = 1 + V - (Lrank-1)*no_rshapes(ntaxa-alpha);
    }
  subrank[0]=alpha;
  subrank[1]=ntaxa-alpha;
  subrank[2]=Lrank;
  subrank[3]=Rrank;
}


struct node *invcrank_original(int ntaxa, int intnod, int ternod, int rk, int alpha, int beta, int gamma, int rL, int rM, int rR)
{
  // rk: LLC rank

  static int i_nod, t_nod;
  i_nod = intnod; t_nod = ternod;
  
  struct node *urnd = talloc();
  if(ntaxa<=1)
    {
      urnd->id = t_nod;
      urnd->terminal = 1;
      urnd->rank = 1;      
    }
  else
    {
      urnd->id = i_nod;
      urnd->terminal = 0;  
      urnd->rank = rk;
      
      int sub_sizeRank[4];
      
      
      if(alpha>1)
	i_nod++;
      else
	t_nod++;
      rRanks(alpha, rL, sub_sizeRank);
      urnd->listnodes[0] = invcrank_original(alpha,i_nod,t_nod,rL,sub_sizeRank[0],sub_sizeRank[1],0,sub_sizeRank[2],sub_sizeRank[3],0);
      //printf("Left Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",urnd->listnodes[0]->id, urnd->listnodes[0]->terminal,i_nod,t_nod);
      
      if(beta>1)
	i_nod++;
      else
	t_nod++;
      rRanks(beta, rM, sub_sizeRank);
      urnd->listnodes[1] = invcrank_original(beta,i_nod,t_nod,rM,sub_sizeRank[0],sub_sizeRank[1],0,sub_sizeRank[2],sub_sizeRank[3],0);
      //printf("Middle Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",urnd->listnodes[1]->id, urnd->listnodes[1]->terminal,i_nod,t_nod);

      if(gamma==0)
	{
	  urnd->listnodes[2] = NULL;
	}
      else
	{
	  if(gamma>1)
	    i_nod++;
	  else
	    t_nod++;
	  rRanks(gamma, rR, sub_sizeRank);
	  urnd->listnodes[2] = invcrank_original(gamma,i_nod,t_nod,rR,sub_sizeRank[0],sub_sizeRank[1],0,sub_sizeRank[2],sub_sizeRank[3],0);
	}	
      
    }

  return urnd;
}



struct node *invcrank_Root(int ntaxa, int intnod, int ternod, int alpha, int beta, int gamma, int rL, int rM, int rR)
{
  static int i_nod, t_nod;
  i_nod = intnod; t_nod = ternod;
  
  struct node *urnd = talloc();

  if(ntaxa<=1)
    {
      urnd->id = t_nod;
      urnd->terminal = 1;
    }
  else
    {      
      urnd->id = i_nod;
      urnd->terminal = 0;  

      int sub_sizeRank[4];
      
      
      if(alpha>1)
	i_nod++;
      else
	t_nod++;

      rRanks(alpha, rL, sub_sizeRank);
      printf("alpha=%d\n",alpha);
      urnd->listnodes[0] = invcrank(urnd,urnd->listnodes[0],alpha,i_nod,t_nod,sub_sizeRank[0],sub_sizeRank[1],sub_sizeRank[2],sub_sizeRank[3]);
      printf("Left Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",urnd->listnodes[0]->id, urnd->listnodes[0]->terminal,i_nod,t_nod);
      
      if(beta>1)
	i_nod++;
      else
	t_nod++;
      rRanks(beta, rM, sub_sizeRank);
      urnd->listnodes[1] = invcrank(urnd,urnd->listnodes[1],beta,i_nod,t_nod,sub_sizeRank[0],sub_sizeRank[1],sub_sizeRank[2],sub_sizeRank[3]);
      //printf("Middle Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",urnd->listnodes[1]->id, urnd->listnodes[1]->terminal,i_nod,t_nod);

      if(gamma==0)
	{
	  urnd->listnodes[2] = NULL;
	}
      else
	{
	  struct node *temp3;
	  if(gamma>1)
	    i_nod++;
	  else
	    t_nod++;
	  rRanks(gamma, rR, sub_sizeRank);
	  urnd->listnodes[2] = invcrank(urnd,urnd->listnodes[2],gamma,i_nod,t_nod,sub_sizeRank[0],sub_sizeRank[1],sub_sizeRank[2],sub_sizeRank[3]);
	}	
      
    }

  return urnd;
}


struct node *invcrank(struct node *nd, struct node *childnd, int ntaxa, int intnod, int ternod, int alpha, int beta, int rL, int rM)
{
  if(ntaxa<=1)
    {    
      //printf("ternod=%d\n",ternod);
      nd->listnodes[1]->id = ternod;
      childnd->terminal = 1;  
    }
  else
    {
      childnd->id = intnod;
      childnd->terminal = 0;  

      int sub_sizeRank[4];
      
      
      if(alpha>1)
	intnod++;
      else
	ternod++;

      rRanks(alpha, rL, sub_sizeRank);
      childnd->listnodes[0] = invcrank(childnd,childnd->listnodes[0],alpha,intnod,ternod,sub_sizeRank[0],sub_sizeRank[1],sub_sizeRank[2],sub_sizeRank[3]);
      //printf("Left Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",urnd->listnodes[0]->id, urnd->listnodes[0]->terminal,i_nod,t_nod);
      
      if(beta>1)
	intnod++;
      else
	ternod++;
      rRanks(beta, rM, sub_sizeRank);
      childnd->listnodes[1] = invcrank(childnd,childnd->listnodes[1],beta,intnod,ternod,sub_sizeRank[0],sub_sizeRank[1],sub_sizeRank[2],sub_sizeRank[3]);
      //printf("Middle Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",urnd->listnodes[1]->id, urnd->listnodes[1]->terminal,i_nod,t_nod);

      
      childnd->listnodes[2] = nd;
      
    }

  return childnd;
}



double derv_i1(int i, int z, int alpha, int beta)
{
  // Furnas (1984), page 221
  // derivative of the first cubic equation 

  return 3*pow(i,2) - 6*(z+1)*i + (3*pow(z,2)+6*z+2);
}

double derv_i2(int i, int z, int alpha, int beta)
{
  // Furnas (1984), page 221
  // derivative of the second cubic equation 

  return 3*pow(i,2) - 6*z*i + 3*pow(z,2)-1;
}


int invW1(int z, int w)
{
  int arg1, arg2;
  arg1 = (2*z-1 - sqrt( pow(2*z-1,2) - 8*(w-z+1) ))/2;
  arg2 = (2*z+1 - sqrt( pow(2*z+1,2) - 8*w) )/2;
  return fmax(arg1,arg2);
}


int invW2(int z, int w)
{
  return w - delta2(z) + delta2(z-invW1(z,w)) + invW1(z,w);
}

int U(int alpha, int beta, int ntaxa)
{
  // Furnas(1984), page 219

  int i,j,k;
  int max_index;
  int iseven = ntaxa/2;
  if(ntaxa == iseven*2)
    max_index = iseven-1;
  else
    max_index = iseven;

  int u=0;
  for(i=1;i<=max_index;i++)
    for(j=i;j<=max_index;j++)
      for(k=j;k<=max_index;k++)
	if(i+j+k==ntaxa)
	  if(ntaxa*i+j < ntaxa*alpha+beta)
	    u += F(i,j,k);

  return u;	  
}

int no_urshapes(int ntaxa)
{
  return no_trptr(ntaxa)+no_dbltr(ntaxa);
}

int no_dbltr(int ntaxa)
{
  // Furnas(1984), page 219
  // output: the number of double truncated trees

  int iseven = ntaxa/2;
  if(ntaxa == iseven*2)
    return delta2(no_rshapes(iseven));
  else
    return 0;
}


int no_trptr(int ntaxa)
{
  // Furnas(1984), page 219
  // output: the number of triple truncated trees
  
  int i,j,k;
  
  int n=0;
  int max_index;
  int iseven = ntaxa/2;
  if(ntaxa == iseven*2)
    max_index = iseven-1;
  else
    max_index = iseven;
  for(i=1;i<=max_index;i++)
    for(j=i;j<=max_index;j++)
      for(k=j;k<=max_index;k++)
	if(i+j+k==ntaxa)
	  {
	    n += F(i,j,k);
	    //printf("i=%d\t j=%d\t k=%d\t i+j+k=%d\t F(i,j,k)=%d\t, n=%d\n",i,j,k,i+j+k,F(i,j,k),n);
	  }
  return n;
}

int no_rshapes(int ntaxa)
{
  // input: the number of taxa
  // output: the number of rooted tree shapes on a given ntaxa
  // reference: Furnas(1984) page 231
  
  int ztable[MAX_TAXA];
  int i,j,even_i;
  
  i=ntaxa;
  if(i <= 2)
    ztable[i] = 1;
  else{
    ztable[i] = 0;
    for(j=1;j<=ntaxa/2;j++)
      {
	if(i == j*2) // when ntaxa is even
	  ztable[i] += no_rshapes(j)*(no_rshapes(j)+1)/2;
	else
	  ztable[i] += no_rshapes(j)*no_rshapes(i-j);
      }
  }
  
  return ztable[ntaxa];
}

int delta2(int m)
{
  // Furnas(1984), page 218
  return m*(m+1)/2;
}

int delta3(int m)
{
  // Furnas(1984), page 218
  return m*(m+1)*(m+2)/6;
}

int F(int i, int j, int k)
{
  // Furnas(1984), page 218

  if(i<j)
    {
      if(j<k)
	return no_rshapes(i)*no_rshapes(j)*no_rshapes(k);
      else if(j==k)
	return no_rshapes(i)*delta2(no_rshapes(j));
    }
  else if(i==j)
    {
      if(j<k)
	return delta2(no_rshapes(i))*no_rshapes(k);
      else if(j==k)
	return delta3(no_rshapes(i));
    }    
}


struct node *invrrank(int ntaxa, int rank, int intnod, int ternod)
{
  // input
  //    ntaxa: the number of taxa
  //    nodeID: the index of root or node of which two children nodes will be decided.
  // output
  //    two children nodes
  // reference: Furnas(1984) page 231
  
  int alpha,p,r,V,even_ntaxa,x,Lrank,Rrank,modulo;
  static int int_nod, ter_nod;
  int_nod = intnod; ter_nod = ternod;

  struct node *nd = talloc();
  
  if(ntaxa<=1)
    {
      nd->id = ter_nod;
      nd->terminal = 1;
    }
  else
    {
      nd->id = int_nod;
      nd->terminal = 0;

      // Decide sizes (no. terminal nodes) of left and right subtrees
      // sizes: alpha, ntaxa-alpha
      alpha=0; p=0;
      while(p < rank)
	{
	  alpha++;
	  //printf("Node: id=%d\t terminal=%d\t, int_nod=%d\t ter_nod=%d\n",nd->id,nd->terminal,int_nod,ter_nod);
	  p += no_rshapes(alpha)*no_rshapes(ntaxa-alpha);
	  //printf("Node: id=%d\t terminal=%d\t, int_nod=%d\t ter_nod=%d\n",nd->id,nd->terminal,int_nod,ter_nod);
	}
      
      //printf("Node: id=%d\t terminal=%d\t, int_nod=%d\t ter_nod=%d\n",nd->id,nd->terminal,intnod,ternod);
      // Decide ranks of subtrees based on their sizes
      p -= no_rshapes(alpha)*no_rshapes(ntaxa-alpha);
      V = rank-p-1;
      if(ntaxa == 2*alpha)
	{
	  x = no_rshapes(alpha);
	  Lrank = 1 + fmax( (2*x-1-sqrt(pow(2*x-1,2)-8*(V-x+1)))/2 , (2*x+1-sqrt(pow(2*x+1,2)-8*V))/2 );
	  Rrank = 1 + V - x*(x+1)/2 + (x-(Lrank-1))*(x-(Lrank-1)+1)/2 + (Lrank-1);
	}
      else
	{
	  Lrank = 1 + V/no_rshapes(ntaxa-alpha);
	  Rrank = 1 + V - (Lrank-1)*no_rshapes(ntaxa-alpha);
	}
      //printf("ntaxa=%d\t rank=%d\t alpha=%d\t Lrank=%d\t Rrank=%d\n",ntaxa,rank,alpha,Lrank,Rrank);

      //printf("Node: id=%d\t terminal=%d\t, int_nod=%d\t ter_nod=%d\n",nd->id,nd->terminal,intnod,ternod);
      struct node *temp1, *temp2;

      
      if(alpha>1)
	int_nod++;
      else
	ter_nod++;
      temp1=invrrank(alpha,Lrank,int_nod,ter_nod);
      nd->listnodes[0] = temp1;
      //printf("Left Child: id=%d\t terminal=%d\t int_nod=%d\t ter_nod=%d\n",nd->leftnode->id,nd->leftnode->terminal,intnod,ternod);

      if(ntaxa-alpha>1)
	int_nod++;
      else
	ter_nod++;
      temp2 = invrrank(ntaxa-alpha,Rrank,int_nod,ter_nod);
      nd->listnodes[1] = temp2;
      //printf("Right Child: id=%d\t terminal=%d\n",nd->rightnode->id,nd->rightnode->terminal);
      nd->listnodes[2] = NULL;
      
    }
  return nd;
}


