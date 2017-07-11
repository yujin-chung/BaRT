// Created by Yujin Chung on Nov 29, 2011
//
// modified on July 14th, 2012 by Yujin Chung
//

#include <math.h>
#include <stdio.h>
#include "shape.h"
#include "misc.h"


long double R_joint(struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F)
{ 

  long double res=0.0, subres, term1, term2;
  int size_vec1=getsize_vectorNode(node1), size_vec2=getsize_vectorNode(node2); //printf("size_vec=%d ",size_vec1);
  int i,j, no_child1, no_child2,found, found2;
  
  int sum_D1 = get_sum_vec(D1,size_vec1), sum_M1 = get_sum_vec(M1,size_vec1), sum_D2 = get_sum_vec(D2,size_vec2), sum_M2 = get_sum_vec(M2,size_vec2);
  int pseudo_sum_M1, pseudo_sum_M2;
  //printf(" size_vec1=%d ",size_vec1 );

  if(size_vec1==0 && size_vec2 ==0)
      res =1;
  else if(size_vec1==0)
    {
      subres = 1;
      if(sum_D2+sum_M2 == 0)
	{
	  subres = 1;
	  found = 0;
	  for(j=0;j<size_vec2 && !found;j++)
	    {
	      if(K2[j]== get_no_tips(node2[j]) && K2[j]<= 2)
		{
		      //printf(" [j=%d, K2[j]=%d]\n ",j,K2[j]);
		  subres *=1 ;
		}
	      else if(K2[j]==get_no_tips(node2[j]) && K2[j]>2)
		{
		  subres *=  exp(log_noTrees_onShape(K2[j], node2[j]));
		  //printf("subres = %f ",subres); 
		  //printf("subres = %f ",exp(log_noTrees_onShape(K2[j], node2[j]))); 
		}
	      else
		found = 1;
	    }
	  if(found == 1)
	    res = 0;
	  else 
	    res = subres;
	}
      else
	res = 0;
    }
  else if(size_vec2==0)
    {
      subres = 1;
      if(sum_D1+sum_M1 == 0)
	{
	  for(i=0;i<size_vec1;i++)
	    if(K1[i] != get_no_tips(node1[i]))
	      subres *= 0;
	  if(subres == 0)
	    res = 0;
	  else
	    res = 1;
	}
      else
	res = 0;      
    }
  else if(sum_D1+sum_M1==0 && sum_D2+sum_M2==0)
    {
      ////////////////////////
      ////     Theorem 1    //
      ///////////////////////
      //printf("yes");
      found=0;
      i=0,j=0;
      for(i=0;i<size_vec1 && !found;i++)
	if(K1[i]!=get_no_tips(node1[i]))
	  found = 1;
      if(found)
	{
	  //printf("[CASE 1-1] \n");
	  res = 0;
	  //printf("res = %d ",res); 
	}
      else
	{	  
	  subres=1.0;
	  found = 0;
	  //printf("[case 1-3:: size_vec2 = %d]\n ",size_vec2);
	  
	  for(j=0;j<size_vec2 && !found;j++)
	    {
	      if(K2[j]== get_no_tips(node2[j]) && K2[j]<= 2)
		{
		  //printf(" [j=%d, K2[j]=%d]\n ",j,K2[j]);
		  subres *=1 ;
		}
	      else if(K2[j]==get_no_tips(node2[j]) && K2[j]>2)
		{
		  subres *=  exp(log_noTrees_onShape(K2[j], node2[j]));
		  //printf("subres = %f ",subres); 
		  //printf("subres = %f ",exp(log_noTrees_onShape(K2[j], node2[j]))); 
		}
	      else
		found = 1;
	    }
	  if(found)
	    {
	      //printf("[CASE 1-2]\n");
	      res = 0;
	      //printf("res = %d\n ",res); 
	    }
	  else
	    {
	      //printf("[CASE 1-3] ");
	      //printf("subres=%f  ",subres);
	      res = subres;
	      //printf("res = %f ",res); 
	    }
	}      
    }
  else
    {
      ////////////////////////
      ////     Theorem 2    //
      ///////////////////////
      pseudo_sum_M1 = sum_M1;
      pseudo_sum_M2 = sum_M2;
      if(size_vec1 ==2 & M1[1]==1)
	if(node1[0]->parent->id == node1[1]->parent->id && is_doubleTruncated(node1[0]->parent)) // children of pseudo-root
	      pseudo_sum_M1 --;
      if(size_vec2 ==2 & M2[1]==1)
	if(node2[0]->parent->id == node2[1]->parent->id && is_doubleTruncated(node2[0]->parent)) // children of pseudo-root
	  pseudo_sum_M2 --;
      if(sum_D1+pseudo_sum_M1 !=  sum_D2+pseudo_sum_M2)
	{
	  //printf("[CASE 1-4]\n");
	  res = 0.0;
	}
      else
	{
	  ////////////////////////
	  ////     Theorem 3    //
	  ///////////////////////
	  //printf("here");
	  found = 0;
	  for(i=0;i<size_vec1 && !found;i++)
	    {
	      no_child1 = no_internalchildrenNodes(node1[i]);
	      //printf("K1[i]=%d get_no_tips(node1[i])=%d ",K1[i],get_no_tips(node1[i]));
	      //if((no_child1==0 && D1[i]>0) || (no_child1<D1[i] && K1[i]<D1[i]-no_child1 && M1[i]==0) || (D1[i]>0 && K1[i]>get_no_tips(node1[i])-2) || (D1[i]==0 && M1[i]==0 && K1[i]!=get_no_tips(node1[i])) || (M1[i]==0 && K1[i] < get_no_directtips(node1[i])) )
	      if((no_child1==0 && D1[i]>0) || (D1[i]>0 && K1[i]>get_no_tips(node1[i])-2) || (D1[i]==0 && M1[i]==0 && K1[i]!=get_no_tips(node1[i])) || (M1[i]==0 && K1[i] < get_no_directtips(node1[i])) || D1[i] >get_no_tips(node1[i])-2 )
		found = 1;
	    }     
	  if(found)
	    {
	      //printf("[CASE 1-5] ");
	      res = 0;
	      //printf("res = %f\n ",res); 
	    }
	  else
	    {
	      for(j=0;j<size_vec2 && !found;j++)
		{
		  no_child2 = no_internalchildrenNodes(node2[j]);
		  //if( (no_child2==0 && D2[j]>0) || (no_child2<D2[j] && K2[j]<D2[j]-no_child2 && M2[j]==0) || (D2[j]>0 && K2[j]>get_no_tips(node2[j])-2) || (D2[j]==0 && M2[j]==0 && K2[j]!=get_no_tips(node2[j])) || (M2[j]==0 && K2[j] < get_no_directtips(node2[j])) )
		  if( (no_child2==0 && D2[j]>0)  || (D2[j]>0 && K2[j]>get_no_tips(node2[j])-2) || (D2[j]==0 && M2[j]==0 && K2[j]!=get_no_tips(node2[j])) || (M2[j]==0 && K2[j] < get_no_directtips(node2[j])) || D2[j] >get_no_tips(node2[j])-2 )
		    {
		      //printf("yes");
		      found = 1;
		    }
		}
	      if(found)
		{
		  //printf("[CASE 1-6]");
		  res = 0;
		  //printf("res = %f\n ",res); 
		  
		}
	      else
		{	      
		  if(sum_M1 == 0)
		    {
		      found = -1;
		    }
		  else
		    {
		      found = get_node_toBeNamed(size_vec1,node1,D1,M1); // found an edge which is the element of edge vector in tree T		
		      /*     
		      if(found >= 0)
			printf("[case 1-7 :: found=%d node1->id = %d]\n ",found, node1[found]->id);
		      */
		    }
		  if(sum_M1 == 0 || found < 0)
		    {		  
		      //////////////////////////////////////////
		      //      Theorem 5                       //
		      //   Formular dismantling a node in T   //
		      //////////////////////////////////////////
		      // Dismantle a node in T
		      
		      //res = 10000;
		      //printf("[CASE 1-7] ");
		      res = R_extend_vec1(node1,node2,D1,D2,K1,K2,M1,M2,F,size_F);
		      //printf("[CASE 1-7]res = %f  ",res); 
		    }
		  else if(sum_M1 !=0 && sum_M2 == 0) // sum_M1 = 0
		    {	  
		      //////////////////////////////////////////
		      //      Theorem 6                       //
		      //   Formular dismantling a node in T   //
		      //////////////////////////////////////////
		      // Dismantle a node in T'
		      
		      //printf("[CASE 1-8]\n");
		      res = R_extend_vec2(node1,node2,D1,D2,K1,K2,M1,M2,F,size_F);
		      //printf("[CASE 1-8]res = %f  ",res); 
		    }
		  else //possibly extension and split
		    {
		      // if found >=0, it indicates the first element of the edge vector in tree T
		      // if found2 >=0, it indicates the clade to have the same name with found in tree T'
		      //printf("split_test\n");
		      /*
		      printf("[start split_test T2:: ");
		      for(j=0;j<size_vec2;j++)
			printf("n2->id =%d (%d,%d,%d) ",node2[j]->id,D2[j],K2[j],M2[j]);
		      printf("]\n");
		      */
		      found2 = split_test(size_vec1,size_vec2,found,node1,node2,D1,D2,K1,K2,M1,M2,F,size_F);
		      /*
		      printf("[after split_test T2:: ");
		      for(j=0;j<size_vec2;j++)
			printf("n2->id =%d (%d,%d,%d) ",node2[j]->id,D2[j],K2[j],M2[j]);
		      printf("]\n");
		      */
		      //printf("found =%d  found2 = %d\n",found,found2);
		      if(found2 == -2)
			res = 0;
		      else if(found2 == -1)
			{		    
			  // Dismantle a node in T'
			  
			  //printf("[CASE 1-9]\n");
			  res =R_extend_vec2(node1,node2,D1,D2,K1,K2,M1,M2,F,size_F);
			}
		      else if(found2 >=0 && found >=0)  // Factorize R_joint: found an edge which is the element of edge vector in tree T
			{		  
			  //////////////////////////////////////////
			  //      Theorem 7                       //    	  
			  //////////////////////////////////////////
			  //printf("[CASE 1-10] R_split_root\n");
			  /*
			   printf("[start CASE 1-10 T1:: ");
			   for(j=0;j<size_vec1;j++)
			     printf("n1->id =%d (%d,%d,%d) ",node1[j]->id,D1[j],K1[j],M1[j]);
			   printf("]\n");			   
			   printf("[start CASE 1-10 T2:: ");
			   for(j=0;j<size_vec2;j++)
			     printf("n2->id =%d (%d,%d,%d) ",node2[j]->id,D2[j],K2[j],M2[j]);
			   printf("]\n");
			  //*/
			  //printf("[case 1-10::R_split_root:: size_vec1 = %d, size_vec2 = %d, found1 = %d found2 = %d nodeID = %d, nodeID2 = %d] \n", size_vec1, size_vec2,found,found2,node1[found]->id,node2[found2]->id);
			  res = R_split_root(size_vec1,size_vec2,found,found2,node1,node2,D1,D2,K1,K2,M1,M2,F,size_F);
			  //printf("[CASE 1-10::R_split_root::  nodeID = %d, nodeID2 = %d:: res = %f]  \n",node1[found]->id,node2[found2]->id,res); 
			  
			}
		      else // none of parent edges of the given nodes are the first element of edge vector in tree T
			{
			  found = get_node_toExtend(size_vec1,node1,D1,M1);
			  if(found >= 0) // found a node to extend
			    {
			      res =100000;//printf("[CASE 1-11]");
			    }
			  else // None of nodes are able to extend
			    {
			      //printf("[CASE 1-12]");
			      res =1000000;
			    }
			  
			  //printf("here");
			  
			}
		    }
		}
	    }
	}  
    }
  //  if(res>0) 
  //  printf("R_joint::res = %f\n ",res); 
  return res;  
}
int in_F(int idx, int F[HALF_MAX_TAXA],int size_F)
{
  int i, found=0;
  
  if(size_F==0)
    return 0;
  else{
    for(i=0;i<size_F && !found;i++)
      if(idx == F[i])
	found = 1;
    
    if(found)
      return 1;
    else
      return 0;  
  }
}

int isSibling(struct node *node1,struct node *node2)
{
  if(node1->parent->id == node2->parent->id)
    return 1;
  else
    return 0;
}

int isDescendant(struct node *parent,struct node *candi_child)
{
  if(parent->terminal==1)
    return 0;
  else if(candi_child->parent == NULL)
    return 0;
  else
    {
      if(parent->id == candi_child->parent->id)
	return 1;
      else 
	return isDescendant(parent,candi_child->parent);
    }
}


int collectTips(struct node *vector_nodes[], struct node *node, int size)
{
  int i,res,find;
  struct node *child1 = node->listnodes[0], *child2 = node->listnodes[1];
  
  //printf("node->id=%d\n",node->id);
  find = 0;
  for(i=0;i<size && find==0;i++)
    if(node->id == vector_nodes[i]->id)
      {
	res = 0;
	find = 1;
      }
  if(!find)
    if(child1->terminal == 1 && child2->terminal == 1)
      {
	// this means that the node is an internal and not in vector_nodes[], so we should not collect tips belong to the node.
	res = 0; 
      }
    else
      {
	if(child1->terminal ==1)
	  {
	    res = 1+collectTips(vector_nodes,child2,size);
	  }
	else if(child2->terminal ==1)
	  {
	    res = 1+collectTips(vector_nodes,child1,size);
	  }
	else
	  {
	    res = collectTips(vector_nodes,child1,size)+ collectTips(vector_nodes,child2,size);
	  }
      }
  
  //printf("res = %d\n",res);
  return res;
}


int isMaximalAntichain(struct node *vector_nodes[], struct node *node, int size)
{
  // nodes in vector_nodes should be antichain
  int i,reach,ndirect_tips,total_ntips;
  struct node *givenNodes[size];

  ndirect_tips = collectTips(vector_nodes,node,size);
  total_ntips = get_no_tips(node)-ndirect_tips;
  //printf("node->id=%d ndirect_tips=%d, total_ntips = %d\n",node->id,ndirect_tips,total_ntips);
  if(total_ntips==0)
    return 1;
  else
    {
      for(i=0;i<size;i++)
	if(isDescendant(node,vector_nodes[i]) || node->id == vector_nodes[i]->id)
	  {
	    total_ntips -= get_no_tips(vector_nodes[i]);
      }
      //printf("total_ntips = %d\n",total_ntips);
      if(total_ntips == 0)
	return 1;
      else
	return 0;
    }

  /*
  for(i=0;i<size;i++)
    givenNodes[i] = vector_nodes[i];
  printf("isMaximalAntichain::node->id = %d total_ntips = %d\n",node->id,total_ntips);
  for(i=0;i<size;i++)
    {
      if(givenNodes[i]->id == node->id)
	total_ntips -= get_no_tips(node)-get_no_directtips(node);
      else if(isDescendant(node,givenNodes[i]))
	{
	  printf("isMaximalAntichain::isDescendant:: i=%d node->id = %d, vector_nodes[i]->id=%d\n",i,node->id, vector_nodes[i]->id);
	  total_ntips -= get_no_tips(givenNodes[i]);
	  printf("isMaximalAntichain::isDescendant:: i=%d givenNodes[i]->id=%d total_ntips = %d\n",i,givenNodes[i]->id,total_ntips);

	  reach =0;
	  while(!reach)
	    {	      
	      if(givenNodes[i]->parent->id != node->id)
		{
		  total_ntips -= get_no_directtips(givenNodes[i]->parent);
		  printf("isMaximalAntichain::isDescendant:: i=%d givenNodes[i]->parent->id=%d total_ntips = %d\n",i,givenNodes[i]->parent->id,total_ntips);
		  givenNodes[i] = givenNodes[i]->parent;
		}
	      else
		reach = 1;
	    }
	}
    }
   printf("isMaximalAntichain::atTheEnd::node->id = %d total_ntips = %d\n",node->id,total_ntips);
  if(total_ntips == 0)
    return 1;
  else
    return 0;
  */
}




struct node *find_SymmSibling_anc(struct node *node1, struct node *node2)
{
  // find a node (ancestor of node1 or itself) which is a symmetric sibling of an ancestor node of node2
  
  int ntips1,ntips2;
  struct node *symm;
  
  //printf("Welcome to findMRCA()\n");
  if(node1->parent == NULL)
    symm = NULL;
  else if(node2->parent == NULL)
    symm = NULL;
  else if(node1->parent->id == node2->parent->id)
    {
      if(isSymmetric(node1,node2))
	symm = node1;      
      else
	symm = NULL;
      
    }
  else
    {
      ntips1 = get_no_tips(node1);
      ntips2 = get_no_tips(node2);
      if(ntips1 == ntips2)
	symm = find_SymmSibling_anc(node1->parent,node2->parent);
      else if(ntips1 > ntips2)
	symm= find_SymmSibling_anc(node1,node2->parent);
      else
	symm = find_SymmSibling_anc(node1->parent,node2);
    }
  
  return symm;
}


int split_test(int size1, int size2, int index, struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA], int F[HALF_MAX_TAXA],int size_F)
{
  ////////////////////
  //    Theorem 7   //
  ////////////////////
  // test if there is a node in T' to be named same as the "index"th node in T
  // if res >= 0, then it returns the index of node to be named
  // if res = -1, then it means that there is no node to be named same as the node in T. Thus, we need to  dismantle a node in T'
  // if res = -2, then the corresponding value of function R_joint is zero (the later case is that )

  int i, j, k, found1, found2, res,temp;
  int index2[HALF_MAX_TAXA];
  //struct node *new_node1[HALF_MAX_TAXA],*new_node2[HALF_MAX_TAXA];
  int new_D1[HALF_MAX_TAXA],new_D2[HALF_MAX_TAXA],new_K1[HALF_MAX_TAXA],new_K2[HALF_MAX_TAXA],new_M1[HALF_MAX_TAXA],new_M2[HALF_MAX_TAXA];
  struct node *symm;

  int q;
  /*
  printf("[start inside of split_test T2:: ");
  for(q=0;q<size2;q++)
    printf("n2->id =%d (%d,%d,%d) ",node2[q]->id,D2[q],K2[q],M2[q]);
  printf("]\n");
  */
  //printf("yes\n");
  //printf("split_test:: index = %d ",index);
  //printf("size2 = %d\n ",size2);
  //printf("split_test:: size_F=%d\t",size_F);
  //if(size_F>0)    
  //printf("split_test:: F[0]=%d\n",F[0]);
  found1 = 0; found2 =0;
  for(i=0; i<size2 && !found1;i++)
    {
      //printf("i=%d node2[i]->id = %d\n",i,node2[i]->id);
      if( size_F==0 || !(in_F(i,F,size_F)) )
	{
	  //printf("split_test:: yes! i=%d\n ",i);
	  //printf("[split_test::T1:: i=%d M1[index]=%d], K1[index]=%d, D1[index]=%d, get_no_tips(node1[index])=%d ]\n",i, M1[index], K1[index],D1[index],get_no_tips(node1[index]));
	  //printf("[split_test::T2: i=%d M2[i]=%d], K2[i]=%d, D2[i]=%d, get_no_tips(node2[i])=%d]\n ",i, M2[i], K2[i],D2[i],get_no_tips(node2[i]));
	  if(M1[index]==1 && M2[i]==1 && K1[index]==0 && K2[i]==0 && D1[index]==D2[i] && get_no_tips(node1[index])==get_no_tips(node2[i]))
	    {
	      //printf(" [i=%d yes]\n ",i);
	      if(i==0) // found the edge able to have the same name
		{
		  //printf("[split_test:: matched_case1:: n1 = %d n2=%d]\n",node1[index]->id,node2[i]->id);
		  found1 = 1;
		  res = i;
		}
	      else //if(node2[i]->id != node2[i-1]->id)
		{
		  temp = 1;
		  for(j=0;j<i;j++)
		    {
		      //printf(" [i=%d yes]\n ",i);		  
		      //printf("node2[i]->id = %d\n",node2[i]->id);
		      //printf("[split_test::node2[j]->id=%d, node2[i]->id=%d]\n",node2[j]->id,node2[i]->id);
		      if(isSibling(node2[j],node2[i]))
			{
			  //printf(" [i=%d yes Sibling]\n ",i);
			  if(isSymmetric(node2[j],node2[i]))
			    temp *= 0;	
			}
		      else
			{	  
			  /*
			  printf("[middle1 of split_test T2:: i=%d j=%d :: ",i,j);
			  for(q=0;q<size2;q++)
			    printf("n2->id =%d (%d,%d,%d) ",node2[q]->id,D2[q],K2[q],M2[q]);
			  printf("]\n");
			  */
			  symm=find_SymmSibling_anc(node2[j],node2[i]);
			    
			  /*
			  printf("[middle2 of split_test T2:: i=%d j=%d :: ",i,j);
			  for(q=0;q<size2;q++)
			    printf("n2->id =%d (%d,%d,%d) ",node2[q]->id,D2[q],K2[q],M2[q]);
			  printf("]\n");
			  //*/
			  if(symm != NULL)
			    {			      
			      //printf("symm->id=%d i=%d node2[i]->id = %d\n ",symm->id,i, node2[i]->id);
			      if(isMaximalAntichain(node2,symm,i))
				{
				  temp *= 0;
				  //printf("yes! it is maximal antichain");
				}
			    }		  
			  /*
			  printf("[middle3 of split_test T2:: i=%d j=%d :: ",i,j);
			  for(q=0;q<size2;q++)
			    printf("n2->id =%d (%d,%d,%d) ",node2[q]->id,D2[q],K2[q],M2[q]);
			  printf("]\n");
			  */
			}		      
		    }
		  if(temp == 1)
		    {
		      //printf("yyyes!!\n");
		      found1 = 1;
		      res = i;
		    }		  
		}
	    }

	  if(!found1)
	    {
	      if(get_no_tips(node1[index])<get_no_tips(node2[i]) && M2[i]==0) // found the node containing a node to be dismantled
		if(i==0)
		  found2 = 1;
		else if(node2[i]->id != node2[i-1]->id)
		  found2 = 1;	  
	    }
	}
    }
  
  if(found1 ==0 && found2 == 1)
    res = -1;
  else if(found1 == 0)
    res = -2;

  //printf(" R_split_root::res=%f ",res);
  return res;
}


long double R_split_root(int size1, int size2, int index1, int index2,struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F)
{
  ////////////////////////////
  // Theorem 7
  ////////////////////////
  int i, j, k, found;
  long double res=0, sub_res1, sub_res2,sub_res3;
  struct node *new_node1[HALF_MAX_TAXA],*new_node2[HALF_MAX_TAXA];
  int new_D1[HALF_MAX_TAXA],new_D2[HALF_MAX_TAXA],new_K1[HALF_MAX_TAXA],new_K2[HALF_MAX_TAXA],new_M1[HALF_MAX_TAXA],new_M2[HALF_MAX_TAXA], new_F[HALF_MAX_TAXA];
  int new_sizeF;
  res=0;
  //printf("yes");
  //printf("[R_split_root:: index1 = %d\t index2 = %d]\n ",index1,index2);
  //printf("size2 = %d\n ",size2);
  /*
  printf("[start R_split_root T1:: ");
  for(j=0;j<size1;j++)
    printf("n1->id =%d (%d,%d,%d) ",node1[j]->id,D1[j],K1[j],M1[j]);
  printf("]\n");

  printf("[start R_split_root T2:: ");
  for(j=0;j<size2;j++)
    printf("n2->id =%d (%d,%d,%d) ",node2[j]->id,D2[j],K2[j],M2[j]);
  printf("]\n");
  //  */
  
  for(i=0;i<size_F;i++)
    new_F[i] = F[i];  
  new_sizeF = 0;
  new_node1[0] = node1[index1]; new_node1[1] = NULL; //new_nodel1[0]->parent = NULL; 
  new_node2[0] = node2[index2]; new_node2[1] = NULL; //new_nodel2[0]->parent = NULL;  
  new_D1[0] = D1[index1]; new_D2[0] = D2[index2];
  new_M1[0] = 0; new_M2[0] = 0;
  sub_res1 = 0;  
  //  printf("[R_split_root::factor1::node1->id=%d new_D1[0]=%d new_M1[0]=%d node2->id=%d new_D2[0]=%d new_M2[0]=%d ]\n",new_node1[0]->id, new_D1[0], new_M1[0],new_node2[0]->id, new_D2[0],new_M2[0]);
  for(k=0;k<=get_no_tips(node2[index2]);k++)
    {	  
      new_K1[0] = k; new_K2[0] = k;
      sub_res1 += R_joint(new_node1, new_node2, new_D1, new_D2, new_K1, new_K2, new_M1, new_M2,new_F,new_sizeF);
    }
  //printf("R_split_root:: sub_res1 = %f\n",sub_res1);

  // construct new node1, D1, M1 and K1 which does not include the "index"th node, d, k and m.
  //printf(" size1 = %d ",size1);
  if(size1>1)
    {
      for(j=0;j<size1;j++)
	{
	  if(j<index1)
	    {
	      new_node1[j] = node1[j]; 
	      new_D1[j] = D1[j]; new_M1[j] = M1[j]; new_K1[j] = K1[j]; 
	    }
	  else if(j>index1)
	    {
	      new_node1[j-1] = node1[j];
	      new_D1[j-1] = D1[j]; new_M1[j-1] = M1[j]; new_K1[j-1] = K1[j]; 
	    }
	}
    }
  new_node1[size1-1]=NULL;
  /*
  if(size1>1)
    for(i=0;i<size1-1;i++)
      printf("R_split_root::sub_res2:: new_node1->ID = %d new_D1 = %d new_K1 = %d new_M1 = %d ",new_node1[i]->id,new_D1[i],new_K1[i],new_M1[i]);
  else
    printf("R_split_root::sub_res2:: size1 = 0 ");
  printf("\n");
  */
  // construct new node2, D2, M2 and K2 which does not include the "i"th node, d, k and m.
  if(size2>1)
    {
      for(j=0;j<size2;j++)
	{
	  if(j<index2)
	    {
	      new_node2[j] = node2[j];
	      new_D2[j] = D2[j]; new_M2[j] = M2[j]; new_K2[j] = K2[j]; 
	    }
	  else if(j>index2)
	    {
	      new_node2[j-1] = node2[j];
	      new_D2[j-1] = D2[j]; new_M2[j-1] = M2[j]; new_K2[j-1] = K2[j]; 
	    }
	}
    }
  new_node2[size2-1]=NULL; 
  /*
  if(size2>1)
    for(i=0;i<size2-1;i++)
      printf("R_split_root::sub_res2:: new_node2->ID = %d new_D2 = %d new_K2 = %d new_M2 = %d ",new_node2[i]->id,new_D2[i],new_K2[i],new_M2[i]);
  else
    printf("R_split_root::sub_res2:: size2 = 0 ");
  printf("\n");
  */
  new_sizeF = 0;
  sub_res2 = R_joint(new_node1, new_node2, new_D1, new_D2, new_K1, new_K2, new_M1, new_M2,new_F,new_sizeF);
  //printf("R_split_root:: sub_res2 = %f\n",sub_res2);

  new_F[size_F]  = index2;
  new_sizeF = size_F+1;
  sub_res3 = R_joint(node1, node2, D1, D2, K1, K2, M1, M2, new_F, new_sizeF);
  //printf("R_split_root:: sub_res3 = %f\n",sub_res3);
  
  res = sub_res1*sub_res2 + sub_res3;
  //printf("R_split_root::res = sub_res1*sub_res2+sub_res3 %f\n",res);

  /*
  for(i=0; i<size2;i++)
    {
      if( size_F==0 || !(in_F(i,F,size_F)) )
	{
	  if(M1[index]==1 && M2[i]==1 && K1[index]==0 && K2[i]==0 && D1[index]==D2[i] && get_no_tips(node1[index])==get_no_tips(node2[i]))
	    {
	      //printf("yes");
	      if(i==0) // found the edge able to have the same name
		{
		  //printf("i=%d ",i);
		  //printf("[split-1] i=%d\t",i);
		  res += subR_split(size1, size2, index, i, node1, node2, D1, D2, K1, K2, M1, M2,F,size_F);
		  F[new_sizeF] = i;
		  new_sizeF++;
		  //printf("R_split_root::res = %f\n ",res);     
		}
	      else if(node2[i]->id != node2[i-1]->id)
		{
		  //printf("node2[i]->id = %d node2[i-1]->id = %d ",node2[i]->id, node2[i-1]->id);
		  if(node2[i]->listnodes[0]->id != node2[i-1]->id && node2[i]->listnodes[1]->id != node2[i-1]->id)
		    {
		      //printf("[split-2] i=%d\t",i);
		      //printf("node2[i]->id = %d node2[i-1]->id = %d ",node2[i]->id, node2[i-1]->id);
		      // not symmetric siblings because nodes are children of the root. That is, node2 contains siblings only.
		      res += subR_split(size1, size2, index, i, node1, node2, D1, D2, K1, K2, M1, M2,F,size_F);
		      F[new_sizeF] = i;
		      new_sizeF++;	
		    }
		}
	    }      
	}
      // else	
//	 printf("i=%d no\n",i); 
    }
  //printf(" R_split_root::res=%f ",res);

  if(new_sizeF != size_F)
    res += R_joint(node1, node2, D1, D2, K1, K2, M1, M2, F, new_sizeF);
  */

  return res;
}


int get_node_toExtend(int size, struct node *node[HALF_MAX_TAXA],int D[HALF_MAX_TAXA],int  M[HALF_MAX_TAXA])
{
  int res = 0, i, found = 0, sum_d;
  
  if(M[0] == 0 & D[0] >0)
    res = 0;
  else
    {
      sum_d = 0;
      for(i=1;i<size && !found; i++)
	{
	  if(M[i-1]==0)
	    sum_d += D[i-1];
	  if(M[i] == 0 && D[i] >0 && sum_d ==0)
	    found = 1;	    
	}
      if(found)
	res = i;
      else
	res = -i; // None of the given nodes are able to extend.
    }
  
  return res;
}

int get_node_toBeNamed(int size, struct node *node[HALF_MAX_TAXA],int  D[HALF_MAX_TAXA],int  M[HALF_MAX_TAXA])
{
  ///////////////////////
  // Theorem 5, 7
  //////////////////////
  // return index of node whose parent edge is the first element of edge vector.
  // returned value = -1 if none of parent edges of the given nodes could be the first element.
  int i, res, found, stop;
  int sum_d_m;
  
  if(M[0]==1)
    res = 0;
  else
    {
      sum_d_m = 0;
      found = 0; stop = 0;
      for(i=1;i<size && !found && !stop;i++)
	{
	  sum_d_m += D[i-1]+M[i-1];
	  if(sum_d_m > 0)
	    stop = 1;
	  else if(M[i]==1 && sum_d_m == 0)
	    found = 1;
	}
      if(found)
	res = i-1;
      else if(stop)
	res = -1;
      else
	res = -1; // none of parent edges of the given nodes could be the first element.
    }
  //printf("get_node_toBeNamed::res=%d\n",res);
  return res;
}

int is_doubleTruncated(struct node *node)
{
  // The given node is the root of a tree and it has "two" children that are internal nodes
 
  // double truncated if #total tips = #tips below the first child node + #tips below the 2nd child node
  if(get_no_tips(node)==get_no_tips(node->listnodes[0])+get_no_tips(node->listnodes[1]) && node->parent ==NULL)
    return 1; 
  else
    return 0;
   

}

long double R_extend_vec1(struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F)
{
  // To extend vector node1 along tree T when node1 contains the root only.

  int idx; // the index of the node to be dismantled
  int i,no_child,size_vec=getsize_vectorNode(node1),sum_d_m, sum_k, found, idx_firstChild;
  int sub_d[3], sub_m[3], sub_k[3], target_D, target_M, target_K; 
  struct node *new_node[HALF_MAX_TAXA], *sub_node[3];
  struct node *target_node;
  int new_D[HALF_MAX_TAXA],new_K[HALF_MAX_TAXA], new_M[HALF_MAX_TAXA];
  long double res=0.0;
  
  found = 0;
  for(i=0;i<size_vec && !found;i++)
    {
      if(D1[i]>0 && M1[i] == 0)
	{
	  found = 1;
	  idx = i;
	}
    }
  target_node = node1[idx];
  target_D = D1[idx]; target_M = M1[idx]; target_K = K1[idx];
  //printf("[R_extend_vec1:: idx = %d, size_vec = %d] ",idx,size_vec);
  //printf("M1[1] = %d ",M1[1]);

  // D1[0] > 0
  no_child = no_internalchildrenNodes(target_node); // no_child = 1,2,3
  // printf("no_child = %d\t",no_child);

  if(no_child ==1) // when no_child=1
    {  
      found = 0;
      for(idx_firstChild=0;idx_firstChild<3 && !found; idx_firstChild++)
	if(!(target_node->listnodes[idx_firstChild]->terminal))
	  found = 1;
      if(found)
	{
	  sub_node[0] = target_node->listnodes[idx_firstChild-1];
	  sub_node[0]->parent = target_node;
	}
      //printf("idx_firstChid=%d",idx_firstChild);
      for(sub_m[0]=0;sub_m[0]<=1;sub_m[0]++)
	{
	  sub_d[0] = target_D-sub_m[0];
	  sub_m[0] = sub_m[0];
	  sub_k[0] = target_K-no_directTips(target_node);
	  if((sub_m[0] == 1 && sub_k[0]==0) || (sub_m[0] ==0 && sub_k[0]>=0))
	    {
	      for(i=0;i<size_vec;i++)
		{
		  if(i<=idx-1)
		    {
		      new_node[i] = node1[i];
		      new_D[i] = D1[i]; new_K[i] = K1[i]; new_M[i] = M1[i];
		    }
		  else if(i==idx)
		    {
		      new_node[i] = sub_node[0];
		      new_D[i] = sub_d[0]; new_K[i] = sub_k[0]; new_M[i] = sub_m[0];
		    }
		  else
		    {
		      new_node[i] = node1[i];
		      new_D[i] = D1[i]; new_K[i] = K1[i]; new_M[i] = M1[i];
		    }
		}
	      
	      new_node[size_vec] = NULL;
	      /*
	      printf("[CASE 2-1 T1: ");
	      for(i=0;i<size_vec;i++)
		printf("n1->id=%d (%d %d %d)  ",new_node[i]->id,new_D[i],new_K[i],new_M[i]);
	      printf("]\n");
	      */
	      //printf(" [CASE 2-1] ");
	      //printf("[CASE 2-1 T1::d0=%d m0=%d k0=%d]\n",sub_d[0],sub_m[0],sub_k[0]);
	      res += R_joint(new_node,node2,new_D,D2,new_K,K2,new_M,M2,F,size_F);
	    }
	}
    }  
  else // no_child = 2,3
    { 
      found = 0;
      for(idx_firstChild=0;idx_firstChild<3 && !found ;idx_firstChild++)
	if(target_node->listnodes[idx_firstChild]->terminal==0)
	  {
	    //printf("idx_firstChild = %d ",idx_firstChild);
	    found = 1;
	  }
      //printf("idx_firstChild = %d ",idx_firstChild);
      
      if(found)
	{
	  sub_node[0] = target_node->listnodes[idx_firstChild-1];
	  sub_node[0]->parent = target_node;
	  sub_node[1] = target_node->listnodes[idx_firstChild];
	  sub_node[1]->parent = target_node;
	  //printf(" index_firstChild = %d ",idx_firstChild);
	  //printf(" #tips(new_node[0])=%d #tips(new_node[1])=%d  ",get_no_tips(new_node[0]),get_no_tips(new_node[1]) );
	}
      for(sub_d[0]=target_D;sub_d[0]>=0;sub_d[0]--)
	for(sub_m[0]=0;sub_m[0]<=1;sub_m[0]++)
	  if(sub_m[0]+sub_d[0]<=target_D)
	    {	    
	      for(sub_k[0]=0;sub_k[0]<=target_K-no_directTips(target_node);sub_k[0]++)
		if((sub_m[0] == 1 && sub_k[0]==0) || (sub_m[0] ==0 && sub_k[0]>=0))
		  {
		    sum_d_m = sub_d[0]+sub_m[0];
		    sum_k = sub_k[0];
		    
		    if(no_child == 2) // When no_child = 2
		      {	      
			sub_node[2] = NULL;
			for(sub_m[1]=0;sub_m[1]<=1;sub_m[1]++)
			  {
			    found = 0;
			    if(is_doubleTruncated(node1[0]))
			      {
				// the given node is pseudo-node. Thereby, sub_m[1] from the root is always 0
				if(sub_m[1] !=sub_m[0] && node1[0]->parent == NULL) 
				  found = 1;
			      }
			    if(is_doubleTruncated(node1[0]) && node1[0]->parent == NULL)
			      {
				sub_d[1] = target_D-sum_d_m; 
				if(sub_m[1]==1)
				  sub_k[1]=0;
				else
				  sub_k[1] = target_K-no_directTips(target_node)-sub_k[0];
			      }
			    else
			      {
				sub_d[1] =  target_D-sum_d_m-sub_m[1]; 
				sub_k[1] = target_K-no_directTips(target_node)-sub_k[0];
			      }
			      
			    //res = 1111111;
			    //printf("K1[0]=%d #directTips=%d  ",K1[0],no_directTips(node1[0]));
			    //printf("[T1::d0=%d m0=%d k0=%d]",new_D[0],new_M[0],new_K[0]);
			    //printf("new_D[1]=%d\n", new_D[1]);			    
			    //printf("yes");
			    sub_k[1] = target_K-no_directTips(target_node)-sub_k[0];
			    if(found == 0 && sub_d[1] >= 0 && (sub_m[1]==0 || sub_m[1]==1) && ((sub_m[1] == 1 && sub_k[1]==0) || (sub_m[1] ==0 && sub_k[1]>=0)))
			      {		
				//printf("yes");
				/*
				if(identicalShape(sub_node[0],sub_node[1]))			    
				  sub_node[1] = sub_node[0];
				*/
				//printf("size_vec=%d\n",size_vec);
				for(i=0;i<size_vec+1;i++)
				  {
				    if(i<=idx-1)
				      {
					new_node[i] = node1[i];
					new_D[i] = D1[i]; new_K[i] = K1[i]; new_M[i] = M1[i];
				      }
				    else if(i>=idx && i<=idx+1)
				      {
					new_node[i] = sub_node[i-idx];
					new_D[i] = sub_d[i-idx]; new_K[i] = sub_k[i-idx]; new_M[i] = sub_m[i-idx];
				      }
				    else
				      {
					//printf("yes ");
					new_node[i] = node1[i-1];
					new_D[i] = D1[i-1]; new_K[i] = K1[i-1]; new_M[i] = M1[i-1];
				      }
				  }
				new_node[size_vec+1] = NULL;
				/*
				printf("[CASE 2-2 T1: ");
				for(i=0;i<size_vec+1;i++)
				  printf("n1->id=%d (%d %d %d)  ",new_node[i]->id,new_D[i],new_K[i],new_M[i]);
				printf("]\n");
				*/
				//printf("K1[0]=%d #tips=%d sum_k=%d sub_k[0]=%d sub_k[1]=%d\n ",K1[0],no_directTips(node1[0]),sum_k,sub_k[0],sub_k[1]);
				//printf("[CASE 2-2]");
				//printf("[d0=%d m0=%d d1=%d m1=%d]",sub_d[0],sub_m[0],sub_d[1],sub_m[1]);
				//printf("[CASE 2-2 T1::n1=%d d0=%d m0=%d k0=%d n2=%d d1=%d m1=%d k1=%d]\n",new_node[0]->id,new_D[0],new_M[0],new_K[0],new_node[1]->id,new_D[1],new_M[1],new_K[1]);
				res += R_joint(new_node,node2,new_D,D2,new_K,K2,new_M,M2,F,size_F);
			      }
			  }
		      }
		    else  // When no_child =3
		      {		      
			//printf("no_child = %d\t",no_child);
			sub_node[2] = target_node->listnodes[2];
			sub_node[2]->parent = target_node;
			
			//printf("sub_d[0]=%d sub_m[0]=%d\t",sub_d[0],sub_m[0]);
			//printf("sum_d_m=%d\t",sum_d_m);
			for(sub_d[1]=target_D-sum_d_m;sub_d[1]>=0;sub_d[1]--)
			  for(sub_m[1]=0;sub_m[1]<=1;sub_m[1]++)
			    if(target_D >= sub_d[0]+sub_m[0]+sub_d[1]+sub_m[1])
			      for(sub_k[1]=0;sub_k[1]<=target_K-no_directTips(target_node)-sub_k[0];sub_k[1]++)
				if((sub_m[1] == 1 && sub_k[1]==0) || (sub_m[1] ==0 && sub_k[1]>=0))
				  {
				    //printf("yes");
				    //printf("sub_d[0]=%d sub_d[1]=%d sub_m[0]=%d sub_m[1]=%d sub_k[0]=%d sub_k[1]=%d\n",sub_d[0],sub_d[1],sub_m[0],sub_m[1],sub_k[0],sub_k[1]);
				    sum_k = sub_k[0]+sub_k[1];
				    sum_d_m = sub_d[0]+sub_m[0]+sub_d[1]+sub_m[1];
				    
				    for(sub_m[2]=0;sub_m[2]<=1;sub_m[2]++)
				      {
					sub_d[2]=target_D-sum_d_m - sub_m[2];
					//printf("sub_d[2]=%d sub_m[2]=%d\n",sub_d[2],sub_m[2]);
					if(sub_d[2]>=0) // M1[0]=0
					  {		 
					    // printf("no_directTips(node1[0])=%d\t",no_directTips(node1[0]));
					    sub_k[2] = target_K-no_directTips(target_node)-sub_k[0] - sub_k[1];
					    //printf("sub_d[2]=%d sub_m[2]=%d sub_k[2]=%d\n",sub_d[2],sub_m[2], sub_k[2]);
					    if((sub_m[2] == 1 && sub_k[2]==0) || (sub_m[2] ==0 && sub_k[2]>=0))
					      {
						//printf("yes");
						sum_k = sub_k[0]+sub_k[1]+sub_k[2];
						//if(is_sameShape(new_node[1],new_node[2]))			    
						//  new_node[2] = new_node[1];

						for(i=0;i<size_vec+2;i++)
						  {
						    if(i<=idx-1)
						      {
							new_node[i] = node1[i];
							new_D[i] = D1[i]; new_K[i] = K1[i]; new_M[i] = M1[i];
						      }
						    else if(i>=idx || i<=idx+2)
						      {
							new_node[i] = sub_node[i-idx];
							new_D[i] = sub_d[i-idx]; new_K[i] = sub_k[i-idx]; new_M[i] = sub_m[i-idx];
						      }
						    else
						      {
							new_node[i] = node1[i-1];
							new_D[i] = D1[i-1]; new_K[i] = K1[i-1]; new_M[i] = M1[i-1];
						      }
						  }
						new_node[size_vec+2] = NULL;
						/*
						printf("[CASE 2-3 T1: ");
						for(i=0;i<size_vec+2;i++)
						  printf("n1->id=%d d0=%d m0=%d k0=%d  ",new_node[i]->id,new_D[i],new_M[i],new_K[i]);
						printf("]\n");			
						*/
						//printf(" [CASE 2-3] ");
						//printf("[CASE 2-3 T1::d0=%d m0=%d k0=%d d1=%d m1=%d k1=%d d2=%d m2=%d k2=%d]\n",new_D[0],new_M[0],new_K[0],new_D[1],new_M[1],new_K[1],new_D[2],new_M[2],new_K[2]);
						res += R_joint(new_node,node2,new_D,D2,new_K,K2,new_M,M2,F,size_F);
					      }
					  }
				      } 
				  }
		      }
		  }
	      
	    }
    }
  
  return res;
}

struct node *findMRCA(struct node *node1, struct node *node2)
{
  // node1 and node2 should be in the same tree
  
  int ntips1,ntips2;
  struct node *mrca;
  
  //printf("Welcome to findMRCA()\n");
  if(node1->parent == NULL)
    mrca = node1;
  else if(node2->parent == NULL)
    mrca = node2;
  else if(node1->parent->id == node2->id)
    mrca = node1->parent;
  else if(node1->id == node2->parent->id)
    mrca = node2->parent;
  else if(node1->parent->id == node2->parent->id)
    mrca = node2->parent;
  else
    {
      ntips1 = get_no_tips(node1);
      ntips2 = get_no_tips(node2);
      if(ntips1 == ntips2)
	mrca = findMRCA(node1->parent,node2->parent);
      else if(ntips1 > ntips2)
	mrca = findMRCA(node1,node2->parent);
      else
	mrca = findMRCA(node1->parent,node2);
    }
  return mrca;
}

int isSymmetric(struct node *node1, struct node *node2)
{
  struct node *mrca = findMRCA(node1,node2);
  int reachMRCA = 0, symmetric=1;
 
  while(!reachMRCA && symmetric)
    {
      if(identicalShape(node1,node2))
	{
	  if(node1->parent->id == mrca->id)
	    if(node2->parent->id == mrca->id)
	      reachMRCA = 1;
	    else
	      symmetric = 0;
	  else if(node2->parent->id == mrca->id)
	    symmetric = 0;
	  else
	    {
	      node1 = node1->parent;
	      node2 = node2->parent;
	    }
	}
      else
	symmetric = 0;
    }
  return symmetric;
}

long double R_extend_vec2(struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F)
{
  // To extend vector node1 along tree T' when node2 contains the root only.

  int idx; // the index of the node to be dismantled
  int found, i;
  int no_child,size_vec=getsize_vectorNode(node2),sum_d_m, sum_k,temp_k,idx_firstChild, isSameShape, symmetric;
  int sub_d[3], sub_m[3], sub_k[3]; 
  struct node *sub_node[3],*new_node[HALF_MAX_TAXA];
  int new_D[HALF_MAX_TAXA],new_K[HALF_MAX_TAXA], new_M[HALF_MAX_TAXA];
  int target_D,target_M,target_K;
  struct node *target_node;
  long double res=0.0, sub_res=0.0;
  long double term1,term1_1,term1_2,term1_3,term2;
  /*
  for(i=0;i<size_vec;i++)
    printf("start R_extend_vec2:: n2->id =%d (%d,%d,%d) ",node2[i]->id,D2[i],K2[i],M2[i]);
  printf("\n");
  */
  found = 0;
  for(i=0;i<size_vec && !found;i++)
    {
      if(D2[i]>0 && M2[i] == 0)
	{
	  found = 1;
	  idx = i;
	}
    }
  //  printf("[R_extend_vec2:: idx = %d, size_vec=%d]\n ",idx,size_vec);
  if(found == 0)
    {
      //printf("no nodes to be dismantled");
      res = 0;
    }
  else
    {
      target_node = node2[idx];
      target_D = D2[idx]; target_M = M2[idx]; target_K = K2[idx];
      
      // D2[0] > 0
      no_child = no_internalchildrenNodes(target_node); // no_child = 1,2,3
      //printf("start R_extend_vec2:: no_child = %d\n",no_child);
      if(no_child ==1) // no_child=1
	{ 
	  found = 0;
	  for(idx_firstChild=0;idx_firstChild<3 && !found; idx_firstChild++)
	    if(!(target_node->listnodes[idx_firstChild]->terminal))
	      found = 1;
	  if(found)
	    {
	      sub_node[0] = target_node->listnodes[idx_firstChild-1];
	      sub_node[0]->parent = target_node;
	    }
	  for(sub_m[0]=0;sub_m[0]<=1;sub_m[0]++)
	    {
	      sub_d[0] = target_D-sub_m[0];
	      sub_m[0] = sub_m[0];
	      sub_k[0] = target_K-no_directTips(target_node);
	      if((sub_m[0] == 1 && sub_k[0]==0) || (sub_m[0] ==0 && sub_k[0]>=0))
		{
		  if(sub_m[0] == 1)
		    {
		      term1 = 0;//log_factorial(target_K);
		    }
		  else
		    term1 = log_factorial(target_K)-log_factorial(sub_k[0]);
		  for(i=0;i<size_vec;i++)
		    {
		      if(i<=idx-1)
			{
			  new_node[i] = node2[i];
			  new_D[i] = D2[i]; new_K[i] = K2[i]; new_M[i] = M2[i];
			}
		      else if(i==idx)
			{
			  new_node[i] = sub_node[0];
			  new_D[i] = sub_d[0]; new_K[i] = sub_k[0]; new_M[i] = sub_m[0];
			}
		      else
			{
			  new_node[i] = node2[i];
			  new_D[i] = D2[i]; new_K[i] = K2[i]; new_M[i] = M2[i];
			}
		    }
		  new_node[size_vec] = NULL;
		  //printf(" [CASE 3-1] ");
		  //printf("[Case 3-1 T2::d0=%d m0=%d k0=%d]\n",sub_d[0],sub_m[0],sub_k[0]);
		  /*
		  printf("[CASE 3-1 T2::");
		  for(i=0;i<size_vec;i++)
		    printf("n2->id=%d d0=%d m0=%d k0=%d  ",new_node[i]->id,new_D[i],new_M[i],new_K[i]);
		  printf("]\n");
		  */
		  sub_res = log(R_joint(node1,new_node,D1,new_D,K1,new_K,M1,new_M,F,size_F))+term1;
		  res += exp(sub_res);
		  /*
		  if(res>0)
		    {
		      printf("[CASE 3-1 T2::");
		      printf("n2->id=%d d0=%d m0=%d k0=%d  ",new_node[0]->id,new_D[0],new_M[0],new_K[0]);
		      printf("]\n");
		     // for(i=0;i<size_vec;i++)
		//	printRoot(new_node[i]);
		      //printf("[case 3-1: term1 = %f] ",term1);
		      printf("[case 3-1: sub_res = %f\t res=%f]\n",sub_res,res);
		    }
		  */
		}
	    }
	}	
      else  // no_child = 2,3
	{      
	  found = 0;
	  for(idx_firstChild=0;idx_firstChild<3 && !found ;idx_firstChild++)
	    if(target_node->listnodes[idx_firstChild]->terminal==0)
	      {
		//printf("idx_firstChild = %d ",idx_firstChild);
		found = 1;
	      }
	  //printf("idx_firstChild = %d ",idx_firstChild);
	  
	  if(found)
	    {
	      sub_node[0] = target_node->listnodes[idx_firstChild-1];
	      sub_node[0]->parent = target_node;
	      sub_node[1] = target_node->listnodes[idx_firstChild];
	      sub_node[1]->parent = target_node;
	      //printf(" index_firstChild = %d ",idx_firstChild);
	      //printf(" #tips(new_node[0])=%d #tips(new_node[1])=%d  ",get_no_tips(new_node[0]),get_no_tips(new_node[1]) );
	    }
	  for(sub_d[0]=target_D;sub_d[0]>=0;sub_d[0]--)
	    for(sub_m[0]=0;sub_m[0]<=1;sub_m[0]++)
	      {
		//printf("sub_d[0]=%d sub_m[0]=%d\t",sub_d[0],sub_m[0]);
		if(sub_m[0]+sub_d[0]<=target_D)
		  {	   
		    //printf("yes\n");
		    for(sub_k[0]=0;sub_k[0]<=target_K-no_directTips(target_node);sub_k[0]++)
		      if((sub_m[0] == 1 && sub_k[0]==0) || (sub_m[0] ==0 && sub_k[0]>=0))
			{
			  //printf("yes\n");
			  sum_d_m = sub_d[0]+sub_m[0];
			  sum_k = sub_k[0];
			  //printf("new_K[0]=%d sub_k[0]=%d\t",new_K[0],sub_k[0]);
			  if(sub_m[0] == 0)
			    term1_1 = log_factorial(sub_k[0]);
			  else 
			    term1_1 = 0;
			  //printf("[CASE 3-2::term1 = %d] \t ",term1); 
			  
			  if(no_child ==2) // the given node is pseudo-node. Thereby, sub_m[1] is always 0
			    {
			      term2 = 1;
			      //printf("yes\n");
			      sub_m[1] = 0; 
			      for(sub_m[1]=0;sub_m[1]<=1;sub_m[1]++)
				{
				  if(is_doubleTruncated(target_node) && target_node->parent==NULL)
				    temp_k=0;
				  else
				    temp_k=sub_m[1];
				  if(temp_k == sub_m[1])
				    {
				      if(is_doubleTruncated(target_node) && target_node->parent==NULL && sub_m[0]==1)
					{
					  sub_m[1]=1;
					  sub_d[1] =  target_D-sum_d_m;					  
					  sub_k[1] = 0;
					}
				      else
					{
					  sub_d[1] =  target_D-sum_d_m-sub_m[1];				      
					  sub_k[1] = target_K-no_directTips(target_node)-sub_k[0];
					}
				      symmetric = isSymmetric(sub_node[0],sub_node[1]);
				      if(!(symmetric && sum_d_m==0 && sub_d[1]+sub_m[1]>0) )
					{
					  //isSameShape = identicalShape(sub_node[0],sub_node[1]);				  
					  //if(!(isSameShape && sum_d_m==0 && sub_d[1]+sub_m[1]>0) )
					  /*
					  if(identicalShape(sub_node[0],sub_node[1]))
					    {
					      sub_node[1] = sub_node[0];				  
					    }
					  */
					  if(sub_d[1] >= 0 && (sub_m[1]==0 || sub_m[1]==1) && ((sub_m[1] == 1 && sub_k[1]==0) || (sub_m[1] ==0 && sub_k[1]>=0)))
					    {
					      sub_node[2] = NULL;		 			  
					      if(sum_d_m+sub_d[1]+sub_m[1] == 0 && identicalShape(sub_node[0],sub_node[1]))
						term2 = 2;
					      else
						term2 = 0;
					      if(sub_m[1] == 0)
						term1_2 = log_factorial(sub_k[1]);
					      else
						term1_2 = 0;
					      term1 = log_factorial(target_K)-term1_1 - term1_2;
					      for(i=0;i<size_vec+1;i++)
						{
						  if(i<=idx-1)
						    {
						      new_node[i] = node2[i];
						      new_D[i] = D2[i]; new_K[i] = K2[i]; new_M[i] = M2[i];
						    }
						  else if(i>=idx && i<=idx+1)
						    {
						      new_node[i] = sub_node[i-idx];
						      new_D[i] = sub_d[i-idx]; new_K[i] = sub_k[i-idx]; new_M[i] = sub_m[i-idx];
						    }
						  else
						    {
						      new_node[i] = node2[i-1];
						      new_D[i] = D2[i-1]; new_K[i] = K2[i-1]; new_M[i] = M2[i-1];
						    }
						}
					      new_node[size_vec+1] = NULL;
					      
					      /*
					      printf("[CASE 3-2 T2::");
					      for(i=0;i<size_vec+1;i++)
						printf("n2->id=%d (%d %d %d)  ",new_node[i]->id,new_D[i],new_K[i],new_M[i]);
					      printf("]\n");
					      */

					      //printf("[CASE 3-2::term1 = %d] \n ",term1); 
					      
					      //printf("K1[0]=%d #tips=%d sum_k=%d sub_k[0]=%d sub_k[1]=%d\n ",K1[0],no_directTips(node2[0]),sum_k,sub_k[0],sub_k[1]);
					      //printf("[CASE 3-2]");
					      //printf("[d0=%d m0=%d d1=%d m1=%d]",sub_d[0],sub_m[0],sub_d[1],sub_m[1]);
					      //printf("[Case 3-2 T2::n1=%d d0=%d m0=%d k0=%d n2=%d d1=%d m1=%d k1=%d]\n",new_node[0]->id,new_D[0],new_M[0],new_K[0],new_node[1]->id,new_D[1],new_M[1],new_K[1]);
					      //printf("[CASE 3-2::term1 = %f term2 = %f] \n ",term1,term2); 
					      sub_res = log(R_joint(node1,new_node,D1,new_D,K1,new_K,M1,new_M,F,size_F))+term1-log_factorial(term2);
					      res += exp(sub_res);
					      //printf("[CASE 3-2::term1 = %f term2 = %f] \n ",term1,term2); 
					      //printf("[CASE 3-2::res = %f] \n ",res); 
					      
					      /*
					      if(res>0)
						{					      
						  printf("[CASE 3-2 T2::");
						  for(i=0;i<size_vec+1;i++)
						    printf("n2->id=%d (%d %d %d)  ",new_node[i]->id,new_D[i],new_K[i],new_M[i]);
						  printf("]\n");
						  //printf("[CASE 3-2::term1 = %f term2 = %f] \n ",term1,term2); 
						  printf("[CASE 3-2::res = %f] \n ",res); 
						}
					      */
					    }
					}
				    }
				}
			    }
			  else // no_child = 3
			    {      
			      // printf("yes\n");
			      //printRoot(node2[0]);
			      term2 = 1;
			      sub_node[2] = target_node->listnodes[2];
			      sub_node[2]->parent =target_node;
			      sum_d_m = sub_d[0]+sub_m[0];
			      
			      if(sub_m[0] == 0)
				term1_1 = log_factorial(sub_k[0]);
			      else 
				term1_1 = 0;
			      for(sub_d[1]=target_D-sum_d_m;sub_d[1]>=0;sub_d[1]--)
				for(sub_m[1]=0;sub_m[1]<=1;sub_m[1]++)
				  for(sub_k[1]=0;sub_k[1]<=target_K-no_directTips(target_node)-sub_k[0];sub_k[1]++)
				    {
				      if((sub_m[1] == 1 && sub_k[1]==0) || (sub_m[1] ==0 && sub_k[1]>=0))
					{	
					  isSameShape = identicalShape(sub_node[0],sub_node[1]);				  
					  if(!(isSameShape && sub_d[0]+sub_m[0]==0 && sub_d[1]+sub_m[1]>0) )
					    {
					      //printf("yes\n");
					      sum_d_m = sub_d[0]+sub_m[0]+sub_d[1]+sub_m[1];
					      /*
					      if(identicalShape(sub_node[0],sub_node[1]))
						{
						  //printf("yes\n");
						  sub_node[1] = sub_node[0];				  
						}
					      */
					      if(sub_m[1] == 0)
						term1_2 = log_factorial(sub_k[1]);	
					      else
						term1_2 = 0;
					      if(sum_d_m == 0 && isSameShape)
						term2 = 2;
					      else
						term2 = 0;				  
					      sub_k[2] = target_K-no_directTips(target_node)-sub_k[0] - sub_k[1];
					      //printf("new_K[0]=%d sub_k[0]=%d sub_k[1]=%d sub_k[2]=%d ",new_K[0],sub_k[0],sub_k[1],sub_k[2]);
					      //printf("no_directTips(node2[0])=%d ",no_directTips(node2[0]));
					      for(sub_d[2]=target_D-sum_d_m;sub_d[2]>=0;sub_d[2]--)
						for(sub_m[2]=0;sub_m[2]<=1;sub_m[2]++)
						  {	
						    if((sub_m[2] == 1 && sub_k[2]==0) || (sub_m[2] ==0 && sub_k[2]>=0))
						      {	
							isSameShape = identicalShape(sub_node[1],sub_node[2]);				  
							if(!(isSameShape && sub_d[1]+sub_m[1]==0 && sub_d[2]+sub_m[2]>0) )
							  {				
							    //printf("yes\n");
							    if(sub_m[2] == 0)
							      term1_3 = log_factorial(sub_k[2]);	
							    else
							      term1_3 = 0;
							    //printf("yes\n");
							    sum_d_m = sub_d[0]+sub_m[0]+sub_d[1]+sub_m[1]+sub_d[2]+sub_m[2];
							    //printf("yes\n");
							    //printf("new_node[0]->id=%d\n",new_node[0]->id);
							    //printf("new_node[1]->id=%d\n",new_node[1]->id);
							    //printf("new_node[2]->id=%d\n",new_node[2]->id);
							    //printf("new_node[0]->id=%d new_node[1]->id=%d new_node[2]->id=%d\n",new_node[0]->id,new_node[1]->id,new_node[2]->id);
							    //printf("yes\n");
							    if(term2 == 2)
							      {
								if(sum_d_m == 0 && identicalShape(sub_node[1],sub_node[2]))
								  term2 = 3;
							      }
							    else
							      {
								if(sum_d_m-sub_d[0]-sub_m[0]==0 &&  identicalShape(sub_node[1],sub_node[2]))
								  term2 = 2;
							      }
							    //printf("sum_d_m=%d\n",sum_d_m);
							    if(sum_d_m == target_D) // M2[0]=0
							      {
								//printf("yes\n");
								term1 = log_factorial(target_K)-term1_1-term1_2-term1_3;
								for(i=0;i<size_vec+2;i++)
								  {
								    if(i<=idx-1)
								      {
									new_node[i] = node2[i];
									new_D[i] = D2[i]; new_K[i] = K2[i]; new_M[i] = M2[i];
								      }
								    else if(i>=idx && i<=idx+2)
								      {
									new_node[i] = sub_node[i-idx];
									new_D[i] = sub_d[i-idx]; new_K[i] = sub_k[i-idx]; new_M[i] = sub_m[i-idx];
								      }
								    else
								      {
									new_node[i] = node2[i-1];
									new_D[i] = D2[i-1]; new_K[i] = K2[i-1]; new_M[i] = M2[i-1];
								      }
								  }
								new_node[size_vec+2] = NULL;
								/*
								printf("[CASE 3-3 T2::");
								for(i=0;i<size_vec+2;i++)
								  printf("n2->id=%d d0=%d m0=%d k0=%d  ",new_node[i]->id,new_D[i],new_M[i],new_K[i]);
								printf("]\n");
								*/
								//printf("[CASE 3-3]");
								//printf("[CASE 3-3 T2::term1=%f term2=%f]",term1,term2);
								//printf("[CASE 3-3 T2::d0=%d m0=%d k0=%d d1=%d m1=%d k1=%d d2=%d m2=%d k2=%d] => R_joint\n",new_D[0],new_M[0],new_K[0],new_D[1],new_M[1],new_K[1],new_D[2],new_M[2],new_K[2]);
								sub_res =  log(R_joint(node1,new_node,D1,new_D,K1,new_K,M1,new_M,F,size_F))+term1-log_factorial(term2);
								res += exp(sub_res);
								/*
								if(res>0)
								  {
								    printf("[CASE 3-3 T2::");
								    for(i=0;i<size_vec+2;i++)
								      printf("n2->id=%d d0=%d m0=%d k0=%d  ",new_node[i]->id,new_D[i],new_M[i],new_K[i]);
								    printf("]\n");
								    //printf(" [case 3-3: term1=%f term2=%f] ",term1,term2);
								    //printf(" [case 3-3 R_joint=%f] ",R_joint(node1,new_node,D1,new_D,K1,new_K,M1,new_M,F,size_F)*term1/factorial(term2));
								    printf("[case 3-3: sub_res = %f\t res=%f]\n",sub_res,res);
								  }
								*/
								
							      }
							  }
						      }
						  }
					    }
					}
				    }
			    } 
			}	    
		  }
	      }    
	}  
    }
  return res;
}


/*
int is_sameShape(struct node *nd1, struct node *nd2)
{
  if(nd1->rank == nd2->rank && get_no_tips(nd1) == get_no_tips(nd2))
    return 1;
  else
    return 0;
}
*/

int no_directTips(struct node *nd)
{
  int i; 
  int n=0;
  
  if(nd->terminal ==1)
    n=1;
  else
    for(i=0;i<3;i++)
      if(nd->listnodes[i] != NULL)
	if(nd->listnodes[i]->terminal == 1)
	  n++;
  
  return n;
}

int no_internalchildrenNodes(struct node *nd)
{
  int i; 
  int n=0;
  
  if(nd->terminal ==0)
    for(i=0;i<3;i++)
      if(nd->listnodes[i] != NULL)
	if(nd->listnodes[i]->terminal == 0)
	  n++;
  
  return n;
}


int get_sum_vec(int vec[HALF_MAX_TAXA], int size)
{
  int i, res=0;

  for(i=0;i<size;i++)
    {
      res += vec[i]; //printf("i=%d\t vec[i]=%d\t",i,vec[i]); 
    }

  return res;
}

int get_no_tips(struct node *nd)
{
  // returns the number of leaves that are descendants of the given node

  int i, res=0;
  
  if(nd->terminal==1)
    res++;
  else
    {
      for(i=0;i<3;i++)
	{
	  if(nd->listnodes[i] != NULL)
	    {
	      if(nd->listnodes[i]->terminal ==1)
		res++;
	      else
		res += get_no_tips(nd->listnodes[i]);
	    }
	}
    }
  
  return res;
}

int get_no_directtips(struct node *nd)
{
  // returns the number of leaves that are directly connected to the given node
  
  int i, res=0;
  
  if(nd->terminal==1)
    res++;
  else
    {
      for(i=0;i<3;i++)
	{
	  if(nd->listnodes[i] != NULL)
	    {
	      if(nd->listnodes[i]->terminal ==1)
		res++;
	    }
	}
    }
  
  return res;
}


int getsize_vectorNode(struct node *node[HALF_MAX_TAXA])
{
  int i=0;
  
  // printf("i=%d\n",i); ---- will be deleted
  while(node[i]!=NULL)
    {
      i++;
      /* will be deleted
      printf("i=%d\n",i);
      printf("res=%d\n",res);
      */
    }

  /* will be deleted
  for(i=0;i<HALF_MAX_TAXA;i++)
    {
      if(node[i]==NULL)
	printf("i=%d\t node is NULL\n",i);
      else
	printf("i=%d\t node is NOT NULL\n",i);	
    }
  */
  return i;
}



long double r_s_joint(struct node *nd1, struct node *nd2, int d, int ntaxa)
{
  // returns the sum of R_joint over k from 0 to Ntaxa, for a given fixed d
  
  struct node *vecNode1[HALF_MAX_TAXA], *vecNode2[HALF_MAX_TAXA];
  int D1[HALF_MAX_TAXA]={0}, D2[HALF_MAX_TAXA]={0}, K1[HALF_MAX_TAXA]={0}, K2[HALF_MAX_TAXA]={0}, M1[HALF_MAX_TAXA]={0}, M2[HALF_MAX_TAXA]={0}, F[HALF_MAX_TAXA]={0};
  int size_F = 0;
  

  vecNode1[0]=nd1; vecNode1[1]=NULL;
  vecNode2[0]=nd2; vecNode2[1]=NULL;
  D1[0]=D2[0]=d;
  M1[0]=M2[0]=0;
  int k; 
  long double r=0, temp;
  for(k=0;k<=ntaxa;k++)
    {      
      K1[0]=K2[0]=k;
      temp =R_joint(vecNode1,vecNode2,D1,D2,K1,K2,M1,M2,F,size_F);
      r += temp;
      //printf("d=%d\t k=%d\t nd1->rank=%d\t nd2->rank=%d\t Rjoint=%f\n",d,k,nd1->rank,nd2->rank, R_joint(vecNode1,vecNode2,D1,D2,K1,K2,M1,M2));
    }
  return r;
}

/*
void gen_r_joint(double gen_r_s_joint[], int ntaxa, struct node *nd1, struct node *nd2)
{
  // generating function for r_s_joint
  int d;
  int n_intnod = ntaxa-3;
  for(d=0;d<=n_intnod;d++)
    {
      gen_r_s_joint[d] = r_s_joint(nd1,nd2,d,ntaxa);
      //printf("in RFdistribution: s=%d\t gen_r_s=%d\n",s, gen_r_s[s]);
    }
}
*/



long double RFjointdistribution(struct node *nd1, struct node *nd2, int ntaxa, int d)
{
  int n_intnod = ntaxa-3;
  //double gen_r_s_joint[n_intnod+1];
  int s;
  long double q_s, temp;
  int iseven;
  
  
  //gen_r_joint(gen_r_s_joint,ntaxa,nd1,nd2);

  q_s = r_s_joint(nd1,nd2,d,ntaxa);
  if(d<n_intnod)
    {
      for(s=d+1;s<=n_intnod;s++)
	{
	  iseven = (s-d)/2;
	  if(2*iseven == s-d)
	    {
	      temp = log(r_s_joint(nd1,nd2,s,ntaxa))+log_factorial(s)-log_factorial(d)-log_factorial(s-d);
	      q_s += exp(temp);
	      //q_s += gen_r_s_joint[s]*factorial(s)/(factorial(d)*factorial(s-d));
	      //printf("even\t q_s=%d\n",factorial(s-d));
	    }
	  else
	    {
	      temp = log(r_s_joint(nd1,nd2,s,ntaxa))+log_factorial(s)-log_factorial(d)-log_factorial(s-d);
	      q_s -= exp(temp);
	      //q_s -= gen_r_s_joint[s]*factorial(s)/(factorial(d)*factorial(s-d));
	      //printf("odd\n");
	    }
	  //printf("in RFdistribution: d=%d\t s=%d\t q_s=%Le\n",d,s,q_s);
	}
    }
  //printf("end: q_s=%Le\n",q_s);
  return q_s;  
}


/*----------------------------------------------------
option = 1 only if the joint distribution is saved
        2 only if the computation time is calculated and the joint distribution is saved
----------------------------------------------------- */
void fullJointDistribution(int ntaxa, int nshapes, long double*** propDistance_twoTrees, int option)
{
  //printf("in fullJointDistribution_test(): read ntaxa = %d\t nshapes = %d\n",ntaxa,nshapes);
  
  long double Z1 = treeSpaceSize(ntaxa);
  int st_sizeRank[nshapes][6];
  FILE *fopen(), *out_jointDist; // FILE declarations
  if(option == 2)
    {
      printf("Start computing the sizes and ranks of subtrees of each tree shape:\n");
      print_cpu_time();
    }

  subtree_sizeRank(ntaxa, st_sizeRank); 
  if(option == 2)
    {
      printf("Done:\n");
      print_cpu_time();
    }

  struct node *nd1, *nd2;
  int intnod = ntaxa+1,ternod=0;
  long double freq;
  char filename[10000];
  int i,j,d;
  //printf("start computing joint distribution...\n");
 if(option >= 1)
   {
     sprintf(filename,"JointDistribution_%dtaxa.txt",ntaxa);
     out_jointDist = fopen(filename,"w");// "JointDistribution","w");
     fprintf(out_jointDist,"P(Distance,Rank2|Rank1)\n\t\tDistance\nRank1\tRank2");
     for(i=ntaxa-3;i>=0;i--)
       fprintf(out_jointDist,"\t%d",i);
     fprintf(out_jointDist,"\n");
     // fprintf(out_jointDist,"%d\t%d\t%d\n",nshapes,nshapes,ntaxa-2);
     fclose(out_jointDist);
   }	  
 if(option==2)
   printf("Start computing the joint distribution:\n");
  for(i=1;i<=nshapes;i++)
    {
      nd1 = invcrank_original(ntaxa,intnod,ternod,i,st_sizeRank[i-1][0],st_sizeRank[i-1][1],st_sizeRank[i-1][2],st_sizeRank[i-1][3],st_sizeRank[i-1][4],st_sizeRank[i-1][5]);
      for(j=1;j<=nshapes;j++)
	{
	  nd2 = invcrank_original(ntaxa,intnod,ternod,j,st_sizeRank[j-1][0],st_sizeRank[j-1][1],st_sizeRank[j-1][2],st_sizeRank[j-1][3],st_sizeRank[j-1][4],st_sizeRank[j-1][5]);       
	  /*
	  if(isprint)
	    {
	      //printRoot(nd1);
	      //printRoot(nd2);
	    }	  
	  */
	  for(d=0;d<=ntaxa-3;d++)
	    {		 
	      freq = RFjointdistribution(nd1,nd2,ntaxa,d);
	      //printf("i=%d j=%d d=%d yes\n",i,j,d);
	      if(freq<0 && freq >-1)
		freq *= -1;
	      if(option==1 || option==2)
		{
		  out_jointDist = fopen(filename,"a");		 
		  //printf("%d\t %d\t %d\t",i,j,ntaxa-3-d); 
		  if(d==0)
		    fprintf(out_jointDist,"%d\t%d\t%Le",i,j,freq/Z1);
		  else
		    fprintf(out_jointDist,"\t%Le",freq/Z1);
		  fclose(out_jointDist);
		}
	      //printf("start\n");
	      propDistance_twoTrees[i-1][j-1][ntaxa-3-d]= (long double) freq/Z1;
	     
	      //printf("end\n");
	      
	      //printf("in fullJointDistribution_test(): the proportion of tree with shape %d, %d distant away from tree with shape %d has been calculated.\n",j,d,i);
	    }
	  if(option==2)
	    {
	      printf("Done: Rank of Shape1 - %d, Rank of shape2 - %d\n",i,j);		  
	      print_cpu_time();
	    }	  
	  if(option>=1)
	    {
	      out_jointDist = fopen(filename,"a");	
	      fprintf(out_jointDist,"\n");
	      fclose(out_jointDist);
	    }
	}
    }  
}



void print_fullJointDistribution(int ntaxa, int nshapes, long double*** propDistance_twoTrees)
{
  fullJointDistribution(ntaxa, nshapes, propDistance_twoTrees,1);
}

void runningTime_fullJointDistribution(int ntaxa, int nshpaes, long double*** propDistance_twoTrees)
{
  fullJointDistribution(ntaxa, nshpaes, propDistance_twoTrees,2);
}


void eachfullJointDistribution(int ntaxa, int nshapes, long double** propDistance_twoTrees, int option, int shapeID)
{
  //printf("in fullJointDistribution_test(): read ntaxa = %d\t nshapes = %d\n",ntaxa,nshapes);
  long double Z1 = treeSpaceSize(ntaxa);
  int st_sizeRank[nshapes][6];
  FILE *fopen(), *out_jointDist; // FILE declarations
  if(option == 2)
    {
      printf("Start computing the sizes and ranks of subtrees of each tree shape:\n");
      print_cpu_time();
    }

  subtree_sizeRank(ntaxa, st_sizeRank); 
  if(option == 2)
    {
      printf("Done:\n");
      print_cpu_time();
    }

  struct node *nd1, *nd2;
  int intnod = ntaxa+1,ternod=0;
  long double freq;
  int i,j,d;
  char filename[10000];
  sprintf(filename,"JointDistribution_%dtaxa_Rank1_%d.txt",ntaxa,shapeID);
  //printf("start computing joint distribution...\n");
  if(option >= 1)
    {
      out_jointDist = fopen(filename,"w");      
      fprintf(out_jointDist,"P(Distance,Rank2|Rank1=%d)\n\tDistance\nRank2",shapeID);
      for(i=ntaxa-3;i>=0;i--)
	fprintf(out_jointDist,"\t%d",i);
      fprintf(out_jointDist,"\n");
      //     out_jointDist = fopen("JointDistribution","w");
      // fprintf(out_jointDist,"%d\t%d\t%d\n",nshapes,nshapes,ntaxa-2);
      fclose(out_jointDist);
    }	  
  if(option==2)
    printf("Start computing the joint distribution:\n");
  //for(i=1;i<=nshapes;i++)
  i = shapeID;
  {
    nd1 = invcrank_original(ntaxa,intnod,ternod,i,st_sizeRank[i-1][0],st_sizeRank[i-1][1],st_sizeRank[i-1][2],st_sizeRank[i-1][3],st_sizeRank[i-1][4],st_sizeRank[i-1][5]);
    for(j=1;j<=nshapes;j++)
      {
	nd2 = invcrank_original(ntaxa,intnod,ternod,j,st_sizeRank[j-1][0],st_sizeRank[j-1][1],st_sizeRank[j-1][2],st_sizeRank[j-1][3],st_sizeRank[j-1][4],st_sizeRank[j-1][5]);       
	/*
	  if(isprint)
	  {
	  //printRoot(nd1);
	  //printRoot(nd2);
	  }	  
	*/
	for(d=0;d<=ntaxa-3;d++)
	  {		 
	    freq = RFjointdistribution(nd1,nd2,ntaxa,d);
	    //printf("i=%d j=%d d=%d yes\n",i,j,d);
	    if(freq<0 && freq >-1)
	      freq *= -1;
	    if(option==1 || option==2)
	      {
		out_jointDist = fopen(filename,"a");
		//out_jointDist = fopen("JointDistribution","a");		 
		//printf("%d\t %d\t %d\t",i,j,ntaxa-3-d); 
		if(d==0)		  
		  fprintf(out_jointDist,"%d\t%Le",j,freq/Z1);
		else
		  fprintf(out_jointDist,"\t%Le",freq/Z1);
		  fclose(out_jointDist);
		}
	      //printf("start\n");
	      propDistance_twoTrees[j-1][ntaxa-3-d]= (long double) freq/Z1;
	     
	      //printf("end\n");
	      
	      //printf("in fullJointDistribution_test(): the proportion of tree with shape %d, %d distant away from tree with shape %d has been calculated.\n",j,d,i);
	    }
	  if(option==2)
	    {
	      printf("Done: Rank of Shape1 - %d, Rank of shape2 - %d\n",i,j);		  
	      print_cpu_time();
	    }	  
	  if(option>=1)
	    {
	      out_jointDist = fopen(filename,"a");
	      //out_jointDist = fopen("JointDistribution","a");	
	      fprintf(out_jointDist,"\n");
	      fclose(out_jointDist);
	    }
	}
    }  
}


void runningTime_eachfullJointDistribution(int ntaxa, int nshpaes, long double** propDistance_twoTrees, int shapeID)
{
  eachfullJointDistribution(ntaxa, nshpaes, propDistance_twoTrees,2,shapeID);
}


/*


void fullJointDistribution(int ntaxa, int nshapes, long double propDistance_twoTrees[][NO_SHAPES][NO_TAXA], int isprint)
{
  printf("ntaxa = %d\t nshapes = %d\n",ntaxa,nshapes);
  long double Z1 = treeSpaceSize(ntaxa);
  int st_sizeRank[nshapes][6];
  subtree_sizeRank(ntaxa, st_sizeRank); 
  
  struct node *nd1, *nd2;
  int intnod = ntaxa+1,ternod=0;
  long double freq;
  int i,j,d;
  
  for(i=1;i<=nshapes;i++)
    {
      nd1 = invcrank_original(ntaxa,intnod,ternod,i,st_sizeRank[i-1][0],st_sizeRank[i-1][1],st_sizeRank[i-1][2],st_sizeRank[i-1][3],st_sizeRank[i-1][4],st_sizeRank[i-1][5]);
      for(j=1;j<=nshapes;j++)
	{
	  nd2 = invcrank_original(ntaxa,intnod,ternod,j,st_sizeRank[j-1][0],st_sizeRank[j-1][1],st_sizeRank[j-1][2],st_sizeRank[j-1][3],st_sizeRank[j-1][4],st_sizeRank[j-1][5]);
	  if(isprint)
	    {
	      printRoot(nd1);
	      printRoot(nd2);
	    }
	  for(d=0;d<=ntaxa-3;d++)
	    {		 
	      freq = RFjointdistribution(nd1,nd2,ntaxa,d);
	      //printf("i=%d j=%d d=%d yes\n",i,j,d);
	      if(freq<0 && freq >-1)
		freq *= -1;
	      if(isprint)
		{		 
		  printf("%d\t %d\t %d\t",i,j,ntaxa-3-d); 
		  printf("Freq=%Le Prob=%Le\n",freq,freq/Z1);
		}
	      else
		{
		  //printf("start\n");
		  propDistance_twoTrees[i-1][j-1][ntaxa-3-d]= (long double) freq/Z1;
		  //printf("end\n");
		}
	      
	    }	  
	  if(isprint)
	    printf("\n");
	}
    }  
}
*/
