/*
  Created Feb. 7, 2011 by Yujin Chung
  updated on June 15, 2012 by Yujin Chung
*/

/* file name: shape.h

   brief: 
          is a header file for shape.c
	  includes structure of node
   
 */


#define HALF_MAX_TAXA 16/2
#define MAX_ARR 3000

struct internal_node
{
  int id;

  struct internal_node *nodelist[3];
  
};

struct node
{
  int id;
  int terminal; // 1 if termnimal node; 0 otherwise (internal node or root node).
  int rank; // LLC rank of subtree
  struct node *listnodes[3]; // pointers to two children nodes
  struct node *parent;
};



/* Global Functions */

// from shape.c

struct node *invcrank_original(int ntaxa, int intnod, int ternod, int rk, int alpha, int beta, int gamma, int rL, int rM, int rR);
struct node *invrrank(int ntaxa, int rank, int intnod, int ternod);
void printRoot(struct node *node);
void printRoot_tips(struct node *node);
void printNode_tips(struct node *node, struct node *childnode);
void subtree_sizeRank(int ntaxa, int sizeRank[][6]);
int no_rshapes(int ntaxa);
int no_urshapes(int ntaxa);
struct node *invcrank(struct node *nd, struct node *childnd,int ntaxa, int intnod, int ternod, int alpha, int beta, int rL, int rM);
struct node *invcrank_Root(int ntaxa, int intnod, int ternod, int alpha, int beta, int gamma, int rL, int rM, int rR);
struct node *reroot(struct node *root);
struct node *noNULL_node(struct node *root);
struct node *pseudo_rooted_Root(struct node *root);
struct node *LLC2urtree(struct node *root);
void printAllShapes(int ntaxa);
int identicalShape(struct node *tr1, struct node *tr2);

// from RFdistribution.c

int no_intnod(struct node *nd, struct node *childnode);
int no_intnod_root(struct node *root);
int no_intchildrenNodes(struct node *root, struct node *nd);
int no_childrenNodes_forRooted(struct node *nd);
double R(struct node *root, struct node *nd, int s, int k);
double r_s(struct node *nd, int s);
double RFdistribution(struct node *nd, int ntaxa, int d);
double log_noTrees_onShape(int ntaxa, struct node *tr);
int factorial(int x);
double log_factorial(int x);


// from RFjointDistribution.c
int no_internalchildrenNodes(struct node *nd);
int no_directTips(struct node *nd);
int getsize_vectorNode(struct node *node[HALF_MAX_TAXA]);
int get_no_tips(struct node *nd);
int get_no_directtips(struct node *nd);
int get_sum_vec(int vec[HALF_MAX_TAXA], int size);
//int is_sameShape(struct node *nd1, struct node *nd2);
int get_node_toExtend(int size, struct node *node[HALF_MAX_TAXA],int D[HALF_MAX_TAXA],int  M[HALF_MAX_TAXA]);
int get_node_toBeNamed(int size, struct node *node[HALF_MAX_TAXA],int  D[HALF_MAX_TAXA],int  M[HALF_MAX_TAXA]);
long double R_split_root(int size1, int size2, int index, int index2, struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F);
long double R_joint(struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F);
long double r_s_joint(struct node *nd1, struct node *nd2, int d, int ntaxa);
long double RFjointdistribution(struct node *nd1, struct node *nd2, int ntaxa, int d);
long double R_extend_vec2(struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F);
long double R_extend_vec1(struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA],int F[HALF_MAX_TAXA],int size_F);
int split_test(int size1, int size2, int index, struct node *node1[HALF_MAX_TAXA],struct node *node2[HALF_MAX_TAXA],int D1[HALF_MAX_TAXA],int D2[HALF_MAX_TAXA],int K1[HALF_MAX_TAXA],int K2[HALF_MAX_TAXA],int M1[HALF_MAX_TAXA],int M2[HALF_MAX_TAXA], int F[HALF_MAX_TAXA],int size_F);
int isSymmetric(struct node *node1, struct node *node2);
struct node *findMRCA(struct node *node1, struct node *node2);
int is_doubleTruncated(struct node *node);
void print_fullJointDistribution(int ntaxa, int nshpaes, long double*** propDistance_twoTrees);
void fullJointDistribution(int ntaxa, int nshapes, long double*** propDistance_twoTrees, int isprint);
void runningTime_fullJointDistribution(int ntaxa, int nshpaes, long double*** propDistance_twoTrees);
void runningTime_eachfullJointDistribution(int ntaxa, int nshpaes, long double** propDistance_twoTrees, int shapeID);
void eachfullJointDistribution(int ntaxa, int nshapes, long double** propDistance_twoTrees, int option, int shapeID);


/* normalizingConstant.c */
long double logZ_L(int ntaxa, int nshapes, int nsites, long double b, long double*** propDistance_twoTrees);
long double log_normalizingConstant_indep(long double*** propDistance_twoTrees, int ntaxa, unsigned int nsites, int nshape, long double b);
long double log_normalizingConstant_Normalapprox(int nsites,int ntaxa, int nshape, long double b, long double*** propDistance_twoTrees);
long double log_normalizingConstant_Normalapprox2(int nsites,int ntaxa, int nshape, long double b, long double*** propDistance_twoTrees);
long double expected_dist2(int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double var_dist2(int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double cov_dist2(int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double treeSpaceSize(int ntaxa);
long double expected_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double var_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double cov_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double mu_L(int nsites, int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double var_L(int nsites, int ntaxa, int nshape,long double*** propDistance_twoTrees);
long double corr_dist(int ntaxa, int nshape,long double*** propDistance_twoTrees);
