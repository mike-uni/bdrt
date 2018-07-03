// PROPOSAL CLASS
# pragma once
# include <iostream>
# include "move_class.h"
# include "tree_class.h"
# include "node_class.h"

using namespace std;

class Prop {
	public:

	// Constructors
	Prop() {};

	// Declarations

	void gprop(Tree * tptr, Node * fptr, float const pmove,
					const vector<int> leaf_nodes,
					const double ALPH, const double BET,
					double in_rule, const int DIM, vector<string> names,
					const MatrixXd& data, int iter,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet);
	void cprop(Tree * tptr, Node * fptr, float const pmove,
				const vector<int> int_nodes, const int iter,
				const MatrixXd& data, const vector<string> names);
	void sprop(Tree * tptr, Node * fptr, float const pmove,
				const vector<int> int_nodes, const int iter);
	void pprop(Tree * tptr, Node * fptr, float const pmove,
					const vector<int> leaf_nodes,
					const double ALPH, const double BET,
				 	const int DIM, vector<string> names,
					const MatrixXd& data, int iter,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet);

};

void Prop::gprop(Tree * tptr, Node * fptr, float const pmove,
					const vector<int> leaf_nodes,
					const double ALPH, const double BET,
					double in_rule, const int DIM, 
					const vector<string> names,
					const MatrixXd& data, int iter,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet) {

	int num_leaves = leaf_nodes.size();
	int rand_leaf = leaf_nodes[rand()%num_leaves];
	
	Node * rleaf;
	rleaf = tptr->find_node(tptr->root, rand_leaf);

	// Grow the tree
	Move m;	
	m.grow(rleaf, ALPH, BET, in_rule, DIM, names, data, 
			iter, fptr->nnode, WMAT, VMAT, WI, VI, FMAT, HMAT, HVH, 
			VH, WFT, WF, FWF, yi, Vdet, Wdet);

	// CALCULATE PROPOSALS 
	float pleaf;
	pleaf = 1/float(num_leaves);
	float pint;
	pint = 1/float(num_leaves + 1);

	// Proposal (T*, T)
	double prop_gf;
	prop_gf = pmove*pleaf*rleaf->prule*rleaf->psplit*
				pow(1- rleaf->left->psplit, 2);
	tptr->prop_new = prop_gf;

	// Proposal (T, T*)
	double prop_gr;
	prop_gr = pmove*pint*(1 - rleaf->psplit);
	tptr->prop_current = prop_gr;
	//cout << "tptr->prop_new is: " << tptr->prop_new << endl;
	//cout << "tptr->prop_current is: " << tptr->prop_current << endl;

	// Proposal prior p(T*)
	tptr->find_prior(tptr->root);		

	// Proposal posterior p(T*)
	tptr->find_post(tptr->root);		
}
	
void Prop::cprop(Tree * tptr, Node * fptr, float const pmove,
				const vector<int> int_nodes, const int iter,
				const MatrixXd& data, const vector<string> names) {

	int num_ints;
	int rand_int;
	if (fptr->nnode == 1) {
		num_ints = 1;
		rand_int = 1;
	}
	else {
		num_ints = int_nodes.size();
		rand_int = int_nodes[rand()%num_ints];
	}
	Node * rleaf;
	rleaf = tptr->find_node(tptr->root, rand_int);

	// Change random node partition
	Move m;
	m.change(rleaf, data, names, iter);

	// CALCULATE PROPOSALS	
	float pcint;	
	pcint = 1/float(num_ints);

	// Proposal (T*, T) and (T, T*)
	double propch = pmove*rleaf->prule*pcint;
	tptr->prop_new = propch;
	tptr->prop_current = propch;
	//cout << "tptr->prop_new is: " << tptr->prop_new << endl;
	//cout << "tptr->prop_current is: " << tptr->prop_current << endl;
	
	// Proposal prior p(T*)
	tptr->find_prior(tptr->root);
	
	// Proposal posterior p(T*)
	tptr->find_post(tptr->root);		
}

void Prop::sprop(Tree * tptr, Node * fptr, float const pmove,
				const vector<int> int_nodes, const int iter) {

	int num_ints;
	int rand_int;
	if (fptr->nnode == 1) {
		num_ints = 1;
		rand_int = 1;
	}
	else {
		num_ints = int_nodes.size();
		rand_int = int_nodes[rand()%num_ints];
	}
	
	Node * rleaf;
	rleaf = tptr->find_node(tptr->root, rand_int);

	// Swap internal node subtrees
	Move m;
	m.mswap(rleaf);

	// CALCULATE PROPOSALS

	// Parent Child Pairs
	float ppair;
	if (rleaf->nnode == 1) { 
		ppair = 1;
	}
	else {
		ppair = 1/float(num_ints-1);
	}

	// Proposal (T*, T) and (T, T*)
	double propsw = pmove*ppair*rleaf->prule;
	tptr->prop_new = propsw;
	tptr->prop_current = propsw;
	//cout << "tptr->prop_new is: " << tptr->prop_new << endl;
	//cout << "tptr->prop_current is: " << tptr->prop_current << endl;

	// Proposal prior p(T*)
	tptr->find_prior(tptr->root);

	// Proposal posterior p(T*)
	tptr->find_post(tptr->root);		
}

void Prop::pprop(Tree * tptr, Node * fptr, float const pmove,
					const vector<int> int_nodes,
					const double ALPH, const double BET,
					const int DIM, 
					const vector<string> names,
					const MatrixXd& data, int iter,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet) {

	int const num_nodes = int_nodes.size();
	int const rand_node = int_nodes[rand()%num_nodes];
	
	Node * rleaf;
	rleaf = tptr->find_node(tptr->root, rand_node);
		
	float pleaf = 1/float(num_nodes+1);	
	float pint = 1/float(num_nodes);
				
	// Proposal (T*, T)
	double ppropf = pmove*pint*(1 - rleaf->psplit);
	tptr->prop_new = ppropf;

	// Proposal (T, T*)
	double ppropr = pmove*pleaf*rleaf->psplit*rleaf->prule*
						pow(1 - rleaf->psplit, 2);
	tptr->prop_current = ppropr;
	//cout << "tptr->prop_new is: " << tptr->prop_new << endl;
	//cout << "tptr->prop_current is: " << tptr->prop_current << endl;

	// Prune
	Move m;
	m.prune(rleaf, DIM, data, fptr->nnode, WMAT, VMAT, WI, VI, FMAT, 
			HMAT, HVH, VH, WFT, WF, FWF, yi, Vdet, Wdet);
				
	// Proposal prior p(T*)
	tptr->find_prior(tptr->root);

	// Proposal posterior p(T*)
	tptr->find_post(tptr->root);		
}
