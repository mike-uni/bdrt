// 	FUNCTIONS FOR EVOLVING THE TREE

# pragma once
# include <iostream>
# include <Eigen/Dense>
# include "priorf.h"
# include "node_class.h"
# include "read_data.h"
# include "posterior_class.h"
# include <algorithm>

using namespace std;
using namespace Eigen;

class Move : public Node {
	public:

	// Constructors
	Move() {};

	// Declarations
	void grow(Node * leaf, const double ALPH, const double BET,
					double in_rule, const int DIM, vector<string> names,
					const MatrixXd& data, int iter, int found_leaf,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet);
	void prune(Node * leaf, const int DIM,
					const MatrixXd& data, int found_leaf,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet);
	void change(Node * leaf, const MatrixXd& data, vector<string> names,
					int iter);
	void mswap(Node * leaf);
};

void Move:: grow(Node * leaf, const double ALPH, const double BET,
					double in_rule, const int DIM, vector<string> names,
					const MatrixXd& data, int iter, int found_leaf,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet) {
 	Lpost l; 
	leaf->isleaf = false;
	
	MatrixXd data1 = data;
	MatrixXd data2 = data;
		
	pair<string, double> lpair;
	lpair = rand_pair(data1, names, data.cols(), data.rows(), iter);
	pair<string, double> rpair;
	rpair = rand_pair(data2, names, data.cols(), data.rows(), iter);

	Node * templ, * tempr;
//	templ = leaf->left;
//	tempr = leaf->right;

	templ = new Node(2*(leaf->nnode), DIM);
	templ->pred = lpair.first;
	templ->split = lpair.second;
	templ->dim = leaf->dim;
	templ->depth = leaf->depth+1;
	templ->psplit = psplitf(ALPH, BET, templ->depth);
	templ->prule = in_rule;
	templ->parent = leaf;
	templ->state = leaf->state;
	templ->svar = leaf->svar;
	templ->atmat = leaf->atmat;
	templ->aimat = leaf->aimat;
	templ->cvec = leaf->cvec;
	templ->dvec = leaf->dvec;
	templ->ppost = leaf->ppost;
	templ->lpost = leaf->lpost;

	tempr = new Node(2*(leaf->nnode)+1, DIM);
	tempr->pred = rpair.first;
	tempr->split = rpair.second;
	tempr->depth = leaf->depth+1;
	tempr->psplit = psplitf(ALPH, BET, tempr->depth);
	tempr->prule = in_rule;
	tempr->parent = leaf;
	tempr->state = leaf->state;
	tempr->svar = leaf->svar;
	tempr->atmat = leaf->atmat;
	tempr->aimat = leaf->aimat;
	tempr->cvec = leaf->cvec;
	tempr->dvec = leaf->dvec;
	tempr->ppost = leaf->ppost;
	tempr->lpost = leaf->lpost;

	if (leaf->nnode == found_leaf) {
		int ind;
		ind = get_pred_index(leaf->pred, names);
		double ssplit = data(iter, ind);
		if (leaf->split < ssplit) {
			// lpostupleft, lpostnupright
			l.lpostup(templ->atmat, templ->aimat, templ->cvec, templ->dvec,
			templ->ppost, templ->lpost, yi, VI, FWF, WI, WF,
			WFT, HVH, VH, Vdet, Wdet);
			l.lpostnup(tempr->atmat, tempr->aimat, tempr->cvec, 
			tempr->dvec, tempr->ppost, tempr->lpost, 
			FWF, WI, WF, WFT, Wdet);
		}
		else {
			// lpostupright, lpostnupleft
			l.lpostup(tempr->atmat, tempr->aimat, tempr->cvec, tempr->dvec,
			tempr->ppost, tempr->lpost, yi, VI, FWF, WI, WF,
			WFT, HVH, VH, Vdet, Wdet);
			l.lpostnup(templ->atmat, templ->aimat, templ->cvec, 
			templ->dvec, templ->ppost, templ->lpost, 
			FWF, WI, WF, WFT, Wdet);
		}
	}
	else {
		// lpostnupboth
		l.lpostnup(templ->atmat, templ->aimat, templ->cvec, templ->dvec, 
				templ->ppost, templ->lpost, FWF, WI, WF, WFT, Wdet);
		l.lpostnup(tempr->atmat, tempr->aimat, tempr->cvec, tempr->dvec, 
				tempr->ppost, tempr->lpost, FWF, WI, WF, WFT, Wdet);
	 }
	leaf->left = templ;
	leaf->right = tempr;
}

void Move:: prune(Node * leaf, const int DIM,
					const MatrixXd& data, int found_leaf,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet) {
	Lpost l;

	VectorXd left_state(DIM);
	left_state = leaf->left->state;
	MatrixXd left_svar(DIM, DIM);
	left_svar = leaf->left->svar;
	VectorXd right_state(DIM);
	right_state = leaf->right->state;
	MatrixXd right_svar(DIM, DIM);
	right_svar = leaf->right->svar;
	
	if (leaf->left->nnode == found_leaf) {
		l.lpostup(leaf->atmat, leaf->aimat, leaf->cvec, leaf->dvec,
		leaf->ppost, leaf->lpost, yi, VI, FWF, WI, WF,
		WFT, HVH, VH, Vdet, Wdet);
	}
	else if (leaf->right->nnode == found_leaf) {
		l.lpostup(leaf->atmat, leaf->aimat, leaf->cvec, leaf->dvec,
		leaf->ppost, leaf->lpost, yi, VI, FWF, WI, WF,
		WFT, HVH, VH, Vdet, Wdet);
	}
	else {
		l.lpostnup(leaf->atmat, leaf->aimat, leaf->cvec, leaf->dvec, 
				leaf->ppost, leaf->lpost, FWF, WI, WF, WFT, Wdet);
	}	
	
	leaf->isleaf = true;
	leaf->state = 0.5*(left_state + right_state);
	leaf->svar = 0.5*(left_svar + right_svar);
	delete leaf->left;
	leaf->left = nullptr;
	delete leaf->right;
	leaf->right = nullptr;
}

void Move:: change(Node * leaf, const MatrixXd& data, 
		vector<string> names, int iter) {
	pair<string, double> newpair;
	newpair = rand_pair(data, names, data.cols(), data.rows(), iter);
	leaf->pred = newpair.first;
	leaf->split = newpair.second; 
}

void Move:: mswap(Node * leaf) {
	//Node * templ = new Node(leaf->left->nnode, leaf->dim);
	//Node * tempr = new Node(leaf->right->nnode, leaf->dim);
	//Node * tempp = new Node(leaf->nnode, leaf->dim);
	
	//templ = leaf->left;
	//templ->aimat = leaf->right->aimat;
	//tempr = leaf->right;
	//tempr->aimat = leaf->left->aimat;
	//tempp = leaf;

	int rand_int = rand()%2;

	if (leaf->nnode == 1) {
		Node * templ;
		Node * tempr;

		templ = leaf->left;
		tempr = leaf->right;
		swap(templ, tempr);
		leaf->left = templ;
		leaf->right = tempr;
	}	
	else {
		Node * tempp;
		Node * templ;
		Node * tempr;

		templ = leaf->left;
		tempr = leaf->right;
		tempp = leaf;

		if (rand_int == 0) {
			swap(tempp, templ);
			leaf->left = tempp;
			leaf = templ;
			}	
		swap(tempp, tempr);
		leaf->right = tempp;
		leaf = tempr;
	}
}
