// PROPOSAL CLASS

# pragma once
# include <iostream>
# include <cmath>
# include <vector>
# include <algorithm>
# include <map>
# include <set>
# include "move.h"
# include "tree.h"
# include "read.h"
# include "dist.h"

using namespace std;

class Prop {
	public:
	Move M;	
	Tree T;
	Dist Di;
	Dataline D;	

	Dist::random_selector<> selector;

	// Constructors
	Prop() {};

	// Declarations
	void gprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, 
			const double& ALPH, const double& BET, const int& DIM,
			const int& npred, double& qratio, const vector<string>& names, 
			map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
			const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
	void pprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, 
			const int& DIM, double& qratio);
	void cprop(map<int, Node*>& tree, const vector<int>& int_nodes,
			double& qratio, const vector<string>& names, 
			map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
			const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
	void sprop(map<int, Node*>& tree, const vector<int> int_nodes,
		double& qratio);
	double gprop_ratio(const pair<int, Node*>& left, const pair<int, Node*>& right, 
		const int& nleaves, const double& psplit, const int& pairs); 
	double pprop_ratio(const pair<int, Node*>& left, const pair<int, Node*>& right, 
		const int& nleaves, const double& psplit, const int& pairs);
	bool check_siblings(const map<int, Node*>& tree, const int& node);
};

double Prop::gprop_ratio(const pair<int, Node*>& left, const pair<int, Node*>& right, 
		const int& nleaves, const double& psplit, const int& npairs) {
	double lprule, rprule, lratio;
	lprule = 1/(double)left.second->prule;
	rprule = 1/(double)right.second->prule;
	cout << "lprule " << lprule << " rprule " << rprule << endl;
	if(nleaves == 1) {
		lratio = 1;
	cout << "lratio " << lratio << endl;
	}
	else {
		lratio = (npairs/double(nleaves-npairs))*((2*nleaves-1)/double(nleaves));
	cout << "lratio " << lratio << endl;
	}

	double qratio = lratio*lprule*rprule*((1-psplit)/(double)psplit);
	cout << "qratio " << qratio << endl;
	return qratio;
}

double Prop::pprop_ratio(const pair<int, Node*>& left, const pair<int, Node*>& right, 
		const int& nleaves, const double& psplit, const int& npairs) {
	double lprule, rprule;
	lprule = left.second->prule;
	rprule = right.second->prule;

	double lratio = (nleaves/double(2*nleaves+3))*((nleaves - npairs)/double(npairs));

	double qratio = lratio*lprule*rprule*psplit/(double)(1 - psplit);
	return qratio;
}

void Prop::gprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, 
		const double& ALPH, const double& BET, const int& DIM,
		const int& npred, double& qratio, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {

	int num_leaves = leaf_nodes.size();
	int rand_leaf = selector(leaf_nodes);
	
	// Grow the tree
	M.grow(tree, rand_leaf, ALPH, BET, DIM, names, alltypes, allcats, dbounds, ibounds);

	// CALCULATE PROPOSALS
	map<int, Node*>::iterator it;
	it = tree.find(rand_leaf);
	pair<int, Node*> randl = *it;
	it = tree.find(2*rand_leaf);
	pair<int, Node*> nleft = *it;
	it = tree.find(2*rand_leaf+1);
	pair<int, Node*> nright = *it;

	vector<int> parent_pairs = T.parent_pairs(leaf_nodes);

	int npairs = parent_pairs.size();
	
	qratio = gprop_ratio(nleft, nright, num_leaves, 
			randl.second->psplit, npairs);

}

void Prop::pprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, 
	const int& DIM, double& qratio) {

	int const num_leaves = leaf_nodes.size();
	vector<int> parent_pairs = T.parent_pairs(leaf_nodes);

	int npairs = parent_pairs.size();

	int rand_parent = selector(parent_pairs);

	// CALCULATE PROPOSALS
	map<int, Node*>::iterator it;
	it = tree.find(rand_parent);
	pair<int, Node*> randl = *it;
	it = tree.find(2*rand_parent);
	pair<int, Node*> nleft = *it;
	it = tree.find(2*rand_parent+1);
	pair<int, Node*> nright = *it;

	qratio = pprop_ratio(nleft, nright, num_leaves, 
			randl.second->psplit, npairs);

	// Prune
	Move M;
	M.prune(tree, rand_parent, DIM); 
}

void Prop::cprop(map<int, Node*>& tree, const vector<int>& int_nodes,
		double& qratio, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	int rand_int;
	if (int_nodes.empty()) {
		rand_int = 1;
	}
	else {
		rand_int = selector(int_nodes);
	}
	
	// Change random node partition
	M.change(tree, rand_int, names, alltypes, allcats, dbounds, ibounds);

	// CALCULATE PROPOSALS
	qratio = 1;	
}

bool Prop::check_siblings(const map<int, Node*>& tree, const int& node) {
	map<int, Node*>::const_iterator lit;
	map<int, Node*>::const_iterator rit;
	lit = tree.find(2*node);
	rit = tree.find(2*node+1);
	if(lit->second->ntype == Dataline::STRING_T &
			rit->second->ntype == Dataline::STRING_T) {
		if(lit->second->pred == rit->second->pred &
				lit->second->ssplit == rit->second->ssplit) {
			return true;}
	}
	else if(lit->second->ntype == Dataline::INT_T &
			rit->second->ntype == Dataline::INT_T) {
		if(lit->second->pred == rit->second->pred &
				lit->second->dsplit == rit->second->dsplit) {
			return true;}
	}
	else if(lit->second->ntype == Dataline::DOUBLE_T &
			rit->second->ntype == Dataline::DOUBLE_T) {
		if(lit->second->pred == rit->second->pred &
				lit->second->dsplit == rit->second->dsplit) {
			return true;}
	}
	else {return false;}
}


void Prop::sprop(map<int, Node*>& tree, const vector<int> int_nodes,
		double& qratio) {
	int num_ints, rand_int1, rand_int2, rand_int3, rand_swap;
	num_ints = int_nodes.size();
	rand_swap = Di.runif(0,1);
	
	while(true) {
		rand_int1 = selector(int_nodes);
		bool match = check_siblings(tree, rand_int1);
		if (rand_int1 == 1) {
			if (match) { // if equal, swap parent with both children
				rand_int2 = 2;
				rand_int3 = 3;
				M.mswap(tree, rand_int1, rand_int2);
				M.mswap(tree, rand_int1, rand_int3);
				break;
			}
			else if (rand_swap == 0) {
				rand_int2 = 2; 
				M.mswap(tree, rand_int1, rand_int2);
				break;
			} //swap parent/lchild
			else if (rand_swap == 1) {
				rand_int2 = 3;
				M.mswap(tree, rand_int1, rand_int2);
				break;
			} //swap parent/rchild
		}
		else {
			Node * left = tree[2*rand_int1]; //left child
			Node * right = tree[2*rand_int1+1]; //right child
			if (match) {
				if (!left->isleaf & !right->isleaf) {
					rand_int2 = 2*rand_int1;
					rand_int3 = 2*rand_int1+1;
					M.mswap(tree, rand_int1, rand_int2);
					M.mswap(tree, rand_int1, rand_int3);
					break;	
				}
			}
			else if (rand_swap == 0 & !left->isleaf) {
				rand_int2 = 2*rand_int1; 
				M.mswap(tree, rand_int1, rand_int2);
				break;
			} 
			else if (rand_swap == 1 & !right->isleaf) {
				rand_int2 = 2*rand_int1+1; 
				M.mswap(tree, rand_int1, rand_int2);
				break;
			} 
		}
	}

	// CALCULATE PROPOSALS
	qratio = 1;
}
