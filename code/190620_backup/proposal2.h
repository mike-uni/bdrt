// PROPOSAL2 CLASS

# pragma once
# include <iostream>
# include <cmath>
# include <vector>
# include <algorithm>
# include <map>
# include <set>
# include "move2.h"
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
	void mgprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, const int nlevels, const int levelt,
		const double& ALPH, const double& BET, const int& MDIM, const int& NDIM, const int ngrow,
		const int& npred, const vector<string>& names, map<string, Dataline::type> alltypes, 
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds, const map<string, 
		vector<int>>& ibounds, const double& ntemp, const double& ctemp, string sttype, 
		double& lfratio, double& lrratio); 
	void mpprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<int>& leaf_nodes,
			const int& DIM, const int prunesize, const double& ntemp, 
			const double& ctemp, string sstype, double& lfratio, double& lrratio);
	void cprop(map<int, Node*>& tree, const vector<int>& int_nodes,
			const vector<string>& names, map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
			const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds, 
			double& lfratio, double& lrratio); 
	void mcprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<string>& names, const int nchange, 
			map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
			const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds,
		       	double& lfratio, double& lrratio); 
	void shiftprop(map<int, Node*>& tree, const vector<int>& int_nodes,
			double& lfratio, double& lrratio);
	void sprop(map<int, Node*>& tree, const vector<int> int_nodes,
			double& lfratio, double& lrratio);
	bool check_siblings(const map<int, Node*>& tree, const int& node);
};

void Prop::mgprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, const int nlevels, const int tlevel,
		const double& ALPH, const double& BET, const int& MDIM, const int& NDIM, const int ngrow,
		const int& npred, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds, 
		const double& ntemp, const double& ctemp, string sttype, double& lfratio, double& lrratio) {

	size_t num_leaves = leaf_nodes.size();
	int rand_leaf = selector(leaf_nodes);
	vector<int> parent_pairs = T.parent_pairs(leaf_nodes);
	size_t npairs = parent_pairs.size();

	map<int, Node*>::iterator it;
	it = tree.find(rand_leaf);

	// Grow the tree
	M.mgrow(tree, rand_leaf, ngrow, ALPH, BET, MDIM, NDIM, names, alltypes, allcats, dbounds, ibounds);

	// CALCULATE PROPOSALS
	map<int, Node*> subtree;
	T.find_subtree(tree, rand_leaf, subtree);
	double growprob = 1.0;
	double pruneprob = 1.0;
	double pgchooseleaf = 1;
	double ppchooseleaf = 1;
	int nint;
	double pgsubtree = 1;
	double ppsubtree = 1;
	vector<int> subint_nodes;
	T.find_int_nodes(subtree, subint_nodes, nint); 
	map<int, Node*>::iterator sit;
	for(auto i : subint_nodes) {
		sit = subtree.find(i);
		growprob = pgsubtree*growprob*(sit->second->psplit);
		pruneprob = ppsubtree*pruneprob*(1-sit->second->psplit);
		pgsubtree = (pgsubtree+1)/double(pgsubtree+2);
		ppsubtree = 1/double(ppsubtree+2);
	}
	if (num_leaves == 1) {
		pgchooseleaf = 1;
		ppchooseleaf = 1;
	}
	else {
		pgchooseleaf = (num_leaves)/double(2*num_leaves-1);
		ppchooseleaf = npairs/double(2*num_leaves-1);
	}

	lfratio = (1/double(ntemp))*log(pruneprob*ppchooseleaf);
	lrratio = (1/double(ctemp))*log(growprob*pgchooseleaf);
	cout << "inside grow proposal: " << sit->second->nnode << " " << (1/double(ntemp))  << " " << log(pruneprob*ppchooseleaf)
		       	  << " " << (1/double(ctemp))  << " " << log(growprob*pgchooseleaf) << endl;
}

void Prop::mpprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<int>& leaf_nodes, const int& DIM, 
	const int prunesize, const double& ntemp, const double& ctemp, string sttype, double& lfratio, double& lrratio) {

	Move M;
	size_t nintnodes = int_nodes.size();
	size_t num_leaves = leaf_nodes.size();
	int rand_int;
	double growprob = 1.0;
	double pruneprob = 1.0;
	double pgchooseleaf = 1;
	double ppchooseleaf = 1;
	map<int, Node*>::iterator sit;

	if(prunesize == 1) {
		vector<int> parent_pairs = T.parent_pairs(leaf_nodes);
		size_t npairs = parent_pairs.size();
		rand_int = selector(parent_pairs);

		ppchooseleaf = npairs/double(2*num_leaves-1);
		pgchooseleaf = num_leaves/double(2*num_leaves-1);
		sit = tree.find(rand_int);;
		pruneprob = pruneprob*(1-sit->second->psplit);
		growprob = growprob*(sit->second->psplit);

		lfratio = (1/double(ntemp))*log(growprob*pgchooseleaf);
		lrratio = (1/double(ctemp))*log(pruneprob*ppchooseleaf);
		// Prune
		M.prune(tree, rand_int, DIM); 
	}
	else if(prunesize == 2) {
		int rep;
		if (num_leaves > 3) {rep = 2;}
		else {rep = 1;}
		for(int i = 0; i < rep; ++i) {
			vector<int> parent_pairs = T.parent_pairs(leaf_nodes);
			size_t npairs = parent_pairs.size();
			rand_int = selector(parent_pairs);
			ppchooseleaf = ppchooseleaf*npairs/double(2*num_leaves-1);
			pgchooseleaf = pgchooseleaf*num_leaves/double(2*num_leaves-1);
			sit = tree.find(rand_int);;
			pruneprob = pruneprob*(1-sit->second->psplit);
			growprob = growprob*(sit->second->psplit);
			// Prune
			M.prune(tree, rand_int, DIM);
			num_leaves = num_leaves-1;
		}	
		lfratio = (1/double(ntemp))*log(growprob*pgchooseleaf);
		lrratio = (1/double(ctemp))*log(pruneprob*ppchooseleaf);
	}
	else {
		rand_int = selector(int_nodes);

		// CALCULATE PROPOSALS
		map<int, Node*> subtree;
		T.find_subtree(tree, rand_int, subtree);
		int nint;
		double pgsubtree = 1;
		double ppsubtree = 1;
		vector<int> subint_nodes;
		T.find_int_nodes(subtree, subint_nodes, nint); 
		for(auto i : subint_nodes) {
			sit = subtree.find(i);;
			growprob = pgsubtree*growprob*(sit->second->psplit);
			pruneprob = ppsubtree*pruneprob*(1-sit->second->psplit);
			pgsubtree = (pgsubtree+1)/double(pgsubtree+2);
			ppsubtree = 1/double(ppsubtree+2);
		}
		if (num_leaves == 1) {
			pgchooseleaf = 1;
			ppchooseleaf = 1;
		}
		else {
			pgchooseleaf = (num_leaves)/double(2*num_leaves-1);
			ppchooseleaf = (num_leaves-1)/double(2*num_leaves-1);
		}
		lfratio = (1/double(ntemp))*log(growprob*pgchooseleaf);
		lrratio = (1/double(ctemp))*log(pruneprob*ppchooseleaf);
		// Prune
		Move M;
		M.prune(tree, rand_int, DIM);
	}	
	cout << "inside prune proposal: " << (1/double(ntemp))  << " " << log(pruneprob)
		       	  << " " << (1/double(ctemp))  << " " << log(growprob);
}

void Prop::cprop(map<int, Node*>& tree, const vector<int>& int_nodes,
		const vector<string>& names, map<string, Dataline::type> alltypes, const map<string, 
		string>& allcats, const map<string, vector<double>>& dbounds, const map<string, 
		vector<int>>& ibounds, double& lfratio, double& lrratio) {
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
	lfratio = 0;
	lrratio = 0;
}

void Prop::mcprop(map<int, Node*>& tree, const vector<int>& int_nodes,
		const vector<string>& names, const int nchange, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds, 
		double& lfratio, double& lrratio) {
	
	// Change random node partition
	M.mchange(tree, int_nodes, names, nchange, alltypes, allcats, dbounds, ibounds);

	// CALCULATE PROPOSALS
	lfratio = 0;
	lrratio = 0;
}

void Prop::shiftprop(map<int, Node*>& tree, const vector<int>& int_nodes, double& lfratio, double& lrratio) {
	int rand_int;
	if (int_nodes.empty()) {
		rand_int = 1;
	}
	else {
		rand_int = selector(int_nodes);
	}
	
	// Change random node partition
	M.shift(tree, rand_int);

	// CALCULATE PROPOSALS
	lfratio = 0;
	lrratio = 0;
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


void Prop::sprop(map<int, Node*>& tree, const vector<int> int_nodes, double& lfratio, double& lrratio) {
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
	lfratio = 0;
	lrratio = 0;
}
