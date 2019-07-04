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
	void grwprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, const double& ALPH, const double& BET, 
		const int& MDIM, const int& NDIM, const vector<string>& names, map<string, Dataline::type> alltypes, 
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds, 
		const map<string, vector<int>>& ibounds, double& lfratio, double& lrratio, const double pgmove, 
		const double ppmove); 
	void prnprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<int>& leaf_nodes, 
		const int& MDIM, double& lfratio, double& lrratio, const double pgmove, const double ppmove);
	void chprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, const map<string, 
		vector<double>>& dbounds, const map<string, vector<int>>& ibounds, double& lfratio, double& lrratio); 
	void swprop(map<int, Node*>& tree, const vector<int> int_nodes, double& lfratio, double& lrratio); 
	void shftprop(map<int, Node*>& tree, const vector<int>& int_nodes, double& lfratio, double& lrratio);
	void mgprop(map<int, Node*>& tree, double& lfratio, double& lrratio, const double pfmove, const double prmove, 
		const double ALPH, const double BET, const double pnotup);
	void tswprop(map<int, Node*>& tree, double& lfratio, double& lrratio);
	void prngrwprop(map<int, Node*>& tree, const double& ALPH, const double& BET, const int& MDIM, const int& NDIM,
		const vector<int> int_nodes, const vector<int> leaf_nodes, const vector<string>& names, map<string, 
		Dataline::type> alltypes, const map<string, string>& allcats, const map<string, vector<double>>& dbounds, 
		const map<string, vector<int>>& ibounds, double& lfratio, double& lrratio);
	bool check_siblings(map<int, Node*>& tree, const int& node);
	double threshprob(map<int, Node*>& subtree);
	double allthreshprob(map<int, Node*>& subtree);
};

void Prop::grwprop(map<int, Node*>& tree, const vector<int>& leaf_nodes, const double& ALPH, const double& BET, 
		const int& MDIM, const int& NDIM, const vector<string>& names, map<string, Dataline::type> alltypes, 
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds, 
		const map<string, vector<int>>& ibounds, double& lfratio, double& lrratio, const double pgmove, 
		const double ppmove) {

	size_t num_leaves = leaf_nodes.size();
	int rand_leaf = selector(leaf_nodes);
	vector<int> parent_pairs = T.parent_pairs(leaf_nodes);
	size_t npairs = parent_pairs.size();
	
	// Grow the tree
	M.grow(tree, rand_leaf, ALPH, BET, MDIM, NDIM, names, alltypes, allcats, dbounds, ibounds);

	// CALCULATE PROPOSALS
	double growprob = 1.0;
	double pruneprob = 1.0;
	double pgchooseleaf = 1;
	double ppchooseleaf = 1;
	map<int, Node*>::iterator sit;
	sit = tree.find(rand_leaf);;
	growprob = pgmove*(sit->second->psplit)*(sit->second->prule);
	pruneprob = ppmove*(1-sit->second->psplit);
	
	if (num_leaves == 1) {
		pgchooseleaf = 1;
		ppchooseleaf = 1;
	}
	else {
		pgchooseleaf = (num_leaves)/double(2*num_leaves-1);
		ppchooseleaf = npairs/double(2*num_leaves-1);
	}

	lfratio = log(pruneprob*ppchooseleaf);
	lrratio = log(growprob*pgchooseleaf);
}

void Prop::prnprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<int>& leaf_nodes, 
	const int& MDIM, double& lfratio, double& lrratio, const double pgmove, const double ppmove) {

	size_t nintnodes = int_nodes.size();
	size_t num_leaves = leaf_nodes.size();
	int rand_int;
	double growprob = 1.0;
	double pruneprob = 1.0;
	double pgchooseleaf = 1;
	double ppchooseleaf = 1;
	map<int, Node*>::iterator sit;

	// Get parent nodes of leaves that are end pairs
	vector<int> parent_pairs = T.parent_pairs(leaf_nodes);
	size_t npairs = parent_pairs.size();
	rand_int = selector(parent_pairs);

	ppchooseleaf = npairs/double(2*num_leaves-1);
	pgchooseleaf = num_leaves/double(2*num_leaves-1);
	sit = tree.find(rand_int);;
	pruneprob = ppmove*pruneprob*(1-sit->second->psplit);
	growprob = pgmove*growprob*(sit->second->psplit)*(sit->second->prule);

	lfratio = log(growprob*pgchooseleaf);
	lrratio = log(pruneprob*ppchooseleaf);

	// Prune
	M.prune(tree, rand_int, MDIM); 
}

void Prop::chprop(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds, 
		double& lfratio, double& lrratio) {
	double pfthresh = 1.0; 
	double prthresh = 1.0; 
	map<int, Node*> subtree;
	int rand_int;
	if (int_nodes.empty()) {
		rand_int = 1;
	}
	else {
		rand_int = selector(int_nodes);
	}
	
	// Change random node partition
	T.find_subtree(tree, rand_int, subtree);
	prthresh = threshprob(subtree);
	subtree.clear();	

	M.change(tree, rand_int, names, alltypes, allcats, dbounds, ibounds);
	T.find_subtree(tree, rand_int, subtree);
	pfthresh = threshprob(subtree);

	// Calculate proposals
	lfratio = log(pfthresh);
	lrratio = log(prthresh);
}

bool Prop::check_siblings(map<int, Node*>& tree, const int& node) {
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

void Prop::swprop(map<int, Node*>& tree, const vector<int> int_nodes, double& lfratio, double& lrratio) {
	int rand_int, rand_swap;
	double pfthresh = 1.0; 
	double prthresh = 1.0; 
	map<int, Node*> subtree;

	cout << "Inside swprop int_nodes: ";
       	for(auto s: int_nodes) {cout << s << " ";}
	cout << endl;

	vector<int> parent_child = T.parent_child(int_nodes);
	cout << "Inside swprop parent_child: ";
       	for(auto s: parent_child) {cout << s << " ";}
	cout << endl;
	rand_int = selector(parent_child); // returns random parent node of parent_child pair

	T.find_subtree(tree, rand_int, subtree);
	prthresh = threshprob(subtree);
	subtree.clear();	

	if(find(int_nodes.begin(), int_nodes.end(), 2*rand_int) != int_nodes.end() & find(int_nodes.begin(), int_nodes.end(), 2*rand_int+1) != int_nodes.end()) {
		bool match = check_siblings(tree, rand_int);
		if(match) {
			M.mswap(tree, 2*rand_int, 2*rand_int+1);
		}
		else {
			rand_swap = Di.runif(0,1);
			if(rand_swap == 0) {
				M.mswap(tree, rand_int, 2*rand_int);
			}
			else{
				M.mswap(tree, rand_int, 2*rand_int+1);
			}
		}

	}
	else if(find(int_nodes.begin(), int_nodes.end(), 2*rand_int) != int_nodes.end()) {
		M.mswap(tree, rand_int, 2*rand_int);
	}
	else {
		M.mswap(tree, rand_int, 2*rand_int+1);
	}

	T.find_subtree(tree, rand_int, subtree);
	pfthresh = threshprob(subtree);

	// Calculate proposals
	lfratio = log(pfthresh);
	lrratio = log(prthresh);
}

void Prop::shftprop(map<int, Node*>& tree, const vector<int>& int_nodes, double& lfratio, double& lrratio) {
	int rand_int;
	double pfthresh = 1.0; 
	double prthresh = 1.0; 
	map<int, Node*> subtree;

	if (int_nodes.empty()) {
		rand_int = 1;
	}
	else {
		rand_int = selector(int_nodes);
	}
	
	// Change random node partition
	T.find_subtree(tree, rand_int, subtree);
	prthresh = allthreshprob(subtree);
	subtree.clear();	

	M.shift(tree, rand_int);
	T.find_subtree(tree, rand_int, subtree);
	pfthresh = allthreshprob(subtree);

	// Calculate proposals
	lfratio = log(pfthresh);
	lrratio = log(prthresh);
}

// Calculates Threshold probability for predictors the match root of subtree
double Prop:: threshprob(map<int, Node*>& subtree) {
	double prule = 1.0;
	int mnode;
	map<int, Node*>::const_iterator it;
	it = subtree.begin();
	while(true) {
		mnode = M.match_pred(subtree, it->second->nnode);
		if(mnode != 0) {
			it = subtree.find(mnode);
			if(it->second->isleaf == false) {
				if(it->second->ntype == Dataline::INT_T) {
					if (it->second->max == it->second->min) {
						prule = prule;
					}
					else {
						prule = prule*(1/double(it->second->max - it->second->min));
					}
				}
				else if(it->second->ntype == Dataline::DOUBLE_T) {
						prule = prule*((it->second->dsplit - it->second->min)/
								double(it->second->max - it->second->min));
				}
				else {
					prule = prule*(1/double(it->second->ncats));
				}
			}
		}
		else { break;}
	}
	return(prule);
}

// Calculates Threshold probability for changing every threshold in subtree
double Prop:: allthreshprob(map<int, Node*>& subtree) {
	double prule = 1.0;
	int mnode;
	map<int, Node*>::const_iterator it;
	it = subtree.begin();
	for(; it != subtree.end(); ++it) {
		if(it->second->isleaf == false) {
			if(it->second->ntype == Dataline::INT_T) {
				if (it->second->max == it->second->min) {
					prule = prule;
				}
				else {
					prule = prule*(1/double(it->second->max - it->second->min));
				}
			}
			else if(it->second->ntype == Dataline::DOUBLE_T) {
					prule = prule*((it->second->dsplit - it->second->min)/
							double(it->second->max - it->second->min));
			}
			else {
				prule = prule*(1/double(it->second->ncats));
			}
		}
	}
	return(prule);
}

void Prop::mgprop(map<int, Node*>& tree, double& lfratio, double& lrratio, const double pfmove, 
		const double prmove, const double ALPH, const double BET, const double pnotup) {
	map<int, Node*>::const_iterator it;
	size_t numnodes = tree.size();
	int numleaves;
	int rand_int;
	int numnotup = 0;
	double pnleaf;
	double pfthresh = 1.0;
	double prthresh = 1.0;
	double pint;
	double forwardprob = 1.0;
	double reverseprob = 1.0;
	vector<int> notup;
	vector<int> leaf_nodes;
	map<int, Node*> subtree;

	// Get leaves that can be pruned 
	T.find_leaf_nup(tree, notup, numnotup);
	if(numnotup > 0) {
		pnleaf = numnotup/double(numnodes);	
		rand_int = selector(notup);
		T.find_subtree(tree, rand_int/2, subtree);
		prthresh = threshprob(subtree);
		subtree.clear();

		M.merge(tree, rand_int, ALPH, BET);
		T.find_subtree(tree, rand_int/2, subtree);
		pfthresh = threshprob(subtree);
		T.find_leaf_nodes(tree, leaf_nodes, numleaves);
		it = tree.find(rand_int/2);
		pint = (numleaves-1)/numnodes;
	
		// Calculate proposals
		forwardprob = pfmove*pnotup*pnleaf*(1-it->second->psplit)*pfthresh;
		reverseprob = prmove*pint*(it->second->psplit)*(it->second->prule)*prthresh;
	}
	else {
		forwardprob = 1.0;
		reverseprob = 1.0;
	}

	lfratio = log(forwardprob);
	lrratio = log(reverseprob);
}

void Prop:: tswprop(map<int, Node*>& tree, double& lfratio, double& lrratio) {
	int rand_int;
	int numnotup = 0;
	double pfthresh = 1.0;
	double prthresh = 1.0;
	vector<int> notup;
	map<int, Node*> subtree;

	// Get leaves that can be pruned 
	T.find_leaf_nup(tree, notup, numnotup);
	if(numnotup != 0) {
		rand_int = selector(notup);
		T.find_subtree(tree, rand_int/2, subtree);
		prthresh = threshprob(subtree);
		subtree.clear();

		M.target_swap(tree, rand_int);
		T.find_subtree(tree, rand_int/2, subtree);
		pfthresh = threshprob(subtree);
	}

	// Calculate proposals
	lfratio = log(pfthresh);
	lrratio = log(prthresh);
}

void Prop::prngrwprop(map<int, Node*>& tree, const double& ALPH, const double& BET, const int& MDIM, const int& NDIM,
		const vector<int> int_nodes, const vector<int> leaf_nodes, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds,
		double& lfratio, double& lrratio) {
	int rand_int;
	double pfthresh;
	double prthresh;
	map<int, Node*>::iterator sit;

	// Get parent nodes of leaves that are end pairs
	vector<int> parent_pairs = T.parent_pairs(leaf_nodes);
	rand_int = selector(parent_pairs);
	sit = tree.find(rand_int);
	prthresh = sit->second->prule;

	// Prune
	M.prune(tree, rand_int, MDIM); 
	
	vector<int> upleaf_nodes;
	int numleaves;
	T.find_leaf_nodes(tree, upleaf_nodes, numleaves);
	int rand_leaf = selector(upleaf_nodes);

	// Grow the tree
	M.grow(tree, rand_leaf, ALPH, BET, MDIM, NDIM, names, alltypes, allcats, dbounds, ibounds);

	// Calculate proposals
	sit = tree.find(rand_leaf);
	pfthresh = sit->second->prule;

	lfratio = log(pfthresh);
	lrratio = log(prthresh);
}
