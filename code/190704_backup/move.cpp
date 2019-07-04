// MOVE CLASS

#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include "prior.h"
#include "read.h"
#include "node.h"
#include "tree.h"
#include "dist.h"

#include "boost/variant.hpp"

using namespace std;
using namespace Eigen;

class Move {
	public:
	Tree T;
	Dist Di;
	Dataline D;
	Dist::random_selector<> selector; 

	// Constructors
	Move() {};

	// Declarations
	void grow(map<int, Node*>& tree, const int node, const double ALPH, 
		const double BET, const int MDIM, const int NDIM, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
	vector<pair<string, int>>::const_iterator prev_it(const vector<pair<string, int>>& prevpreds, 
		const string& pred);
	int select_split(map<int, Node*>& tree, const int node, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds);
	void prune(map<int, Node*>& tree, const int& node, const int& MDIM);
	void merge(map<int, Node*>& tree, const int& node, const double& ALPH, const double& BET); 
	void change(map<int, Node*>& tree, const int& node, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds);
	void path_update(map<int, Node*>& tree, const pair<string, int>& predpair);
	void mswap(map<int, Node*>& tree, int node1, int node2);
	void target_swap(map<int, Node*>& tree, const int& node);
	void swaprules(map<int, Node*>& tree, int node1, int node2);
	int match_pred(map<int, Node*>& tree, const int& node);
	int find_branch(map<int, Node*>& tree, const int& node, const int& mnode);
	void update_boundaries(map<int, Node*>& tree, const int& node);
	void shift(map<int, Node*>& tree, const int node);
	void shift_split(map<int, Node*>& subtree, const int node);
};

void Move:: grow(map<int, Node*>& tree, const int node, const double ALPH, 
	const double BET, const int MDIM, const int NDIM, const vector<string>& names, 
	map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
	const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
		
	map<int, Node*>::const_iterator it;
	it = tree.find(node);

	int left = 2*node;
	int right = 2*node+1;

	int numpred = names.size();

	tree.insert(it, pair<int, Node*>(left, new Node(left, MDIM, NDIM)));
	tree.insert(it, pair<int, Node*>(right, new Node(right, MDIM, NDIM)));
	map<int, Node*>::const_iterator lit;
	lit = tree.find(left);
	map<int, Node*>::const_iterator rit;
	rit = tree.find(right);

	int rightsplit, leftsplit;
	rightsplit = select_split(tree, right, names, alltypes, allcats, dbounds, ibounds);
	leftsplit = select_split(tree, left, names, alltypes, allcats, dbounds, ibounds);
	if (rightsplit*leftsplit == 1) {

		//New left node
		lit->second->depth = it->second->depth+1;
		lit->second->psplit = psplitf(ALPH, BET, lit->second->depth);
		double lprule;
		if(lit->second->ntype == Dataline::INT_T) {
			if (lit->second->max == lit->second->min) {
				lprule = (1/double(numpred));
			}
			else {
				lprule = (1/double(abs(lit->second->max - lit->second->min)))*(1/double(numpred));
			}
		}
		else if(lit->second->ntype == Dataline::DOUBLE_T) {
				lprule = ((lit->second->dsplit - lit->second->min)/(double)(lit->second->max - lit->second->min))*
					(1/double(numpred));
		}
		else {
			lprule = (1/double(lit->second->ncats))*(1/double(numpred));
		}
		lit->second->prule = lprule;

		lit->second->state = it->second->state;
		lit->second->svar = it->second->svar;
		lit->second->rmat = it->second->rmat;
		lit->second->atmat = it->second->atmat;
		lit->second->aimat = it->second->aimat;
		lit->second->adsqp = it->second->adsqp;
		lit->second->adetp = it->second->adetp;
		lit->second->HTVYp = it->second->HTVYp;
		lit->second->dvec = it->second->dvec;
		lit->second->lpost = it->second->lpost;
		lit->second->isleaf = true;

		//New right node
		rit->second->depth = it->second->depth+1;
		rit->second->psplit = psplitf(ALPH, BET, rit->second->depth);
		double rprule;
		if( rit->second->ntype == Dataline::INT_T) {
			if (rit->second->max == rit->second->min) {
				rprule = (1/double(numpred));
			}
			else {
				rprule = (1/double(rit->second->max - rit->second->min))*(1/double(numpred));
			}
		}
		else if(rit->second->ntype == Dataline::DOUBLE_T) {
				rprule = ((rit->second->dsplit - rit->second->min)/(double)(rit->second->max - rit->second->min))*
					(1/double(numpred));
		}
		else {
			rprule = (1/double(rit->second->ncats))*(1/double(numpred));
		}
		rit->second->prule = rprule;

		rit->second->state = it->second->state;
		rit->second->svar = it->second->svar;
		rit->second->rmat = it->second->rmat;
		rit->second->atmat = it->second->atmat;
		rit->second->aimat = it->second->aimat;
		rit->second->adsqp = it->second->adsqp;
		rit->second->adetp = it->second->adetp;
		rit->second->HTVYp = it->second->HTVYp;
		rit->second->dvec = it->second->dvec;
		rit->second->lpost = it->second->lpost;
		rit->second->isleaf = true;

		it->second->isleaf = false;
	}
	else {
		it = tree.find(left);
		delete it->second;
		tree.erase(it);
		it = tree.find(right);
		delete it->second;
		tree.erase(it);
	}
}

vector<pair<string, int>>::const_iterator Move:: prev_it(const vector<pair<string, int>>& prevpreds, 
		const string& pred) {
	vector<pair<string, int>>::const_reverse_iterator rit;
	for(rit = prevpreds.rbegin(); rit != prevpreds.rend(); ++rit) {
		if(rit->first == pred) {return (rit+1).base();}
	}
	return prevpreds.end();
}


int Move::select_split(map<int, Node*>& tree, const int node, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	string npred = Move::selector(names);
	Dataline::type ntype;
	map<string, Dataline::type>::const_iterator tit;
	tit = alltypes.find(npred);
	ntype = tit->second;

	map<int, Node*>::const_iterator lit;
	lit = tree.find(node);
	//lit->second->pathpreds = it->second->pathpreds;

	//if(node%2 == 0){lit->second->isleft = true;}
	//else {lit->second->isright = true;}

	if(ntype == Dataline::STRING_T) {
		lit->second->pred = npred;
		lit->second->ntype = ntype;
		vector<pair<string, int>>::const_iterator pit;
		pit = prev_it(lit->second->pathpreds, npred);
		if(pit != lit->second->pathpreds.end()) { //if pred previously chosen
			map<int, Node*>::const_iterator npit;
			npit = tree.find(pit->second);
			Node * prevnode = npit->second;
			if(lit->second->isright) {
				if(prevnode->rcats == "") {
					return 0;
				}
				else {
					pair<string, int> newpair = make_pair(npred, node);
					lit->second->pathpreds.push_back(newpair);
					lit->second->cats = prevnode->rcats;
					lit->second->ncats = prevnode->rcats.length();
					char nsplit = Di.rchar(lit->second->cats);
					lit->second->ssplit = nsplit;
					pair<string, string> lrcats = Di.part_cats(nsplit, lit->second->cats);
					lit->second->lcats = lrcats.first;
					lit->second->rcats = lrcats.second;
					lit->second->nlcats = lrcats.first.length();
					lit->second->nrcats = lrcats.second.length();
					lit->second->min = 0;
					lit->second->max = 0;
					lit->second->dsplit = 0;
					lit->second->bwidth = 0;
					return 1;
				}
			}
			else {
				if(prevnode->lcats == "") {
					return 0;
				} 
				else {
					pair<string, int> newpair = make_pair(npred, node);
					lit->second->pathpreds.push_back(newpair);
					lit->second->cats = prevnode->lcats;
					lit->second->ncats = prevnode->lcats.length();
					char nsplit = Di.rchar(lit->second->cats);
					lit->second->ssplit = nsplit;
					pair<string, string> lrcats = Di.part_cats(nsplit, lit->second->cats);
					lit->second->lcats = lrcats.first;
					lit->second->rcats = lrcats.second;
					lit->second->nlcats = lrcats.first.length();
					lit->second->nrcats = lrcats.second.length();
					lit->second->min = 0;
					lit->second->max = 0;
					lit->second->dsplit = 0;
					lit->second->bwidth = 0;
					return 1;
				}
			}
		}
		else {	
			pair<string, int> newpair = make_pair(npred, node);
			lit->second->pathpreds.push_back(newpair);
			map<string, string>::const_iterator sit;
			sit = allcats.find(npred);
			lit->second->cats = sit->second;
			lit->second->ncats = sit->second.length();
			char nsplit = Di.rchar(lit->second->cats);
			lit->second->ssplit = nsplit;
			pair<string, string> lrcats = Di.part_cats(nsplit, lit->second->cats);
			lit->second->lcats = lrcats.first;
			lit->second->rcats = lrcats.second;
			lit->second->nlcats = lrcats.first.length();
			lit->second->nrcats = lrcats.second.length();
			lit->second->min = 0;
			lit->second->max = 0;
			lit->second->dsplit = 0;
			lit->second->bwidth = 0;
			return 1;
		}
	}
	else if(ntype == Dataline::DOUBLE_T) {
		lit->second->pred = npred;
		lit->second->ntype = ntype;
		vector<pair<string, int>>::const_iterator pit;
		pit = prev_it(lit->second->pathpreds, npred);
		if(pit !=  lit->second->pathpreds.end()) { //if pred previously chosen
			map<int, Node*>::const_iterator npit;
			npit = tree.find(pit->second);
			Node * prevnode = npit->second;
			pair<string, int> newpair = make_pair(npred, node);
			lit->second->pathpreds.push_back(newpair);
			lit->second->bwidth = prevnode->bwidth;
			if(lit->second->isright) {
				lit->second->min = prevnode->dsplit;
				lit->second->max = prevnode->max;
			}
			else {
				lit->second->min = prevnode->min;
				lit->second->max = prevnode->dsplit;
			}
			lit->second->dsplit = Di.cunif(lit->second->min, lit->second->max);
			lit->second->cats = "";
			lit->second->lcats = "";
			lit->second->rcats = "";
			lit->second->ssplit = "";
		}
		else {	
			pair<string, int> newpair = make_pair(npred, node);
			lit->second->pathpreds.push_back(newpair);
			map<string, vector<double>>::const_iterator dit;
			dit = dbounds.find(npred);
			lit->second->min = dit->second[0];
			lit->second->max = dit->second[1];
			lit->second->bwidth = abs(lit->second->max - lit->second->min);
			lit->second->dsplit = Di.cunif(lit->second->min, lit->second->max);
			lit->second->cats = "";
			lit->second->lcats = "";
			lit->second->rcats = "";
			lit->second->ssplit = "";
		}
		return 1;
	}
	else {
		lit->second->pred = npred;
		lit->second->ntype = ntype;
		vector<pair<string, int>>::const_iterator pit;
		pit = prev_it(lit->second->pathpreds, npred);
		if(pit !=  lit->second->pathpreds.end()) { //if pred previously chosen
			map<int, Node*>::const_iterator npit;
			npit = tree.find(pit->second);
			Node * prevnode = npit->second;
			if((prevnode->max - prevnode->min) == 0) {
				return 0;
			}
			else {
				pair<string, int> newpair = make_pair(npred, node);
				lit->second->pathpreds.push_back(newpair);
				if(lit->second->isright) {
					lit->second->min = prevnode->dsplit;
					lit->second->max = prevnode->max;
				}
				else {
					lit->second->min = prevnode->min;
					lit->second->max = prevnode->dsplit;
				}
				lit->second->dsplit = Di.runif(lit->second->min, lit->second->max);
				lit->second->cats = "";
				lit->second->lcats = "";
				lit->second->rcats = "";
				lit->second->ssplit = "";
				return 1;
			}
		}
		else {	
			pair<string, int> newpair = make_pair(npred, node);
			lit->second->pathpreds.push_back(newpair);
			map<string, vector<int>>::const_iterator iit;
			iit = ibounds.find(npred);
			lit->second->min = iit->second[0];
			lit->second->max = iit->second[1];
			lit->second->dsplit = Di.runif(lit->second->min, lit->second->max);
			lit->second->cats = "";
			lit->second->lcats = "";
			lit->second->rcats = "";
			lit->second->ssplit = "";
			return 1;
		}
	}
}


void Move:: prune(map<int, Node*>& tree, const int& node, const int& MDIM) {
	VectorXd substate(MDIM);
	substate = VectorXd::Zero(MDIM);
	MatrixXd subvar(MDIM, MDIM);
	subvar = MatrixXd::Zero(MDIM, MDIM);
	double post = 0.0;
	double nl = 1.0;
	
	map<int, Node*> subtree;
	T.find_subtree(tree, node, subtree);
	T.get_subtree_info(subtree, node, substate, subvar, post, nl);

	map<int, Node*>::const_iterator it;
	it = tree.find(node);
		
	it->second->isleaf = true;
	it->second->state = substate;
	it->second->svar = subvar;
	it->second->lpost = post;
	
	T.del_subtree(tree, node, subtree);	
}

void Move:: merge(map<int, Node*>& tree, const int& node, const double& ALPH, const double& BET) {
	int parent = node/2;
	int sib;
	if(node%2 == 0) {sib = node+1;}
	else{sib = node-1;}
	cout << "Inside merge node, parent, sib: " << node << " " << parent << " " << sib << endl;

	cout << "Inside merge tree.display at start: " << endl;
	T.display(tree);
	cout << endl;

	// Copy current tree and then clear tree
	map<int, Node*> copytree;
	for(auto s : tree) {
		copytree[s.first] = new Node;
		*copytree[s.first] = *tree[s.first];
	}
	tree.clear();

	// Get vector of all nodes
	vector<int> allnodes;
	T.find_nodes(copytree, allnodes);

	// Get subtree nodes
	vector<int> subtreenodes;
	T.find_subtree_nodes(copytree, parent, subtreenodes); 
	sort(subtreenodes.begin(), subtreenodes.end());

	// Remove subtree nodes from all nodes
	vector<int> diff;
	set_difference(allnodes.begin(), allnodes.end(), subtreenodes.begin(), subtreenodes.end(), 
			inserter(diff, diff.begin()));
	
	// Remove parent and node
	subtreenodes.erase(find(subtreenodes.begin(), subtreenodes.end(), node));
	subtreenodes.erase(find(subtreenodes.begin(), subtreenodes.end(), parent));

	// Modify subtree node numbers to replace sib, node, parent and children
	// Place new node as key in tree and copy contents across
	int nodediff = sib - parent;
	int newnode;
	double floorlog;
	vector<int> newnodes;
	for(int s : subtreenodes) {
		floorlog = (floor(log2(s))-floor(log2(sib)));
		newnode = s-nodediff*pow(2, floorlog);
		tree[newnode] = new Node;
		*tree[newnode] = *copytree[s];
		newnodes.push_back(newnode);
	}

	cout << "Inside merge, allnodes, subtreenodes, diff, newnodes: " << endl;
	for(auto s : allnodes){cout << s << " ";}
	cout << endl;
	for(auto s : subtreenodes){cout << s << " ";}
	cout << endl;
	for(auto s : diff){cout << s << " ";}
	cout << endl;
	for(auto s : newnodes){cout << s << " ";}
	cout << endl;

	// Inserts the rest of the nodes
	for(int s : diff) {
		tree[s] = new Node;
		*tree[s] = *copytree[s];
	}

	// Deletes copytree
	for(auto s : copytree) {
		delete s.second;
		copytree.erase(s.first);
	}
	
	// Update the new subtree
	map<int, Node*>::const_iterator it;
	for(int s : newnodes) {
		it = tree.find(s);
		it->second->nnode = s;
		it->second->depth -= 1;
		it->second->psplit = psplitf(ALPH, BET, it->second->depth);
		if(s == 1) {
			it->second->isleft = false;
			it->second->isright = false;
		}
		if(s%2 == 0) {
			it->second->isleft = true;
			it->second->isright = false;
		}
		else {
			it->second->isleft = false;
			it->second->isright = true;
		}
	}
	update_boundaries(tree, parent);

	cout << "Inside merge tree.display after: " << endl;
	T.display(tree);
	cout << endl;
}
		
void Move:: change(map<int, Node*>& tree, const int& node, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	map<int, Node*>::const_iterator it;
	it = tree.find(node);
	int numpred = names.size();
	string cpred = it->second->pred;

	if(!it->second->pathpreds.empty()) {
		it->second->pathpreds.pop_back();
	}

	int nsplit = select_split(tree, node, names, alltypes, allcats, dbounds, ibounds);
	if(nsplit == 1) {
		pair<string, int> newpair = make_pair(it->second->pred, node);
		path_update(tree, newpair);
		double prule;
		if(it->second->ntype == Dataline::INT_T) {
			if (it->second->max == it->second->min) {
				prule = (1/double(numpred));
			}
			else {
				prule = (1/double(abs(it->second->max - it->second->min)))*(1/double(numpred));
			}
		}
		else if(it->second->ntype == Dataline::DOUBLE_T) {
				prule = ((it->second->dsplit - it->second->min)/(double)(it->second->max - it->second->min))*
					(1/double(numpred));
		}
		else {
			prule = (1/double(it->second->ncats))*(1/double(numpred));
		}
		it->second->prule = prule;
		update_boundaries(tree, node);
	}
	else {
		pair<string, int> oldpair = make_pair(cpred, node);
		it->second->pathpreds.push_back(oldpair);
	}
}

void Move:: path_update(map<int, Node*>& tree, const pair<string, int>& predpair) {
	map<int, Node*>::const_iterator it;
	for(it = tree.find(predpair.second); it != tree.end(); ++it) {
		for(int i = 0; i < it->second->pathpreds.size(); ++i) {
			if(it->second->pathpreds[i].second == predpair.second) {
				it->second->pathpreds[i].first = predpair.first;
			}
		}
	}
}
		
void Move:: mswap(map<int, Node*>& tree, int node1, int node2) {
	map<int, Node*>::const_iterator it1 = tree.find(node1);
	map<int, Node*>::const_iterator it2 = tree.find(node2);

	swaprules(tree, node1, node2);
	pair<string, int> newpair1 = make_pair(it1->second->pred, node1);
	pair<string, int> newpair2 = make_pair(it2->second->pred, node2);
	path_update(tree, newpair1);
	path_update(tree, newpair2);

}

void Move:: target_swap(map<int, Node*>& tree, const int& node) {
	cout << "Inside target_swap: " << node << " " << node/2 << endl;
	int parent = node/2;
	mswap(tree, node, parent);
}

void Move:: swaprules(map<int, Node*>& tree, int node1, int node2) {
	map<int, Node*>::const_iterator it1 = tree.find(node1);
	map<int, Node*>::const_iterator it2 = tree.find(node2);

	pair<string, int> n1pair = make_pair(it1->second->pred, node1);
	pair<string, int> n2pair = make_pair(it2->second->pred, node2);
	it1->second->pathpreds.pop_back();
	it2->second->pathpreds.pop_back();

	cout << "Inside swaprules isleaf node1, node2: " << it1->second->isleaf << " " << it2->second->isleaf << endl;

	string nsplit;
	string npred;
	double ndsplit;
	double nprule;

	if(it1->second->ntype == Dataline::STRING_T) {
		nsplit = it1->second->ssplit;
		npred = it1->second->pred;
		nprule = it1->second->prule;
	}
	else {
		ndsplit = it1->second->dsplit;
		npred = it1->second->pred;
		nprule = it1->second->prule;
	}

	if(it2->second->ntype == Dataline::STRING_T && it1->second->ntype == Dataline::STRING_T) {
		it1->second->ssplit = it2->second->ssplit;
		it1->second->pred = it2->second->pred;
		it1->second->prule = it2->second->prule;
		it2->second->ssplit = nsplit;
		it2->second->pred = npred;
		it2->second->prule = nprule;
	}
	else if(it2->second->ntype == Dataline::STRING_T && it1->second->ntype == Dataline::DOUBLE_T) {
		it1->second->dsplit = it2->second->dsplit;
		it1->second->pred = it2->second->pred;
		it1->second->prule = it2->second->prule;
		it2->second->ssplit = nsplit;
		it2->second->pred = npred;
		it2->second->prule = nprule;
	}
	else if(it2->second->ntype == Dataline::STRING_T && it1->second->ntype == Dataline::INT_T) {
		it1->second->dsplit = it2->second->dsplit;
		it1->second->pred = it2->second->pred;
		it1->second->prule = it2->second->prule;
		it2->second->ssplit = nsplit;
		it2->second->pred = npred;
		it2->second->prule = nprule;
	}
	else if(it1->second->ntype == Dataline::STRING_T && it2->second->ntype == Dataline::DOUBLE_T) {
		it1->second->ssplit = it2->second->ssplit;
		it1->second->pred = it2->second->pred;
		it1->second->prule = it2->second->prule;
		it2->second->dsplit = ndsplit;
		it2->second->pred = npred;
		it2->second->prule = nprule;
	}
	else if(it1->second->ntype == Dataline::STRING_T && it2->second->ntype == Dataline::INT_T) {
		it1->second->ssplit = it2->second->ssplit;
		it1->second->pred = it2->second->pred;
		it1->second->prule = it2->second->prule;
		it2->second->dsplit = ndsplit;
		it2->second->pred = npred;
		it2->second->prule = nprule;
	}
	else {
		it1->second->dsplit = it2->second->dsplit;
		it1->second->pred = it2->second->pred;
		it1->second->prule = it2->second->prule;
		it2->second->dsplit = ndsplit;
		it2->second->pred = npred;
		it2->second->prule = nprule;
	}


	vector<pair<string, int>> pathpredst = it1->second->pathpreds;


	it1->second->pathpreds = it2->second->pathpreds;
	it1->second->pathpreds.push_back(n1pair);
	

	it2->second->pathpreds = pathpredst;
	it2->second->pathpreds.push_back(n2pair);

	if(node1 == node2+1 || node2 == node1+1) {
		update_boundaries(tree, node1);
		update_boundaries(tree, node2);
	}
	else if(node1 == 2*node2) {
		update_boundaries(tree, node1);
	}
	else if(node2 == 2*node1) {
		update_boundaries(tree, node2);
	}
	else if(node1 == 2*node2+1) {
		update_boundaries(tree, node1);
	}
	else if(node2 == 2*node1+1) {
		update_boundaries(tree, node2);
	}
}

int Move::match_pred(map<int, Node*>& tree, const int& node) { // returns first instance of matching predictors
	map<int, Node*>::const_iterator it;
	map<int, Node*>::const_iterator lit;
	it = tree.find(node);
	string pred = it->second->pred;
	lit = tree.find(node);
	++lit;
	if(lit == tree.end()) {return 0;}
	else {
		for(; lit != tree.end(); ++lit) {
			if(lit->second->pred == pred) {return(lit->second->nnode);}
			else {return 0;}
		}
	}
}

int Move::find_branch(map<int, Node*>& tree, const int& node, const int& mnode) {
	map<int, Node*>::const_iterator it;
	it = tree.find(node);
	bool isleft = T.is_node(tree, 2*node); // checks if mnode exists on left
	bool isright = T.is_node(tree, 2*node+1); // checks if right child node exists
	if(isleft) {
		if(2*node == mnode) {return(1);}
		else{find_branch(tree, 2*node, mnode); }
	}
	if(isright) {
		if(2*node+1 == mnode) {return(0);}
		else{find_branch(tree, 2*node+1, mnode); }
	}
}

void Move::update_boundaries(map<int, Node*>& tree, const int& node) {
	map<int, Node*>::const_iterator it;
	map<int, Node*>::const_iterator lit;
	map<int, Node*>::const_iterator rit;
	it = tree.find(node);
	string cpred = it->second->pred;
	Dataline::type ntype = it->second->ntype;
	int mnode = match_pred(tree, node); // Returns the first match of predictors in child tree
	cout << "Inside update boundaries node, mnode: " << node << " " << mnode << endl;

	if(mnode != 0) {
		int branch = find_branch(tree, node, mnode);
		if(branch == 1) {
			lit = tree.find(mnode);
			if(ntype == Dataline::STRING_T) {
				lit->second->cats = it->second->lcats;
				lit->second->ncats = it->second->nlcats;
				char nsplit = Di.rchar(lit->second->cats);
				lit->second->ssplit = nsplit;
				pair<string, string> lrcats = Di.part_cats(nsplit, lit->second->cats);
				lit->second->lcats = lrcats.first;
				lit->second->rcats = lrcats.second;
				lit->second->nlcats = lrcats.first.length();
				lit->second->nrcats = lrcats.second.length();
				update_boundaries(tree, mnode);
			}
			else if(ntype == Dataline::DOUBLE_T) {
				lit->second->min = it->second->min;
				lit->second->max = it->second->dsplit;
				lit->second->dsplit = Di.cunif(lit->second->min, lit->second->max);
				update_boundaries(tree, mnode);
			}
			else {
				lit->second->min = it->second->min;
				lit->second->max = it->second->dsplit;
				lit->second->dsplit = Di.runif(lit->second->min, lit->second->max);
				update_boundaries(tree, mnode);
			}
		}
		else {
			rit = tree.find(mnode);
			if(ntype == Dataline::STRING_T) {
				rit->second->cats = it->second->lcats;
				rit->second->ncats = it->second->nlcats;
				char nsplit = Di.rchar(rit->second->cats);
				rit->second->ssplit = nsplit;
				pair<string, string> lrcats = Di.part_cats(nsplit, rit->second->cats);
				rit->second->lcats = lrcats.first;
				rit->second->rcats = lrcats.second;
				rit->second->nlcats = lrcats.first.length();
				rit->second->nrcats = lrcats.second.length();
				update_boundaries(tree, mnode);
			}
			else if(ntype == Dataline::DOUBLE_T) {
				rit->second->min = it->second->min;
				rit->second->max = it->second->dsplit;
				rit->second->dsplit = Di.cunif(rit->second->min, rit->second->max);
				update_boundaries(tree, mnode);
			}
			else {
				rit->second->min = it->second->min;
				rit->second->max = it->second->dsplit;
				rit->second->dsplit = Di.runif(rit->second->min, rit->second->max);
				update_boundaries(tree, mnode);
			}
		}
	}
}

void Move::shift(map<int, Node*>& tree, const int node) {
	map<int, Node*> subtree;
	T.find_subtree(tree, node, subtree);
	cout << "Inside shift subtree: ";
	T.display(subtree);
	cout << endl;
	shift_split(subtree, node);
	T.del_subtree(tree, node, subtree);
	T.copyin_subtree(tree, subtree);
}

void Move::shift_split(map<int, Node*>& subtree, const int node) {
	Dataline::type ntype;
	map<int, Node*>::const_iterator it;
	it = subtree.find(node);

	for(it; it != subtree.end(); ++it) {
		ntype = it->second->ntype;
		string cpred = it->second->pred;
		int cnode = it->second->nnode;
		if (ntype == Dataline::STRING_T) {
			char nsplit = Di.rchar(it->second->cats);
			it->second->ssplit = nsplit;
			pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
			it->second->lcats = lrcats.first;
			it->second->rcats = lrcats.second;
			it->second->nlcats = lrcats.first.length();
			it->second->nrcats = lrcats.second.length();
			update_boundaries(subtree, cnode);
		}
		else if (ntype == Dataline::DOUBLE_T) {
			double ndsplit = Di.cunif(it->second->min, it->second->max);
			cout << "Inside shift_split cnode, ndsplit: " << cnode << " " << ndsplit << endl;
			it->second->dsplit = ndsplit;
			update_boundaries(subtree, cnode);
		}
		else {
			it->second->dsplit = Di.runif(it->second->min, it->second->max);
			update_boundaries(subtree, cnode);
		}
	}
}

/*
int main() {
	ifstream tin("testdata.txt");
	Dataline Data(9);
	Data.make_names();
	string typestring;
	getline(tin, typestring);
	string boundstring;
	getline(tin, boundstring);
	Data.read_alltypes(typestring);
	Data.read_bounds(boundstring);
	
	Tree T(1,1);
	T.D = Data;
	T.init_tree(T.tree);

	Move M;
	M.D = Data;
	M.grow(T.tree, 1, 0.95, 1, 0.5, 1); 
//	T.display(T.tree);
	M.grow(T.tree, 2, 0.95, 1, 0.5, 1); 
//	T.display(T.tree);
	M.grow(T.tree, 4, 0.95, 1, 0.5, 1); 
//	T.display(T.tree);
	M.change(T.tree, 4); 
	M.grow(T.tree, 8, 0.95, 1, 0.5, 1); 
//	T.display(T.tree);
	M.change(T.tree, 8); 
	M.grow(T.tree, 16, 0.95, 1, 0.5, 1); 
	M.change(T.tree, 16); 

	T.display(T.tree);
	cout << endl;

	int fnode;
	string newstring;
	getline(tin, newstring);
	Data.dline = Data.read_line(newstring);
	T.find_uleaf(T.tree, Data.dline, fnode);
	for(auto s : Data.dline) {
		cout << s.first << " " << Data.get_value(s.first, Data.dline) << " ";
	}
	cout << endl;
	cout << "found: " << fnode << endl;

	M.prune(T.tree, 16, 1);
	M.mswap(T.tree, 2, 4);

	T.display(T.tree);
	cout << endl;

	// OLD MERGE
	// Get node - the one to be removed
	map<int, Node*>::const_iterator it;
	it = tree.find(node);

	// Get parent
	map<int, Node*>::const_iterator pit;
	pit = tree.find(parent);
	cout << "Inside merge pit->nnode: " << pit->second->nnode << endl;

	// Get subtree nodes
	vector<int> subtreenodes;
	T.find_subtree_nodes(tree, parent, subtreenodes); 
	cout << "Inside merge subtreenodes before: " << endl;
	for(auto s : subtreenodes){cout << s << " ";}
	cout << endl;

	// Get sibling - the one to be kept
	map<int, Node*>::const_iterator sit;
	sit = tree.find(sib);
	cout << "Inside merge sit->nnode: " << sit->second->nnode << endl;

	// Delete the parent and node from the tree
	delete it->second;
	delete pit->second;
	tree.erase(it);
	tree.erase(pit);

	// Move sibling to parent
	tree[parent] = sit->second;
	tree.erase(sit);
	cout << "Inside merge tree.display after move sibilng to parent: " << endl;
	T.display(tree);
	cout << endl;
	cout << "Inside merge subtreenodes after: " << endl;
	for(auto s : subtreenodes){cout << s << " ";}
	cout << endl;

	// Update new parent node (old left or right leaf)
	it = tree.find(parent);
	it->second->depth = it->second->depth-1;
	it->second->psplit = psplitf(ALPH, BET, it->second->depth);
	it->second->isleaf = false;

	// Renumber subtree nodes and keys
	size_t nsub = subtreenodes.size();
	if(nsub > 1) {
		for(int i = 0; i < nsub; ++i) {
			int cnode = subtreenodes[i];
			it = tree.find(cnode);
			if(cnode%2 == 0) {
				tree[2*parent] = it->second;
				sit = tree.find(cnode/2);
				sit->second->nnode = cnode/2;
				sit->second->depth = sit->second->depth-1;
				sit->second->psplit = psplitf(ALPH, BET, sit->second->depth);
			}
			else{
				tree[(2*parent+1] = it->second;
				sit = tree.find((cnode/2)+1);
				sit->second->nnode = (cnode/2)+1;
				sit->second->depth = sit->second->depth-1;
				sit->second->psplit = psplitf(ALPH, BET, sit->second->depth);
			}
			tree.erase(it);

int Move::shift_split(const map<int, Node*> subtree, const int& node) {
	map<int, Node*>::const_iterator it;
	it = subtree.find(node);

	// Search subtree by type and for each type get pair (pred, node)
	vector<pair<string, int>> stringpreds;
	vector<pair<string, int>> doublepreds;
	vector<pair<string, int>> intpreds;
	Dataline::type ntype;

	for(it; it != subtree.end(); ++it) {
		ntype = it->second->ntype;
		string cpred = it->second->pred;
		int cnode = it->second->nnode;
		pair<string, int> newpair = make_pair(cpred, cnode);
		if (ntype == Dataline::STRING_T) {
			stringpreds.push_back(newpair);
		}
		else if (ntype == Dataline::DOUBLE_T) {
			doublepreds.push_back(newpair);
		}
		else {
			intpreds.push_back(newpair);
		}
	}

	// For each vector of types
	string nspred;
	int nsnode;
	vector<pair<string, int>>::const_iterator sit;
	for(sit = stringpreds.begin(); sit != stringpreds.end(); ++sit) {
		string cpred = sit->first;
		int cnode = sit->second;
		if(cpred != nspred) {
			it = subtree.find(cnode);
			if(it->second->cats == "") {
				return 0;
			}
			char nsplit = Di.rchar(it->second->cats);
			it->second->ssplit = nsplit;
			pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
			it->second->lcats = lrcats.first;
			it->second->rcats = lrcats.second;
			it->second->nlcats = lrcats.first.length();
			it->second->nrcats = lrcats.second.length();
			nspred = cpred;
			nsnode = cnode;
			//return 1;
		}
		else {
			it = subtree.find(cnode);
			map<int, Node*>::const_iterator pit;
			pit = subtree.find(nsnode);
			Node * prevnode = pit->second;
			if(it->second->isright) {
				if(prevnode->rcats == "") {
					return 0;
				}
				else {
					it->second->cats = prevnode->rcats;
					it->second->ncats = prevnode->rcats.length();
					char nsplit = Di.rchar(it->second->cats);
					it->second->ssplit = nsplit;
					pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
					it->second->lcats = lrcats.first;
					it->second->rcats = lrcats.second;
					it->second->nlcats = lrcats.first.length();
					it->second->nrcats = lrcats.second.length();
					//return 1;
				}
			}
			else {
				if(prevnode->lcats == "") {
					return 0;
				} 
				else {
					it->second->cats = prevnode->lcats;
					it->second->ncats = prevnode->lcats.length();
					char nsplit = Di.rchar(it->second->cats);
					it->second->ssplit = nsplit;
					pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
					it->second->lcats = lrcats.first;
					it->second->rcats = lrcats.second;
					it->second->nlcats = lrcats.first.length();
					it->second->nrcats = lrcats.second.length();
					//return 1;
				}
			}
		}
	}

	string ndpred;
	int ndnode;
	vector<pair<string, int>>::const_iterator dit;
	for(dit = doublepreds.begin(); dit != doublepreds.end(); ++dit) {
		string cpred = dit->first;
		int cnode = dit->second;
		if(cpred != ndpred) {
			it = subtree.find(cnode);
			it->second->dsplit = Di.cunif(it->second->min, it->second->max);
			ndpred = cpred;
			ndnode = cnode;
			//return 1;
		}
		else {
			it = subtree.find(cnode);
			map<int, Node*>::const_iterator pit;
			pit = subtree.find(ndnode);
			Node * prevnode = pit->second;
			if(it->second->isright) {
				it->second->min = prevnode->dsplit;
				it->second->max = prevnode->max;
			}
			else {
				it->second->min = prevnode->min;
				it->second->max = prevnode->dsplit;
			}
			it->second->dsplit = Di.cunif(it->second->min, it->second->max);
			//return 1;
		}
	}

	string nipred;
	int ninode;
	vector<pair<string, int>>::const_iterator iit;
	for(iit = doublepreds.begin(); iit != doublepreds.end(); ++iit) {
		string cpred = dit->first;
		int cnode = dit->second;
		if(cpred != nipred) {
			it = subtree.find(cnode);
			it->second->dsplit = Di.runif(it->second->min, it->second->max);
			nipred = cpred;
			ninode = cnode;
			//return 1;
		}
		else {
			it = subtree.find(cnode);
			map<int, Node*>::const_iterator pit;
			pit = subtree.find(ndnode);
			Node * prevnode = pit->second;
			if((prevnode->max - prevnode->min) == 0) {
				return 0;
			}
			else {
				if(it->second->isright) {
					it->second->min = prevnode->dsplit;
					it->second->max = prevnode->max;
				}
				else {
					it->second->min = prevnode->min;
					it->second->max = prevnode->dsplit;
				}
				it->second->dsplit = Di.runif(it->second->min, it->second->max);
				//return 1;
			}
		}
	}
	return 1;
}
}*/
