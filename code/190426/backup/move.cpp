// MOVE CLASS

#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include "priorf.h"
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
		const double BET, const int MDIM, const vector<string>& names, 
		map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
	void prune(map<int, Node*>& tree, const int& node, const int& MDIM);
	void get_subtree_info(map<int, Node*>& subtree, int node,
		Ref<VectorXd> svec, Ref<MatrixXd> svar, double& npost, double& nleaf);
	void del_subtree(map<int, Node*>& tree, int node, map<int, Node*>& subtree);
	void change(map<int, Node*>& tree, const int& node, const vector<string>& names, 
			const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
			const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds);
	void swapkeep(map<int, Node*>& tree, int node1, int node2);
	void mswap(map<int, Node*>& tree, int node1, int node2);
	vector<pair<string, int>>::const_iterator prev_it(const vector<pair<string, int>>& prevpreds, 
		const string& pred);
	void select_split(const string& pred, const map<int, Node*> tree, const int& node,
			const map<int, Node*>::const_iterator& it, const Dataline::type& intype,
			const vector<string>& names, const map<string, Dataline::type>& alltypes, 
			const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
			const map<string, vector<int>>& ibounds);
	void path_update(map<int, Node*> tree, const pair<string, int>& predpair);
};

void Move:: grow(map<int, Node*>& tree, const int node, const double ALPH, 
	const double BET, const int MDIM, const vector<string>& names, 
	map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
	const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
		
	map<int, Node*>::const_iterator it;
	it = tree.find(node);

	int left = 2*node;
	int right = 2*node+1;

	string lpred = Move::selector(names);
	string rpred = Move::selector(names);
	int numpred = names.size();

	Dataline::type ltype;
	Dataline::type rtype;
	map<string, Dataline::type>::const_iterator ltit;
	map<string, Dataline::type>::const_iterator rtit;
	ltit = alltypes.find(lpred);
	rtit = alltypes.find(rpred);
	ltype = ltit->second;
	rtype = rtit->second;
	
	tree.insert(it, pair<int, Node*>(left, new Node(left, MDIM)));
	tree.insert(it, pair<int, Node*>(right, new Node(right, MDIM)));
	map<int, Node*>::const_iterator lit;
	lit = tree.find(left);
	map<int, Node*>::const_iterator rit;
	rit = tree.find(right);


	//New left node
	select_split(lpred, tree, left, it, ltype, names, alltypes, allcats, dbounds, ibounds);
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
			lprule = (abs(lit->second->max - lit->second->min)/(double)lit->second->bwidth)*(1/double(numpred));
	}
	else {
		lprule = (1/double(lit->second->ncats))*(1/double(numpred));
	}
	lit->second->prule = lprule;

	lit->second->state = it->second->state;
	lit->second->svar = it->second->svar;
	lit->second->atmat = it->second->atmat;
	lit->second->aimat = it->second->aimat;
	lit->second->cvec = it->second->cvec;
	lit->second->dvec = it->second->dvec;
	lit->second->ppost = (it->second->ppost)+log(0.5);
	lit->second->lpost = 0.0;
	lit->second->isleaf = true;

	//New right node
	select_split(rpred, tree, right, it, rtype, names, alltypes, allcats, dbounds, ibounds);
	rit->second->depth = it->second->depth+1;
	rit->second->psplit = psplitf(ALPH, BET, rit->second->depth);
	double rprule;
	if( rit->second->ntype == Dataline::INT_T) {
		if (rit->second->max == rit->second->min) {
			rprule = (1/double(numpred));
		}
		else {
			rprule = (1/double(abs(rit->second->max - rit->second->min)))*(1/double(numpred));
		}
	}
	else if(rit->second->ntype == Dataline::DOUBLE_T) {
			rprule = (abs(rit->second->max - rit->second->min)/(double)rit->second->bwidth)*(1/double(numpred));
	}
	else {
		rprule = (1/double(rit->second->ncats))*(1/double(numpred));
	}
	rit->second->prule = rprule;

	rit->second->state = it->second->state;
	rit->second->svar = it->second->svar;
	rit->second->atmat = it->second->atmat;
	rit->second->aimat = it->second->aimat;
	rit->second->cvec = it->second->cvec;
	rit->second->dvec = it->second->dvec;
	rit->second->ppost = (it->second->ppost)+log(0.5);
	rit->second->lpost = 0.0;
	rit->second->isleaf = true;

	it->second->isleaf = false;
}

vector<pair<string, int>>::const_iterator Move:: prev_it(const vector<pair<string, int>>& prevpreds, 
		const string& pred) {
	vector<pair<string, int>>::const_reverse_iterator rit;
	for(rit = prevpreds.rbegin(); rit != prevpreds.rend(); ++rit) {
		if(rit->first == pred) {return (rit+1).base();}
	}
	return prevpreds.end();
}


void Move::select_split(const string& pred, const map<int, Node*> tree, const int& node,
		const map<int, Node*>::const_iterator& it, const Dataline::type& intype,
		const vector<string>& names, const map<string, Dataline::type>& alltypes, 
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
		const map<string, vector<int>>& ibounds) {
	string npred = pred;
	Dataline::type ntype;
	ntype = intype;
	map<int, Node*>::const_iterator lit;
	lit = tree.find(node);
	lit->second->pathpreds = it->second->pathpreds;

	if(node%2 == 0){lit->second->isleft = true;}
	else {lit->second->isright = true;}

	bool exists = true;
	while(exists) {
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
						npred = Move::selector(names);
						map<string, Dataline::type>::const_iterator rtit;
						rtit = alltypes.find(npred);
						ntype = rtit->second;
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
						exists = false;
					}
				}
				else {
					if(prevnode->lcats == "") {
						npred = Move::selector(names);
						map<string, Dataline::type>::const_iterator rtit;
						rtit = alltypes.find(npred);
						ntype = rtit->second;
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
						exists = false;
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
				exists = false;
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
			exists = false;
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
					npred = Move::selector(names);
					map<string, Dataline::type>::const_iterator rtit;
					rtit = alltypes.find(npred);
					ntype = rtit->second;
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
					exists = false;
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
				exists = false;
			}
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
	it->second->ppost = post;
	
	T.del_subtree(tree, node, subtree);	
}
		
void Move:: change(map<int, Node*>& tree, const int& node, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	map<int, Node*>::const_iterator it;
	it = tree.find(node);
	int numpred = names.size();

	if(!it->second->pathpreds.empty()) {
		it->second->pathpreds.pop_back();
	}

	string npred = Move::selector(names);
	it->second->pred = npred;

	Dataline::type ntype;
	map<string, Dataline::type>::const_iterator ntit;
	ntit = alltypes.find(npred);
	ntype = ntit->second;

	select_split(npred, tree, node, it, ntype, names, alltypes, allcats, dbounds, ibounds);
	it = tree.find(node);
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
			prule = (abs(it->second->max - it->second->min)/(double)it->second->bwidth)*(1/double(numpred));
	}
	else {
		prule = (1/double(it->second->ncats))*(1/double(numpred));
	}
	it->second->prule = prule;
}

void Move:: path_update(map<int, Node*> tree, const pair<string, int>& predpair) {
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

	swapkeep(tree, node1, node2);
	swap(*it1->second, *it2->second);

	pair<string, int> newpair1 = make_pair(it1->second->pred, node1);
	pair<string, int> newpair2 = make_pair(it2->second->pred, node2);
	path_update(tree, newpair1);
	path_update(tree, newpair2);

}

void Move:: swapkeep(map<int, Node*>& tree, int node1, int node2) {
	map<int, Node*>::const_iterator it1 = tree.find(node1);
	map<int, Node*>::const_iterator it2 = tree.find(node2);

	pair<string, int> n1pair = make_pair(it1->second->pred, node2);
	pair<string, int> n2pair = make_pair(it2->second->pred, node1);
	it1->second->pathpreds.pop_back();
	it2->second->pathpreds.pop_back();

	int nnodet = it1->second->nnode;
	int deptht = it1->second->depth;
	double psplitt = it1->second->psplit;
	bool isleaft = it1->second->isleaf;
	bool isleftt = it1->second->isleft;
	bool isrightt = it1->second->isright;
	vector<pair<string, int>> pathpredst = it1->second->pathpreds;


	it1->second->nnode = it2->second->nnode;
	it1->second->depth = it2->second->depth;
	it1->second->psplit = it2->second->psplit;
	it1->second->isleaf = it2->second->isleaf;
	it1->second->isleft = it2->second->isleft;
	it1->second->isright = it2->second->isright;
	it1->second->pathpreds = it2->second->pathpreds;
	it1->second->pathpreds.push_back(n1pair);
	

	it2->second->nnode = nnodet;
	it2->second->depth = deptht;
	it2->second->psplit = psplitt;
	it2->second->isleaf = isleaft;
	it2->second->isleft = isleftt;
	it2->second->isright = isrightt;
	it2->second->pathpreds = pathpredst;
	it2->second->pathpreds.push_back(n2pair);
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
}*/
