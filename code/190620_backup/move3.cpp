// MOVE2 CLASS

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
	void prune(map<int, Node*>& tree, const int& node, const int& MDIM);
	void change(map<int, Node*>& tree, const int& node, const vector<string>& names, 
			const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
			const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds);
	void mchange(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<string>& names, const int nchange, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds);
	void shift(const map<int, Node*> tree, const int& node);
	void swapkeep(map<int, Node*>& tree, int node1, int node2);
	void mswap(map<int, Node*>& tree, int node1, int node2);
	vector<pair<string, int>>::const_iterator prev_it(const vector<pair<string, int>>& prevpreds, 
		const string& pred);
	int select_split(const string& pred, const map<int, Node*> tree, const int& node,
			const Dataline::type& intype, const vector<string>& names, const map<string, Dataline::type>& alltypes, 
			const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
			const map<string, vector<int>>& ibounds);
	int shift_split(const map<int, Node*> subtree, const int& node);
	void path_update(map<int, Node*> tree, const pair<string, int>& predpair);
};

void Move:: grow(map<int, Node*>& tree, const int node, const double ALPH, 
	const double BET, const int MDIM, const int NDIM, const vector<string>& names, 
	map<string, Dataline::type> alltypes, const map<string, string>& allcats, 
	const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	//cout << "Inside grow: " << endl;
		
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

	tree.insert(it, pair<int, Node*>(left, new Node(left, MDIM, NDIM)));
	tree.insert(it, pair<int, Node*>(right, new Node(right, MDIM, NDIM)));
	map<int, Node*>::const_iterator lit;
	lit = tree.find(left);
	lit->second->pathpreds = it->second->pathpreds;
	map<int, Node*>::const_iterator rit;
	rit = tree.find(right);
	rit->second->pathpreds = it->second->pathpreds;

	int rightsplit, leftsplit;
	rightsplit = select_split(rpred, tree, right, rtype, names, alltypes, allcats, dbounds, ibounds);
	leftsplit = select_split(lpred, tree, left, ltype, names, alltypes, allcats, dbounds, ibounds);
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
		tree.erase(it);
		it = tree.find(right);
		tree.erase(it);
	}
}

void Move:: mchange(map<int, Node*>& tree, const vector<int>& int_nodes, const vector<string>& names, const int nchange, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	//cout << "Inside mchange: " << endl;
	for(int i = 0; i < nchange; ++i) {
		int nnode = selector(int_nodes);
	//	cout << "Inside mchange nnode: " << nnode << endl;
		change(tree, nnode, names, alltypes, allcats, dbounds, ibounds);
	}
}

void Move:: change(map<int, Node*>& tree, const int& node, const vector<string>& names,
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, 
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	//cout << "Inside change: " << endl;
	map<int, Node*>::const_iterator it;
	it = tree.find(node);
	int numpred = names.size();
	string cpred = it->second->pred;
		
	if(!it->second->pathpreds.empty()) {
		it->second->pathpreds.pop_back();
	}

	string npred = Move::selector(names);
	it->second->pred = npred;

	Dataline::type ntype;
	map<string, Dataline::type>::const_iterator ntit;
	ntit = alltypes.find(npred);
	ntype = ntit->second;

	int nsplit;
	nsplit = select_split(npred, tree, node, ntype, names, alltypes, allcats, dbounds, ibounds);
	if(nsplit == 1) {
		pair<string, int> newpair = make_pair(npred, node);
		path_update(tree, newpair);
		double prule;
		if(it->second->ntype == Dataline::INT_T) {
			if (it->second->max == it->second->min) {
				prule = (1/double(numpred));
			}
			else {
				prule = (1/double(it->second->max - it->second->min))*(1/double(numpred));
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
	}
	else {
		pair<string, int> oldpair = make_pair(cpred, node);
		it->second->pathpreds.push_back(oldpair);
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


int Move::select_split(const string& pred, const map<int, Node*> tree, const int& node,
		const Dataline::type& intype, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, 
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
		const map<string, vector<int>>& ibounds) {
	string npred = pred;
	//cout << "Inside select split npred: "<< npred << endl;
	Dataline::type ntype;
	ntype = intype;
	map<int, Node*>::const_iterator lit;
	lit = tree.find(node);
	//cout << "Inside select split: " << endl;
	//for(auto s : lit->second->pathpreds) {
	//       cout << s.first << " " << s.second  << " ";
	//}
	//cout << endl;	

	if(node%2 == 0){lit->second->isleft = true;}
	else {lit->second->isright = true;}

	if(ntype == Dataline::STRING_T) {
		lit->second->pred = npred;
		lit->second->ntype = ntype;
		vector<pair<string, int>>::const_iterator pit;
		pit = prev_it(lit->second->pathpreds, npred);
		//cout << "Testing pit: " << endl;
		if(pit != lit->second->pathpreds.end()) { //if pred previously chosen
			map<int, Node*>::const_iterator npit;
			npit = tree.find(pit->second);
			Node * prevnode = npit->second;
			//cout << "Testing prevnode nnode: " << prevnode->nnode << endl;
			if(lit->second->isright) {
				if(prevnode->rcats == "") {
					return 0;
				}
				else {
					pair<string, int> newpair = make_pair(npred, node);
					lit->second->pathpreds.push_back(newpair);
					lit->second->cats = prevnode->rcats;
					lit->second->ncats = prevnode->rcats.length();
					//cout << "Inside select split right cats: "<< lit->second->cats << endl;
					char nsplit = Di.rchar(lit->second->cats);
					lit->second->ssplit = nsplit;
					//cout << "Inside select split right nsplit "<< nsplit << endl;
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
					//cout << "Inside select split left cats: "<< lit->second->cats << endl;
					char nsplit = Di.rchar(lit->second->cats);
					lit->second->ssplit = nsplit;
					//cout << "Inside select split left nsplit "<< nsplit << endl;
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

void Move::shift(const map<int, Node*> tree, const int& node) {
	//cout << "Inside shift: " << endl;
	map<int, Node*> subtree;
	T.find_subtree(tree, node, subtree);
	int shift = shift_split(subtree, node);
}

int Move::shift_split(const map<int, Node*> subtree, const int& node) {
	//cout << "Inside shift_split: " << endl;
	map<int, Node*>::const_iterator it;
	it = subtree.find(node);

	// Get all previous preds by type including nodes
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
		//	cout << "Inside shift_split Di.rcar1: " << endl;
			char nsplit = Di.rchar(it->second->cats);
			it->second->ssplit = nsplit;
		//	cout << "Inside shift_split Di.part_cats: " << endl;
			pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
			it->second->lcats = lrcats.first;
			it->second->rcats = lrcats.second;
			it->second->nlcats = lrcats.first.length();
			it->second->nrcats = lrcats.second.length();
			nspred = cpred;
			nsnode = cnode;
			return 1;
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
					//cout << "Inside shift_split Di.rcar2: " << endl;
					char nsplit = Di.rchar(it->second->cats);
					it->second->ssplit = nsplit;
					//cout << "Inside shift_split Di.part_cats: " << endl;
					pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
					it->second->lcats = lrcats.first;
					it->second->rcats = lrcats.second;
					it->second->nlcats = lrcats.first.length();
					it->second->nrcats = lrcats.second.length();
					return 1;
				}
			}
			else {
				if(prevnode->lcats == "") {
					return 0;
				} 
				else {
					it->second->cats = prevnode->lcats;
					it->second->ncats = prevnode->lcats.length();
					//cout << "Inside shift_split Di.rcar3: " << endl;
					char nsplit = Di.rchar(it->second->cats);
					it->second->ssplit = nsplit;
					//cout << "Inside shift_split Di.part_cats: " << endl;
					pair<string, string> lrcats = Di.part_cats(nsplit, it->second->cats);
					it->second->lcats = lrcats.first;
					it->second->rcats = lrcats.second;
					it->second->nlcats = lrcats.first.length();
					it->second->nrcats = lrcats.second.length();
					return 1;
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
			//cout << "Inside shift_split Di.cunif: " << endl;
			it->second->dsplit = Di.cunif(it->second->min, it->second->max);
			ndpred = cpred;
			ndnode = cnode;
			return 1;
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
			//cout << "Inside shift_split Di.cunif: " << endl;
			it->second->dsplit = Di.cunif(it->second->min, it->second->max);
			return 1;
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
			//cout << "Inside shift_split Di.runif: " << endl;
			it->second->dsplit = Di.runif(it->second->min, it->second->max);
			nipred = cpred;
			ninode = cnode;
			return 1;
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
				//cout << "Inside shift_split Di.runif: " << endl;
				it->second->dsplit = Di.runif(it->second->min, it->second->max);
				return 1;
			}
		}
	}
}

void Move:: prune(map<int, Node*>& tree, const int& node, const int& MDIM) {
	//cout << "Inside prune: " << endl;
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

void Move:: merge(map<int, Node*>& tree, const int& node, const int& MDIM) {
	map<int, Node*>::const_iterator nit;
	map<int, Node*>::const_iterator sit;
	int parent = node/2;
	bool leftparent = true;
	bool rightparent = true;
	VectorXd substate(MDIM);
	substate = VectorXd::Zero(MDIM);
	MatrixXd subvar(MDIM, MDIM);
	subvar = MatrixXd::Zero(MDIM, MDIM);
	if (node%2 == 0) {
		sit = tree.find(node+1);
		if (sit->second->isleaf) {
			prune(tree, node, MDIM);
		}
		else {
			nit = tree.find(node);
			delete(nit->second);
			delete(sit->second);
			tree.erase(node);
			tree.erase(node+1);
			sit = tree.find(parent);
			++sit;
			for(sit; sit != tree.end(); ++sit) {
				if(sit->second->isleft) {
					sit->second->nnode = 2*leftparent;
					if(sit->second->isleaf = false){
						leftparent = true;
					}
				}
				else {
					sit->second->nnode = 2*rightparent+1;
					if(sit->second->isleaf = false){
						rightparent = 2*rightparent+1;
					}

			}
		}
	}
	else {
		tit = tree.find(node-1);
		if (tit->second->isleaf) {
			prune(tree, node, MDIM);
		}
		else {
		}
	}
}

		



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

void Move::growprune(

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
	//cout << "Inside swap: " << endl;
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
