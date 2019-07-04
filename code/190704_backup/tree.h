//TREE CLASS

# pragma once
# include <iostream>
# include <string>
# include <Eigen/Dense>
# include <vector>
# include <algorithm>
# include <typeinfo>
# include <limits>
# include <set>
# include "node.h"
# include "prior.h"
# include "read.h"
# include "dist.h"

# include "boost/variant.hpp"

using namespace std;
using namespace Eigen;


class Tree {
	public:
		Dataline Data;
		Dist Di;
		map<int, Node*> tree;
		int tname = 0;
		int threadid = 0;
		int numleaves = 0;
		int move = 0;
		int level = 0;
		int nlevel = 0;
		int accept = 0;
		int found = 0;
		double fkprior = 0.0;
		double prior = 1.0;
		double lfratio = 0.0;
		double lrratio = 0.0;
		double logtreepost = 0.0;
		double slpostn = 0.0;
		double qstar = 0.0;
		double wdet = 0;
		double vdet = 0;
		MatrixXd WT;
		MatrixXd WI;
		MatrixXd FTWI;
		MatrixXd WIF;
		MatrixXd FTWIF;
	
		// Constructors
		Tree() {};	
		Tree(const int& nnode, const int& MDIMS, const int& NDIMS) : WT(MDIMS, MDIMS), WI(MDIMS, MDIMS),
      			FTWI(MDIMS, MDIMS), WIF(MDIMS, MDIMS), FTWIF(MDIMS, MDIMS){
			tree[1] = new Node(nnode, MDIMS, NDIMS);
		};

		void init_tree(map<int, Node*>& tree, const int& DATALENGTH, const vector<string>& names,
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats,
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
		void find_uleaf(const map<int, Node*>& tree, Dataline::line fline, int& found);
		void display(const map<int, Node*>& tree);
		void find_nodes(const map<int, Node*>& tree, vector<int>& outvec);
		void find_leaf_nodes(const map<int, Node*>& tree, vector<int>& outvec, int& nleaves);
		void find_leaf_nup(const map<int, Node*>& tree, vector<int>& outvec, int& nleaves); 
		void find_int_nodes(const map<int, Node*>& tree, vector<int>& outvec, int& intnodes);
		void find_subtree(const map<int, Node*>& tree, const int& node, 
			map<int, Node*>& subtree);
		void find_subtree_nodes(const map<int, Node*>& tree, const int& node, vector<int>& subnodes);
		void find_prior(const map<int, Node*>& tree); 
		void find_post(const map<int, Node*>& tree); 
		void get_subtree_info(const map<int, Node*>& subtree, const int& node,
			Ref<VectorXd> svec, Ref<MatrixXd> svar, double& npost, double& nleaf);
		void del_subtree(map<int, Node*>& tree, const int& node, map<int, Node*>& subtree);
		void copyin_subtree(map<int, Node*>& tree, const map<int, Node*>& subtree);
		bool is_node(const map<int, Node*>& tree, const int& node);
		vector<int> parent_pairs(const vector<int>& leaves); 
		vector<int> parent_child(const vector<int>& int_nodes); 
};

void Tree:: init_tree(map<int, Node*>& tree, const int& DATALENGTH, const vector<string>& names,
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats,
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) {
	bool init = true;
	//cout <<  Data.p << endl;
	int randint = Di.runif(0, DATALENGTH-1);
	//for(auto t : Data.allcats){cout << t.second << " ";} 
	//cout << endl;
	string pred = names[randint];
	tree[1]->pred = pred;
	pair<string, int> predpair = make_pair(pred, 1);
	tree[1]->pathpreds.push_back(predpair);
	map<string, Dataline::type>::const_iterator tit;
	tit = alltypes.find(pred);
	//cout << tit->first << " " << tit->second << endl;
	while(init) {
		if (tit->second == Dataline::STRING_T) {
			map<string, string>::const_iterator sit;
			sit = allcats.find(pred);
			tree[1]->ntype = tit->second;
			tree[1]->cats = sit->second;
			tree[1]->ncats = sit->second.length();
			char nsplit = Di.rchar(sit->second);
			tree[1]->ssplit = nsplit; 
			pair<string, string> lrcats = 
				Di.part_cats(nsplit, sit->second);
			tree[1]->lcats = lrcats.first;
			tree[1]->rcats = lrcats.second;
			tree[1]->nlcats = lrcats.first.length();
			tree[1]->nrcats = lrcats.second.length();
			tree[1]->prule = (1/(double)(tree[1]->ncats))*(1/(double)DATALENGTH);
			init = false;
		}
		else if (tit->second == Dataline::DOUBLE_T) {
			map<string, vector<double>>::const_iterator dit;
			dit = dbounds.find(pred);
			tree[1]->ntype = tit->second;
			tree[1]->min = dit->second[0];
			tree[1]->max = dit->second[1];
			tree[1]->bwidth = abs(tree[1]->max - tree[1]->min);
			tree[1]->dsplit = Di.cunif(tree[1]->min, tree[1]->max);
			tree[1]->prule = ((tree[1]->dsplit - tree[1]->min)/(double)(tree[1]->max - tree[1]->min))*
				(1/(double)DATALENGTH);
			init = false;
		}
		else {
			map<string, vector<int>>::const_iterator iit;
			iit = ibounds.find(pred);
			tree[1]->ntype = tit->second;
			tree[1]->min = iit->second[0];
			tree[1]->max = iit->second[1];
			tree[1]->dsplit = Di.runif(tree[1]->min, tree[1]->max);
			if(tree[1]->min == tree[1]->max) {
				tree[1]->prule = 1/(double)DATALENGTH;
			}
	 		else {		
				tree[1]->prule = (1/(double)(tree[1]->max - tree[1]->min))*(1/(double)DATALENGTH);
			}
			cout << tree[1]->prule << endl;
			init = false;
		}
	}
}

void Tree:: find_uleaf(const map<int, Node*>& tree, Dataline::line fline,
		int& found) {
	map<int, Node*>::const_iterator it;
	Dataline::line_it lit;

	for (it = tree.begin(); it != tree.end();) {
		string lpred = it->second->pred;
		lit = fline.find(lpred);

		if(it->second->isleaf) {
			found = it->second->nnode;
			it = tree.end();
			continue;
		}
		if(it->second->ntype == Dataline::STRING_T) {
			string nsplit = Data.get_string(lpred, fline);
			if (it->second->lcats.find(nsplit)!=std::string::npos) {
				it = tree.find(2*(it->second->nnode));
			}
			else {
				it = tree.find(2*(it->second->nnode)+1);
			}
		}
		else if(it->second->ntype == Dataline::DOUBLE_T) {
			double tsplit = it->second->dsplit;
			double nsplit = Data.get_double(lpred, fline);
			if (nsplit < tsplit) {
				it = tree.find(2*(it->second->nnode));
			}
			else {
				it = tree.find(2*(it->second->nnode)+1);
			}
		}
		else {
			int tsplit = it->second->dsplit;
			int nsplit = Data.get_int(lpred, fline);
			if (nsplit < tsplit) {
				it = tree.find(2*(it->second->nnode));
			}
			else {
				it = tree.find(2*(it->second->nnode)+1);
			}
		}
	}
	cout << "Inside uleaf found: " << found << endl;
}
void Tree:: find_nodes(const map<int, Node*>& tree, vector<int>& outvec) {
	map<int, Node*>::const_iterator it;

	for (it = tree.begin(); it != tree.end(); ++it) {
		outvec.push_back(it->second->nnode);
	}
}

bool Tree:: is_node(const map<int, Node*>& tree, const int& node) {
	if(tree.find(node) != tree.end()) {return(true);}
	else {return(false);}
}

void Tree:: find_leaf_nodes(const map<int, Node*>& tree, vector<int>& outvec, int& nleaves) {
	map<int, Node*>::const_iterator it;

	for (it = tree.begin(); it != tree.end(); ++it) {
		if (it->second->isleaf) {
			outvec.push_back(it->second->nnode);
			++nleaves;
		}
	}
}

void Tree:: find_leaf_nup(const map<int, Node*>& tree, vector<int>& outvec, int& nleaves) {
	map<int, Node*>::const_iterator it;

	for (it = tree.begin(); it != tree.end(); ++it) {
		if (it->second->isleaf) {
			if(it->second->isup == false) {
				outvec.push_back(it->second->nnode);
				++nleaves;
			}
		}
	}
}

void Tree:: find_int_nodes(const map<int, Node*>& tree, vector<int>& outvec, int& intnodes) {
	map<int, Node*>::const_iterator it;

	for (it = tree.begin(); it != tree.end(); ++it) {
		if (!it->second->isleaf) {
			outvec.push_back(it->second->nnode);
			++intnodes;
		}
	}
}

void Tree:: find_subtree(const map<int, Node*>& tree, const int& node, map<int, Node*>& subtree) {
	map<int, Node*>::const_iterator it;
	it = tree.find(node);

	if (it != tree.end()) {
		subtree[it->first] = new Node;
		*subtree[it->first] = *it->second;
		find_subtree(tree, 2*node, subtree);
		find_subtree(tree, 2*node+1, subtree);
	}
}

void Tree:: copyin_subtree(map<int, Node*>& tree, const map<int, Node*>& subtree) {
	map<int, Node*>::const_iterator it;
	it = subtree.begin();

	while(it != subtree.end()) {
		tree[it->first] = new Node;
		*tree[it->first] = *it->second;
		++it;
	}
}

void Tree:: find_subtree_nodes(const map<int, Node*>& tree, const int& node, vector<int>& subnodes) {
	map<int, Node*>::const_iterator it;
	it = tree.find(node);

	if (it != tree.end()) {
		subnodes.push_back(node);
		find_subtree_nodes(tree, 2*node, subnodes);
		find_subtree_nodes(tree, 2*node+1, subnodes);
	}
}

void Tree:: find_prior(const map<int, Node*>& tree) {
	map<int, Node*>::const_iterator it;
	size_t nleaves = tree.size();
	if (nleaves == 1) {
		it = tree.begin();
		prior = prior*(it->second->psplit)*(it->second->prule);
	}
	else {
		for (it = tree.begin(); it != tree.end(); ++it) {
			if (!it->second->isleaf) {
				prior = prior*(it->second->psplit)*(it->second->prule);
			}
			else {
				prior = prior*(1-it->second->psplit);
			}
		}
	}
}

void Tree:: find_post(const map<int, Node*>& tree) {
	map<int, Node*>::const_iterator it;

	for (it = tree.begin(); it != tree.end(); ++it) {
		if (it->second->isleaf) {
			logtreepost += it->second->lpost;
		}
	}
}

void Tree:: get_subtree_info(const map<int, Node*>& subtree, const int& node,
	Ref<VectorXd> svec, Ref<MatrixXd> svar, double& npost, double& nleaf) {

	map<int, Node*>::const_iterator it;

	for (it = subtree.begin(); it != subtree.end(); ++it) {
		if (it->second->isleaf == true) {
			svec = (svec + it->second->state)/(nleaf);
			svar = (svar + it->second->svar)/(nleaf);
			npost = npost + it->second->lpost;
			nleaf += 1;
		}
		else {
			svec = svec;
			svar = svar;
			npost = npost;
		}
	}
}

void Tree:: del_subtree(map<int, Node*>& tree, const int& node, map<int, Node*>& subtree) {
	map<int, Node*>::iterator it;
	it = subtree.find(node);
	delete it->second;
	subtree.erase(it);
	
	for (it = subtree.begin(); it != subtree.end(); ++it) {
		map<int, Node*>::const_iterator it1;
		it1 = tree.find(it->first);
		delete it1->second;
		tree.erase(it1);
	}
}

vector<int> Tree:: parent_pairs(const vector<int>& leaves) {
	set<int>lparents;
	set<int>rparents;
	vector<int> parent_pairs;
	if(leaves.size() == 1) {
		parent_pairs.push_back(1);
	}
	else {	
		for(int i = 0; i<leaves.size(); ++i) {
			if(leaves[i]%2 == 0) {
				lparents.insert(leaves[i]/int(2));
			}
			else {
				rparents.insert(leaves[i]/int(2));
			}
		}
		set_intersection(lparents.begin(), lparents.end(),
				rparents.begin(), rparents.end(),
				back_inserter(parent_pairs));
	}
	return parent_pairs;
}

vector<int> Tree:: parent_child(const vector<int>& int_nodes) {
	vector<int> parent_child;
	parent_child.push_back(1);
	if(int_nodes.size() > 2) {
		for(int i = 1; i < int_nodes.size(); ++i) {
			int cnode = int_nodes[i];
			if(find(int_nodes.begin(), int_nodes.end(), 2*cnode) != int_nodes.end()) {
				parent_child.push_back(cnode);
			}
			if(find(int_nodes.begin(), int_nodes.end(), 2*cnode+1) != int_nodes.end()) {
				parent_child.push_back(cnode);
			}
		}
	}
	auto last = unique(parent_child.begin(), parent_child.end());
	parent_child.erase(last, parent_child.end());
	return parent_child;
}

void Tree:: display(const map<int, Node*>& tree) {
	map<int, Node*>::const_iterator it;

	for (it = tree.begin(); it != tree.end(); ++it) {
		cout << 
		it->first << " " <<
		it->second->nnode << " " <<
		it->second->pred << " " <<
		//it->second->min << " " <<
		//it->second->max << " " <<
		//it->second->lcats << " " <<
		//it->second->rcats << " " <<
		it->second->dsplit << " " <<
		it->second->ssplit << " "; //<< endl;
		//for(auto s : it->second->pathpreds) {
		//	cout << s.first << " " << s.second << " ";
		//}
		//cout << endl;
		//it->second->psplit << " " <<
		//it->second->prule << " " <<
		//it->second->depth << " " <<
		//it->second->ppost << " " <<
		//it->second->lpost << " ";
		//if (it->second->isleaf) {cout << 1 << endl;}
		//else {cout << 0 << endl;}
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
	T.display(T.tree);
}

	for (int i = 2; i < 4; ++i) {
		string newline;
		Dataline::dline = Data.read_line(getline(tin, newline));
		int randint = Di.runif(0, Data.p);
		string pred = Data.names(randint);
		val = Data.dline.find(pred);
		T.tree[i] = new Node(i, 1, Data.ncats, 0, 3);
		T.tree[i]->pred = val->first;
		if(val->second.first == Data.LDOUBLE_T) {
			T.tree[i]->dsplit = val->second.first;
		}
		else if (val->second.first == Data.STRING_T)
		{
			T.tree[i]->ssplit = val->second.first;
		}
		T.tree[i]->lpost = 1;
		T.tree[i]->depth = T.tree[i]->nnode/2;
		T.tree[i]->psplit = psplitf(0.95, 1, T.tree[i]->depth);
		T.tree[i]->prule = prule(5, 5);
		if ((T.tree[i]->nnode/2 + 1) == i) {
			T.tree[i]->isleaf = false;
		}
	}

	T.display(T.tree);i
}

	VectorXd data(5);
	data << 1,3,5,6,3;
	int found = T.find_leaf(T.tree, data, mynames);
	cout << "found is: " << found << endl;

	VectorXd data1(5);
	data1 << 3,6,5,6,3;
	found = T.find_leaf(T.tree, data1, mynames);
	cout << "found is: " << found << endl;

	T.find_prior(T.tree);
	cout << "prior is: " << T.prior << endl;

	T.find_post(T.tree);
	cout << "post is: " << T.tree_post << endl;
	
	map<int, Node*> subtree;
	T.find_subtree(T.tree, 2, subtree);
	cout << subtree[2]->nnode << endl;
	T.display(subtree);

	vector<int> leaf_nodes;
	T.find_leaf_nodes(T.tree, leaf_nodes);
	int nleafs = leaf_nodes.size(); 
	cout << "nleafs: " << nleafs << endl;

map<string, double> Tree:: get_pred_path(map<int, Node*> tree, int node) {
	map<string, double> predpath;
	map<int, Node*>::const_iterator it;
	map<int, Node*>::const_iterator nit = tree.find(node);
	string nodepred =  nit->second->pred;
	double nodesplit = nit->second->split;

	for (it = tree.begin(); it != tree.end();) {
		predpath.insert(pair<string, double>(it->second->pred, 
			it->second->split));
		if (nodesplit < it->second->split) { 
			it = tree.find(2*(it->second->nnode));
		}
		else {
			it = tree.find(2*(it->second->nnode)+1);}
	}
	return(predpath);
}
}*/

