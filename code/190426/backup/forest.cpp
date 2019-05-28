// FOREST CLASS

# pragma once
# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <algorithm>
# include <vector>
# include <time.h>
# include "read.h"
# include "priorf.h"
# include "kalman.h"
# include "posterior.h"
# include "node.h"
# include "tree.h"
# include "move.h"
# include "proposal.h"
# include "dist.h"

using namespace std;
using namespace Eigen;

class Forest {
	public:
	map<int, Tree*> forest;
	Dist Di;
	Tree T;
	Lpost L;
	Prop P;
	Dataline D;
	
	//Constructor
	Forest() {};
	
	// Declarations
	void init_tree(map<int, Node*>& tree, Ref<MatrixXd> ATMAT, 
		Ref<MatrixXd> AIMAT, const double& ALPHA, const double& BETA, 
		Ref<VectorXd> CVEC, Ref<VectorXd> DVEC, const double& ppost,
		const double& lpost, const VectorXd& MU, const MatrixXd& SIG, 
		const int& DATALENGTH, const vector<string>& names,
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats,
		const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
	void kupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
		const int& NDIM, const int& MDIM, const VectorXd& yi, const MatrixXd& WMAT, 
		const MatrixXd& VMAT, const MatrixXd& FMAT, const MatrixXd& HMAT,
		const MatrixXd& GMAT, const VectorXd& input); 
	void lupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
		const VectorXd& yi, const MatrixXd& WI, const MatrixXd& VI, const MatrixXd& HVIH, 
		const MatrixXd& VIH, const MatrixXd& WIFT, const MatrixXd& WIF, const MatrixXd& FWIF, 
		const MatrixXd& UPGWI, const MatrixXd& GUWI, const double& Vdet, const double& Wdet, 
		double& sumlpostnup); 
	void proposal(Tree * tptr, map<int, Node*>& tree, const vector<int>& leaf_nodes, 
		const vector<int>& int_nodes, const double& ALPH, const double& BET, const int& DIM,
		const int& npred, const vector<string>& names, const map<string, Dataline::type>& alltypes,
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
		const map<string, vector<int>> ibounds); 
	void ch_tree(Tree * ltree, const double& ltreeprior, const double& ltreepost,
		const double& copyprior, const double& copypost, int& newtree);
	void nodeout(map<int, Node*> tree, ofstream& outfile, const int& k, const int& iter);
	void treeout(Tree * tptr, ofstream& outfile, const int& iter);
};
	
void Forest::init_tree(map<int, Node*>& tree, Ref<MatrixXd> ATMAT, 
	Ref<MatrixXd> AIMAT, const double& ALPHA, const double& BETA, 
	Ref<VectorXd> CVEC, Ref<VectorXd> DVEC, const double& ppost,
	const double& lpost, const VectorXd& MU, const MatrixXd& SIG, 
	const int& DATALENGTH, const vector<string>& names,
	const map<string, Dataline::type>& alltypes, const map<string, string>& allcats,
	const map<string, vector<double>>& dbounds, const map<string, vector<int>>& ibounds) { 

	T.init_tree(tree, DATALENGTH, names, alltypes, allcats, dbounds, ibounds);
	map<int, Node*>::const_iterator it;
	it = tree.find(1);
	it->second->psplit = psplitf(ALPHA, BETA, it->second->depth);
	it->second->ppost = ppost;
	it->second->lpost = lpost;

	it->second->state = MU;
	it->second->svar = SIG;
	it->second->atmat = ATMAT;
	it->second->aimat = AIMAT;
	it->second->cvec = CVEC;
	it->second->dvec = DVEC;
}

void Forest::kupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
	const int& NDIM, const int& MDIM, const VectorXd& yi, const MatrixXd& WMAT, 
	const MatrixXd& VMAT, const MatrixXd& FMAT, const MatrixXd& HMAT,
	const MatrixXd& GMAT, const VectorXd& input) {

	Kalman K(NDIM, MDIM);

	vector<int>::iterator vit;
	map<int, Node*>::const_iterator mit;
	vector<int> tnodes = nodes;

	if (tree.size() > 1) {
		tnodes.erase(remove(tnodes.begin(), tnodes.end(), node), tnodes.end());

		for (vit = tnodes.begin(); vit != tnodes.end(); ++vit) {
			mit = tree.find(*vit);
			K.nzhat(mit->second->state, FMAT);
			K.nsighat(mit->second->svar, WMAT);
		}
	}
	mit = tree.find(node);
	// Predict
	K.zvec(K.z_vec, FMAT, GMAT, mit->second->state, input);
	K.rmat(K.r_mat, FMAT, WMAT, mit->second->svar);
	// Update
	K.yvec(K.y_vec, yi, HMAT, K.z_vec);
	K.smat(K.s_mat, HMAT, K.r_mat, VMAT);
	K.kalg(K.kal_g, K.r_mat, FMAT, K.s_mat);
	K.zhat(mit->second->state, K.z_vec, K.kal_g, K.y_vec);
	K.sighat(mit->second->svar, K.r_mat, HMAT, K.kal_g);
}

void Forest::lupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
	const VectorXd& yi, const MatrixXd& WI, const MatrixXd& VI, const MatrixXd& HVIH, 
	const MatrixXd& VIH, const MatrixXd& WIFT, const MatrixXd& WIF, const MatrixXd& FWIF, 
	const MatrixXd& UPGWI, const MatrixXd& GUWI, const double& Vdet, const double& Wdet, 
	double& sumlpostnup) {
	
	vector<int> tnodes = nodes;
	vector<int>::iterator vit;
	map<int, Node*>::const_iterator mit;

	if (tree.size() > 1) {
		tnodes.erase(remove(tnodes.begin(), tnodes.end(), node),
						tnodes.end());

		for (vit = tnodes.begin(); vit != tnodes.end(); ++vit) {
			mit = tree.find(*vit);
			L.lpostnup(mit->second->atmat, mit->second->aimat, 
				mit->second->cvec, mit->second->dvec, 
				mit->second->ppost, mit->second->lpost, FWIF, WI,
				WIF, WIFT, GUWI, UPGWI, Wdet);
			sumlpostnup = sumlpostnup + mit->second->lpost;
		}
	}
	mit = tree.find(node);
	L.lpostup(mit->second->atmat, mit->second->aimat, mit->second->cvec, 
		mit->second->dvec, mit->second->ppost, mit->second->lpost, 
		yi, VI, FWIF, WI, WIF, WIFT, HVIH, VIH, UPGWI, GUWI, Vdet, Wdet);
}

void Forest::proposal(Tree * tptr, map<int, Node*>& tree, const vector<int>& leaf_nodes, 
	const vector<int>& int_nodes, const double& ALPH, const double& BET, const int& DIM,
	const int& npred, const vector<string>& names, const map<string, Dataline::type>& alltypes,
	const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
	const map<string, vector<int>> ibounds) {

	int numleaves = leaf_nodes.size();

	int move;
	if (numleaves == 1) {
		move = Di.runif(1, 2);
	}
	else {
		move = Di.runif(1, 4);
	}
	tptr->move = move; 
	switch (move) {
		case 1: {
		P.gprop(tree, leaf_nodes, ALPH, BET, DIM, npred, tptr->qratio, names,
				alltypes, allcats, dbounds, ibounds);
		break;
		}
		case 2: {
		P.cprop(tree, int_nodes, tptr->qratio, names, alltypes, allcats, dbounds,
				ibounds);
		break;
		}
		case 3: {
		P.sprop(tree, int_nodes, tptr->qratio);
		break;
		}
		case 4: {
		P.pprop(tree, leaf_nodes, DIM, tptr->qratio);
		break;
		}
	};
}

void Forest::ch_tree(Tree * ltree, const double& ltreeprior, const double& ltreepost,
	const double& copyprior, const double& copypost, int& newtree) {
	double qcur = log(ltreeprior)+ltreepost;
	double qnew = log(copyprior)+copypost;
	double prop_ratio = qnew - qcur + log(ltree->qratio);
	double u_init = log(Di.cunif(0, 1));
	if (u_init < prop_ratio) {
		newtree = 1;
	}
	else {
		newtree = 0;
	}
}

void Forest:: nodeout(map<int, Node*> tree, ofstream& outfile,
		const int& k, const int& iter)  {

	map<int, Node*>::iterator it;
	IOFormat linemat(StreamPrecision, DontAlignCols, ",", ",", "", "", "", "");

	for (it = tree.begin(); it != tree.end(); ++it) {
		outfile << iter << " " << k << " "
			<< it->second->nnode << " " 
			<< it->second->pred << " ";
			if (it->second->ssplit != "") { 
				outfile << it->second->ssplit << " ";
			}
			else {
				outfile << it->second->dsplit << " ";
			}
			outfile << it->second->depth << " "
			<< it->second->psplit << " "
			<< it->second->prule << " "
			<< it->second->state.format(linemat) << " "
			<< it->second->svar.format(linemat) << " "
			<< it->second->atmat.format(linemat) << " "
			<< it->second->aimat.format(linemat) << " "
			<< it->second->cvec.format(linemat) << " "
			<< it->second->dvec.format(linemat) << " "
			<< it->second->ppost << " "
			<< it->second->lpost << " ";
		if (it->second->isleaf == true) {
			outfile << 1 << "\n"; 
		}
		else { 
			outfile << 0 << "\n"; 
		}
	}
}

void Forest:: treeout(Tree * tptr, ofstream& outfile, const int& iter)  {
	outfile << iter << " " 
		<< tptr->tname << " " 
		<< tptr->threadid << " " 
		<< tptr->move << " " 
		<< tptr->accept << " " 
		<< tptr->numleaves << " " 
		<< tptr->numendleaves << " " 
		<< tptr->found << " " 
		<< tptr->cfound << " " 
		<< tptr->prior << " " 
		<< tptr->qratio << " " 
		<< tptr->logtreepost << "\n";
}	
