// FOREST2 CLASS

# pragma once
# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <algorithm>
# include <vector>
# include <time.h>
# include "read.h"
# include "prior.h"
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
	void init_tree(map<int, Node*>& tree, const double& ALPHA, const double& BETA, 
		const VectorXd& MU, const MatrixXd& SIG, const int& DATALENGTH, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, const map<string, 
		vector<double>>& dbounds, const map<string, vector<int>>& ibounds); 
	void lupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
		const VectorXd& yi, const MatrixXd& WI, const MatrixXd& VI, const MatrixXd& HVIH, 
		const MatrixXd& HTVI, const MatrixXd& FTWI, const MatrixXd& WIF, const MatrixXd& FTWIF, 
		const VectorXd& GU, const double& Vdet, const double& Wdet, double& sumlpostnup, const int UPITER);
	void kupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
		const int& NDIM, const int& MDIM, const VectorXd& yi, const MatrixXd& WMAT, 
		const MatrixXd& VMAT, const MatrixXd& FMAT, const MatrixXd& HMAT,
		const MatrixXd& GMAT, const VectorXd& input); 
	void proposal(Tree * tptr, map<int, Node*>& tree, const vector<int>& leaf_nodes, 
		const vector<int>& int_nodes, const double& ALPH, const double& BET, const int& MDIM, const int& NDIM,
		const int& npred, const vector<string>& names, const map<string, Dataline::type>& alltypes,
		const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
		const map<string, vector<int>> ibounds, const int numlevels, const int TREESIZE, size_t treesize, 
		const double pgmove, const double ppmove, const double pnotup, const double pfmove, const double prmove);
	void ch_tree(Tree * ltree, const double& ltreeprior, const double& ltreepost,
		const double& copyprior, const double& copypost, int& newtree, const double ltreezk,
		const double copyzk); 
	void ch_temp(Tree * ltree, const double& ltreeprior, const double& ltreepost,
		const vector<double>& temp, const vector<double>& fkprior);
	void nodeout(map<int, Node*> tree, ofstream& outfile, const int& k, const int& iter);
	void treeout(Tree * tptr, ofstream& outfile, const int& iter);
	double logzk(const map<int, Node*> tree);
};
	
void Forest::init_tree(map<int, Node*>& tree, const double& ALPHA, const double& BETA,
		const VectorXd& MU, const MatrixXd& SIG, const int& DATALENGTH, const vector<string>& names, 
		const map<string, Dataline::type>& alltypes, const map<string, string>& allcats, const map<string, 
		vector<double>>& dbounds, const map<string, vector<int>>& ibounds) { 

	T.init_tree(tree, DATALENGTH, names, alltypes, allcats, dbounds, ibounds);
	map<int, Node*>::const_iterator it;
	it = tree.find(1);
	it->second->psplit = psplitf(ALPHA, BETA, it->second->depth);
	it->second->state = MU;
	it->second->zpred = MU;
	it->second->ypred = VectorXd::Zero(it->second->ndim);
	it->second->svar = SIG;
	it->second->rmat = SIG;
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
			mit->second->zpred = mit->second->state;
			mit->second->ypred = HMAT*mit->second->zpred;
			mit->second->rmat = mit->second->svar;
		}
	}
	mit = tree.find(node);
	// Predict
	K.zvec(mit->second->zpred, FMAT, GMAT, mit->second->state, input);
	K.rmat(mit->second->rmat, FMAT, WMAT, mit->second->svar);
	// Update
	K.yvec(mit->second->ypred, yi, HMAT, mit->second->zpred);
	K.smat(K.s_mat, HMAT, mit->second->rmat, VMAT);
	K.kalg(K.kal_g, mit->second->rmat, HMAT, K.s_mat);
	K.zhat(mit->second->state, mit->second->zpred, K.kal_g, mit->second->ypred);
	K.sighat(mit->second->svar, mit->second->rmat, HMAT, K.kal_g);
}

void Forest::lupdate(const map<int, Node*>& tree, const int& node, const vector<int>& nodes, 
	const VectorXd& yi, const MatrixXd& WI, const MatrixXd& VI, const MatrixXd& HVIH, 
	const MatrixXd& HTVI, const MatrixXd& FTWI, const MatrixXd& WIF, const MatrixXd& FTWIF, 
	const VectorXd& GU, const double& Vdet, const double& Wdet, double& sumlpostnup,
	const int UPITER) {
	
	vector<int> tnodes = nodes;
	vector<int>::iterator vit;
	map<int, Node*>::const_iterator mit;

	if (tree.size() > 1) {
		tnodes.erase(remove(tnodes.begin(), tnodes.end(), node),
	tnodes.end());

		cout << "Inside lupdate not found: ";	
		for (vit = tnodes.begin(); vit != tnodes.end(); ++vit) {
			cout << *vit << " ";
			mit = tree.find(*vit);
			L.lpostnup(mit->second->atmat, mit->second->aimat, 
				mit->second->dvec, mit->second->lpost, mit->second->adetp,
				mit->second->adsqp, mit->second->HTVYp, FTWIF, WI,
				WIF, FTWI, GU, Wdet);
			sumlpostnup += mit->second->lpost;
			++mit->second->nup;
			if(mit->second->nup > UPITER) {
				mit->second->isup = false;
			}
		}
	}
	cout << endl;	
	mit = tree.find(node);
	cout << "Inside lupdate found: " << node << endl;	
	L.lpostup(mit->second->atmat, mit->second->aimat, mit->second->dvec, mit->second->lpost,
		       mit->second->adetp, mit->second->adsqp, mit->second->HTVYp, yi, VI, FTWIF, 
		       WI, WIF, FTWI, HVIH, HTVI, GU, Vdet, Wdet);
}

void Forest::proposal(Tree * tptr, map<int, Node*>& tree, const vector<int>& leaf_nodes, 
	const vector<int>& int_nodes, const double& ALPH, const double& BET, const int& MDIM, const int& NDIM,
	const int& npred, const vector<string>& names, const map<string, Dataline::type>& alltypes,
	const map<string, string>& allcats, const map<string, vector<double>>& dbounds,
	const map<string, vector<int>> ibounds, const int numlevels,
	const int TREESIZE, size_t treesize, const double pgmove, const double ppmove, const double pnotup,
	const double pfmove, const double prmove) {

	size_t numleaves = leaf_nodes.size();
	int move;
	int ngrow;
	int clevel = tptr->level;
	int nlevel;
	int levelprob;

	// Propose new levels
	if(clevel == 1) {
		nlevel = 2;
		tptr->nlevel = nlevel;
	}
	else if (clevel == numlevels) {
		nlevel = numlevels-1;
		tptr->nlevel = nlevel;
	}
	else {
		int ch_level = Di.runif(0,1);
		if (ch_level == 1) {
			nlevel = clevel+1;
			tptr->nlevel = nlevel;
	}
		else {
			nlevel = clevel-1;
			tptr->nlevel = nlevel;
		}
	}

	if (treesize < TREESIZE) {
		P.grwprop(tree, leaf_nodes, ALPH, BET, MDIM, NDIM, names, alltypes, allcats, 
				dbounds, ibounds, tptr->lfratio, tptr->lrratio, pgmove, ppmove);
		tptr->move = 1; 
	}
	else {
		move = Di.runif(1, 8);
		switch (move) {
			case 1: {
				P.shftprop(tree, int_nodes, tptr->lfratio, tptr->lrratio);
				tptr->move = 8; 
				break;
			}
			case 2: {
				P.chprop(tree, int_nodes, names, alltypes, allcats, dbounds,
					ibounds, tptr->lfratio, tptr->lrratio);
				tptr->move = 3; 
				break;
			}
			case 3: {
				if(treesize > 2) {
					P.swprop(tree, int_nodes, tptr->lfratio, tptr->lrratio);
				}
				tptr->move = 4; 
				break;
			}
			case 4: {
				P.mgprop(tree, tptr->lfratio, tptr->lrratio, pfmove, prmove, ALPH, BET, pnotup);
				tptr->move = 5; 
				break;
			}
			case 5: {
				if(treesize > 1) {
					P.prngrwprop(tree, ALPH, BET, MDIM, NDIM, int_nodes, leaf_nodes, names, alltypes, allcats, 
						dbounds, ibounds, tptr->lfratio, tptr->lrratio);
				}
				tptr->move = 7; 
				break;
			}
			case 6: {
				P.grwprop(tree, leaf_nodes, ALPH, BET, MDIM, NDIM, names, alltypes, allcats, 
						dbounds, ibounds, tptr->lfratio, tptr->lrratio, pgmove, ppmove);
				tptr->move = 1; 
				break;
			}
			case 7: {
				if(treesize > 1) {
					P.prnprop(tree, int_nodes, leaf_nodes, MDIM, tptr->lfratio, tptr->lrratio,
							pgmove, ppmove);
				}
				tptr->move = 2; 
				break;
			}
			case 8: {
				if(treesize > 1) {
					P.tswprop(tree, tptr->lfratio, tptr->lrratio);
				}
				tptr->move = 6; 
				break;
			}
		};
	}
}

double Forest::logzk(const map<int, Node*> tree) {
	map<int, Node*>::const_iterator it;
	double lzk = 0.0;
	double dzk = 1.0;
	it = tree.begin();
	for(; it != tree.end(); ++it) {
		dzk = Di.dmvnorm(it->second->zpred, it->second->zpred, it->second->rmat);
		lzk = lzk + log(dzk);
	}
	return(lzk);
}

void Forest::ch_tree(Tree * ltree, const double& ltreeprior, const double& ltreepost,
	const double& copyprior, const double& copypost, int& newtree, const double ltreezk,
	const double copyzk) {
	double qcur = log(ltreeprior)+ltreezk+ltreepost;
	double qnew = log(copyprior)+copyzk+copypost;
	double prop_ratio = (qnew + ltree->lfratio) - (qcur + ltree->lrratio);
	double u_init = log(Di.cunif(0, 1));
	if (u_init < prop_ratio) {
		newtree = 1;
	}
	else {
		newtree = 0;
	}
}

void Forest::ch_temp(Tree * ltree, const double& ltreeprior, const double& ltreepost,
	const vector<double>& temp, const vector<double>& fkprior) {
	int lprobratio;
	int clevel = ltree->level;
	int nlevel = ltree->nlevel;
	double ctemp = temp[clevel-1];
	double ntemp = temp[nlevel-1];
	int maxlevel = temp.size();
	if(clevel == maxlevel || clevel == 1) {
		lprobratio = 2;
	}
	else {
		lprobratio = 1;
	}

	double qcur = (1/double(ctemp))*(log(ltreeprior)+ltreepost) + fkprior[clevel-1];
	double qnew = (1/double(ntemp))*(log(ltreeprior)+ltreepost) + fkprior[nlevel-1];
	double prop_ratio = (qnew - qcur) + log(lprobratio);
	double u_init = log(Di.cunif(0, 1));
	if (u_init < prop_ratio) {
		ltree->level = nlevel;
		ltree->fkprior = fkprior[nlevel-1];
	}
	else {
		ltree->level = clevel;
		ltree->fkprior = fkprior[clevel-1];
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
			<< it->second->ypred.format(linemat) << " "
			<< it->second->zpred.format(linemat) << " "
			<< it->second->svar.format(linemat) << " "
			<< it->second->atmat.format(linemat) << " "
			<< it->second->aimat.format(linemat) << " "
			<< it->second->dvec.format(linemat) << " "
			<< it->second->HTVYp.format(linemat) << " "
			<< it->second->adetp << " "
			<< it->second->adsqp << " "
			<< it->second->lpost << " ";
		if (it->second->isleaf == true) {
			outfile << 1 << " "; 
		}
		else { 
			outfile << 0 << " "; 
		}
		if (it->second->isright == true) {
			outfile << 1 << " "; 
		}
		else { 
			outfile << 0 << " "; 
		}
		if (it->second->isleft == true) {
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
		<< tptr->found << " " 
		<< tptr->level << " " 
		<< tptr->fkprior << " " 
		<< tptr->prior << " " 
		<< tptr->lfratio << " " 
		<< tptr->lrratio << " " 
		<< tptr->logtreepost << " "
		<< tptr->slpostn << " " 
		<< tptr->qstar << "\n";
}	
