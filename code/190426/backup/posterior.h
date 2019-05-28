// FUNCTIONS FOR THE POSTERIOR CALCULATIONS 

# pragma once
# include <iostream>
# include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Lpost {
	public:

	// Constructors
	Lpost() {};

	//Declarations

	void AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
		const MatrixXd& WIF, const MatrixXd& WIFT, 
		const MatrixXd& WI, const MatrixXd& HVIH); 
	void AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
		const MatrixXd& WIF, const MatrixXd& WIFT, 
		const MatrixXd& WI);
	void AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,
		const MatrixXd& FWIF);
	void Ci(Ref<VectorXd> CVec, Ref<MatrixXd> AImatI, 
		const MatrixXd& UPGWI,
		Ref<VectorXd> DVec, const MatrixXd& WIFT);
	void Di(Ref<VectorXd> DVec, Ref<VectorXd> CVec, 
		const MatrixXd& GUWI,
		const MatrixXd& VIH, const VectorXd& Yi);
	void Dni(Ref<VectorXd> DVec, Ref<VectorXd> CVec, 
		const MatrixXd& GUWI);
	void lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
		Ref<VectorXd> CVec, Ref<VectorXd> DVec,
		double &fppost, double &flpost, 
		const VectorXd& yi, const MatrixXd& VTI,
		const MatrixXd& FWIF, const MatrixXd& WI, 
		const MatrixXd& WIF, const MatrixXd& WIFT, 
		const MatrixXd& HVIH, const MatrixXd& VIH,
		const MatrixXd& UPGWI, const MatrixXd& GUWI,
		const double Vdet, const double Wdet);
	void lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
		Ref<VectorXd> CVec, Ref<VectorXd> DVec,
		double &fppost, double &flpost, 
		const MatrixXd& FWIF, const MatrixXd& WI, 
		const MatrixXd& WIF, const MatrixXd& WIFT, 
		const MatrixXd& GUWI, const MatrixXd& UPGWI,
		const double Wdet); 
};

void Lpost:: AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
	const MatrixXd& WIF, const MatrixXd& WIFT, 
	const MatrixXd& WI, const MatrixXd& HVIH) {
	ATmat = HVIH + WI - (WIFT*AImatI*WIF);
}

void Lpost:: AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
	const MatrixXd& WIF, const MatrixXd& WIFT, 
	const MatrixXd& WI) {
	ATmat = WI - (WIFT*AImatI*WIF);
}

void Lpost:: AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,
	const MatrixXd& FWIF) {
//	cout << ATmat << endl;
//	cout << AImat << endl;
//	cout << FWIF << endl;
	AImat = ATmat+FWIF;
}

void Lpost:: Ci(Ref<VectorXd> CVec, Ref<MatrixXd> AImatI, 
	const MatrixXd& UPGWI,
	Ref<VectorXd> DVec, const MatrixXd& WIFT) {
	CVec = WIFT*AImatI*DVec + UPGWI;
}

void Lpost:: Di(Ref<VectorXd> DVec, Ref<VectorXd> CVec, 
	const MatrixXd& GUWI,
	const MatrixXd& VIH, const VectorXd& Yi) {
	DVec = CVec - GUWI + VIH.transpose()*Yi;
}

void Lpost:: Dni(Ref<VectorXd> DVec, Ref<VectorXd> CVec, 
	const MatrixXd& GUWI) {
	DVec = CVec - GUWI;
}

void Lpost:: lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
	Ref<VectorXd> CVec, Ref<VectorXd> DVec,
	double &fppost, double &flpost, 
	const VectorXd& yi, const MatrixXd& VTI,
	const MatrixXd& FWIF, const MatrixXd& WI, 
	const MatrixXd& WIF, const MatrixXd& WIFT, 
	const MatrixXd& HVIH, const MatrixXd& VIH,
	const MatrixXd& UPGWI, const MatrixXd& GUWI,
	const double Vdet, const double Wdet) {

	double adetp;
	double adetc;
	double adsqp;
	double adsqc;
	double ysq;
	MatrixXd AImatI;
	
	// Update the previous ATmat into the latest AImat: AImat = ATmat + FIWF
	AI_up(ATmat, AImat, FWIF);

	// Get the inverse and previous-updated inputs 
	AImatI = AImat.inverse();
	adsqp = DVec.transpose()*AImatI*DVec;
	adetp = log(AImat.determinant());
	// Update ppost
	fppost = fppost+0.5*(adetp+adsqp);  

	// Update ATmat using AImatI
	AT_up(ATmat, AImatI, WIF, WIFT, WI, HVIH); // Update AT with prev AI

	// Update Cvec and Dvec using AImatI and pDvec
	Ci(CVec, AImatI, UPGWI, DVec, WIFT); // Update current CVec
	Di(DVec, CVec, GUWI, VIH, yi); // Update current DVec

	// Calculate current values for lpost
	ysq = yi.transpose()*VTI*yi;
	adsqc = DVec.transpose()*ATmat.inverse()*DVec;
	adetc = log(ATmat.determinant());

	//flpost = get(lpost);
    	flpost = fppost + (double)0.5*(adetc-log(Wdet)- log(Vdet))-(double)0.5*(ysq - adsqc);
}

void Lpost:: lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
	Ref<VectorXd> CVec, Ref<VectorXd> DVec,
	double &fppost, double &flpost, 
	const MatrixXd& FWIF, const MatrixXd& WI, 
	const MatrixXd& WIF, const MatrixXd& WIFT, 
	const MatrixXd& GUWI, const MatrixXd& UPGWI,
	const double Wdet) {

	double adetp;
	double adetc;
	double adsqp;
	double adsqc;
	MatrixXd AImatI;

	// Update the previous ATmat into the latest AImat: AImat = ATmat + FIWF
	AI_up(ATmat, AImat, FWIF);
	
	// Using the pAImat and pDvec, create new values for ppost
	AImatI = AImat.inverse();
	adsqp = (DVec.transpose()*AImatI*DVec);
	adetp = log(AImat.determinant());

	// Update ppost
	fppost = fppost+0.5*(adetp+adsqp);  

	// Update ATmat using previous AImatI
	AT_nup(ATmat, AImatI, WIF, WIFT, WI); // Update AT with prev AI

	// Update Cvec and Dvec using AImat and pDvec
	Ci(CVec, AImatI, UPGWI, DVec, WIFT); // Update current CVec
	Dni(DVec, CVec, GUWI); // Update current DVec

	// Calculate current values for lpost
	adsqc = (DVec.transpose()*ATmat.inverse()*DVec);
	adetc = log(ATmat.determinant());

    	flpost = fppost + (double)0.5*(adetc-log(Wdet))+(double)0.5*(adsqc);
}
