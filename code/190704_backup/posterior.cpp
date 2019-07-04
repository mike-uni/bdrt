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
	void AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, const MatrixXd& WIF, const MatrixXd& FTWI, 
		const MatrixXd& WI, const MatrixXd& HVIH, const VectorXd& GU); 
	void AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
		const MatrixXd& WIF, const MatrixXd& FTWI, const MatrixXd& WI, const VectorXd& GU);
	void AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,const MatrixXd& FTWIF);
	void Di(Ref<VectorXd> DVec, const MatrixXd& AImatI, const MatrixXd& FTWI, const MatrixXd& WI, 
		const VectorXd& GU, const MatrixXd& WIF, Ref<VectorXd> HTVY);
	void lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, Ref<VectorXd> DVec, double& lpost, 
		double& adetp, double& adsqp, Ref<VectorXd> HTVYp, const VectorXd& yi, const MatrixXd& VTI, 
		const MatrixXd& FTWIF, const MatrixXd& WI, const MatrixXd& WIF, const MatrixXd& FTWI, 
		const MatrixXd& HVIH, const MatrixXd& HTVI, const VectorXd& GU, const double Vdet, const double Wdet);
	void lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, Ref<VectorXd> DVec, double& lpost, 
		double& adetp, double& adsqp, Ref<VectorXd> HTVYp, const MatrixXd& FTWIF, 
		const MatrixXd& WI, const MatrixXd& WIF, const MatrixXd& FTWI,
		const VectorXd& GU, const double Wdet); 
};

void Lpost:: AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
	const MatrixXd& WIF, const MatrixXd& FTWI, 
	const MatrixXd& WI, const MatrixXd& HVIH, const VectorXd& GU) {
	//cout << " Inside ATup: " << endl;
/*	cout << ATmat << endl;
	cout << AImatI << endl;
	cout << WIF << endl;
	cout << WIFT << endl;
	cout << WI << endl;
*/	ATmat = HVIH + WI - WIF*AImatI*FTWI;
	//cout << (WIF*AImatI*FTWI) << endl;
}

void Lpost:: AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, const MatrixXd& WIF, const MatrixXd& FTWI, 
	const MatrixXd& WI, const VectorXd& GU) {
	ATmat = WI - WIF*AImatI*FTWI;
}

void Lpost:: AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, const MatrixXd& FTWIF) {
/*	cout << " Inside AIup: " << endl;
	cout << ATmat << endl;
	cout << AImat << endl;
	cout << FWIF << endl;
*/	AImat = ATmat+FTWIF;
}

void Lpost:: Di(Ref<VectorXd> DVec, const MatrixXd& AImatI, const MatrixXd& FTWI, const MatrixXd& WI, const VectorXd& GU, 
		const MatrixXd& WIF, Ref<VectorXd> HTVY) {
	int mdim = WIF.rows();
	DVec = HTVY+AImatI*FTWI*DVec+WI*GU-WIF*AImatI*FTWI*GU;
}

void Lpost:: lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, Ref<VectorXd> DVec, double& lpost, 
		double& adetp, double& adsqp, Ref<VectorXd> HTVYp, const VectorXd& yi, const MatrixXd& VTI, 
		const MatrixXd& FTWIF, const MatrixXd& WI, const MatrixXd& WIF, const MatrixXd& FTWI, 
		const MatrixXd& HVIH, const MatrixXd& HTVI, const VectorXd& GU, 
		const double Vdet, const double Wdet) {
	//cout << "Inside lpostup: " << "TBC" << endl;
	double adetup;
	double adetc;
	double adsqup;
	double adsqc;
	double ysq;
	double daigu;
	double ugwgu;
	MatrixXd AImatI;
	Ref<VectorXd> Dup = DVec;

	// Subtract previous terms from previous posterior
	lpost = lpost - (adsqp + adetp);
	
	// Update the previous ATmat into the latest AImat: AImat = ATmat + FTWIF
	AI_up(ATmat, AImat, FTWIF);

	// Create updated form of previous terms and add them back to previous posterior
	AImatI = AImat.inverse();
	adetup = 0.5*log((2*M_PI*AImat).determinant());
	Di(Dup, AImatI, FTWI, WI, GU, WIF, HTVYp); // Update temporary Dup with updated AImatI
	adsqup = 0.5*Dup.transpose()*AImatI*Dup;

	lpost = lpost + adsqup + adetup;

	// Update ATmat using newly updated AImatI
	AT_up(ATmat, AImatI, WIF, FTWI, WI, HVIH, GU);
	//cout << ATmat << endl;
	//cout << ATmat.determinant() << endl;

	// Update Dvec using AImatI and new HTVY
	VectorXd HTVY = HTVI*yi;
	HTVYp = HTVY;
	Di(DVec, AImatI, FTWI, WI, GU, WIF, HTVY); // Update temporary Dvec

	// Calculate current values for lpost
	ysq = yi.transpose()*VTI*yi;
	daigu = 2*Dup.transpose()*AImatI*FTWI*GU;
	ugwgu = GU.transpose()*(WI-WIF*AImatI*FTWI)*GU;
	adsqc = 0.5*DVec.transpose()*ATmat.inverse()*DVec;
	adetc = 0.5*(log((2*M_PI*ATmat).determinant()));
	adsqp = adsqc;
	adetp = adetc;

	//flpost = get(lpost);
    	lpost = lpost-0.5*(log(Wdet)+log(Vdet)+daigu+ugwgu+ysq)+(adetc + adsqc);
	//cout << "Inside lpostup: " << endl;
	//cout << adetc << " " << log(Wdet) << " " << log(Vdet) << " " << ysq << " " << adsqc << endl;
}

void Lpost:: lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, Ref<VectorXd> DVec, double& lpost, 
		double& adetp, double& adsqp, Ref<VectorXd> HTVYp, const MatrixXd& FTWIF, 
		const MatrixXd& WI, const MatrixXd& WIF, const MatrixXd& FTWI,
		const VectorXd& GU, const double Wdet) {

	double adetup;
	double adetc;
	double adsqup;
	double adsqc;
	double daigu;
	double ugwgu;
	MatrixXd AImatI;
	Ref<VectorXd> Dup = DVec;

	// Subtract previous terms from previous posterior
	lpost = lpost - (adsqp + adetp);
	
	// Update the previous ATmat into the latest AImat: AImat = ATmat + FTWIF
	AI_up(ATmat, AImat, FTWIF);

	// Create updated form of previous terms and add them back to previous posterior
	AImatI = AImat.inverse();
	adetup = 0.5*log((2*M_PI*AImat).determinant());
	Di(Dup, AImatI, FTWI, WI, GU, WIF, HTVYp); // Update temporary Dup with updated AImatI
	adsqup = 0.5*Dup.transpose()*AImatI*Dup;

	lpost = lpost + adsqup + adetup;

	// Update ATmat using newly updated AImatI
	AT_nup(ATmat, AImatI, WIF, FTWI, WI, GU);
	//cout << ATmat << endl;
	//cout << ATmat.determinant() << endl;

	// Update Dvec using AImatI and new HTVY
	int mdim = HTVYp.rows();
	VectorXd HTVY = VectorXd::Zero(mdim);
	HTVYp.setZero();
	Di(DVec, AImatI, FTWI, WI, GU, WIF, HTVY); // Update DVec

	// Calculate current values for lpost
	daigu = 2*Dup.transpose()*AImatI*FTWI*GU;
	ugwgu = GU.transpose()*(WI-WIF*AImatI*FTWI)*GU;
	adsqc = 0.5*DVec.transpose()*ATmat.inverse()*DVec;
	adetc = 0.5*(log((2*M_PI*ATmat).determinant()));
	adsqp = adsqc;
	adetp = adetc;

	//flpost = get(lpost);
    	lpost = lpost-0.5*(log(Wdet)+daigu+ugwgu)+(adetc + adsqc);
	//cout << "Inside lpostup: " << endl;
	//cout << adetc << " " << log(Wdet) << " " << log(Vdet) << " " << ysq << " " << adsqc << endl;
}
