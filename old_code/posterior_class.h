// FUNCTIONS FOR THE POSTERIOR CALCULATIONS 

# pragma once
# include <iostream>
# include <Eigen/Dense>
# include "node_class.h"

using namespace std;
using namespace Eigen;

class Lpost : public Node {
	public:
	// Constructors
	Lpost() {};

double get(double value) {return value;};
void set(double &value, double item) {item = value;}; 

	//Declarations

void AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
					const MatrixXd& WIF, const MatrixXd& WIFT, 
					const MatrixXd& WI, const MatrixXd& HVIH);
void AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
					const MatrixXd& WIF, const MatrixXd& WIFT, 
					const MatrixXd& WI);
void AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,
					const MatrixXd& FWIF);
void Ci(Ref<VectorXd> CVec, Ref<MatrixXd> AImat, 
				Ref<VectorXd> DVec, const MatrixXd& WIFT);
void Di(Ref<VectorXd> DVec, Ref<VectorXd> CVec,
				const MatrixXd& VH, const VectorXd& Yi);
void Dni(Ref<VectorXd> DVec, Ref<VectorXd> CVec); 
void lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const VectorXd& yi, const MatrixXd& VTI,
						const MatrixXd& FWIF, const MatrixXd& WI, 
						const MatrixXd& WIF, const MatrixXd& WIFT, 
						const MatrixXd& HVIH, const MatrixXd& VIH, 
						const double Vdet, const double Wdet);
void lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const MatrixXd& FWIF, const MatrixXd& WI, 
						const MatrixXd& WIF, const MatrixXd& WIFT, 
						const double Wdet); 
};

void Lpost:: AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
					const MatrixXd& WIF, const MatrixXd& WIFT, 
					const MatrixXd& WI, const MatrixXd& HVIH) {
    ATmat = HVIH + WI - (0.25)*(WIFT*AImatI*WIF);
}

void Lpost:: AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImatI, 
					const MatrixXd& WIF, const MatrixXd& WIFT, 
					const MatrixXd& WI) {
    ATmat = WI - (0.25)*(WIFT*AImatI*WIF);
}

void Lpost:: AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,
					const MatrixXd& FWIF) {
	AImat = ATmat+FWIF;
}

void Lpost:: Ci(Ref<VectorXd> CVec, Ref<MatrixXd> AImatI, 
				Ref<VectorXd> DVec, const MatrixXd& WIFT) {
    CVec = 0.25*WIFT*AImatI*DVec;
}

void Lpost:: Di(Ref<VectorXd> DVec, Ref<VectorXd> CVec,
				const MatrixXd& VIH, const VectorXd& Yi) {
    DVec = CVec + VIH*Yi;
}

void Lpost:: Dni(Ref<VectorXd> DVec, Ref<VectorXd> CVec) {
    DVec = CVec;
}

void Lpost:: lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const VectorXd& yi, const MatrixXd& VTI,
						const MatrixXd& FWIF, const MatrixXd& WI, 
						const MatrixXd& WIF, const MatrixXd& WIFT, 
						const MatrixXd& HVIH, const MatrixXd& VIH, 
						const double Vdet, const double Wdet) {
	double adetp;
	double adetc;
	double adsqp;
	double adsqc;
	double ysq;
	MatrixXd AImatI;

	if (AImat.isZero()) { 
		AImat = ATmat; // The initial case
		adsqp = 0.0;
		adetp = 0.0;
		cout << "Inside lpostup AI_mat before is: " << AImat << endl;
		AImatI = AImat.inverse();
		cout << "Inside lpostup AI_mat after is: " << AImat << endl;
	}
	else {
		cout << "Inside lpostup AI_mat before is: " << AImat << endl;
		AI_up(ATmat, AImat, FWIF); // Update previous At to Ai
		cout << "Inside lpostup AI_mat after is: " << AImat << endl;
		AImatI = AImat.inverse();
		adsqp = (DVec.transpose()*AImatI*DVec);
		adetp = log(AImat.determinant());
	}

	ysq = (yi.transpose()*VTI*yi);

	fppost = fppost-0.5*(log(Wdet)+log(Vdet)+ysq)+0.5*(adetp+adsqp);  

	cout << "Inside lpostup AT_mat before is: " << ATmat << endl;
	AT_up(ATmat, AImatI, WIF, WIFT, WI, HVIH); // Update AT with prev AI
	cout << "Inside lpostup AT_mat after is: " << ATmat << endl;
	Ci(CVec, AImatI, DVec, WIFT); // Update current Ci vector

	Di(DVec, CVec, VIH, yi); // Update current Di vector

	adsqc = (DVec.transpose()*ATmat.inverse()*DVec);
	adetc = log(ATmat.determinant());

	cout << "Inside lpostup flpost before is: " << flpost << endl;
	cout << "fppost = " << fppost << " adetp = " << adetp << " adetc = " << adetc 
			<< " adsqp = " << adsqp << " adsqc = " << adsqc << endl;	

	//flpost = get(lpost);
    flpost = fppost + 0.5*(adetc + (0.25)*adsqc);
	cout << "Inside lpostup flpost after is: " << flpost << endl;
/*	cout << "log(2*M_PI*(adetp*adetc)) is: "<< log(2*M_PI*(adetp*adetc)) 
		<< endl;
	cout << "adsqp+adsqc is: " << adsqp+adsqc << endl;
	//set(flpost, lpost);

	cout << "fppost after is: " << fppost << endl;
	cout << "flpost after is: " << flpost << endl;
*/}

void Lpost:: lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const MatrixXd& FWIF, const MatrixXd& WI, 
						const MatrixXd& WIF, const MatrixXd& WIFT, 
						const double Wdet) {
	double adetp;
	double adetc;
	double adsqp;
	double adsqc;
	MatrixXd AImatI;

	if (AImat.isZero()) { 
		AImat = ATmat; // The initial case
		adsqp = 0.0;
		adetp = 0.0;
		cout << "Inside lpostup AI_mat before is: " << AImat << endl;
		AImatI = AImat.inverse();
		cout << "Inside lpostup AI_mat after is: " << AImat << endl;
	}
	else {
		cout << "Inside lpostup AI_mat before is: " << AImat << endl;
		AI_up(ATmat, AImat, FWIF); // Update previous At to Ai
		cout << "Inside lpostup AI_mat after is: " << AImat << endl;
		AImatI = AImat.inverse();
		adsqp = (DVec.transpose()*AImatI*DVec);
		adetp = log(AImat.determinant());
	}

	fppost = fppost-0.5*(log(Wdet))+0.5*(adetp+adsqp);  

	cout << "Inside lpostnup AT_mat before is: " << ATmat << endl;
	AT_nup(ATmat, AImatI, WIF, WIFT, WI); // Update At
	cout << "Inside lpostnup AT_mat after is: " << ATmat << endl;
	Ci(CVec, AImat, DVec, WIFT); // Update current Ci vector

	Dni(DVec, CVec); // Update current Di vector

	adsqc = (DVec.transpose()*ATmat.inverse()*DVec);
	adetc = log(ATmat.determinant());

	cout << "Inside lpostnup flpost before is: " << flpost << endl;
	cout << "fppost = " << fppost << " adetp = " << adetp << " adetc = " << adetc 
			<< " adsqp = " << adsqp << " adsqc = " << adsqc << endl;	

    flpost = fppost + 0.5*(adetc + (0.25)*adsqc);
	cout << "Inside lpostnup flpost after is: " << flpost << endl;
	//set(flpost, lpost);
}

/*
int main() {
	int const DIMS = 2;
	Node a(1, DIMS);
	Lpost l;
	
	VectorXd MU(DIMS);
	MU << 0,0;
	VectorXd y_i(DIMS);
	y_i << 1, 1;
	Matrix2d F_MAT;
	F_MAT << 1,0,0,1;
	Matrix2d H_MAT;
	H_MAT << 1,0,0,1;
	Matrix2d W_MAT;
	W_MAT << 1,0,0,1; 
	Matrix2d V_MAT;
	V_MAT << 1,0,0,1;
	Matrix2d WI;
	WI = W_MAT.inverse();
	Matrix2d VI;
	VI = V_MAT.inverse();
	Matrix2d WIFT;
	WIFT = WI*F_MAT.transpose();
	Matrix2d WIF;
	WIF = WI*F_MAT;
	Matrix2d FWIF;
	FWIF = F_MAT.transpose()*WI*F_MAT;
	Matrix2d VH;
	VH = VI*H_MAT;
	Matrix2d HVIH;
	HVIH = H_MAT.transpose()*VI*H_MAT;
	double Vid;
	Vid = VI.determinant();
	double WId;
	WId = WI.determinant();

	Matrix2d W0;
	W0 << 0.9,0,0,0.9;
	double lW0det = -log(2*M_PI*W0.determinant());
	double W0sq = MU.transpose()*W0.inverse()*MU;
	
	VectorXd D0(DIMS);
	D0 = MU.transpose()*W0;
	Matrix2d A0;
	A0 = W0.inverse()+FWIF;
	double adsq0 = D0.transpose()*A0.inverse()*D0;
	
	a.aimat.setZero();
	cout << "AImat init is:\n" << a.aimat << endl;
	a.atmat = A0;
	cout << "ATmat init is:\n" << a.atmat << endl;
	a.cvec = D0;
	cout << "Cvec init is:\n" << a.cvec << endl;
	a.dvec = D0;
	cout << "DVec init is:\n" << a.dvec << endl;
	a.ppost = lW0det + W0sq + adsq0;
	a.lpost = 1;
	cout << "ppostp init is:" << a.ppost << endl;

	for(int i = 0; i < 3; i++) { 
		l.lpostup(a.atmat, a.aimat, a.cvec, a.dvec, a.ppost,
					a.lpost, y_i, VI, FWIF, WI, WIF, WIFT, 
					HVIH, VH, Vid, WId);

		cout << "AImat " << i+1 << " is:\n" << a.aimat << endl;
		cout << "ATmat " << i+1 << " is:\n" << a.atmat << endl;
		cout << "Cvec " << i+1 << " is:\n" << a.cvec << endl;
		cout << "DVec " << i+1 << " is:\n" << a.dvec << endl;
		cout << "ppost " << i+1 << " is: " << a.ppost << endl;
		cout << "lpost " << i+1 << " is: " << a.lpost << endl;
		//cout << "a.leaf_post " << i+1 << " is: " << a.leaf_post << endl;
		cout << endl;


		l.lpostnup(a.atmat, a.atmat, a.cvec, a.dvec, a.ppost, 
					a.lpost, FWIF, WI, WIF, WIFT, WId);

		cout << "postpp 2 is:" << a.ppost << endl;
		cout << "a.lpost 2 is: " << a.lpost << endl;
	}
}*/
