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
					const MatrixXd& WF, const MatrixXd& WFT, 
					const MatrixXd& WMAT, const MatrixXd& HVH);
void AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
					const MatrixXd& WF, const MatrixXd& WFT, 
					const MatrixXd& WMAT);
void AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,
					const MatrixXd& FWF);
void Ci(Ref<VectorXd> CVec, Ref<MatrixXd> AImat, 
				Ref<VectorXd> DVec, const MatrixXd& WFT);
void Di(Ref<VectorXd> DVec, Ref<VectorXd> CVec,
				const MatrixXd& VH, const VectorXd& Yi);
void lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const VectorXd& yi, const MatrixXd& VTI,
						const MatrixXd& Fwf, const MatrixXd& Wi, 
						const MatrixXd& Wf, const MatrixXd& Wft, 
						const MatrixXd& Hvh, const MatrixXd& Vh, 
						const double Vdet, const double Wdet);
void lpostnup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const MatrixXd& Fwf, const MatrixXd& Wi, 
						const MatrixXd& Wf, const MatrixXd& Wft, 
						const double Wdet); 
};

void Lpost:: AT_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
					const MatrixXd& WF, const MatrixXd& WFT, 
					const MatrixXd& WMAT, const MatrixXd& HVH) {
    ATmat = -(WFT*AImat*WF) + HVH + WMAT;
}

void Lpost:: AT_nup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
					const MatrixXd& WF, const MatrixXd& WFT, 
					const MatrixXd& WMAT) {
    ATmat = -(WFT*AImat*WF) + WMAT;
}

void Lpost:: AI_up(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat,
					const MatrixXd& FWF) {
	AImat = ATmat+FWF;
}

void Lpost:: Ci(Ref<VectorXd> CVec, Ref<MatrixXd> AImat, 
				Ref<VectorXd> DVec, const MatrixXd& WFT) {
    CVec.noalias() = WFT*AImat*DVec;
}

void Lpost:: Di(Ref<VectorXd> DVec, Ref<VectorXd> CVec,
				const MatrixXd& VH, const VectorXd& Yi) {
    DVec = CVec + VH*Yi;
}

void Lpost:: lpostup(Ref<MatrixXd> ATmat, Ref<MatrixXd> AImat, 
						Ref<VectorXd> CVec, Ref<VectorXd> DVec,
						double &fppost, double &flpost, 
						const VectorXd& yi, const MatrixXd& VTI,
						const MatrixXd& Fwf, const MatrixXd& Wi, 
						const MatrixXd& Wf, const MatrixXd& Wft, 
						const MatrixXd& Hvh, const MatrixXd& Vh, 
						const double Vdet, const double Wdet) {
	double adetp;
	double adetc;
	double adsqp;
	double adsqc;
	double ysq;
	
	bool isempty = AImat.isZero(0);
	if (isempty) {
		AImat = ATmat;
		adsqp = DVec.transpose()*AImat.inverse()*DVec;
		adetp = AImat.determinant();
	}
	else {
		AI_up(ATmat, AImat, Fwf); // Update previous At to Ai
		adsqp = DVec.transpose()*AImat.inverse()*DVec;
		adetp = AImat.determinant();
	}

	AT_up(ATmat, AImat, Wf, Wft, Wi, Hvh); // Update At
	Ci(CVec, AImat, DVec, Wft); // Update current Ci vector
	Di(DVec, CVec, Vh, yi); // Update current Di vector
 	
	ysq = yi.transpose()*VTI*yi;
	adsqc = DVec.transpose()*ATmat.inverse()*DVec;
	adetc = ATmat.determinant();

/*	cout << "AImat is: " << AImat << endl;
	cout << "ATmat is: " << ATmat << endl;
	cout << "Cvec  is: " << CVec << endl;
	cout << "DVec  is: " << DVec << endl;
	cout << "adetp  is: " << adetp << endl;
	cout << "adetc  is: " << adetc << endl;
	cout << "adsqp  is: " << adsqp << endl;
	cout << "adsqc  is: " << adsqc << endl;
	cout << "fppost is: " << fppost << endl;
	cout << "flpost is: " << flpost << endl;
*/	
	//fppost = get(ppost);
	fppost = fppost-log(2*M_PI*Vdet)-log(2*M_PI*Wdet)-0.5*ysq;
	//set(fppost, ppost);

	//flpost = get(lpost);
    flpost = fppost+0.5*(log(2*M_PI*(adetp*adetc)))+0.5*(adsqp+adsqc);
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
						const MatrixXd& Fwf, const MatrixXd& Wi, 
						const MatrixXd& Wf, const MatrixXd& Wft, 
						const double Wdet) {
	double adetp;
	double adetc;
	double adsqp;
	double adsqc;
	double ysq;

	bool isempty = AImat.isZero(0);
	if (isempty) {
		AImat = ATmat;
		adsqp = DVec.transpose()*AImat.inverse()*DVec;
		adetp = AImat.determinant();
	}
	else {
		AI_up(ATmat, AImat, Fwf); // Update previous At to Ai
		adsqp = DVec.transpose()*AImat.inverse()*DVec;
		adetp = AImat.determinant();
	}

	AT_nup(ATmat, AImat, Wf, Wft, Wi); // Update At
	Ci(CVec, AImat, DVec, Wft); // Update current Ci vector
	DVec = CVec;
 	
	adsqc = DVec.transpose()*ATmat.inverse()*DVec;
	adetc = ATmat.determinant();

	//fppost = get(ppost);
	fppost = fppost-log(2*M_PI*(Wdet));
	//set(fppost, ppost);

	//flpost = get(lpost);
    flpost = fppost+0.5*(log(2*M_PI*(adetp*adetc)))+0.5*(adsqp+adsqc);
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
	Matrix2d WFT;
	WFT = WI*F_MAT.transpose();
	Matrix2d WF;
	WF = WI*F_MAT;
	Matrix2d FWF;
	FWF = F_MAT.transpose()*WI*F_MAT;
	Matrix2d VH;
	VH = VI*H_MAT;
	Matrix2d HVH;
	HVH = H_MAT.transpose()*VI*H_MAT;
	double Vid;
	Vid = VI.determinant();
	double Wid;
	Wid = WI.determinant();

	Matrix2d W0;
	W0 << 0.9,0,0,0.9;
	double lW0det = -log(2*M_PI*W0.determinant());
	double W0sq = MU.transpose()*W0.inverse()*MU;
	
	VectorXd D0(DIMS);
	D0 = MU.transpose()*W0;
	Matrix2d A0;
	A0 = W0.inverse()+FWF;
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
					a.lpost, y_i, VI, FWF, WI, WF, WFT, 
					HVH, VH, Vid, Wid);

		cout << "AImat " << i+1 << " is:\n" << a.aimat << endl;
		cout << "ATmat " << i+1 << " is:\n" << a.atmat << endl;
		cout << "Cvec " << i+1 << " is:\n" << a.cvec << endl;
		cout << "DVec " << i+1 << " is:\n" << a.dvec << endl;
		cout << "ppost " << i+1 << " is: " << a.ppost << endl;
		cout << "lpost " << i+1 << " is: " << a.lpost << endl;
		//cout << "a.leaf_post " << i+1 << " is: " << a.leaf_post << endl;
		cout << endl;


		l.lpostnup(a.atmat, a.atmat, a.cvec, a.dvec, a.ppost, 
					a.lpost, FWF, WI, WF, WFT, Wid);

		cout << "postpp 2 is:" << a.ppost << endl;
		cout << "a.lpost 2 is: " << a.lpost << endl;
	}
}*/
