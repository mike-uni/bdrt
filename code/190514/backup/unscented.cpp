// THE UNSCENTED KALMAN FILTER

# pragma once
# include <iostream>
# include <Eigen/Dense>
# include <functional>

using namespace std;
using namespace Eigen;

class Unkal {
	public:
		double a = 0.05;
		double b = 2.0;
		int L;
		double kap;
		double lam = a*a*(L+kap)-L;
		double w0m = lam/(double)(L+lam);
		double w0c = w0m+1-a*a+b;
		double wim = 1/(double)(2*(L+lam));
		double gam = sqrt(lam + L); 
		double alpha = 10;
		double beta = 0.2;
		double gamma = 0.1;

		//MatrixXd sigpmat;
		//MatrixXd sigpup;
		//VectorXd zpred;
		//MatrixXd spred;
	
	// Constructors
	Unkal() {};
	Unkal(const int& MDIMS) : L(MDIMS), kap(3 - MDIMS){};
	//Unkal(const int& MDIMS) : sigpmat(MDIMS, 2*MDIMS+1), sigpup(MDIMS, 2*MDIMS+1), 
	//zpred(MDIMS), spred(MDIMS, MDIMS), L(MDIMS), kap(3 - MDIMS){};

	//Declarations
	void mackey(Ref<MatrixXd>sigpup, const MatrixXd& sigpmat, const double& alpha, const double& beta, 
			const double& gamma); 
	void sig_pmat(Ref<MatrixXd> sigpmat, const VectorXd& zprev, const double& gamma, const MatrixXd& LSMAT);
	void z_pred(Ref<VectorXd>zpred, const MatrixXd& sigpup, const double& w0m, const double& wim);
	void s_pred(Ref<MatrixXd>spred, const MatrixXd& sigpup, const VectorXd& zpred, const double& wim,
			const double& w0c, const MatrixXd& LW);
	void p_mat(Ref<MatrixXd> PMAT, const MatrixXd& XU, const MatrixXd& YU, const VectorXd& zpred,
			const VectorXd& ypred, const double& wim, const double& w0c); 
	void k_mat(Ref<MatrixXd> KMAT, const MatrixXd& PMAT, const MatrixXd& SY); 
	void z_hat(Ref<VectorXd>zhat, const MatrixXd& KMAT, const VectorXd& zpred, const VectorXd& yi, 
			const VectorXd& ypred); 
	void u_mat(Ref<MatrixXd> UMAT, const MatrixXd& KMAT, const MatrixXd& SY); 
	void s_hat(Ref<MatrixXd> SU, const MatrixXd& SX, const MatrixXd& UMAT); 
};

void Unkal::sig_pmat(Ref<MatrixXd> sigpmat, const VectorXd& zprev, const double& gamma, const MatrixXd& LSMAT) {
	sigpmat << zprev, (gamma*LSMAT).colwise() + zprev, (-gamma*LSMAT).colwise() + zprev;
}

void Unkal::mackey(Ref<MatrixXd>sigpup, const MatrixXd& sigpmat, const double& alpha, const double& beta, 
		const double& gamma) {
	MatrixXd FMAT(2,2);	
	FMAT << (1-gamma), 0, 0, beta;
	VectorXd temp(2);
	double ztemp;
	for (int i = 0; i < 5; ++i) {
		ztemp = sigpmat(1,i);
		ztemp = ztemp/(double)(1+pow(ztemp, alpha));
		temp << sigpmat(0,i), ztemp; 
		sigpup.col(i) << FMAT*temp;
	}
}	

void Unkal::z_pred(Ref<VectorXd>zpred, const MatrixXd& sigpup, const double& w0m, const double& wim) {
	size_t DIMS = sigpup.rows();
	zpred = w0m*(sigpup.col(0));
	zpred = ((wim*sigpup.block(0, 1, DIMS, 2*DIMS)).rowwise().sum()).colwise() + zpred;
}

void Unkal::s_pred(Ref<MatrixXd>tspred, const MatrixXd& sigpup, const VectorXd& tzpred, const double& wim,
		const double& w0c, const MatrixXd& LW) {
	size_t DIMS = sigpup.rows();
	MatrixXd temp(DIMS, 2*DIMS + DIMS);
	temp << sqrt(wim)*(sigpup.block(0, 1, DIMS, 2*DIMS).colwise() - tzpred), LW;
	//cout << "temp " << endl;
	//cout << temp << endl; 
	//MatrixXd temp1;
	temp = temp.transpose().colPivHouseholderQr().matrixR().triangularView<Upper>();
	//cout << "temp " << endl;
	//cout << temp << endl;
	//MatrixXd temp2;
	//temp2 = temp.transpose().colPivHouseholderQr().matrixQ();//.triangularView<Upper>();
	//cout << "temp2 " << endl;
	//cout << temp2 << endl;
	MatrixXd tochol(DIMS, DIMS);
	tochol = temp.transpose()*temp;
	//cout << "tochol " << endl;
	//cout << tochol << endl; 
	LLT<MatrixXd> chSPRED(tochol);
	chSPRED = chSPRED.rankUpdate((sigpup.col(0) - tzpred), w0c);
	//cout << "temp " << endl;
	//cout << temp << endl; 
	//MatrixXd LSPRED = chSPRED.matrixL();
	tspred = chSPRED.matrixL(); //temp.triangularView<Lower>();//LSPRED;
	//cout << "spred " << endl;
	//cout << spred << endl; 
}

void Unkal::p_mat(Ref<MatrixXd> PMAT, const MatrixXd& XU, const MatrixXd& YU, const VectorXd& zpred,
		const VectorXd& ypred, const double& wim, const double& w0c) {
	size_t MDIM = XU.rows();
	size_t NDIM = YU.rows();
	PMAT = (wim*(XU.block(0,1,MDIM,2*MDIM).colwise()-zpred)*(YU.block(0,1,NDIM,2*MDIM).colwise()-ypred).transpose())+
	(w0c*((XU.col(0)-zpred)*(YU.col(0)-ypred).transpose()));
}

void Unkal::k_mat(Ref<MatrixXd> KMAT, const MatrixXd& PMAT, const MatrixXd& SY) {
	KMAT = (PMAT.transpose().colPivHouseholderQr().solve(SY.transpose())).transpose().colPivHouseholderQr().solve(SY);
}

void Unkal::z_hat(Ref<VectorXd>zhat, const MatrixXd& KMAT, const VectorXd& zpred, const VectorXd& yi, const VectorXd& ypred) {
	zhat = zpred + KMAT*(yi - ypred);
}

void Unkal::u_mat(Ref<MatrixXd> UMAT, const MatrixXd& KMAT, const MatrixXd& SY) {
	//MatrixXd UMAT1;
	UMAT = KMAT*SY;
	//cout << UMAT1 << endl;
}

void Unkal::s_hat(Ref<MatrixXd> SU, const MatrixXd& SX, const MatrixXd& UMAT) {
	LDLT<MatrixXd> chS(SX);
	chS = chS.rankUpdate(UMAT, -1);
	SU = chS.matrixL();
}
/*
int main() {
	Unkal U(2);
	
	// Initialise z0 and S0
	MatrixXd S0(U.L,U.L); 
	S0 << 1,0.5,0.5,1;
	LLT<MatrixXd> chSMAT(S0);
	MatrixXd LSMAT = chSMAT.matrixL();
	VectorXd zhat(2);
	zhat << 1,1;
	
	// Sigma point calculation and update
	MatrixXd sigpmat(U.L, 2*U.L+1);
	U.sig_pmat(sigpmat, zhat, U.gam, LSMAT);
	
	MatrixXd sigpup(U.L, 2*U.L+1);
	U.mackey(sigpup, sigpmat, U.alpha, U.beta, U.gamma);
	
	VectorXd zpred(U.L);
	U.z_pred(zpred, sigpup, U.w0m, U.wim);

	MatrixXd W(2,2);
	W << 1,0,0,1;
	LLT<MatrixXd> chW(W);
	MatrixXd LW = chW.matrixL();

	MatrixXd SX(U.L, U.L);
	U.s_pred(SX, sigpup, zpred, U.wim, U.w0c, LW);

	MatrixXd XU(U.L, 2*U.L+1);
	U.sig_pmat(XU, zpred, gam, SX);

	cout << "sigpmat " << endl;
	cout << sigpmat << endl;
	cout << "sigpup " << endl;
	cout << sigpup << endl;
	cout << "zpred " << endl;
	cout << zpred << endl;
	cout << "SX " << endl;
	cout << SX << endl;
	cout << "XU " << endl;
	cout << XU << endl;
	
	// Y pred
	int NDIM = 1;
	VectorXd H(U.L);
	H << 1,0;
	MatrixXd V(1,1);
	V << 0.5;	
	
	MatrixXd YU(NDIM, 2*U.L+1);
	YU = H.transpose()*XU;
	cout << "YU " << endl;
	cout << YU << endl;

	VectorXd ypred(NDIM);
	U.z_pred(ypred, YU, U.w0m, U.wim);
	cout << "ypred " << endl;
	cout << ypred << endl;

	// Measurement Update
	MatrixXd SY(NDIM, NDIM);
	U.s_pred(SY, YU, ypred, U.wim, U.w0c, V);
	cout << "SY " << endl;
	cout << SY << endl;

	MatrixXd PMAT(U.L, NDIM);
	U.p_mat(PMAT, XU, YU, zpred, ypred, U.wim, U.w0c); 
	cout << "PMAT " << endl;
	cout << PMAT << endl;
	
	MatrixXd KMAT(U.L, NDIM);	
	U.k_mat(KMAT, PMAT, SY);
	cout << "KMAT " << endl;
	cout << KMAT << endl;	

	VectorXd yi(1);
	yi << 0.5;

	VectorXd zup(U.L);
	U.z_hat(zup, KMAT, zpred, yi, ypred); 
	cout << "zup " << endl;
	cout << zup << endl;

	MatrixXd UMAT(U.L, NDIM);
	U.u_mat(UMAT, KMAT, SY); 
	cout << "UMAT " << endl;
	cout << UMAT << endl;
	
	MatrixXd SUPMAT(U.L, U.L);
	U.s_hat(SUPMAT, SX, UMAT);
	cout << "SUPMAT " << endl;
	cout << SUPMAT << endl;
}*/
