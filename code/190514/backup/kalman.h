// FUNCTIONS FOR THE KALMAN FILTER

# pragma once
# include <iostream>
# include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Kalman {
	public:
		VectorXd z_vec;
		MatrixXd r_mat;
		VectorXd y_vec;
		MatrixXd s_mat;
		MatrixXd kal_g;
	
	// Constructors
	Kalman() {};
	Kalman(const int& NDIMS, const int& MDIMS) : z_vec(MDIMS), r_mat(MDIMS, MDIMS), 
	y_vec(NDIMS), s_mat(NDIMS, NDIMS), kal_g(MDIMS, MDIMS) {};
	//Declarations
	void zvec(Ref<VectorXd> z_predict, const MatrixXd& FMAT, const MatrixXd& GMAT,
		const VectorXd& zprev, const VectorXd& uvec); 
	void rmat(Ref<MatrixXd> r_predict, const MatrixXd& FMAT, 
					const MatrixXd& WMAT, const MatrixXd& sigmat);
	void yvec(Ref<VectorXd> y_update, const VectorXd& yi, 
						const MatrixXd& HMAT, const VectorXd& zvec);
	void smat(Ref<MatrixXd> s_update, const MatrixXd& HMAT, 
					const MatrixXd& rmat, const MatrixXd& VMAT);
	void kalg(Ref<MatrixXd> Kalg_update, const MatrixXd& RMAT, 
					const MatrixXd& FMAT, const MatrixXd& smat);
	void zhat(Ref<VectorXd> state_update, const VectorXd& zprev, 
					const MatrixXd& kalg, const VectorXd& yvec);
	void sighat(Ref<MatrixXd> sigma_update, const MatrixXd& rmat, 
						const MatrixXd& HMAT, const MatrixXd& kalg);
	void nzhat(Ref<VectorXd> state_update, const MatrixXd& FMAT);
	void nsighat(Ref<MatrixXd> sigma_update, const MatrixXd& WMAT);
};

// PREDICT

// A priori state estimate (with input)
void Kalman::zvec(Ref<VectorXd> z_predict, const MatrixXd& FMAT, const MatrixXd& GMAT,
		const VectorXd& zprev, const VectorXd& uvec) {
	z_predict = FMAT*zprev + GMAT*uvec;
} 

// A priori covariance estimate matrix
void Kalman::rmat(Ref<MatrixXd> r_predict, const MatrixXd& FMAT, const MatrixXd& WMAT, 
		const MatrixXd& sigmat) {
	r_predict = FMAT*sigmat*FMAT.transpose();
	r_predict += WMAT;
}

// UPDATE

// Innovation or measurement residual
void Kalman::yvec(Ref<VectorXd> y_update, const VectorXd& yi, 
						const MatrixXd& HMAT, const VectorXd& zvec) {
	y_update = yi - (HMAT*zvec);
}
// Innovation or residual covariance
void Kalman::smat(Ref<MatrixXd> s_update, const MatrixXd& HMAT, 
					const MatrixXd& rmat, const MatrixXd& VMAT) {
	s_update.noalias() = HMAT*rmat*HMAT.transpose();
	s_update += VMAT;
}

// Kalman Gain
void Kalman::kalg(Ref<MatrixXd> Kalg_update, const MatrixXd& RMAT, 
					const MatrixXd& FMAT, const MatrixXd& smat) {
	Kalg_update.noalias() = RMAT*FMAT.transpose()*smat.inverse();
}

// A posteriori state estimate
void Kalman:: zhat(Ref<VectorXd> state_update, const VectorXd& zprev, 
					const MatrixXd& kalg, const VectorXd& yvec) {
	state_update = zprev + kalg*yvec;
}

// A posteriori estimate covariance
void Kalman::sighat(Ref<MatrixXd> sigma_update, const MatrixXd& rmat, 
						const MatrixXd& FMAT, const MatrixXd& kalg) {
	sigma_update.noalias() = rmat - kalg*FMAT*rmat;
}

// NOT UPDATE

void Kalman::nzhat(Ref<VectorXd> state_update, const MatrixXd& FMAT) {
                state_update = FMAT*state_update;
}

void Kalman::nsighat(Ref<MatrixXd> sigma_update, 
						const MatrixXd& WMAT) {
                sigma_update += WMAT;
}
/*
int main() {
	int const DIMS = 2;
	//Node a(1);
	Kalman k(DIMS);
	
	VectorXd state(DIMS);	
	state << 0, 0;
	k.state = state;
	VectorXd y_i(DIMS);
	y_i << 1, 1;
	MatrixXd F_MAT(DIMS, DIMS);
	F_MAT << 1,0,0,1;
	Matrix2d H_MAT;
	H_MAT << 1,0,0,1;
	MatrixXd W_MAT(DIMS, DIMS);
	W_MAT << 1,0,0,1; 
	Matrix2d V_MAT;
	V_MAT << 1,0,0,1;
	MatrixXd svar(DIMS, DIMS);
	svar << 1,0,0,1;
	k.svar = svar;
	//VectorXd ztvec(DIMS);
	k.zvec(k.z_vec, F_MAT, k.state); 
	//k.zvec = ztvec;
	cout << "The predicted state is: " << endl
		<< k.z_vec << endl;
	//VectorXd ytvec(DIMS);
	k.yvec(k.y_vec, y_i, H_MAT, k.z_vec);
	//k.yvec = ytvec;
	cout << "The measurement residual is: " << endl 
		<< k.y_vec << endl;

	//MatrixXd rtmat(DIMS, DIMS);
	k.rmat(k.r_mat, F_MAT, W_MAT, k.svar);
	//k.r_mat = rtmat;
	cout << "The predicted state covariance is: " << endl 
	<< k.r_mat << endl;
	//MatrixXd stmat(DIMS, DIMS);
	k.smat(k.s_mat, H_MAT, k.r_mat, V_MAT);
	//k.smat = stmat;
	cout << "The residual covariance is: "<< endl
		<< k.s_mat << endl;
	//MatrixXd tkalg(DIMS, DIMS); 
	k.kalg(k.kal_g, k.r_mat, H_MAT, k.s_mat);
	//k.kalg_g = tkalg;
	cout << "The Kalman gain is: " << endl 
		<< k.kal_g << endl;
	//VectorXd tstate(DIMS);
	k.zhat(k.state, k.z_vec, k.kal_g, k.y_vec);
	//kstate = tstate;
	cout << "The state estimate is: " << endl 
		<< k.state << endl;
	//MatrixXd tsvar(DIMS, DIMS);
	k.sighat(k.svar, k.r_mat, F_MAT, k.kal_g);
	//k.svar = tsvar;
	cout << "The state covariance estimate is: " << endl
		<< k.svar << endl;
}*/
