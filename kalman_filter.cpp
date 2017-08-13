#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft_; 
  Ft_ = F_.transpose();
  P_ = F_ * P_ * Ft_  ;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  //cout << "Entering update function " << endl;
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_ ;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // new state
  x_ = x_ + K * y;

  //size of P_ matrix;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);

  P_ = (I - K * H_) * P_ ;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //cout << "Starting UpdateEKF" << endl;
  float x_tmp = x_(0);
  float y_tmp = x_(1);
  float vx_tmp = x_(2);
  float vy_tmp = x_(3);

  float rho = sqrt(x_tmp*x_tmp + y_tmp*y_tmp);
  float theta = atan2(y_tmp,x_tmp);
  float ro_dot =  ( x_tmp * vx_tmp + y_tmp * vy_tmp)/ rho;

  //cout << "debug1" << endl;
  VectorXd Zpred = VectorXd(3);
  Zpred << rho, theta, ro_dot; 

  VectorXd y = Zpred - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_ ;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // new state
  x_ = x_ + K * y;

  //size of P_ matrix;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);

  //cout << "debug2" << endl;
  //cout << K.size() << endl;

  //cout << H_.size() << endl;
  //cout << P_.size() << endl;
  //cout << I.size() << endl;

  P_ = (I - K * H_) * P_; 
  cout << "Done UpdateEKF" << endl;
}
