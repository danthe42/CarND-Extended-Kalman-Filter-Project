#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in ) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();  
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_* x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si; 

  // new estimate
  x_ = x_ + ( K * y );
  size_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity( x_size, x_size );
  P_ = ( I - K * H_ ) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // In case of radar we have the following measurement values:
  //  rho - radial distance from origin
  //  phi - angle (measured counter-clockwise, the forward direction is 0 degrees.) 
  //  rho_dot - range rate (velocity in the direction of phi, "away" from us)

  // calculate h(x')
  float rho = sqrt( x_[0] * x_[0] + x_[1] * x_[1] );

  float rho_dot = 0.0;
  float phi = 0.0;

  if (rho < 0.0001 || 
        (fabs(x_[0])<0.00001 && fabs(x_[1])<0.00001))
  {
    // It's and error. Normally shouldn't happen because the radar's destination can never be so close to the sensor.
    // We need to handle these cases: 
    //  a) we need to divide with c1 later, it would cause division by zero exception 
    //  b) or atan2(y,x) is "undefined" at the origin, where both values are zero. Some platform/compiler handle this case somehow, others just give back NaN or throw exception... 
    cout <<  "UpdateEKF: Division by Zero. Plase check the radar sensor. This can indicate hardware failure, or just that the sensor's surface is dirty.";    
    // in this case use the default 0 for rho_dot and phi, it's still better than application crash
  } 
  else 
  {
    phi = atan2( x_[1], x_[0] );
    rho_dot = ( x_[0]*x_[2] + x_[1]*x_[3]) / rho;
  }

  VectorXd h_func(3);
  h_func <<                      rho              ,
                                 phi              ,
                               rho_dot            ;
  VectorXd y = z - h_func;

  // put y[1] in the proper range
  while (y[1]<-M_PI) y[1] += 2*M_PI;
  while (y[1]>M_PI) y[1] -= 2*M_PI;

  // Here, H_ is already replaced with Hj, the jacobian matrix.
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si; 

  // new estimate
  x_ = x_ + ( K * y );
  size_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity( x_size, x_size );
  P_ = ( I - K * H_ ) * P_;
}
