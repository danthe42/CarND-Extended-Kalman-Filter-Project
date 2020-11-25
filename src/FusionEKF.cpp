#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  noise_ax = 9;
  noise_ay = 9;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    previous_timestamp_ = measurement_pack.timestamp_;            

    MatrixXd f = MatrixXd( 4, 4 );
    f <<  1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;
               
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      int rho = measurement_pack.raw_measurements_[0];
      int phi = measurement_pack.raw_measurements_[1];
      int rho_dot = measurement_pack.raw_measurements_[2];

      VectorXd x = VectorXd(4);
      x <<  rho * cos(phi),
            rho * sin(phi),
            rho_dot * cos(phi),
            rho_dot * sin(phi);

      // state covariance matrix. With radar we are as certain in velocity as in position, so all those probabilities are 1.  
      MatrixXd p = MatrixXd( 4, 4 );
      p << 1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1, 0,
           0, 0, 0, 1;     

      ekf_.Init( x, p, f, ekf_.H_, R_radar_, ekf_.Q_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      VectorXd x = VectorXd(4);
      x <<  measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1],
            0,
            0;
      
      // state covariance matrix. At the beginning, we have a big uncertainity of vx and vy (valocity), so those values are big initial numbers.  
      MatrixXd p = MatrixXd( 4, 4 );
      p << 1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1000, 0,
           0, 0, 0, 1000;     

      ekf_.Init( x, p, f, H_laser_, R_laser_, ekf_.Q_);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float dt = ( measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;            

  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;

  // Update the state transition matrix F
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // calculate the Process covariance matrix
  ekf_.Q_ = MatrixXd( 4, 4 );
  ekf_.Q_ <<  dt4/4*noise_ax ,  0               , dt3/2*noise_ax , 0 ,
              0              ,  dt4/4*noise_ay  , 0              , dt3/2*noise_ay , 
              dt3/2*noise_ax ,  0               , dt2*noise_ax   , 0 ,
              0              ,  dt3/2*noise_ay  , 0              , dt2*noise_ay;


  /**
   * Prediction
   */

  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    ekf_.H_ = tools.CalculateJacobian( ekf_.x_ );
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
