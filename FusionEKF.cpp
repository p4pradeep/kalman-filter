#include "FusionEKF.h"
#include "tools.h"
#include "math.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  //cout << "Inside FusionEKF";
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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  //cout << "Intialize all the matrices"; 
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;


  noise_ax = 9;
  noise_ay = 9;

  VectorXd x = VectorXd(4);

  MatrixXd P = MatrixXd(4,4);
  P << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1000, 0,
       0, 0, 0, 1000;

  MatrixXd F = MatrixXd(4,4);
  F << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1, 0,
       0, 0, 0, 1;

  MatrixXd H = MatrixXd(2,4);

  MatrixXd R = MatrixXd(2,2);

  MatrixXd Q = MatrixXd(4,4);

  ekf_.Init(x, P, F, H, R, Q);

  //cout << "End" << endl; 
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  //cout << "inside contructor" << endl;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: Intialize the first measurment" << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	// ekf_.x_ = measurement_pack.raw_measurements_;
	//cout << ekf_.x_;
	float Ro = measurement_pack.raw_measurements_[0];
	float theta = measurement_pack.raw_measurements_[1];
	
	ekf_.x_[0] = Ro * cos ( theta );
	ekf_.x_[1] = Ro * sin ( theta );
        ekf_.x_[2] = 0;
        ekf_.x_[3] = 0;
	cout << "ekf_.x = ekf_.x" << endl;	
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	ekf_.x_[0] = measurement_pack.raw_measurements_[0];
	ekf_.x_[1] = measurement_pack.raw_measurements_[1];
        ekf_.x_[2] = 0;
        ekf_.x_[3] = 0;
    }
//
    ekf_.F_ = MatrixXd::Identity(4,4);
    previous_timestamp_ = measurement_pack.timestamp_;
    

    // done initializing, no need to predict or update
    is_initialized_ = true;
    
    //cout << "End of measurement intialization step" << endl ;


    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */


  //cout << "Intialize matrices for Predict step" << endl; 
  // State covariance matric


  //cout << measurement_pack.timestamp_ << endl; 
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
 
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  //cout << ekf_.F_ << endl;
  //cout << noise_ax << endl;
  //cout << noise_ay << endl;
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dt_4 * noise_ax/4, 0, dt_3 * noise_ax/2, 0,
             0, dt_4 * noise_ay/4, 0, dt_3 * noise_ay/2,
             dt_3 * noise_ax/2, 0, dt_2 * noise_ax, 0,
             0, dt_3 * noise_ay/2, 0, dt_2 * noise_ay;

  //cout << ekf_.F_.size() << endl;
  //cout << ekf_.x_.size() << endl;

  //cout << "Start Predict step" << endl; 
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  //cout << "Start update step" << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updatesj;
    //cout << "Perform Radar measurement update" << endl;
    MatrixXd Hj = tools.CalculateJacobian(ekf_.x_);
    //cout << "Calculated Jacobin" << endl;
    ekf_.H_ = Hj;
    ekf_.R_ = R_radar_;
    //cout << "Update measurement" << endl;
    //ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    //cout << "Perform Laser measurement update" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
