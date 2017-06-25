#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */


  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  //lambda_ = 3 - n_aug_;
  lambda_ = 3 - n_x_;


  //create predicted sigma points matrix in state space
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_  + 1);
  

  ///* time when the state is true, in us
  time_us_ = 0.0;


  // the current NIS for radar
  NIS_radar_ = 0.0;


  // the current NIS for laser
  NIS_laser_ = 0.0;

}



UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //skip prediction/update if sensor type is ignored
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
    (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/

    if (!is_initialized_) {

      // first measurement
      cout << "Unscented Kalman filter initializing" << endl;
   
      // initial measurement
      x_ = VectorXd::Zero(5);
      //x_ << 1,1,1,1,0.1;

      // initial state covariance matrix P
      P_ << 0.15, 0, 0, 0, 0,
            0, 0.15, 0, 0, 0,
            0,    0, 1, 0, 0,
            0,    0, 0, 1, 0,
            0,    0, 0, 0, 1;

      // init timestamp
      time_us_ = meas_package.timestamp_;

      if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) 
      {
        /**
        Initialize state.
        */
        x_(0) = meas_package.raw_measurements_(0); // px
        x_(1) = meas_package.raw_measurements_(1); // py
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) 
      {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */

        float rho     = meas_package.raw_measurements_(0);
        float phi     = meas_package.raw_measurements_(1);
        float rho_dot = meas_package.raw_measurements_(2);
      
        // normalize phi (y(1)) between -pi and pi
        //while (phi < -M_PI) phi += 2 * M_PI;     // phi < - pi => phi+2*pi < pi
        //while (phi > M_PI)  phi -= 2 * M_PI;    // phi > pi   => phi-2*pi> -pi

        x_(0) = rho * cos(phi); //px
        x_(1) = rho * sin(phi); //py
      }
      // done initializing, no need to predict or update
      is_initialized_ = true;
    
      cout << "x_initial: " << x_ << endl;
      cout << "P_initial: " << P_ << endl;
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

    //compute the time elapsed between the current and previous measurements
    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; //delta_t - expressed in seconds
    time_us_ = meas_package.timestamp_;

    //Call the Kalman Filter predict() function
    Prediction(delta_t);


    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
    TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
    */

    //VectorXd z = measurement_pack.raw_measurements_;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {
      // Laser updates
      UpdateLidar(meas_package);
      cout << "x_update_lidar= " << x_ << endl;
      cout << "P_update_lidar= " << P_ << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) 
    {
      // Radar updates
      UpdateRadar(meas_package);
      cout << "x_update_radar= " << x_ << endl;
      cout << "P_update_radar= " << P_ << endl;
    } 
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  /*****************************************************************************
   *  Generate sigma points 
   ****************************************************************************/

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P (Standard Cholesky decomposition (LL^T) of a matrix)
  MatrixXd A = P_.llt().matrixL(); // P = A*A^T or A = sqrt(P)

  //define spreading parameter for non-augmented
  lambda_ = 3 - n_x_;

  //set sigma points as columns of matrix Xsig
  Xsig.col(0) = x_;  
  
  for (int i = 0; i < n_x_; ++i)
  {
    Xsig.col(i+1)      =     x_ + sqrt(lambda_ + n_x_) * A.col(i) ;
    Xsig.col(n_x_+1+i) = x_ - sqrt(lambda_ + n_x_) * A.col(i) ;
  }
  //Xsig.block(0, 1, n_x,n_x)      = x.replicate(1,n_x) + sqrt(lambda + n_x) * A ;
  //Xsig.block(0, n_x +2, n_x,n_x) = x.replicate(1,n_x) - sqrt(lambda + n_x) * A ;
  

  /*****************************************************************************
   *  Augmented sigma points 
   ****************************************************************************/

  //define spreading parameter for augmented
  lambda_ = 3 - n_aug_;


  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);


  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);


  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


 //create augmented mean state
  x_aug = VectorXd ::Zero(n_aug_);
  x_aug.head(n_x_) = x_;
  
  //create augmented covariance matrix
  P_aug = MatrixXd ::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  //P_aug(n_aug_ - 2, n_aug_ - 2) = std_a_*std_a_;
  //P_aug(n_aug_ - 1, n_aug_ - 1) = std_yawdd_*std_yawdd_;
  
  MatrixXd Q = MatrixXd(2, 2);
  Q = MatrixXd ::Zero(2, 2);
  Q(0,0) = std_a_ * std_a_;
  Q(1,1) = std_yawdd_ * std_yawdd_;
  P_aug.bottomRightCorner(2,2) = Q;

  //create square root matrix P_aug
  MatrixXd A_aug = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1)          = x_aug + sqrt(n_aug_ + lambda_) * A_aug.col(i);  
    Xsig_aug.col(n_aug_ + 1 + i) = x_aug - sqrt(n_aug_ + lambda_) * A_aug.col(i);
  }


  /*****************************************************************************
   *  Predict sigma points 
   ****************************************************************************/

  //avoid division by zero
  //write predicted sigma points into right column
  //predict sigma points  
  for (int i =0; i < 2 * n_aug_ + 1; ++i){      
      double px =          Xsig_aug(0,i);      
      double py =          Xsig_aug(1,i);      
      double v =           Xsig_aug(2,i);      
      double yawn =        Xsig_aug(3,i);      
      double yawn_dot =    Xsig_aug(4,i);      
      double nu_a =        Xsig_aug(5,i);      
      double nu_yawn_dot = Xsig_aug(6,i);            
      
      //predicted state values      
      double px_pred, py_pred;            
      
      // before adding noise      
      if (fabs(yawn_dot) > 0.001)
      {          
          px_pred = px + v/yawn_dot*(sin(yawn + delta_t*yawn_dot) - sin(yawn));          
          py_pred = py + v/yawn_dot*(-cos(yawn + delta_t*yawn_dot) + cos(yawn));      
          
      }      
      else
      {          
          px_pred = px + v * cos(yawn) * delta_t;          
          py_pred = py + v * sin(yawn) * delta_t;      
          
      }      

      double v_pred = v + 0;      
      double yawn_pred = yawn + yawn_dot * delta_t;      
      double yawn_dot_pred = yawn_dot + 0;            
      
      // adding noise      
      px_pred = px_pred + 1/2*delta_t*delta_t*cos(yawn)*nu_a;      
      py_pred = py_pred + 1/2*delta_t*delta_t*sin(yawn)*nu_a;      
      //px_pred = px_pred + 1/2*pow(delta_t,2)*cos(yawn);      
      //py_pred = py_pred + 1/2*pow(delta_t,2)*sin(yawn);      
      v_pred = v_pred + delta_t*nu_a;      
      yawn_pred = yawn_pred + 1/2*pow(delta_t,2)*nu_yawn_dot;      
      yawn_dot_pred = yawn_dot_pred  + delta_t*nu_yawn_dot;            
      
      // predicted sigma points      
      Xsig_pred_(0,i) = px_pred;      
      Xsig_pred_(1,i) = py_pred;      
      Xsig_pred_(2,i) = v_pred;      
      Xsig_pred_(3,i) = yawn_pred;      
      Xsig_pred_(4,i) = yawn_dot_pred;  
      
  }  


  /*****************************************************************************
   *  Predict Mean and Covariance
   ****************************************************************************/

  //set vector for weights
  double weight_0 = lambda_/(lambda_ + n_aug_);

  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) 
  {  
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
 
  //predict state mean (5x15 * 15x1 = 5x1)
  x_ =VectorXd::Zero(n_x_);
  x_ = Xsig_pred_ * weights_;
  //for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  //{  //iterate over sigma points
  //  x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  //}
  
  //predict state covariance matrix  (5x5)
  P_ = MatrixXd::Zero(n_x_, n_x_);
  
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    VectorXd x_diff =   Xsig_pred_.col(i) - x_;  // (5x1 - 5x1 = 5x1)
    
    // x_diff(3) : yawn angle
    // normalize yawn angle into [-pi, pi], 
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI; //yawn angle > pi => yawn angle-2*pi > -pi
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI; //yawn angle < -pi => yawn angle+2*pi < pi

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose(); //(5x5 +  5x1 * 1x5 = 5x5)
  }
}





/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // extract measurement z
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure px, py
  int n_z = 2;

 
  /*****************************************************************************
   *  Predict LIDAR Measurement 
   ****************************************************************************/

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  //transform sigma points into measurement space
  for (int i=0; i< 2 * n_aug_ + 1; i++)
  {
    double px =  Xsig_pred_(0,i);
    double py =  Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = px ;  //px
    Zsig(1,i) = py ; // py   
  }
  

  //calculate mean predicted measurement (2x15 * 15x1 = 2x1)
  z_pred = VectorXd::Zero(n_z);
  z_pred = Zsig * weights_;
  //for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  //{
  //  z_pred = z_pred + weights_(i) * Zsig.col(i);
  //}



  //calculate measurement covariance matrix S (2x2)
  S = MatrixXd::Zero(n_z, n_z);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++){

    //residual 
    VectorXd z_diff =   Zsig.col(i) - z_pred;
    

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

 //add measurement noise
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R(0,0) = std_laspx_ * std_laspx_;
  R(1,1) = std_laspy_ * std_laspy_;

  S = S + R ;




  /*****************************************************************************
   *  Update LIDAR Measurement -  Mean and Covariance 
   ****************************************************************************/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix (5x1 * 1x2 = 5x2)
  Tc = MatrixXd::Zero(n_x_, n_z);
  
  // Xsig_pred : sigma points in prediction 
  // Zsig : sigma points in measurement 
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd x_diff =   Xsig_pred_.col(i) - x_; // n_x column vector 
    VectorXd z_diff =   Zsig.col(i) - z_pred; // n_z column vector 
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K; (5x2)
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff =   z - z_pred; // n_z column vector
  

  // calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  
  P_ = P_ - K * S * K.transpose();
}







/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // extract measurement z
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;


  /*****************************************************************************
   *  Predict RADAR Measurement 
   ****************************************************************************/

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  //transform sigma points into measurement space
  for (int i=0; i< 2 * n_aug_ + 1; i++)
  {
    double px =  Xsig_pred_(0,i);
    double py =  Xsig_pred_(1,i);
    double v =   Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    
    Zsig(0,i) = sqrt(px*px + py*py) ;  //rho
    Zsig(1,i) = atan2(py,px); // phi
    Zsig(2,i) = (px*cos(yaw)*v + py*sin(yaw)*v)/sqrt(px*px + py*py); //rho_dot
    
  }
  
  //calculate mean predicted measurement (3x15 * 15x1 = 3x1)
  z_pred = VectorXd::Zero(n_z);
  z_pred = Zsig * weights_;
  //for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  //{
  //  z_pred = z_pred + weights_(i) * Zsig.col(i);
  //}


  //calculate measurement covariance matrix S (3x3)
  S = MatrixXd::Zero(n_z, n_z);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff =   Zsig.col(i) - z_pred;
    
    //z_diff(1) : phi
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;//phi > pi => phi - 2*pi > -pi
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;//phi < -pi => phi + 2*pi < pi
    
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

 //add measurement noise
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;  
  R(2,2) = std_radrd_ * std_radrd_;
  
  S = S + R ;




  /*****************************************************************************
   *  Update RADAR Measurement -  Mean and Covariance 
   ****************************************************************************/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix (5x1 * 1x3 = 5x3)
  Tc = MatrixXd::Zero(n_x_, n_z);
  
  // Xsig_pred : sigma points in prediction 
  // Zsig : sigma points in measurement 
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd x_diff =   Xsig_pred_.col(i) - x_; // n_x column vector 
    VectorXd z_diff =   Zsig.col(i) - z_pred; // n_z column vector 
    
    //z_diff(1) : phi into (-pi, pi)
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI; // phi>pi => phi-2*pi> -pi
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;// phi< - pi => phi+2*pi < pi
    
    //x_diff(3) : yawn angle into (-pi, pi)
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;// yawn angle > pi => yawn angle-2*pi> -pi
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;// yawn angle < - pi => yawn angle+2*pi < pi

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K; (5x3)
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff =   z - z_pred; // n_z column vector
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  

  // calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // update state mean and covariance matrix
  //x = x + K*(z - z_pred);
  x_ = x_ + K * z_diff;
  
  P_ = P_ - K * S * K.transpose();
}
