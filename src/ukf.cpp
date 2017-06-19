#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define PI 3.14159265358979323846

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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

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

  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_x_;

  P_ == MatrixXd::Identity(n_x_, n_x_);

  n_sigma_ = 2 * n_aug_ + 1;

  //create predicted sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(n_sigma_);

  //set weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  double weight = 1.0/(2.0*(lambda_+n_aug_));

  for (size_t i=1; i < n_sigma_; ++i) {
  	weights_(i) = weight;
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_) {

	  double px = 0.0;
	  double py = 0.0;
	  double v = 0.0;

	  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	  	// Convert from polar to cartesian coordinates
	    double rho = meas_package.raw_measurements_[0];
	    double phi = meas_package.raw_measurements_[1];
	    double rhod = meas_package.raw_measurements_[2];
	      
	    // Initialize position
	    px = rho * cos(phi);
	    py = rho * sin(phi);

	    double vx = rhod * cos(phi);
	    double vy = rhod * sin(phi);
	    v = sqrt(vx*vx + vy*vy);
	  }
	  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

	    // Initialize position
	    px = meas_package.raw_measurements_[0];
	    py = meas_package.raw_measurements_[1];
	  }

	  x_ << px, py, v, 0.0, 0.0;

		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;

		return;
	}

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0 ;

  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
  	UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
  	UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {

  // Create Sigma Points of current state
  MatrixXd sigma_points;
  GenerateAugmentedSigmaPoints(sigma_points);

  // Predict Sigma Points after timestep dt
  PredictSigmaPoints(sigma_points, dt);

  // Predict new mean and covariance
  PredictMeanAndCovariance();

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred(n_z);
  
  //measurement covariance matrix S
  MatrixXd S(n_z,n_z);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_laspx_*std_laspx_, 0.0,
        0.0, std_laspy_*std_laspy_;
	
  //fill sigma points in measurement space
  //juxt px and py for lidar measurements
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
      Zsig.col(i) << Xsig_pred_(0,i), Xsig_pred_(1,i);
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
      z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  //and cross correlation Tc
  S.fill(0.0);
  Tc.fill(0.0);
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
  		VectorXd z_diff = Zsig.col(i)-z_pred;

      S += weights_(i) * z_diff * z_diff.transpose();
      
      Tc += weights_(i)*(Xsig_pred_.col(i)-x_)*z_diff.transpose();
  }
  
  S += R;
  
  //create vector for incoming lidar measurement
  VectorXd z = VectorXd(n_z);
  z <<	meas_package.raw_measurements_[0],
  			meas_package.raw_measurements_[1];
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ += K * (z-z_pred);
  P_ -= K*S*K.transpose();

  cout << "x:" << endl << x_ << endl;
  cout << "P:" << endl << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred(n_z);
  
  //measurement covariance matrix S
  MatrixXd S(n_z,n_z);

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0.0, 0.0,
        0.0, std_radphi_*std_radphi_, 0.0,
        0.0, 0.0, std_radrd_*std_radrd_;

  //transform sigma points into measurement space
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double psi = Xsig_pred_(3,i);
      
      double rho = sqrt(px*px+py*py);
      double phi = atan2(py,px);
      double rhod = (px*cos(psi)*v+py*sin(psi)*v)/rho;
      
      Zsig.col(i) << rho, phi, rhod;
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
      z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
  		VectorXd z_diff = Zsig.col(i)-z_pred;

			while (z_diff(1)> PI) z_diff(1) -= 2.0*PI;
			while (z_diff(1)<-PI) z_diff(1) += 2.0*PI;

      S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  S += R;

  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<	meas_package.raw_measurements_[0],
  			meas_package.raw_measurements_[1],
  			meas_package.raw_measurements_[2];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (size_t i=0; i < 2 * n_aug_ + 1; ++i) {
  		VectorXd z_diff = Zsig.col(i)-z_pred;

  		while (z_diff(1) > PI) z_diff(1) -= 2.0*PI;
  		while (z_diff(1) <-PI) z_diff(1) += 2.0*PI;

  		VectorXd x_diff = Xsig_pred_.col(i)-x_;

  		while (x_diff(3) > PI) x_diff(3) -= 2.0*PI;
  		while (x_diff(3) <-PI) x_diff(3) += 2.0*PI;

      Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z-z_pred;

  while (z_diff(1) > PI) z_diff(1) -= 2.0*PI;
  while (z_diff(1) <-PI) z_diff(1) += 2.0*PI;

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  cout << "x:" << endl << x_ << endl;
  cout << "P:" << endl << P_ << endl;
}

void UKF::GenerateAugmentedSigmaPoints(MatrixXd& sigma_out) const {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  sigma_out = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug.tail(2) << 0.0, 0.0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) << std_a_*std_a_, 0.0, 0.0, std_yawdd_*std_yawdd_;

  //create square root of augmented covariance matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  sigma_out.col(0) = x_aug;
  
  for (int i = 0; i < n_aug_; ++i) {
      sigma_out.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
      sigma_out.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
  }
}

void UKF::PredictSigmaPoints(const MatrixXd& sigma_in, const double dt) {

 	double dt2_2 =0.5 * dt * dt;

  //predict sigma points
  for (size_t i = 0; i < sigma_in.cols(); ++i) {
      VectorXd sig_in = sigma_in.col(i);
      VectorXd sig_out = VectorXd(n_x_);
      
  		double sin_psi = sin(sig_in(3));
  		double cos_psi = cos(sig_in(3));

      if (fabs(sig_in(4)) > 0.001) {
          sig_out <<    sig_in(2)/sig_in(4) * (sin(sig_in(3)+sig_in(4)*dt) - sin_psi),
                        sig_in(2)/sig_in(4) * (-cos(sig_in(3)+sig_in(4)*dt) + cos_psi),
                        0.0,
                        sig_in(4)*dt,
                        0.0;
      } else {
          sig_out <<    sig_in(2)*cos_psi*dt,
                        sig_in(2)*sin_psi*dt,
                        0.0,
                        0.0,
                        0.0;
      }
      
      sig_out(0) += dt2_2 * cos_psi * sig_in(5);
      sig_out(1) += dt2_2 * sin_psi * sig_in(5);
      sig_out(2) += dt * sig_in(5);
      sig_out(3) += dt2_2 * sig_in(6);
      sig_out(4) += dt * sig_in(6);
      
      sig_out += sig_in.topRows(5);
      
      Xsig_pred_.col(i) = sig_out;
  }
}

void UKF::PredictMeanAndCovariance() {

  //predict state mean
  x_.fill(0.0);
  for (size_t i=0; i < n_sigma_; ++i) {
      x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix  
  P_.fill(0.0);
  for (size_t i=0; i < n_sigma_; ++i) {
  		VectorXd x_diff = Xsig_pred_.col(i)-x_;

  		while(x_diff(3) > PI) x_diff(3) -= 2*PI;
  		while(x_diff(3) < -PI) x_diff(3) += 2*PI;

      P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}