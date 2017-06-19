#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0.0,0.0,0.0,0.0;

  if(estimations.size() != ground_truth.size()) {
    cout << "CalculateRMSE: Size of estimation data and ground truth data differ" << endl;
    return rmse;
  }

  if(estimations.empty()) {
    cout << "CalculateRMSE: Estimation data is empty" << endl;
    return rmse;
  }


  for(size_t i=0; i < estimations.size(); ++i) {
    VectorXd dif = estimations[i]-ground_truth[i];
    rmse += (VectorXd)(dif.array()*dif.array());
  }

  rmse = rmse.array() / estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;
}