#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  for (int i=0; i<estimations.size(); i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    rmse = rmse.array() + (residual.array() * residual.array());
  }

  // Take the mean
  rmse = rmse / estimations.size();
  // Then the square root
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  float px_plus_py = px * px + py * py;

  if (px_plus_py == 0) {
      cout << "Error: Division by zero" << endl;
      return Hj;
  }

  float sqrt_of_px_plus_py = pow(px_plus_py, 0.5);
  float px_plus_py_32 = pow(px_plus_py, 3.0/2.0);

  Hj << px / sqrt_of_px_plus_py, py / sqrt_of_px_plus_py, 0, 0,
        -py / px_plus_py, px / px_plus_py, 0, 0,
        py * (vx * py - vy * px) / px_plus_py_32, px * (vy * px - vx * py) / px_plus_py_32, px / sqrt_of_px_plus_py, py / sqrt_of_px_plus_py;

  return Hj;
}
