#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

/**
 * Calculate the RMSE here.
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rv(4);
  rv << 0,0,0,0;
  if (estimations.size() != ground_truth.size() || ground_truth.size()==0)
  {
    cout << "CalculateRMSE: Something is not consistent here.";
    return rv;
  }

  for (size_t i=0; i<estimations.size(); i++)
  {
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array()*diff.array();
    rv += diff;
  }

  return (rv/estimations.size()).array().sqrt();
}

/**
 * Calculate a Jacobian here.
 * It's necessary for the measurement function for processing the radar measurements
 */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  MatrixXd Hj(3,4);

  // state
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];

  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = c1*c2;
  float c4 = vx*py-vy*px;

  if (fabs(c1) < 0.0001)
  {
    cout <<  "CalculateJacobian: Division by Zero";
    return Hj;
  }

  Hj << px/c2   , py/c1 , 0 , 0 ,
      -py/c1    , px/c2 , 0 , 0 ,
      py*c4/c3  , -px*c4/c3,  px/c2,  py/c2;

  return Hj;     
}
