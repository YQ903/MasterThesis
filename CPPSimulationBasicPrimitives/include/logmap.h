#ifndef LOGMAP_H
#define LOGMAP_H
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <dqrobotics/DQ.h>
#include <dqrobotics/robot_modeling/DQ_SerialManipulator.h>
#include <functional>

using namespace Eigen;
using namespace DQ_robotics;


MatrixXd logmap(const MatrixXd &M_1, const MatrixXd &M_2, const double &step);
MatrixXd expmap(const MatrixXd &M_1, const MatrixXd &M_2);
MatrixXd para_trans(const MatrixXd &M_old, const MatrixXd &M_new, const MatrixXcd &T);

#endif