#ifndef CONVERT_MATRIX2CSV_H
#define CONVERT_MATRIX2CSV_H

#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
using namespace Eigen;

void write2csv (const std::string & path, const MatrixXd& matrix_);

#endif