#include "../include/convert_matrix2csv.h"

void write2csv (const std::string & path, const MatrixXd& matrix_){
    std::ofstream file(path);
    MatrixXd matrix = matrix_.transpose();
    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j < matrix.cols() - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }
        file.close();
        std::cout << "Matrix written to matrix.csv" << std::endl;
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}
