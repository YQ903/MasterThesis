// Not the same as in Matlab
// input arguments are Matrix instead of Tensor
#include "../include/logmap.h"

MatrixXd logmap(const MatrixXd &M_1, const MatrixXd &M_2, const double &step){
    int size = M_1.rows();
    PartialPivLU<MatrixXd> lu(M_2); // WHy different from Matlab mldivide????
    MatrixXd div = lu.solve(M_1);
    // std::cout<<"div: "<<std::endl<<div<<std::endl;
    EigenSolver<MatrixXd> eigensolver(div); // WHy different from Matlab eig()????


    double M_1_det = M_1.determinant();
    double M_2_det = M_2.determinant();
    
    auto EV = (eigensolver.eigenvalues()); // / M_1_det * M_2_det;

    VectorXcd EV_log = EV.array().log();
    MatrixXcd M_log = EV_log.asDiagonal();
    // std::cout<<"eigenvalue log: "<<std::endl<<EV<<std::endl;
    // std::cout<<"eigenvalue log M: "<<std::endl<<M_log<<std::endl;
    // std::cout<<"M_1: "<<M_1<<std::endl;

    MatrixXcd EVec = eigensolver.eigenvectors();
    MatrixXcd M_pow = EVec.inverse();
    MatrixXcd M_temp = M_2 * EVec * M_log * M_pow *  step;
    MatrixXd result(size,size);
    result = M_temp.real();
    // std::cout<<"result: "<<std::endl<<result<<std::endl;
    return result;
}


MatrixXd expmap(const MatrixXd &M_1, const MatrixXd &M_2){
    int size = M_1.rows();
    PartialPivLU<MatrixXd> lu(M_2); // WHy different from Matlab mldivide????
    MatrixXd div = lu.solve(M_1);

    EigenSolver<MatrixXd> eigensolver(div); // WHy different from Matlab eig()????

    double M_1_det = M_1.determinant();
    double M_2_det = M_2.determinant();
    
    auto EV = (eigensolver.eigenvalues());// / M_1_det * M_2_det;

    VectorXcd EV_log = EV.array().exp();
    MatrixXcd M_log = EV_log.asDiagonal();
    // std::cout<<"eigenvalue log: "<<std::endl<<EV<<std::endl;
    // std::cout<<"eigenvalue log M: "<<std::endl<<M_log<<std::endl;
    // std::cout<<"M_1: "<<M_1<<std::endl;

    MatrixXcd EVec = eigensolver.eigenvectors();
    MatrixXcd M_pow = EVec.inverse();
    MatrixXcd M_temp = M_2 * EVec * M_log * M_pow;
    MatrixXd result(size,size);
    result = M_temp.real();
    // std::cout<<"result: "<<std::endl<<result<<std::endl;
    return result;
}



MatrixXd para_trans(const MatrixXd &M_old, const MatrixXd &M_new, const MatrixXcd &T){

    int size = M_old.rows();
    PartialPivLU<MatrixXd> lu(M_new); // WHy different from Matlab mldivide????
    MatrixXd div = lu.solve(M_old);

    // EigenSolver<MatrixXd> eilver(div); 

    // auto EV = (eigensolver.eigenvalues());// / M_1_det * M_2_det;
    // MatrixXd M_temp;
    // M_temp = (div.array().sqrt()).real();


    // // std::cout<<"div: "<<std::endl<<div<<std::endl;

    // MatrixXd Ac_1;
    MatrixXcd T_tilde;
    MatrixXd T_tilde_real;
    MatrixXcd Ac;

    // Ac_1 = (M_new * M_old.inverse());

    // std::cout<<" Ac_1"<<Ac_1<<std::endl;
        

    // // PartialPivLU<MatrixXd> lu(Ac_1); // WHy different from Matlab mldivide????
    // // MatrixXd div = lu.solve(Ac_1);
    // // // std::cout<<"div: "<<std::endl<<div<<std::endl;
    // // SelfAdjointEigenSolver<MatrixXd> sqrt_solver(div); // WHy different from Matlab eig()????

    // // PartialPivLU<MatrixXd> lu(Ac_1); // WHy different from Matlab mldivide????
    // // MatrixXd div = lu.solve(Ac_1);
    // // EigenSolver<MatrixXd> es(A);

    // // std::cout<<"div: "<<std::endl<<div<<std::endl;
    // EigenSolver<MatrixXd> eigensolver(Ac_1); // WHy different from Matlab eig()????
    
    // std::cout<<"q45555555555555555555554444444442"<<std::endl;
        
    // // auto EV_3 = (eigensolver.eigenvalues());

    // VectorXcd EV = eigensolver.eigenvalues();

    // MatrixXcd EVec = eigensolver.eigenvectors();

    // std::cout<<"q666666666666666666666666664442"<<std::endl;


    // std::cout<<" EV"<<EV<<std::endl;

    // std::cout<<" EVec"<<EVec<<std::endl;
        
    // Matrix<double, 6, 1> Ev_2;
    // // Ev_2 << sqrt((EV.real()).row(0).value()),sqrt((EV.real()).row(1).value()),sqrt((EV.real()).row(2).value());
    // double a = sqrt((EV.real()).row(0).value());
    // double b = sqrt((EV.real()).row(1).value());
    // double c = sqrt((EV.real()).row(2).value());
    // double d = sqrt((EV.real()).row(3).value());
    // double e = sqrt((EV.real()).row(4).value());
    // double f = sqrt((EV.real()).row(5).value());
    
    // std::cout<<"q77777777777777777777777777777777777777777777442"<<std::endl;
        
    // Ev_2 << a, b, c, d, e, f;


    // std::cout<<" Ev_2 "<<Ev_2 <<std::endl;

    // // Ev_2 << sqrt(EV.row(0).value()),sqrt(EV.row(1).value()),sqrt(EV.row(2).value());
    // MatrixXcd Ev_m = Ev_2.asDiagonal();

    MatrixXd pro = (M_new * M_old.completeOrthogonalDecomposition().pseudoInverse());


    // Ac = pro.cwiseSqrt();

    EigenSolver<MatrixXd> eigensolver(pro); 

    auto EV = (eigensolver.eigenvalues());// / M_1_det * M_2_det;

    VectorXcd EV_sqrt = (EV.array()).sqrt();
    // VectorXcd EV_sqrt = (EV.array()).sqrt();

    MatrixXcd M_sqrt = EV_sqrt.asDiagonal();

    MatrixXcd EVec = eigensolver.eigenvectors();
    MatrixXcd M_pow = EVec.inverse();

    Ac = EVec * M_sqrt * M_pow;


    MatrixXd Ac_real = Ac.real();
    
    // std::cout<<"8888888888888888888877777777442"<<std::endl;
    
    // // M_diff = logmap(Me_d, Me_ct); // 6x6

    T_tilde = Ac_real * T.real() * Ac_real.transpose();
    T_tilde_real = T_tilde.real();
    return T_tilde_real;
}