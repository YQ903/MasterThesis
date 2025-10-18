// manipulability_sampling.cpp
#include <iostream>
#include <random>
#include <limits>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include "../include/FrankaRobot.h"
#include "../include/geomJac.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

int main() {
    using namespace DQ_robotics;

    constexpr int DOF = 7;
    constexpr int N   = 20000;
    constexpr double EPS = 1e-12;

    FrankaRobot fr_builder;
    auto robot = fr_builder.kinematics();

    // Franka joint limits
    const VectorXd q_max = ((VectorXd(7) <<
        2.8973,  1.7628,  2.8973, -0.0698,  2.8973,  3.7525,  2.8973
    ).finished());
    const VectorXd q_min = ((VectorXd(7) <<
       -2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973
    ).finished());

    std::mt19937 rng(123456);
    auto sample_q = [&](VectorXd& q){
        q.resize(DOF);
        for (int i = 0; i < DOF; ++i) {
            std::uniform_real_distribution<double> dist(q_min(i), q_max(i));
            q(i) = dist(rng);
        }
    };

    double best_scalar  = std::numeric_limits<double>::infinity();
    VectorXd best_q(DOF);
    VectorXd best_singulars(6);

    double worst_scalar = 0.0;
    VectorXd worst_q(DOF);

    // Added for unit-balanced evaluation
    double best_scalar_scaled = std::numeric_limits<double>::infinity();
    double L_opt_at_best = 1.0;
    VectorXd best_singulars_scaled(6);

    for (int iter = 0; iter < N; ++iter) {
        VectorXd q(DOF);
        sample_q(q);

        MatrixXd poseJac = robot.pose_jacobian(q);        // 8×7
        MatrixXd J = geomJac7(robot, poseJac, q, DOF);    // 6×7

        // Unscaled SVD (original)
        Eigen::JacobiSVD<MatrixXd> svd(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
        VectorXd s = svd.singularValues();                // 6×1

        double s_min = std::numeric_limits<double>::infinity();
        double s_max = 0.0;
        for (int k = 0; k < s.size(); ++k) {
            if (s(k) < s_min) s_min = s(k);
            if (s(k) > s_max) s_max = s(k);
        }
        const double scalar = s_max / std::max(s_min, EPS);

        if (scalar < best_scalar) {
            best_scalar   = scalar;
            best_q        = q;
            best_singulars = s;
        }
        if (scalar > worst_scalar) {
            worst_scalar = scalar;
            worst_q      = q;
        }

        // === Unit balancing with characteristic length L* ===
        Eigen::Matrix<double,3,7> Jw = J.topRows<3>();
        Eigen::Matrix<double,3,7> Jv = J.bottomRows<3>();

        double tr_w = (Jw * Jw.transpose()).trace();
        double tr_v = (Jv * Jv.transpose()).trace();
        double L_opt = (tr_w > 1e-16) ? std::sqrt(std::max(tr_v, 0.0) / tr_w) : 1.0;

        Eigen::Matrix<double,6,7> J_scaled = J;
        if (L_opt > 1e-16) {
            J_scaled.bottomRows<3>() = Jv / L_opt;
        }

        Eigen::JacobiSVD<MatrixXd> svd_scaled(J_scaled, Eigen::ComputeThinU | Eigen::ComputeThinV);
        VectorXd s_scaled = svd_scaled.singularValues();

        double smin_sc = std::numeric_limits<double>::infinity();
        double smax_sc = 0.0;
        for (int k = 0; k < s_scaled.size(); ++k) {
            if (s_scaled(k) < smin_sc) smin_sc = s_scaled(k);
            if (s_scaled(k) > smax_sc) smax_sc = s_scaled(k);
        }
        double scalar_scaled = smax_sc / std::max(smin_sc, EPS);

        if (scalar_scaled < best_scalar_scaled) {
            best_scalar_scaled   = scalar_scaled;
            L_opt_at_best        = L_opt;
            best_q               = q;              // keep same variable name
            best_singulars_scaled = s_scaled;
        }
    }

    cout << std::fixed << std::setprecision(6);
    cout << "=== Sampling finished (" << N << " configs) ===\n";
    cout << "[Ground-Truth Manipulability (min condition #, unscaled)] = " << best_scalar << "\n";
    cout << "[Isotropic-closest joint config (rad)]:\n";
    for (int j = 0; j < DOF; ++j) {
        cout << "  q" << (j+1) << " = " << best_q(j) << "\n";
    }
    cout << "[Singular values at best (unscaled, 6)] : ";
    for (int k = 0; k < best_singulars.size(); ++k) {
        cout << best_singulars(k) << (k+1<best_singulars.size() ? ", " : "\n");
    }
    cout << "\n(Optional) Worst condition # observed (unscaled) = " << worst_scalar << "\n";

    cout << "\n--- Unit-balanced result (recommended) ---\n";
    cout << "[Min scaled condition #] = " << best_scalar_scaled << "\n";
    cout << "[Characteristic length L* at best] = " << L_opt_at_best << " (meters)\n";
    //cout << "[Singular values at best (scaled J, 6)] : ";
    //for (int k = 0; k < best_singulars_scaled.size(); ++k) {
    //    cout << best_singulars_scaled(k) << (k+1<best_singulars_scaled.size() ? ", " : "\n");
    //}

    // CSV: unscaled best
    //{
    //    std::ofstream fout("manipulability_best.csv");
    //    fout << std::setprecision(16);
    //    fout << "scalar";
    //    for (int j=0;j<DOF;++j) fout << ",q" << (j+1);
    //    for (int k=0;k<6;++k)  fout << ",s" << (k+1);
    //    fout << "\n";

    //    fout << best_scalar;
    //    for (int j=0;j<DOF;++j) fout << "," << best_q(j);
    //    for (int k=0;k<6;++k)   fout << "," << best_singulars(k);
    //    fout << "\n";
    //}

    // CSV: scaled best
    {
        std::ofstream fout("manipulability_best_scaled.csv");
        fout << std::setprecision(16);
        fout << "scalar_scaled,L_opt";
        for (int j=0;j<DOF;++j) fout << ",q" << (j+1);
        for (int k=0;k<6;++k)  fout << ",s_scaled_" << (k+1);
        fout << "\n";

        fout << best_scalar_scaled << "," << L_opt_at_best;
        for (int j=0;j<DOF;++j) fout << "," << best_q(j);
        for (int k=0;k<6;++k)   fout << "," << best_singulars_scaled(k);
        fout << "\n";
    }

    return 0;
}
