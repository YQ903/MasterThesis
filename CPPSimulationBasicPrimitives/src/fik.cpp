#include <iostream>
#include <array> 
#include <vector>
#include <chrono>
//#include <omp.h>
#include "Eigen/Dense"
using namespace std;
using namespace std::chrono;

#define PI 3.14159265359
#define d1 0.333
#define d3 0.316
#define a4 0.0825
#define a5 0.0825
#define d5 0.384
#define a7 0.088
// dE =   + 0.1034
#define dE 0.2104
// b1 = sqrt(d3*d3 + a4*a4)
#define b1 0.3265918706887849
// b2 = sqrt(d5*d5 + a5*a5)
#define b2 0.39276233271534583
// beta1 = arctan(a4/d3)
#define beta1 0.25537561488738186
// beta2 = arctan(a5/d5)
#define beta2 0.21162680876562978

// toletance for entering in singularity mode
# define SING_TOL 1e-3

// error threshold for swivel angle solver
# define ERR_THRESH 0.01 // this slightly smaller than 1deg
// max number of points in discretisation for swivel angle solver
const unsigned int MAX_N_POINTS = 1000;

// first row ROE:
// ROE[0], ROE[1], ROE[2]
// ROE[3], ROE[4], ROE[5]
// ROE[6], ROE[7], ROE[8]

const array<double,7> q_low = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
const array<double,7> q_up = {2.8973, 1.762, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
const array<double,7> q_mid = {0.0, 0.0, 0.0, -1.5708, 0.0, 1.8675, 0.0};

// Jacobian at home configuration (orientation only)
//const Eigen::Matrix<double, 3, 7> J0_S({{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//                                        {0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0},
//                                        {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0}});
const Eigen::Matrix<double,3,7> J0_S = (Eigen::Matrix<double,3,7>() << 
    0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,   // row 0
    0.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0,   // row 1
    1.0, 0.0, 1.0,  0.0, 1.0,  0.0, -1.0   // row 2
).finished();


Eigen::Matrix3d tmp_R;
Eigen::Matrix<double, 3, 7> tmp_J;
Eigen::Matrix<double, 6, 7> tmp_J_6d;
Eigen::Matrix<double, 3, 7> J_old;
Eigen::Matrix<double, 3, 4> J_old_low;
Eigen::Vector3d s;

void R_axis_angle(const Eigen::Vector3d& s, double theta){
    double x = s[0];
    double y = s[1];
    double z = s[2];
    double ct = cos(theta);
    double st = sin(theta);
    double one_minus_ct = 1-ct;
    tmp_R << ct + x * x * one_minus_ct, x * y * one_minus_ct - z * st, x * z * one_minus_ct + y * st,
         y * x * one_minus_ct + z * st, ct + y * y * one_minus_ct, y * z * one_minus_ct - x * st,
         z * x * one_minus_ct - y * st, z * y * one_minus_ct + x * st, ct + z * z * one_minus_ct;
}

void R_axis_angle(const array<double,3>& s, double theta){
    double x = s[0];
    double y = s[1];
    double z = s[2];
    double ct = cos(theta);
    double st = sin(theta);
    double one_minus_ct = 1-ct;
    tmp_R << ct + x * x * one_minus_ct, x * y * one_minus_ct - z * st, x * z * one_minus_ct + y * st,
         y * x * one_minus_ct + z * st, ct + y * y * one_minus_ct, y * z * one_minus_ct - x * st,
         z * x * one_minus_ct - y * st, z * y * one_minus_ct + x * st, ct + z * z * one_minus_ct;
}

array<double,3> Cross(const Eigen::Vector3d& u, const Eigen::Vector3d& v){
    return array<double,3> {u[1]*v[2]-v[1]*u[2], v[0]*u[2]-u[0]*v[2], u[0]*v[1]-v[0]*u[1]};
}

array<double,3> Cross(const array<double,3>& u, const Eigen::Vector3d& v){
    return array<double,3> {u[1]*v[2]-v[1]*u[2], v[0]*u[2]-u[0]*v[2], u[0]*v[1]-v[0]*u[1]};
}

array<double,3> Cross(const Eigen::Vector3d& u, const array<double,3>& v){
    return array<double,3> {u[1]*v[2]-v[1]*u[2], v[0]*u[2]-u[0]*v[2], u[0]*v[1]-v[0]*u[1]};
}

array<double,3> Cross(const array<double,3>& u, const array<double,3>& v){
    return array<double,3> {u[1]*v[2]-v[1]*u[2], v[0]*u[2]-u[0]*v[2], u[0]*v[1]-v[0]*u[1]};
}

void Cross_(const array<double,3>& u, const Eigen::Vector3d& v, array<double,3>& w){
    w[0] = u[1]*v[2]-v[1]*u[2]; 
    w[1] = v[0]*u[2]-u[0]*v[2]; 
    w[2] = u[0]*v[1]-v[0]*u[1];
}

void Cross_(const array<double, 3>& u, const array<double, 3>& v, array<double, 3>& w) {
    w[0] = u[1] * v[2] - v[1] * u[2];
    w[1] = v[0] * u[2] - u[0] * v[2];
    w[2] = u[0] * v[1] - v[0] * u[1];
}

double Dot(const Eigen::Vector3d& u, const Eigen::Vector3d& v){
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

double Dot(const array<double,3>& u, const array<double,3>& v){
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

double Dot(const array<double,3>& u, const Eigen::Vector3d& v){
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

double Dot(const Eigen::Vector3d& u, const array<double,3>& v){
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

double Norm(const Eigen::Vector3d& u){
    return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

double Norm(const array<double,3>& u){
    return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

Eigen::Matrix3d R_basis(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k){
    Eigen::Matrix3d R;
    R << i[0], j[0], k[0],
         i[1], j[1], k[1],
         i[2], j[2], k[2];
    return R;
}

void J_dir(const array<double,3>& s2, const array<double,3>& s3, const array<double,3>& s4, const array<double,3>& s5, const array<double,3>& s6, const array<double,3>& s7){
    tmp_J << 0, s2[0], s3[0], s4[0], s5[0], s6[0], s7[0],
             0, s2[1], s3[1], s4[1], s5[1], s6[1], s7[1],
             1, s2[2], s3[2], s4[2], s5[2], s6[2], s7[2];
}



// ------------------------------------------------------------
/*
J_6d(const array<double,3>& s2, const array<double,3>& s3, const array<double,3>& s4, const array<double,3>& s5, const array<double,3>& s6, const array<double,3>& s7, const array<double,3>& r4, const array<double,3>& r5, const array<double,3>& rE_W){
    tmp_J_6d.block<3, 7>(0, 0) << 0, s2[0], s3[0], s4[0], s5[0], s6[0], s7[0],
                                  0, s2[1], s3[1], s4[1], s5[1], s6[1], s7[1],
                                  1, s2[2], s3[2], s4[2], s5[2], s6[2], s7[2];
    // position vectors from E
    array<double,3> r1_E = {-rE_W[0], -rE_W[1], d1 - rE_W[2]};
    array<double,3> r2_E = r1_E;
    array<double,3> r3_E = r1_E;
    array<double,3> r4_E = {r4[0] - rE_W[0], r4[1] - rE_W[1], d1 + r4[2] - rE_W[2]};
    array<double,3> r5_E = {r5[0] - rE_W[0], r5[1] - rE_W[1], d1 + r5[2] - rE_W[2]};
    array<double,3> m_E;
    tmp_J_6d.block<3, 1>(3, 0) << r1_E[1], -r1_E[0], 0; // r1_E x (0,0,1)
    Cross_(r2_E, s2, m_E);
    tmp_J_6d.block<3, 1>(3, 1) << m_E[0], m_E[1], m_E[2];
    Cross_(r3_E, s3, m_E);
    tmp_J_6d.block<3, 1>(3, 2) << m_E[0], m_E[1], m_E[2];
    Cross_(r4_E, s4, m_E);
    tmp_J_6d.block<3, 1>(3, 3) << m_E[0], m_E[1], m_E[2];
    Cross_(r5_E, s5, m_E);
    tmp_J_6d.block<3, 1>(3, 4) << m_E[0], m_E[1], m_E[2];
    Cross(r5_E, s6, m_E); // r6 = r5
    tmp_J_6d.block<3, 1>(3, 5) << m_E[0], m_E[1], m_E[2];
    tmp_J_6d.block<3, 1>(3, 6) << 0, 0, 0;
}
*/

void save_J_sol(const array<double,3>& s2,
                const array<double,3>& s3,
                const array<double,3>& s4,
                const array<double,3>& s5,
                const array<double,3>& s6,
                const array<double,3>& s7,
                const array<double,3>& r4,
                const array<double,3>& r5,
                const array<double,3>& rE_W,
                vector<array<array<double, 6>, 7>>& Jsols,
                const int index,
                const bool Jacobian_ee_is_8) {
    array<double,3> r1_E = {-rE_W[0], -rE_W[1], d1 - rE_W[2]};
    array<double,3> r4_E = {r4[0] - rE_W[0], r4[1] - rE_W[1], d1 + r4[2] - rE_W[2]};
    array<double,3> r5_E = {r5[0] - rE_W[0], r5[1] - rE_W[1], d1 + r5[2] - rE_W[2]};
    if(Jacobian_ee_is_8){
        // r_PF_O = r_PE_O + r_EF_O
        //        = r_PE_O - rEW_O + rOW_O + r_EF_O
        //        = r_PO_O - rEW_O + (0,0,d1) + 0.1034*s_7_O
        r1_E = {r1_E[0] + 0.1034*s7[0], r1_E[1] + 0.1034*s7[1], r1_E[2] + 0.1034*s7[2]};
        r4_E = {r4_E[0] + 0.1034*s7[0], r4_E[1] + 0.1034*s7[1], r4_E[2] + 0.1034*s7[2]};
        r5_E = {r5_E[0] + 0.1034*s7[0], r5_E[1] + 0.1034*s7[1], r5_E[2] + 0.1034*s7[2]};
    }
    array<double,3> m_E;
    Jsols[2*index][0] = {0, 0, 1, r1_E[1], -r1_E[0], 0}; // r1_E x (0,0,1)
    Jsols[2*index+1][0] = {0, 0, 1, r1_E[1], -r1_E[0], 0};
    Cross_(r1_E, s2, m_E); // r2_E = r1_E
    Jsols[2*index][1] = {s2[0], s2[1], s2[2], m_E[0], m_E[1], m_E[2]};
    Jsols[2*index+1][1] = {-s2[0], -s2[1], -s2[2], -m_E[0], -m_E[1], -m_E[2]}; // due to wrist
    Cross_(r1_E, s3, m_E); //  r3_E = r1_E
    Jsols[2*index][2] = {s3[0], s3[1], s3[2], m_E[0], m_E[1], m_E[2]};
    Jsols[2*index+1][2] = {s3[0], s3[1], s3[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r4_E, s4, m_E);
    Jsols[2*index][3] = {s4[0], s4[1], s4[2], m_E[0], m_E[1], m_E[2]};
    Jsols[2*index+1][3] = {s4[0], s4[1], s4[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r5_E, s5, m_E);
    Jsols[2*index][4] = {s5[0], s5[1], s5[2], m_E[0], m_E[1], m_E[2]};
    Jsols[2*index+1][4] = {s5[0], s5[1], s5[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r5_E, s6, m_E); // r6 = r5
    Jsols[2*index][5] = {s6[0], s6[1], s6[2], m_E[0], m_E[1], m_E[2]};
    Jsols[2*index+1][5] = {s6[0], s6[1], s6[2], m_E[0], m_E[1], m_E[2]};
    Jsols[2*index][6] = {s7[0], s7[1], s7[2], 0, 0, 0};
    Jsols[2*index+1][6] = {s7[0], s7[1], s7[2], 0, 0, 0};
}

void save_J_sol(const array<double, 3>& s2,
    const array<double, 3>& s3,
    const array<double, 3>& s4,
    const array<double, 3>& s5,
    const array<double, 3>& s6,
    const array<double, 3>& s7,
    const array<double, 3>& r4,
    const array<double, 3>& r5,
    const array<double, 3>& rE_W,
    array<array<array<double, 6>, 7>, 8>& Jsols,
    const int index,
    const bool Jacobian_ee_is_8) {
    array<double, 3> r1_E = { -rE_W[0], -rE_W[1], d1 - rE_W[2] };
    array<double, 3> r4_E = { r4[0] - rE_W[0], r4[1] - rE_W[1], d1 + r4[2] - rE_W[2] };
    array<double, 3> r5_E = { r5[0] - rE_W[0], r5[1] - rE_W[1], d1 + r5[2] - rE_W[2] };
    if (Jacobian_ee_is_8) {
        // r_PF_O = r_PE_O + r_EF_O
        //        = r_PE_O - rEW_O + rOW_O + r_EF_O
        //        = r_PO_O - rEW_O + (0,0,d1) + 0.1034*s_7_O
        r1_E = { r1_E[0] + 0.1034 * s7[0], r1_E[1] + 0.1034 * s7[1], r1_E[2] + 0.1034 * s7[2] };
        r4_E = { r4_E[0] + 0.1034 * s7[0], r4_E[1] + 0.1034 * s7[1], r4_E[2] + 0.1034 * s7[2] };
        r5_E = { r5_E[0] + 0.1034 * s7[0], r5_E[1] + 0.1034 * s7[1], r5_E[2] + 0.1034 * s7[2] };
    }
    array<double, 3> m_E;
    Jsols[2 * index][0] = { 0, 0, 1, r1_E[1], -r1_E[0], 0 }; // r1_E x (0,0,1)
    Jsols[2 * index + 1][0] = { 0, 0, 1, r1_E[1], -r1_E[0], 0 };
    Cross_(r1_E, s2, m_E); // r2_E = r1_E
    Jsols[2 * index][1] = { s2[0], s2[1], s2[2], m_E[0], m_E[1], m_E[2] };
    Jsols[2 * index + 1][1] = { -s2[0], -s2[1], -s2[2], -m_E[0], -m_E[1], -m_E[2] }; // due to wrist
    Cross_(r1_E, s3, m_E); //  r3_E = r1_E
    Jsols[2 * index][2] = { s3[0], s3[1], s3[2], m_E[0], m_E[1], m_E[2] };
    Jsols[2 * index + 1][2] = { s3[0], s3[1], s3[2], m_E[0], m_E[1], m_E[2] };
    Cross_(r4_E, s4, m_E);
    Jsols[2 * index][3] = { s4[0], s4[1], s4[2], m_E[0], m_E[1], m_E[2] };
    Jsols[2 * index + 1][3] = { s4[0], s4[1], s4[2], m_E[0], m_E[1], m_E[2] };
    Cross_(r5_E, s5, m_E);
    Jsols[2 * index][4] = { s5[0], s5[1], s5[2], m_E[0], m_E[1], m_E[2] };
    Jsols[2 * index + 1][4] = { s5[0], s5[1], s5[2], m_E[0], m_E[1], m_E[2] };
    Cross_(r5_E, s6, m_E); // r6 = r5
    Jsols[2 * index][5] = { s6[0], s6[1], s6[2], m_E[0], m_E[1], m_E[2] };
    Jsols[2 * index + 1][5] = { s6[0], s6[1], s6[2], m_E[0], m_E[1], m_E[2] };
    Jsols[2 * index][6] = { s7[0], s7[1], s7[2], 0, 0, 0 };
    Jsols[2 * index + 1][6] = { s7[0], s7[1], s7[2], 0, 0, 0 };
}




void save_J_sol(const array<double, 3>& s2,
                const array<double, 3>& s3,
                const array<double, 3>& s4,
                const array<double, 3>& s5,
                const array<double, 3>& s6,
                const array<double, 3>& s7,
                const array<double, 3>& r4,
                const array<double, 3>& r5,
                const array<double, 3>& rE_W,
                array<array<array<double, 6>, 7>, 8>& Jsols,
                const int index,
                const char Jacobian_ee) {
    array<double, 3> r1_ee, r4_ee, r5_ee;
    if (Jacobian_ee == '6') {
        // r_P6_O = r_PO_O + r_O6_O
        //          r_PO_O - r_6O_O remember r_6O_O = r_5O_O
        r1_ee = { -r5[0], -r5[1], -r5[2] };
        r4_ee = { r4[0] - r5[0], r4[1] - r5[1] , r4[2] - r5[2] };
        r5_ee = { 0, 0 , 0 };
    } else if (Jacobian_ee == '8' || Jacobian_ee == 'F') {
        // r_PF_O = r_PE_O + r_EF_O
        //        = r_PE_O - rEW_O + rOW_O + r_EF_O
        //        = r_PO_O - rEW_O + (0,0,d1) + 0.1034*s_7_O
        r1_ee = { -rE_W[0] + 0.1034 * s7[0], -rE_W[1] + 0.1034 * s7[1], d1 - rE_W[2] + 0.1034 * s7[2] };
        r4_ee = { r4[0] - rE_W[0] + 0.1034 * s7[0], r4[1] - rE_W[1] + 0.1034 * s7[1], d1 + r4[2] - rE_W[2] + 0.1034 * s7[2] };
        r5_ee = { r5[0] - rE_W[0] + 0.1034 * s7[0], r5[1] - rE_W[1] + 0.1034 * s7[1], d1 + r5[2] - rE_W[2] + 0.1034 * s7[2] };
    } else {
        // r_PE = r_PO_O + r_OE_O
        //      = r_PO_O - r_EO_O
        //      = r_PO_O - (r_EW_O + r_WO_O) = r_PO_O - r_EW_O + r_OW_O
        r1_ee = { -rE_W[0], -rE_W[1], d1 - rE_W[2] };
        r4_ee = { r4[0] - rE_W[0], r4[1] - rE_W[1], d1 + r4[2] - rE_W[2] };
        r5_ee = { r5[0] - rE_W[0], r5[1] - rE_W[1], d1 + r5[2] - rE_W[2] };
    }

    array<double, 3> m_ee;
    Jsols[2 * index][0] = { 0, 0, 1, r1_ee[1], -r1_ee[0], 0 }; // r1_ee x (0,0,1)
    Jsols[2 * index + 1][0] = { 0, 0, 1, r1_ee[1], -r1_ee[0], 0 };
    Cross_(r1_ee, s2, m_ee); // r2_ee = r1_ee
    Jsols[2 * index][1] = { s2[0], s2[1], s2[2], m_ee[0], m_ee[1], m_ee[2] };
    Jsols[2 * index + 1][1] = { -s2[0], -s2[1], -s2[2], -m_ee[0], -m_ee[1], -m_ee[2] }; // due to wrist
    Cross_(r1_ee, s3, m_ee); //  r3_ee = r1_ee
    Jsols[2 * index][2] = { s3[0], s3[1], s3[2], m_ee[0], m_ee[1], m_ee[2] };
    Jsols[2 * index + 1][2] = { s3[0], s3[1], s3[2], m_ee[0], m_ee[1], m_ee[2] };
    Cross_(r4_ee, s4, m_ee);
    Jsols[2 * index][3] = { s4[0], s4[1], s4[2], m_ee[0], m_ee[1], m_ee[2] };
    Jsols[2 * index + 1][3] = { s4[0], s4[1], s4[2], m_ee[0], m_ee[1], m_ee[2] };
    Cross_(r5_ee, s5, m_ee);
    Jsols[2 * index][4] = { s5[0], s5[1], s5[2], m_ee[0], m_ee[1], m_ee[2] };
    Jsols[2 * index + 1][4] = { s5[0], s5[1], s5[2], m_ee[0], m_ee[1], m_ee[2] };
    Cross_(r5_ee, s6, m_ee); // r6 = r5
    Jsols[2 * index][5] = { s6[0], s6[1], s6[2], m_ee[0], m_ee[1], m_ee[2] };
    Jsols[2 * index + 1][5] = { s6[0], s6[1], s6[2], m_ee[0], m_ee[1], m_ee[2] };
    if (Jacobian_ee == '6') {
        Jsols[2 * index][6] = { 0, 0, 0, 0, 0, 0 };
        Jsols[2 * index + 1][6] = { 0, 0, 0, 0, 0, 0 };
    }
    else {
        Jsols[2 * index][6] = { s7[0], s7[1], s7[2], 0, 0, 0 };
        Jsols[2 * index + 1][6] = { s7[0], s7[1], s7[2], 0, 0, 0 };
    }
}






void pushback_J_sol(const array<double,3>& s2,
                    const array<double,3>& s3,
                    const array<double,3>& s4,
                    const array<double,3>& s5,
                    const array<double,3>& s6,
                    const array<double,3>& s7,
                    const array<double,3>& r4,
                    const array<double,3>& r5,
                    const array<double,3>& rE_W,
                    vector<array<array<double, 6>, 7>>& Jsols,
                    const bool Jacobian_ee_is_8) {
    array<array<double,6>,7> J;
    array<double,3> r1_E = {-rE_W[0], -rE_W[1], d1 - rE_W[2]};
    array<double,3> r4_E = {r4[0] - rE_W[0], r4[1] - rE_W[1], d1 + r4[2] - rE_W[2]};
    array<double,3> r5_E = {r5[0] - rE_W[0], r5[1] - rE_W[1], d1 + r5[2] - rE_W[2]};
    if(Jacobian_ee_is_8){
        // r_PF_O = r_PE_O + r_EF_O
        //        = r_PE_O - rEW_O + rOW_O + r_EF_O
        //        = r_PO_O - rEW_O + (0,0,d1) + 0.1034*s_7_O
        r1_E = {r1_E[0] + 0.1034*s7[0], r1_E[1] + 0.1034*s7[1], r1_E[2] + 0.1034*s7[2]};
        r4_E = {r4_E[0] + 0.1034*s7[0], r4_E[1] + 0.1034*s7[1], r4_E[2] + 0.1034*s7[2]};
        r5_E = {r5_E[0] + 0.1034*s7[0], r5_E[1] + 0.1034*s7[1], r5_E[2] + 0.1034*s7[2]};
    }
    array<double,3> m_E;
    J[0] = {0, 0, 1, r1_E[1], -r1_E[0], 0}; // r1_E x (0,0,1)
    Cross_(r1_E, s2, m_E); // r2_E = r1_E
    J[1] = {s2[0], s2[1], s2[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r1_E, s3, m_E); //  r3_E = r1_E
    J[2] = {s3[0], s3[1], s3[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r4_E, s4, m_E);
    J[3] = {s4[0], s4[1], s4[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r5_E, s5, m_E);
    J[4] = {s5[0], s5[1], s5[2], m_E[0], m_E[1], m_E[2]};
    Cross_(r5_E, s6, m_E); // r6 = r5
    J[5] = {s6[0], s6[1], s6[2], m_E[0], m_E[1], m_E[2]};
    J[6] = {s7[0], s7[1], s7[2], 0, 0, 0};
    Jsols.push_back(J);
    J[1] = {-s2[0], -s2[1], -s2[2], -J[1][3], -J[1][4], -J[1][5]}; // due to wrist
    Jsols.push_back(J);
}

// ------------------------------------------------------------


double signed_angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& s){
    return atan2(Dot(Cross(v1, v2), s), Dot(v1, v2));
}

double signed_angle(const array<double,3>& v1, const array<double,3>& v2, const array<double,3>& s){
    //return atan2(s[2]*(v1[0]*v2[1] - v1[1]*v2[0]) - s[1]*(v1[0]*v2[2] - v1[2]*v2[0]) + s[0]*(v1[1]*v2[2] - v1[2]*v2[1]), v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    return atan2(Dot(Cross(v1, v2), s), Dot(v1, v2));
}

double signed_d_circ(double theta1, double theta2){
    return atan2(sin(theta2-theta1), cos(theta2-theta1));
}

void check_limits(array<double,7>& q, int n){
    for(int i; i<n; i++){
        q[i] = q_mid[i] + atan2(sin(q[i]-q_mid[i]), cos(q[i]-q_mid[i]));
        if(q[i] < q_low[i] || q[i] > q_up[i]) q[i] = NAN;
    }
}

// whole rotational part
array<double,6> q_from_J(const Eigen::Matrix<double, 3, 7>& J){
    const int dof = 7;
    array<double, dof-1> q;
    J_old = J0_S; // [3,dof]
    for(int i=0; i < dof-1; i++){
        s = J.col(i);
        q[i] = signed_angle(J_old.col(i+1), J.col(i+1), s);
        if(i == dof-2) break;
        R_axis_angle(s, q[i]);
        J_old.block(0, i+1, 3, dof-(i+1)) = tmp_R * J_old.block(0,i+1,3,dof-(i+1));
    }
    return q;
}

// whole Jacobian and q
array<double, 7> J_to_q(const array<array<double,6>,7>& J, const array<array<double, 3>, 3> R, const char ee='E') {
    // J is the transpose of the Jacobian
    // R is the rotation matrix of frame ee
    // ee must be a frame attached to the gripper: "E", "F" or "8".
    Eigen::Matrix<double, 3, 7> Jrot;
    Jrot << J[0][0], J[1][0], J[2][0], J[3][0], J[4][0], J[5][0], J[5][0],
            J[0][1], J[1][1], J[2][1], J[3][1], J[4][1], J[5][1], J[5][1],
            J[0][2], J[1][2], J[2][2], J[3][2], J[4][2], J[5][2], J[5][2];
    const int dof = 7;
    array<double, dof> q;
    array<double, 3> i7, ie, s6, s7;
    s6 = { J[5][0] , J[5][1] , J[5][2] };
    s7 = { J[6][0] , J[6][1] , J[6][2] };
    Cross_(s6, s7, i7);
    ie = { R[0][0], R[1][0], R[2][0] };
    q[6] = signed_angle(i7, ie, s7) + (ee=='E' ? -PI / 4 : 0);
    J_old = J0_S; // [3,dof]
    for (int i = 0; i < dof - 1; i++) {
        s = Jrot.col(i);
        q[i] = signed_angle(J_old.col(i + 1), Jrot.col(i + 1), s);
        if (i == dof - 2) break;
        R_axis_angle(s, q[i]);
        J_old.block(0, i + 1, 3, dof - (i + 1)) = tmp_R * J_old.block(0, i + 1, 3, dof - (i + 1));
    }
    return q;
}

// lower matrix
array<double,3> q_from_low_J(const Eigen::Matrix<double, 3, 7>& J){
    const int dof = 4;
    array<double, dof-1> q;
    J_old_low = J0_S.block(0, 0, 3, dof); // [3,dof]
    for(int i=0; i < dof-1; i++){
        s = J.col(i);
        q[i] = signed_angle(J_old_low.col(i+1), J.col(i+1), s);
        if(i == dof-2) break;
        R_axis_angle(s, q[i]);
        J_old_low.block(0, i+1, 3, dof-(i+1)) =  tmp_R * J_old_low.block(0,i+1,3,dof-(i+1));
    }
    return q;
}

Eigen::Matrix3d R_rpy(const double r, const double p, const double y){
    Eigen::Matrix3d R;
    R << cos(p)*cos(y), cos(y)*sin(p)*sin(r) - cos(r)*sin(y), sin(r)*sin(y) + cos(r)*cos(y)*sin(p),
         cos(p)*sin(y), cos(r)*cos(y) + sin(p)*sin(r)*sin(y), cos(r)*sin(p)*sin(y) - cos(y)*sin(r),
         -sin(p), cos(p)*sin(r),  cos(p)*cos(r);
    return R;
}

Eigen::Matrix3d R_z(const double theta){
    Eigen::Matrix3d R;
    R << cos(theta), -sin(theta), 0,
         sin(theta), cos(theta), 0,
         0, 0, 1;
    return R;
}

Eigen::Matrix4d T_from_R_and_r(const Eigen::Matrix3d& R, const double x, const double y, const double z){
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block(0,0,3,3) = R;
    T(0,3) = x;
    T(1,3) = y;
    T(2,3) = z;
    return T;
}

void column_1s_times_vec(const array<double, 9>& R, const array<double,3>& v, array<double,3>& res){
    res[0] = R[0]*v[0] + R[1]*v[1] + R[2]*v[2]; 
    res[1] = R[3]*v[0] + R[4]*v[1] + R[5]*v[2]; 
    res[2] = R[6]*v[0] + R[7]*v[1] + R[8]*v[2];
}

void rotate_by_axis_angle(const array<double, 3>& s, const double theta, const array<double,3>& v, array<double,3>& res){
    R_axis_angle(s, theta);
    res[0] = tmp_R(0,0)*v[0] + tmp_R(0,1)*v[1] + tmp_R(0,2)*v[2];
    res[1] = tmp_R(1,0)*v[0] + tmp_R(1,1)*v[1] + tmp_R(1,2)*v[2];
    res[2] = tmp_R(2,0)*v[0] + tmp_R(2,1)*v[1] + tmp_R(2,2)*v[2];
}

Eigen::Matrix4d T_rpy(const double r, const double p, const double y, const double px, const double py, const double pz){
    Eigen::Matrix4d T;
    T << cos(p)*cos(y), cos(y)*sin(p)*sin(r) - cos(r)*sin(y), sin(r)*sin(y) + cos(r)*cos(y)*sin(p), px,
         cos(p)*sin(y), cos(r)*cos(y) + sin(p)*sin(r)*sin(y), cos(r)*sin(p)*sin(y) - cos(y)*sin(r), py,
         -sin(p), cos(p)*sin(r),  cos(p)*cos(r), pz,
         0, 0, 0, 1;
    return T;
}

Eigen::Matrix4d T_rot_z(const double theta, const double px = 0.0, const double py = 0.0, const double pz = 0.0){
    Eigen::Matrix4d T;
    T << cos(theta), -sin(theta), 0, px,
         sin(theta), cos(theta), 0, py,
         0, 0, 1, pz,
         0, 0, 0, 1;
    return T;
}

Eigen::Matrix<double, 6, 6> Adj_trans(const double x, const double y, const double z){
    Eigen::Matrix<double, 6, 6> Adj = Eigen::Matrix<double, 6, 6>::Identity();
    Adj.block<3, 3>(3, 0) << 0, -z, y,
                             z, 0, -x,
                             -y, x, 0;
    return Adj;
}

// Eigen::Matrix4d franka_fk5(const array<double,7>& q, const unsigned int ee=6){
//     // this function is not optimised for speed
//     array<Eigen::Matrix4d, 9> Ti;
//     // T01
//     Ti[0] = T_from_R_and_r(R_rpy(0, 0, 0), 0, 0, 0.333) * T_from_R_and_r(R_z(q[0]), 0, 0, 0);
//     // T12
//     Ti[1] = T_from_R_and_r(R_rpy(-PI/2, 0, 0), 0, 0, 0) * T_from_R_and_r(R_z(q[1]), 0, 0, 0);
//     // T23
//     Ti[2] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0, -0.316, 0) * T_from_R_and_r(R_z(q[2]), 0, 0, 0);
//     // T34
//     Ti[3] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0.0825, 0, 0) * T_from_R_and_r(R_z(q[3]), 0, 0, 0);
//     // T45
//     Ti[4] = T_from_R_and_r(R_rpy(-PI/2, 0, 0), -0.0825, 0.384, 0) * T_from_R_and_r(R_z(q[4]), 0, 0, 0);
//     // T56
//     Ti[5] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0, 0, 0) * T_from_R_and_r(R_z(q[5]), 0, 0, 0);
//     // T67
//     Ti[6] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0.088, 0, 0) * T_from_R_and_r(R_z(q[6]), 0, 0, 0);
//     // T78
//     Ti[7] = T_from_R_and_r(R_rpy(0, 0, 0), 0, 0, 0.107);
//     // T8E
//     Ti[8] = T_from_R_and_r(R_z(-PI/4), 0, 0, 0.1034);
//     Eigen::Matrix4d T0ee = Ti[0];
//     for(int i=1; i<ee; i++)
//         T0ee = T0ee*Ti[i];
//     return T0ee;
// }

Eigen::Matrix4d franka_fk(const array<double,7>& q, const unsigned int ee=9){
    // this function is not optimised for speed
    array<Eigen::Matrix4d, 9> Ti;
    // T01
    Ti[0] = T_from_R_and_r(R_rpy(0, 0, 0), 0, 0, 0.333) * T_from_R_and_r(R_z(q[0]), 0, 0, 0);
    // T12
    Ti[1] = T_from_R_and_r(R_rpy(-PI/2, 0, 0), 0, 0, 0) * T_from_R_and_r(R_z(q[1]), 0, 0, 0);
    // T23
    Ti[2] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0, -0.316, 0) * T_from_R_and_r(R_z(q[2]), 0, 0, 0);
    // T34
    Ti[3] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0.0825, 0, 0) * T_from_R_and_r(R_z(q[3]), 0, 0, 0);
    // T45
    Ti[4] = T_from_R_and_r(R_rpy(-PI/2, 0, 0), -0.0825, 0.384, 0) * T_from_R_and_r(R_z(q[4]), 0, 0, 0);
    // T56
    Ti[5] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0, 0, 0) * T_from_R_and_r(R_z(q[5]), 0, 0, 0);
    // T67
    Ti[6] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0.088, 0, 0) * T_from_R_and_r(R_z(q[6]), 0, 0, 0);
    // T78
    Ti[7] = T_from_R_and_r(R_rpy(0, 0, 0), 0, 0, 0.107);
    // T8E
    Ti[8] = T_from_R_and_r(R_z(-PI/4), 0, 0, 0.1034);
    Eigen::Matrix4d T0ee = Ti[0];
    for(int i=1; i<ee; i++)
        T0ee = T0ee*Ti[i];
    return T0ee;
}

array<array<double, 6>,7> partial_Jacobian(const array<double, 7>& q, const unsigned int ee = 9) {
    // returns J^T for a given vector of joint angles, q. The end-effector frame is ee \in {1,...,9}, where 9 is frame E, and 8 is frame F.
    // OUTPUT: J^T \in R^(7,6): array<array<double,6>,7>
    // INPUT: q \in R^7, array<double,7>
    //        ee 
    Eigen::Vector3d s, r, m;
    array<Eigen::Matrix4d, 9> Ti;
    unsigned int cols = ee >= 7 ? 7 : ee;
    Eigen::Matrix<double, 6, 7> J6d;
    Ti[0] = T_rpy(0, 0, 0, 0, 0, 0.333) * T_rot_z(q[0]); // T01
    Ti[1] = T_rpy(-PI / 2, 0, 0, 0, 0, 0) * T_rot_z(q[1]); // T12
    Ti[2] = T_rpy(PI / 2, 0, 0, 0, -0.316, 0) * T_rot_z(q[2]); // T23
    Ti[3] = T_rpy(PI / 2, 0, 0, 0.0825, 0, 0) * T_rot_z(q[3]); // T34
    Ti[4] = T_rpy(-PI / 2, 0, 0, -0.0825, 0.384, 0) * T_rot_z(q[4]); // T45
    Ti[5] = T_rpy(PI / 2, 0, 0, 0, 0, 0) * T_rot_z(q[5]); // T56
    Ti[6] = T_rpy(PI / 2, 0, 0, 0.088, 0, 0) * T_rot_z(q[6]); // T67
    Ti[7] = T_rpy(0, 0, 0, 0, 0, 0.107); // T78
    Ti[8] = T_rot_z(-PI / 4, 0, 0, 0.1034); // T8E
    Eigen::Matrix4d T = Ti[0]; // T01
    J6d.col(0) << 0, 0, 1, 0, 0, 0;
    for (int i = 1; i < ee; i++) {
        T = T * Ti[i]; // T0{i+1} = T0{i}*T{i}{i+1}
        if (i < 7) {
            s = T.block<3, 1>(0, 2);
            r = T.block<3, 1>(0, 3);
            m = r.cross(s);
            J6d.col(i) << s[0], s[1], s[2], m[0], m[1], m[2];
        }
    }
    J6d = Adj_trans(-T(0, 3), -T(1, 3), -T(2, 3)) * J6d;
    array<array<double, 6>,7> Jarr;
    for (int i = 0; i < cols; i++)
        Jarr[i] = { J6d(0,i), J6d(1,i), J6d(2,i), J6d(3,i), J6d(4,i), J6d(5,i) };
    for (int i = cols; i < 7; i++)
        Jarr[i] = { 0, 0, 0, 0, 0, 0 };
    return Jarr;
}

array<array<double,6>,7> Jacobian(const array<double,7>& q, const bool Jacobian_ee_is_8 = false){
    // returns J^T for a given vector of joint angles, q. The default end-effector point is E.
    // OUTPUT: J^T \in R^(7,6): array<array<double,6>,7>
    // INPUT: q \in R^7, array<double,7>
    //        Jacobian_ee_is_8, if true, the end-effector point is O_8 (origin of frames F and 8)
    Eigen::Vector3d s, r, m;
    array<Eigen::Matrix4d, 9> Ti;
    Eigen::Matrix<double, 6, 7> J6d;
    Ti[0] = T_rpy(0, 0, 0, 0, 0, 0.333) * T_rot_z(q[0]); // T01
    Ti[1] = T_rpy(-PI/2, 0, 0, 0, 0, 0) * T_rot_z(q[1]); // T12
    Ti[2] = T_rpy(PI/2, 0, 0, 0, -0.316, 0) * T_rot_z(q[2]); // T23
    Ti[3] = T_rpy(PI/2, 0, 0, 0.0825, 0, 0) * T_rot_z(q[3]); // T34
    Ti[4] = T_rpy(-PI/2, 0, 0, -0.0825, 0.384, 0) * T_rot_z(q[4]); // T45
    Ti[5] = T_rpy(PI/2, 0, 0, 0, 0, 0) * T_rot_z(q[5]); // T56
    Ti[6] = T_rpy(PI/2, 0, 0, 0.088, 0, 0) * T_rot_z(q[6]); // T67
    Ti[7] = T_rpy(0, 0, 0, 0, 0, 0.107); // T78
    Ti[8] = T_rot_z(-PI/4, 0, 0, 0.1034); // T8E
    Eigen::Matrix4d T = Ti[0]; // T01
    J6d.col(0) << 0,0,1,0,0,0;
    for(int i=1; i<7; i++){
        T = T*Ti[i]; // T0{i+1} = T0{i}*T{i}{i+1}
        s = T.block<3, 1>(0, 2);
        r = T.block<3, 1>(0, 3);
        m = r.cross(s);
        J6d.col(i) << s[0], s[1], s[2], m[0], m[1], m[2];
    }
    T = T*Ti[7]; //TO8
    if(!Jacobian_ee_is_8)
        T = T*Ti[8]; //TOE
    J6d = Adj_trans(-T(0,3), -T(1,3), -T(2,3))*J6d;
    array<array<double, 6>, 7> Jarr;
    for(int i=0; i<7; i++)
        Jarr[i] = {J6d(0,i), J6d(1,i), J6d(2,i), J6d(3,i), J6d(4,i), J6d(5,i)};
    return Jarr;
}



array<Eigen::Matrix4d,9> franka_fk_all_frames(const array<double,7>& q){
    array<Eigen::Matrix4d, 9> Ti, T0;
    // T01
    Ti[0] = T_from_R_and_r(R_rpy(0, 0, 0), 0, 0, 0.333) * T_from_R_and_r(R_z(q[0]), 0, 0, 0);
    // T12
    Ti[1] = T_from_R_and_r(R_rpy(-PI/2, 0, 0), 0, 0, 0) * T_from_R_and_r(R_z(q[1]), 0, 0, 0);
    // T23
    Ti[2] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0, -0.316, 0) * T_from_R_and_r(R_z(q[2]), 0, 0, 0);
    // T34
    Ti[3] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0.0825, 0, 0) * T_from_R_and_r(R_z(q[3]), 0, 0, 0);
    // T45
    Ti[4] = T_from_R_and_r(R_rpy(-PI/2, 0, 0), -0.0825, 0.384, 0) * T_from_R_and_r(R_z(q[4]), 0, 0, 0);
    // T56
    Ti[5] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0, 0, 0) * T_from_R_and_r(R_z(q[5]), 0, 0, 0);
    // T67
    Ti[6] = T_from_R_and_r(R_rpy(PI/2, 0, 0), 0.088, 0, 0) * T_from_R_and_r(R_z(q[6]), 0, 0, 0);
    // T78
    Ti[7] = T_from_R_and_r(R_rpy(0, 0, 0), 0, 0, 0.107);
    // T8E
    Ti[8] = T_from_R_and_r(R_z(-PI/4), 0, 0, 0.1034);
    T0[0] = Ti[0];
    for(int i=1; i<9; i++)
        T0[i] = T0[i-1]*Ti[i];
    return T0;
}

vector<array<double,7>> franka_ik_q7(const array<double, 3>& r, const array<double, 9>& ROE, const double q7, const double q1_sing = PI/2) {
    // ROE is in row-first
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    array<double,3> k_E_O = {ROE[2], ROE[5], ROE[8]};
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    Eigen::Vector3d i_6_O = tmp_R * i_E_O;
    array<double,3> s6 = Cross(k_E_O, i_6_O);
    array<double,3> r6 = {r[0] - dE*k_E_O[0] - a7*i_6_O[0], r[1] - dE*k_E_O[1] - a7*i_6_O[1], r[2] - d1 - dE*k_E_O[2] - a7*i_6_O[2]};
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    if(tmp>1){
        if((tmp - 1) * (tmp - 1) < SING_TOL){
            tmp = 1;
        } 
        else{
            cout << "ERROR: unable to assembly kinematic chain";
            return vector<array<double,7>>();
        }
    }
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double,3> k_C_O = {-r6[0]/l, -r6[1]/l, -r6[2]/l};
    array<double,3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
    array<double,3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double,3>,4> s5s;
    double sa2, ca2;
    int n_alphs = 1;
    unsigned int n_sols = 0;
    if (d3+d5 < l && l < b1+b2) n_alphs = 2;
    double v[3];
    for(int i=0; i<n_alphs; i++){
        sa2 = sin(alpha2);
        ca2 = cos(alpha2);
        tmp = -rz * ca2 / (ry * sa2);
        if (tmp*tmp> 1)
        {
            if (tmp*tmp - 1 < SING_TOL)
                tmp = (tmp > 0) ? 1 : -1;
            else
                continue;
        }
        tmp = asin(tmp);
        v[0] = -sa2*cos(tmp);
        v[1] = -sa2*sin(tmp);
        v[2] = -ca2;
        s5s[n_sols] = {i_C_O[0]*v[0] + j_C_O[0]*v[1] + k_C_O[0]*v[2],
                       i_C_O[1]*v[0] + j_C_O[1]*v[1] + k_C_O[1]*v[2],
                       i_C_O[2]*v[0] + j_C_O[2]*v[1] + k_C_O[2]*v[2]};
        tmp = 2*sa2*cos(tmp);
        //s5[n_sols+1] = s5s[n_sols] + (2*sa2*cos(tmp)*i_C_O);
        s5s[n_sols+1] = {s5s[n_sols][0] + tmp*i_C_O[0],
                        s5s[n_sols][1] + tmp*i_C_O[1],
                        s5s[n_sols][2] + tmp*i_C_O[2]};
        n_sols+=2;
        alpha2 = beta2 - actmp;
    }
    vector<array<double,7>> sols(2*n_sols);
    array<double,3> s4, r4, s3, s2, s5;
    array<double,6> sol1;
    array<double,3> sol2;
    for(int i=0; i<n_sols; i++){
        s5 = s5s[i];
        s4 = Cross(s5, r6);
        tmp = Norm(s4);
        s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
        //r4 = r6 - d5 * s5 + a5 * Cross(s5, s4);
        r4 = Cross(s5, s4);
        r4 = {r6[0] - d5*s5[0] + a5*r4[0], r6[1] - d5*s5[1] + a5*r4[1], r6[2] - d5*s5[2] + a5*r4[2]};
        //s3 = R_axis_angle(s4, beta1) * r4;
        R_axis_angle(s4, beta1);
        s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
              tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
              tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
        tmp = Norm(s3);
        s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
        tmp = s3[1]*s3[1] + s3[0]*s3[0];
        if(tmp > SING_TOL) {
            s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
        }
        else{
            s2 = {sin(q1_sing), cos(q1_sing), 0};
        }
        J_dir(s2, s3, s4, s5, s6, k_E_O);
        sol1 = q_from_J(tmp_J);
        tmp_J.col(1) = -1*tmp_J.col(1);
        sol2 = q_from_low_J(tmp_J);
        sols[2*i] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7} ;
        check_limits(sols[2*i], 7);
        sols[2*i+1] = {sol2[0], sol2[1], sol2[2], sols[2*i][3], sols[2*i][4], sols[2*i][5], sols[2*i][6]} ;
        check_limits(sols[2*i+1], 3);
    }
    return sols;
}


unsigned int franka_ik_q7_arr(const array<double, 3>& r,
                              const array<double, 9>& ROE,
                              const double q7,
                              array<array<double, 7>, 8>& qsols,
                              const double q1_sing = PI / 2) {
    // ROE is in row-first
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    array<double, 3> k_E_O = { ROE[2], ROE[5], ROE[8] };
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    Eigen::Vector3d i_6_O = tmp_R * i_E_O;
    //array<double, 3> s6 = Cross(k_E_O, i_6_O);
    array<double, 3> s6;
    Cross_(k_E_O, i_6_O, s6);
    array<double, 3> r6 = { r[0] - dE * k_E_O[0] - a7 * i_6_O[0], r[1] - dE * k_E_O[1] - a7 * i_6_O[1], r[2] - d1 - dE * k_E_O[2] - a7 * i_6_O[2] };
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    if (tmp > 1) {
        if ((tmp - 1) * (tmp - 1) < SING_TOL) {
            tmp = 1;
        }
        else {
            cout << "ERROR: unable to assembly kinematic chain";
            for (int i = 0; i < 8; ++i) {
                fill(qsols[i].begin(), qsols[i].end(), NAN);
            }
            return 0;
        }
    }
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double, 3> k_C_O = { -r6[0] / l, -r6[1] / l, -r6[2] / l };
    //array<double, 3> i_C_O = Cross(k_C_O, s6);
    array<double, 3> i_C_O;
    Cross_(k_C_O, s6, i_C_O);
    tmp = Norm(i_C_O);
    i_C_O = { i_C_O[0] / tmp, i_C_O[1] / tmp, i_C_O[2] / tmp };
    //array<double, 3> j_C_O = Cross(k_C_O, i_C_O);
    array<double, 3> j_C_O;
    Cross_(k_C_O, i_C_O, j_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double, 3>, 4> s5s;
    double sa2, ca2;
    int n_alphs = 1;
    unsigned int n_sols = 0;
    if (d3 + d5 < l && l < b1 + b2) n_alphs = 2;
    double v[3];
    for (int i = 0; i < n_alphs; i++) {
        sa2 = sin(alpha2);
        ca2 = cos(alpha2);
        tmp = -rz * ca2 / (ry * sa2);
        if (tmp * tmp > 1)
            continue;
        tmp = asin(tmp);
        v[0] = -sa2 * cos(tmp);
        v[1] = -sa2 * sin(tmp);
        v[2] = -ca2;
        s5s[n_sols] = { i_C_O[0] * v[0] + j_C_O[0] * v[1] + k_C_O[0] * v[2],
                       i_C_O[1] * v[0] + j_C_O[1] * v[1] + k_C_O[1] * v[2],
                       i_C_O[2] * v[0] + j_C_O[2] * v[1] + k_C_O[2] * v[2] };
        tmp = 2 * sa2 * cos(tmp);
        //s5[n_sols+1] = s5s[n_sols] + (2*sa2*cos(tmp)*i_C_O);
        s5s[n_sols + 1] = { s5s[n_sols][0] + tmp * i_C_O[0],
                        s5s[n_sols][1] + tmp * i_C_O[1],
                        s5s[n_sols][2] + tmp * i_C_O[2] };
        n_sols += 2;
        alpha2 = beta2 - actmp;
    }
    array<double, 6> sol1;
    array<double, 3> sol2;
    array<double, 3> s4, r4, s3, s2, s5;
    for (int i = 0; i < n_sols; i++) {
        s5 = s5s[i];
        //s4 = Cross(s5, r6);
        Cross_(s5, r6, s4);
        tmp = Norm(s4);
        s4 = { s4[0] / tmp, s4[1] / tmp, s4[2] / tmp };
        //r4 = Cross(s5, s4);
        Cross_(s5, s4, r4);
        r4 = { r6[0] - d5 * s5[0] + a5 * r4[0], r6[1] - d5 * s5[1] + a5 * r4[1], r6[2] - d5 * s5[2] + a5 * r4[2] };
        R_axis_angle(s4, beta1);
        s3 = { tmp_R(0,0) * r4[0] + tmp_R(0,1) * r4[1] + tmp_R(0,2) * r4[2],
              tmp_R(1,0) * r4[0] + tmp_R(1,1) * r4[1] + tmp_R(1,2) * r4[2],
              tmp_R(2,0) * r4[0] + tmp_R(2,1) * r4[1] + tmp_R(2,2) * r4[2] };
        tmp = Norm(s3);
        s3 = { s3[0] / tmp, s3[1] / tmp, s3[2] / tmp };
        tmp = s3[1] * s3[1] + s3[0] * s3[0];
        if (tmp > SING_TOL) {
            s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
        }
        else {
            s2 = { sin(q1_sing), cos(q1_sing), 0 };
        }
        J_dir(s2, s3, s4, s5, s6, k_E_O);
        sol1 = q_from_J(tmp_J);
        tmp_J.col(1) = -1 * tmp_J.col(1);
        sol2 = q_from_low_J(tmp_J);
        qsols[2 * i] = { sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7 };
        check_limits(qsols[2 * i], 7);
        qsols[2 * i + 1] = { sol2[0], sol2[1], sol2[2], qsols[2 * i][3], qsols[2 * i][4], qsols[2 * i][5], qsols[2 * i][6] };
        check_limits(qsols[2 * i + 1], 3);
    }
    for (int i = 2 * n_sols; i < 8; i++)
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    return 2 * n_sols;
}




vector<array<double,7>> franka_ik_q4(const array<double, 3>& r, const array<double, 9>& ROE, const double q4, const double q1_sing=0, const double q7_sing=0) {
    array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    // ROE is in row-first
    array<double,3> tmp_v = {r_EO_O[1]*ROE[8]-r_EO_O[2]*ROE[5], 
                            r_EO_O[2]*ROE[2]-r_EO_O[0]*ROE[8],
                            r_EO_O[0]*ROE[5]-r_EO_O[1]*ROE[2]};
    if(tmp_v[0]*tmp_v[0] + tmp_v[1]*tmp_v[1] + tmp_v[2]*tmp_v[2] < SING_TOL)
        return franka_ik_q7(r, ROE, q7_sing, q1_sing);
    array<double,3> r_O7O_O = {r_EO_O[0]-dE*ROE[2], r_EO_O[1]-dE*ROE[5], r_EO_O[2]-dE*ROE[8]};
    array<double,3> r_O7O_E = {ROE[0]*r_O7O_O[0] + ROE[3]*r_O7O_O[1] + ROE[6]*r_O7O_O[2],
                               ROE[1]*r_O7O_O[0] + ROE[4]*r_O7O_O[1] + ROE[7]*r_O7O_O[2],
                               ROE[2]*r_O7O_O[0] + ROE[5]*r_O7O_O[1] + ROE[8]*r_O7O_O[2]};
    double alpha = q4 + beta1 + beta2 - PI;
    double lo2 = b1*b1 + b2*b2 - 2*b1*b2*cos(alpha);
    double lp2 = lo2 - r_O7O_E[2]*r_O7O_E[2];
    if (lp2*lp2 < SING_TOL) lp2 = 0;
    if (lp2 < 0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double gamma2 = beta2 + asin(b1*sin(alpha)/sqrt(lo2));
    double cg2 = cos(gamma2), sg2 = sin(gamma2);
    double Lp2 = r_O7O_E[0]*r_O7O_E[0] + r_O7O_E[1]*r_O7O_E[1], phi = atan2(-r_O7O_E[1], -r_O7O_E[0]);
    double tmp = (Lp2 + a7*a7 - lp2)/(2*sqrt(Lp2)*a7);
    if((tmp-1)*(tmp-1) < SING_TOL)
        tmp = 1.0;
    if(tmp > 1.0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double psi = acos(tmp), ry, rz;
    double q7s[2] = {-phi - psi - 3*PI/4, -phi + psi - 3*PI/4};
    double gammas[2] = {0,0};
    size_t ind = 0;
    array<double,3> s2, s3, s4, s5, s6, r4, r6, i_C_O, j_C_O, k_C_O;
    array<double,6> sol1;
    array<double,3> sol2;
    vector<array<double,7>> sols(0); 
    for (auto q7 : q7s){
        tmp_v = {cos(-q7 + 3*PI/4), sin(-q7 + 3*PI/4), 0};
        s6 = {ROE[0]*tmp_v[0] + ROE[1]*tmp_v[1], ROE[3]*tmp_v[0] + ROE[4]*tmp_v[1], ROE[6]*tmp_v[0] + ROE[7]*tmp_v[1]};
        tmp_v = {-a7*cos(-q7 + PI/4), -a7*sin(-q7 + PI/4), 0};
        r6 = {ROE[0]*tmp_v[0] + ROE[1]*tmp_v[1], ROE[3]*tmp_v[0] + ROE[4]*tmp_v[1], ROE[6]*tmp_v[0] + ROE[7]*tmp_v[1]};
        r6 = {r6[0]+r_O7O_O[0], r6[1]+r_O7O_O[1], r6[2]+r_O7O_O[2]};
        tmp = Norm(r6);
        k_C_O = {-r6[0]/tmp, -r6[1]/tmp, -r6[2]/tmp};
        Cross_(k_C_O, s6, i_C_O);
        tmp = Norm(i_C_O);
        i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
        Cross_(k_C_O, i_C_O, j_C_O);
        ry = s6[0]*j_C_O[0] + s6[1]*j_C_O[1] + s6[2]*j_C_O[2];
        rz = s6[0]*k_C_O[0] + s6[1]*k_C_O[1] + s6[2]*k_C_O[2];
        tmp = -rz*cg2/(ry*sg2);
        if (tmp*tmp > 1) continue;
        tmp = asin(tmp);
        gammas[0] = tmp;
        gammas[1] = PI - tmp;
        for (auto gamma : gammas){
            tmp_v = {-sg2*cos(gamma), -sg2*sin(gamma), -cg2};
            s5 = {i_C_O[0]*tmp_v[0] + j_C_O[0]*tmp_v[1] + k_C_O[0]*tmp_v[2],
                  i_C_O[1]*tmp_v[0] + j_C_O[1]*tmp_v[1] + k_C_O[1]*tmp_v[2],
                  i_C_O[2]*tmp_v[0] + j_C_O[2]*tmp_v[1] + k_C_O[2]*tmp_v[2]};
            Cross_(s5, r6, s4);
            tmp = Norm(s4);
            s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
            Cross_(s5, s4, r4);
            r4 = {r6[0] - d5*s5[0] + a5*r4[0], r6[1] - d5*s5[1] + a5*r4[1], r6[2] - d5*s5[2] + a5*r4[2]};
            R_axis_angle(s4, beta1);
            s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
                  tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
                  tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
            tmp = Norm(s3);
            s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
            tmp = s3[1]*s3[1] + s3[0]*s3[0];
            if(tmp > SING_TOL) 
                s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
            else
                s2 = {sin(q1_sing), cos(q1_sing), 0};
            J_dir(s2, s3, s4, s5, s6, array<double,3>{ROE[2],ROE[5],ROE[8]});
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1*tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            sols.push_back(array<double,7>{sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7});
            check_limits(sols.back(), 7);
            ind = sols.size() - 1;
            sols.push_back(array<double,7>{sol2[0], sol2[1], sol2[2], sols[ind][3], sols[ind][4], sols[ind][5], sols[ind][6]});
            check_limits(sols.back(), 3);
        }
    }
    return sols;
}

unsigned int franka_ik_q4_arr(const array<double, 3>& r,
                              const array<double, 9>& ROE,
                              const double q4,
                              array<array<double,7>,8>& qsols,
                              const double q1_sing=0,
                              const double q7_sing=0) {
    array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    // ROE is in row-first
    array<double,3> tmp_v = {r_EO_O[1]*ROE[8]-r_EO_O[2]*ROE[5],
                            r_EO_O[2]*ROE[2]-r_EO_O[0]*ROE[8],
                            r_EO_O[0]*ROE[5]-r_EO_O[1]*ROE[2]};
    if(tmp_v[0]*tmp_v[0] + tmp_v[1]*tmp_v[1] + tmp_v[2]*tmp_v[2] < SING_TOL)
        return franka_ik_q7_arr(r, ROE, q7_sing, qsols, q1_sing);
    array<double,3> r_O7O_O = {r_EO_O[0]-dE*ROE[2], r_EO_O[1]-dE*ROE[5], r_EO_O[2]-dE*ROE[8]};
    array<double,3> r_O7O_E = {ROE[0]*r_O7O_O[0] + ROE[3]*r_O7O_O[1] + ROE[6]*r_O7O_O[2],
                               ROE[1]*r_O7O_O[0] + ROE[4]*r_O7O_O[1] + ROE[7]*r_O7O_O[2],
                               ROE[2]*r_O7O_O[0] + ROE[5]*r_O7O_O[1] + ROE[8]*r_O7O_O[2]};
    double alpha = q4 + beta1 + beta2 - PI;
    double lo2 = b1*b1 + b2*b2 - 2*b1*b2*cos(alpha);
    double lp2 = lo2 - r_O7O_E[2]*r_O7O_E[2];
    if (lp2*lp2 < SING_TOL) lp2 = 0;
    if (lp2 < 0){
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
        }
        return 0;
    }
    double gamma2 = beta2 + asin(b1*sin(alpha)/sqrt(lo2));
    double cg2 = cos(gamma2), sg2 = sin(gamma2);
    double Lp2 = r_O7O_E[0]*r_O7O_E[0] + r_O7O_E[1]*r_O7O_E[1], phi = atan2(-r_O7O_E[1], -r_O7O_E[0]);
    double tmp = (Lp2 + a7*a7 - lp2)/(2*sqrt(Lp2)*a7);
    if((tmp-1)*(tmp-1) < SING_TOL)
        tmp = 1.0;
    if(tmp > 1.0){
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
                fill(qsols[i].begin(), qsols[i].end(), NAN);
        }
        return 0;
    }
    double psi = acos(tmp), ry, rz;
    double q7s[2] = {-phi - psi - 3*PI/4, -phi + psi - 3*PI/4};
    double gammas[2] = {0,0};
    unsigned int ind = 0;
    array<double,3> s2, s3, s4, s5, s6, r4, r6, i_C_O, j_C_O, k_C_O;
    array<double,6> sol1;
    array<double,3> sol2;
    //vector<array<double,7>> sols(0);
    for (auto q7 : q7s){
        tmp_v = {cos(-q7 + 3*PI/4), sin(-q7 + 3*PI/4), 0};
        s6 = {ROE[0]*tmp_v[0] + ROE[1]*tmp_v[1], ROE[3]*tmp_v[0] + ROE[4]*tmp_v[1], ROE[6]*tmp_v[0] + ROE[7]*tmp_v[1]};
        tmp_v = {-a7*cos(-q7 + PI/4), -a7*sin(-q7 + PI/4), 0};
        r6 = {ROE[0]*tmp_v[0] + ROE[1]*tmp_v[1], ROE[3]*tmp_v[0] + ROE[4]*tmp_v[1], ROE[6]*tmp_v[0] + ROE[7]*tmp_v[1]};
        r6 = {r6[0]+r_O7O_O[0], r6[1]+r_O7O_O[1], r6[2]+r_O7O_O[2]};
        tmp = Norm(r6);
        k_C_O = {-r6[0]/tmp, -r6[1]/tmp, -r6[2]/tmp};
        Cross_(k_C_O, s6, i_C_O);
        tmp = Norm(i_C_O);
        i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
        Cross_(k_C_O, i_C_O, j_C_O);
        ry = s6[0]*j_C_O[0] + s6[1]*j_C_O[1] + s6[2]*j_C_O[2];
        rz = s6[0]*k_C_O[0] + s6[1]*k_C_O[1] + s6[2]*k_C_O[2];
        tmp = -rz*cg2/(ry*sg2);
        if (tmp*tmp > 1) continue;
        tmp = asin(tmp);
        gammas[0] = tmp;
        gammas[1] = PI - tmp;
        for (auto gamma : gammas){
            tmp_v = {-sg2*cos(gamma), -sg2*sin(gamma), -cg2};
            s5 = {i_C_O[0]*tmp_v[0] + j_C_O[0]*tmp_v[1] + k_C_O[0]*tmp_v[2],
                  i_C_O[1]*tmp_v[0] + j_C_O[1]*tmp_v[1] + k_C_O[1]*tmp_v[2],
                  i_C_O[2]*tmp_v[0] + j_C_O[2]*tmp_v[1] + k_C_O[2]*tmp_v[2]};
            Cross_(s5, r6, s4);
            tmp = Norm(s4);
            s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
            Cross_(s5, s4, r4);
            r4 = {r6[0] - d5*s5[0] + a5*r4[0], r6[1] - d5*s5[1] + a5*r4[1], r6[2] - d5*s5[2] + a5*r4[2]};
            R_axis_angle(s4, beta1);
            s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
                  tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
                  tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
            tmp = Norm(s3);
            s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
            tmp = s3[1]*s3[1] + s3[0]*s3[0];
            if(tmp > SING_TOL)
                s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
            else
                s2 = {sin(q1_sing), cos(q1_sing), 0};
            J_dir(s2, s3, s4, s5, s6, array<double,3>{ROE[2],ROE[5],ROE[8]});
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1*tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            //sols.push_back(array<double,7>{sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7});
            qsols[2*ind] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7};
            //check_limits(sols.back(), 7);
            check_limits(qsols[2*ind], 7);
            //ind = sols.size() - 1;
            //sols.push_back(array<double,7>{sol2[0], sol2[1], sol2[2], sols[ind][3], sols[ind][4], sols[ind][5], sols[ind][6]});
            qsols[2*ind+1] = {sol2[0], sol2[1], sol2[2], qsols[2*ind][3], qsols[2*ind][4], qsols[2*ind][5], qsols[2*ind][6]};
            check_limits(qsols[2*ind+1], 3);
            ind++;
        }
    }
    for (int i = 2*ind; i < 8; ++i) {
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    }
    return 2*ind;
}

vector<array<double,7>> franka_ik_q6_parallel(const array<double, 3>& r_EO_O, const array<double, 9>& ROE, const int sgn, const double q1_sing=0) {
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    array<double,3> r_QO_O = {r_EO_O[0] + (-dE + sgn*d5)*s7[0], r_EO_O[1] + (-dE + sgn*d5)*s7[1], r_EO_O[2] + (-dE + sgn*d5)*s7[2]};
    array<double,3> r_OQ_Q = {-ROE[0]*r_QO_O[0] - ROE[3]*r_QO_O[1] - ROE[6]*r_QO_O[2], 
                              -ROE[1]*r_QO_O[0] - ROE[4]*r_QO_O[1] - ROE[7]*r_QO_O[2],
                              -ROE[2]*r_QO_O[0] - ROE[5]*r_QO_O[1] - ROE[8]*r_QO_O[2]};
    double tmp = b1*b1 - r_OQ_Q[2]*r_OQ_Q[2];
    if (tmp*tmp < SING_TOL)
        tmp = 0;
    if(tmp<0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double lp = sqrt(tmp);
    array<double,3> r_OpQ_Q = {r_OQ_Q[0], r_OQ_Q[1], 0};
    double l_OpQ = sqrt(r_OQ_Q[0]*r_OQ_Q[0] + r_OQ_Q[1]*r_OQ_Q[1]);
    double alphas[2], Ls[2];
    double q7;
    Ls[0] = a5+lp, 
    Ls[1] = a5-lp;
    array<double,3> tmp_v, r_6pQ_Q, i_4_Q, r_4Q_Q, s6_Q, r6_Q, s4_Q, s3_Q, s2, s3, s4, s5, s6;
    Eigen::Matrix<double, 3, 4> partial_J_Q, partial_J_O; 
    Eigen::Matrix3d ROQ;
    ROQ << ROE[0], ROE[1], ROE[2],
           ROE[3], ROE[4], ROE[5],
           ROE[6], ROE[7], ROE[8];
    const array<double,3> k{{0,0,1}};
    array<double,3> s5_Q{{0,0,-1.0*sgn}};
    vector<array<double,7>> sols(0); 
    int tmp_sgn;
    size_t ind;
    array<double,6> sol1;
    array<double,3> sol2;
    for(auto L:Ls){
        tmp = (-L*L + a7*a7 + l_OpQ*l_OpQ) / (2*a7*l_OpQ);
        if((tmp-1)*(tmp-1) < SING_TOL)
            tmp = 1;
        else if((tmp+1)*(tmp+1) < SING_TOL)
            tmp = -1;
        if(tmp*tmp>1) 
            continue;
        alphas[0] = acos(tmp);
        alphas[1] = -acos(tmp);
        for(auto alpha:alphas){
            rotate_by_axis_angle(k, alpha, r_OpQ_Q, r_6pQ_Q);
            r_6pQ_Q = {a7*r_6pQ_Q[0]/l_OpQ, a7*r_6pQ_Q[1]/l_OpQ, a7*r_6pQ_Q[2]/l_OpQ};
            i_4_Q = {r_OpQ_Q[0] - r_6pQ_Q[0], r_OpQ_Q[1] - r_6pQ_Q[1], r_OpQ_Q[2] - r_6pQ_Q[2]};
            tmp = Norm(i_4_Q);
            tmp_sgn = L<0 ? -1:1;
            i_4_Q = {tmp_sgn*i_4_Q[0]/tmp, tmp_sgn*i_4_Q[1]/tmp, tmp_sgn*i_4_Q[2]/tmp};
            r_4Q_Q = {r_6pQ_Q[0] + a5*i_4_Q[0], r_6pQ_Q[1] + a5*i_4_Q[1], r_6pQ_Q[2] + a5*i_4_Q[2]};
            Cross_(r_6pQ_Q, k, s6_Q);
            r6_Q = {r_6pQ_Q[0], r_6pQ_Q[1], r_6pQ_Q[2] -sgn*d5};
            Cross_(i_4_Q, s5_Q, s4_Q);
            //s4_Q = {i_4_Q[1]*s5_Q[2], -i_4_Q[0]*s5_Q[2], 0};
            tmp_v = {r_4Q_Q[0] - r_OQ_Q[0], r_4Q_Q[1] - r_OQ_Q[1], r_4Q_Q[2] - r_OQ_Q[2]};
            rotate_by_axis_angle(s4_Q, beta1, tmp_v, s3_Q);
            tmp = Norm(s3_Q);
            //s3_Q = {s3_Q[0]/tmp,s3_Q[1]/tmp,s3_Q[2]/tmp};
            partial_J_Q << s3_Q[0]/tmp, s4_Q[0], s5_Q[0], s6_Q[0],
                           s3_Q[1]/tmp, s4_Q[1], s5_Q[1], s6_Q[1],
                           s3_Q[2]/tmp, s4_Q[2], s5_Q[2], s6_Q[2];
            partial_J_O = ROQ*partial_J_Q;
            s3 = {partial_J_O(0,0), partial_J_O(1,0), partial_J_O(2,0)};
            s4 = {partial_J_O(0,1), partial_J_O(1,1), partial_J_O(2,1)};
            s5 = {partial_J_O(0,2), partial_J_O(1,2), partial_J_O(2,2)};
            s6 = {partial_J_O(0,3), partial_J_O(1,3), partial_J_O(2,3)};
            //k_6_Q = -r_6pQ_Q / norm(r_6pQ_Q)
            //q7 = np.arctan2(cross(k_6_Q, [1, 0, 0])[2], np.dot(k_6_Q, [1, 0, 0])) + PI / 4
            q7 = atan2(r_6pQ_Q[1], -r_6pQ_Q[0]) + PI / 4;
            tmp = s3[1]*s3[1] + s3[0]*s3[0];
            if(tmp > SING_TOL) 
                s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
            else
                s2 = {sin(q1_sing), cos(q1_sing), 0};
            J_dir(s2, s3, s4, s5, s6, s7);
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1*tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            sols.push_back(array<double,7>{sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7});
            check_limits(sols.back(), 7);
            ind = sols.size() - 1;
            sols.push_back(array<double,7>{sol2[0], sol2[1], sol2[2], sols[ind][3], sols[ind][4], sols[ind][5], sols[ind][6]});
            check_limits(sols.back(), 3);
        }
    }
    return sols;
}

unsigned int franka_ik_q6_parallel_arr(const array<double, 3>& r_EO_O,
                                       const array<double, 9>& ROE,
                                       const int sgn,
                                       array<array<double,7>,8>& qsols,
                                       const double q1_sing=0) {
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    array<double,3> r_QO_O = {r_EO_O[0] + (-dE + sgn*d5)*s7[0], r_EO_O[1] + (-dE + sgn*d5)*s7[1], r_EO_O[2] + (-dE + sgn*d5)*s7[2]};
    array<double,3> r_OQ_Q = {-ROE[0]*r_QO_O[0] - ROE[3]*r_QO_O[1] - ROE[6]*r_QO_O[2],
                              -ROE[1]*r_QO_O[0] - ROE[4]*r_QO_O[1] - ROE[7]*r_QO_O[2],
                              -ROE[2]*r_QO_O[0] - ROE[5]*r_QO_O[1] - ROE[8]*r_QO_O[2]};
    double tmp = b1*b1 - r_OQ_Q[2]*r_OQ_Q[2];
    if (tmp*tmp < SING_TOL)
        tmp = 0;
    if(tmp<0){
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
        }
        return 0;
    }
    double lp = sqrt(tmp);
    array<double,3> r_OpQ_Q = {r_OQ_Q[0], r_OQ_Q[1], 0};
    double l_OpQ = sqrt(r_OQ_Q[0]*r_OQ_Q[0] + r_OQ_Q[1]*r_OQ_Q[1]);
    double alphas[2], Ls[2];
    double q7;
    Ls[0] = a5+lp,
    Ls[1] = a5-lp;
    array<double,3> tmp_v, r_6pQ_Q, i_4_Q, r_4Q_Q, s6_Q, r6_Q, s4_Q, s3_Q, s2, s3, s4, s5, s6;
    Eigen::Matrix<double, 3, 4> partial_J_Q, partial_J_O;
    Eigen::Matrix3d ROQ;
    ROQ << ROE[0], ROE[1], ROE[2],
           ROE[3], ROE[4], ROE[5],
           ROE[6], ROE[7], ROE[8];
    const array<double,3> k{{0,0,1}};
    array<double,3> s5_Q{{0,0,-1.0*sgn}};
    //vector<array<double,7>> sols(0);
    int tmp_sgn;
    unsigned int ind=0;
    array<double,6> sol1;
    array<double,3> sol2;
    for(auto L:Ls){
        tmp = (-L*L + a7*a7 + l_OpQ*l_OpQ) / (2*a7*l_OpQ);
        if((tmp-1)*(tmp-1) < SING_TOL)
            tmp = 1;
        else if((tmp+1)*(tmp+1) < SING_TOL)
            tmp = -1;
        if(tmp*tmp>1)
            continue;
        alphas[0] = acos(tmp);
        alphas[1] = -acos(tmp);
        for(auto alpha:alphas){
            rotate_by_axis_angle(k, alpha, r_OpQ_Q, r_6pQ_Q);
            r_6pQ_Q = {a7*r_6pQ_Q[0]/l_OpQ, a7*r_6pQ_Q[1]/l_OpQ, a7*r_6pQ_Q[2]/l_OpQ};
            i_4_Q = {r_OpQ_Q[0] - r_6pQ_Q[0], r_OpQ_Q[1] - r_6pQ_Q[1], r_OpQ_Q[2] - r_6pQ_Q[2]};
            tmp = Norm(i_4_Q);
            tmp_sgn = L<0 ? -1:1;
            i_4_Q = {tmp_sgn*i_4_Q[0]/tmp, tmp_sgn*i_4_Q[1]/tmp, tmp_sgn*i_4_Q[2]/tmp};
            r_4Q_Q = {r_6pQ_Q[0] + a5*i_4_Q[0], r_6pQ_Q[1] + a5*i_4_Q[1], r_6pQ_Q[2] + a5*i_4_Q[2]};
            Cross_(r_6pQ_Q, k, s6_Q);
            r6_Q = {r_6pQ_Q[0], r_6pQ_Q[1], r_6pQ_Q[2] -sgn*d5};
            Cross_(i_4_Q, s5_Q, s4_Q);
            //s4_Q = {i_4_Q[1]*s5_Q[2], -i_4_Q[0]*s5_Q[2], 0};
            tmp_v = {r_4Q_Q[0] - r_OQ_Q[0], r_4Q_Q[1] - r_OQ_Q[1], r_4Q_Q[2] - r_OQ_Q[2]};
            rotate_by_axis_angle(s4_Q, beta1, tmp_v, s3_Q);
            tmp = Norm(s3_Q);
            //s3_Q = {s3_Q[0]/tmp,s3_Q[1]/tmp,s3_Q[2]/tmp};
            partial_J_Q << s3_Q[0]/tmp, s4_Q[0], s5_Q[0], s6_Q[0],
                           s3_Q[1]/tmp, s4_Q[1], s5_Q[1], s6_Q[1],
                           s3_Q[2]/tmp, s4_Q[2], s5_Q[2], s6_Q[2];
            partial_J_O = ROQ*partial_J_Q;
            s3 = {partial_J_O(0,0), partial_J_O(1,0), partial_J_O(2,0)};
            s4 = {partial_J_O(0,1), partial_J_O(1,1), partial_J_O(2,1)};
            s5 = {partial_J_O(0,2), partial_J_O(1,2), partial_J_O(2,2)};
            s6 = {partial_J_O(0,3), partial_J_O(1,3), partial_J_O(2,3)};
            //k_6_Q = -r_6pQ_Q / norm(r_6pQ_Q)
            //q7 = np.arctan2(cross(k_6_Q, [1, 0, 0])[2], np.dot(k_6_Q, [1, 0, 0])) + PI / 4
            q7 = atan2(r_6pQ_Q[1], -r_6pQ_Q[0]) + PI / 4;
            tmp = s3[1]*s3[1] + s3[0]*s3[0];
            if(tmp > SING_TOL)
                s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
            else
                s2 = {sin(q1_sing), cos(q1_sing), 0};
            J_dir(s2, s3, s4, s5, s6, s7);
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1*tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            qsols[2*ind] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7};
            check_limits(qsols[2*ind], 7);
            qsols[2*ind+1] = {sol2[0], sol2[1], sol2[2], qsols[2*ind][3], qsols[2*ind][4], qsols[2*ind][5], qsols[2*ind][6]};
            check_limits(qsols[2*ind+1], 3);
            ind++;
        }
    }
    for (int i = 2*ind; i < 8; ++i) {
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    }
    return 2*ind;
}


vector<array<double,7>> franka_ik_q6(const array<double, 3>& r, const array<double, 9>& ROE, const double q6, const double q1_sing=0, const double q7_sing=0) {
    array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    // ROE is in row-first
    array<double,3> tmp_v = {r_EO_O[1]*ROE[8]-r_EO_O[2]*ROE[5], 
                            r_EO_O[2]*ROE[2]-r_EO_O[0]*ROE[8],
                            r_EO_O[0]*ROE[5]-r_EO_O[1]*ROE[2]};
    if(tmp_v[0]*tmp_v[0] + tmp_v[1]*tmp_v[1] + tmp_v[2]*tmp_v[2] < SING_TOL)
        return franka_ik_q7(r, ROE, q7_sing, q1_sing);
    if(sin(q6)*sin(q6) < SING_TOL)
        return franka_ik_q6_parallel(r_EO_O, ROE, cos(q6)>=0 ? 1:-1, q1_sing);
    // NON-PARALLEL CASE:
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    double gamma1 = PI-q6;
    double cg1 = cos(gamma1);
    double sg1 = sin(gamma1);
    array<double,3> r_O7O_O = {r_EO_O[0]-dE*ROE[2], r_EO_O[1]-dE*ROE[5], r_EO_O[2]-dE*ROE[8]};
    array<double,3> r_PO_O = {r_O7O_O[0] + (a7/tan(gamma1))*s7[0], r_O7O_O[1] + (a7/tan(gamma1))*s7[1], r_O7O_O[2] + (a7/tan(gamma1))*s7[2]};
    double lP = Norm(r_PO_O);
    double lC = a7/sg1;
    double Cx = -(ROE[0]*r_PO_O[0] + ROE[3]*r_PO_O[1] + ROE[6]*r_PO_O[2]);
    double Cy = -(ROE[1]*r_PO_O[0] + ROE[4]*r_PO_O[1] + ROE[7]*r_PO_O[2]);
    double Cz = -(ROE[2]*r_PO_O[0] + ROE[5]*r_PO_O[1] + ROE[8]*r_PO_O[2]);
    double c = sqrt(a5*a5 + (lC+d5)*(lC+d5));
    double tmp = (-b1*b1 + lP*lP + c*c)/(2*lP*c);
    if((tmp-1)*(tmp-1) < SING_TOL)
        tmp = 1.0;
    if(tmp > 1.0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double tau = acos(tmp);
    unsigned int n_gamma_sols = 1;
    if((d3 + d5 + lC < lP) && (lP < b1 + c)) n_gamma_sols = 2;
    double gamma2s[2];
    if(d5 < -lC)
        gamma2s[0] = tau + atan(a5/(d5+lC))+PI;
    else
        gamma2s[0] = tau + atan(a5/(d5+lC));
    if(n_gamma_sols>1)
        gamma2s[1] = gamma2s[0] - 2*tau;
    array<array<double,3>,4> s5s;
    double q7s[4];
    double d, u1, u2;
    unsigned int n_sols = 0;
    for(int i=0; i<n_gamma_sols; i++){
        d = lP*cos(gamma2s[i]);
        tmp = (d + Cz*cg1) / (sqrt(Cx*Cx*sg1*sg1 + Cy*Cy*sg1*sg1));
        if((tmp-1)*(tmp-1) < SING_TOL)
            tmp = 1;
        else if((tmp+1)*(tmp+1) < SING_TOL)
            tmp = -1;
        if(tmp*tmp>1) 
            continue;
        u1 = asin(tmp);
        u2 = atan2(Cx*sg1, Cy*sg1);
        q7s[n_sols] = 5*PI/4 - u1 + u2;
        tmp_v = {-sg1*cos(u1 - u2), -sg1*sin(u1 - u2), cg1};
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
        q7s[n_sols] = PI/4 + u1 + u2;
        tmp_v = {-sg1*cos(PI - u1 - u2), -sg1*sin(PI - u1 - u2), cg1};
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
    }
    array<double,3> s2, s3, s4, s6, r4, r6;
    array<double,6> sol1;
    array<double,3> sol2;
    vector<array<double,7>> sols(2*n_sols); 
    for(int i; i<n_sols; i++){
        r6 = {r_PO_O[0] - lC*s5s[i][0], r_PO_O[1] - lC*s5s[i][1], r_PO_O[2] - lC*s5s[i][2]};
        tmp_v = {r_O7O_O[0]-r6[0], r_O7O_O[1]-r6[1], r_O7O_O[2]-r6[2]};
        Cross_(s7, tmp_v, s6);
        tmp = Norm(s6);
        s6 = {s6[0]/tmp,s6[1]/tmp,s6[2]/tmp};
        Cross_(s5s[i], r6, s4);
        tmp = Norm(s4);
        s4 = {s4[0]/tmp,s4[1]/tmp,s4[2]/tmp};
        Cross_(s5s[i], s4, tmp_v);
        r4 = {r6[0] - d5*s5s[i][0] + a5*tmp_v[0], r6[1] - d5*s5s[i][1] + a5*tmp_v[1], r6[2] - d5*s5s[i][2] + a5*tmp_v[2]};
        rotate_by_axis_angle(s4, beta1, r4, s3);
        tmp = Norm(s3);
        s3 = {s3[0]/tmp,s3[1]/tmp,s3[2]/tmp};
        tmp = s3[1]*s3[1] + s3[0]*s3[0];
        if(tmp > SING_TOL) 
            s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
        else
            s2 = {sin(q1_sing), cos(q1_sing), 0};
        J_dir(s2, s3, s4, s5s[i], s6, s7);
        sol1 = q_from_J(tmp_J);
        tmp_J.col(1) = -1*tmp_J.col(1);
        sol2 = q_from_low_J(tmp_J);
        sols[2*i] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7s[i]};
        check_limits(sols[2*i], 7);
        sols[2*i+1] = {sol2[0], sol2[1], sol2[2], sols[2*i][3], sols[2*i][4], sols[2*i][5], sols[2*i][6]};
        check_limits(sols[2*i+1], 3);
    }
    return sols;
}


unsigned int franka_ik_q6_arr(const array<double, 3>& r,
                              const array<double, 9>& ROE,
                              const double q6,
                              array<array<double,7>,8>& qsols,
                              const double q1_sing=0,
                              const double q7_sing=0) {
    array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    // ROE is in row-first
    array<double,3> tmp_v = {r_EO_O[1]*ROE[8]-r_EO_O[2]*ROE[5],
                            r_EO_O[2]*ROE[2]-r_EO_O[0]*ROE[8],
                            r_EO_O[0]*ROE[5]-r_EO_O[1]*ROE[2]};
    if(tmp_v[0]*tmp_v[0] + tmp_v[1]*tmp_v[1] + tmp_v[2]*tmp_v[2] < SING_TOL)
        return franka_ik_q7_arr(r, ROE, q7_sing, qsols, q1_sing);
    if(sin(q6)*sin(q6) < SING_TOL)
        return franka_ik_q6_parallel_arr(r_EO_O, ROE, cos(q6)>=0 ? 1:-1, qsols, q1_sing);
    // NON-PARALLEL CASE:
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    double gamma1 = PI-q6;
    double cg1 = cos(gamma1);
    double sg1 = sin(gamma1);
    array<double,3> r_O7O_O = {r_EO_O[0]-dE*ROE[2], r_EO_O[1]-dE*ROE[5], r_EO_O[2]-dE*ROE[8]};
    array<double,3> r_PO_O = {r_O7O_O[0] + (a7/tan(gamma1))*s7[0], r_O7O_O[1] + (a7/tan(gamma1))*s7[1], r_O7O_O[2] + (a7/tan(gamma1))*s7[2]};
    double lP = Norm(r_PO_O);
    double lC = a7/sg1;
    double Cx = -(ROE[0]*r_PO_O[0] + ROE[3]*r_PO_O[1] + ROE[6]*r_PO_O[2]);
    double Cy = -(ROE[1]*r_PO_O[0] + ROE[4]*r_PO_O[1] + ROE[7]*r_PO_O[2]);
    double Cz = -(ROE[2]*r_PO_O[0] + ROE[5]*r_PO_O[1] + ROE[8]*r_PO_O[2]);
    double c = sqrt(a5*a5 + (lC+d5)*(lC+d5));
    double tmp = (-b1*b1 + lP*lP + c*c)/(2*lP*c);
    if((tmp-1)*(tmp-1) < SING_TOL)
        tmp = 1.0;
    if(tmp > 1.0){
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
        }
        return 0;
    }
    double tau = acos(tmp);
    unsigned int n_gamma_sols = 1;
    if((d3 + d5 + lC < lP) && (lP < b1 + c)) n_gamma_sols = 2;
    double gamma2s[2];
    if(d5 < -lC)
        gamma2s[0] = tau + atan(a5/(d5+lC))+PI;
    else
        gamma2s[0] = tau + atan(a5/(d5+lC));
    if(n_gamma_sols>1)
        gamma2s[1] = gamma2s[0] - 2*tau;
    array<array<double,3>,4> s5s;
    double q7s[4];
    double d, u1, u2;
    unsigned int n_sols = 0;
    for(int i=0; i<n_gamma_sols; i++){
        d = lP*cos(gamma2s[i]);
        tmp = (d + Cz*cg1) / (sqrt(Cx*Cx*sg1*sg1 + Cy*Cy*sg1*sg1));
        if((tmp-1)*(tmp-1) < SING_TOL)
            tmp = 1;
        else if((tmp+1)*(tmp+1) < SING_TOL)
            tmp = -1;
        if(tmp*tmp>1)
            continue;
        u1 = asin(tmp);
        u2 = atan2(Cx*sg1, Cy*sg1);
        q7s[n_sols] = 5*PI/4 - u1 + u2;
        tmp_v = {-sg1*cos(u1 - u2), -sg1*sin(u1 - u2), cg1};
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
        q7s[n_sols] = PI/4 + u1 + u2;
        tmp_v = {-sg1*cos(PI - u1 - u2), -sg1*sin(PI - u1 - u2), cg1};
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
    }
    array<double,3> s2, s3, s4, s6, r4, r6;
    array<double,6> sol1;
    array<double,3> sol2;
    //vector<array<double,7>> sols(2*n_sols);
    for(int i; i<n_sols; i++){
        r6 = {r_PO_O[0] - lC*s5s[i][0], r_PO_O[1] - lC*s5s[i][1], r_PO_O[2] - lC*s5s[i][2]};
        tmp_v = {r_O7O_O[0]-r6[0], r_O7O_O[1]-r6[1], r_O7O_O[2]-r6[2]};
        Cross_(s7, tmp_v, s6);
        tmp = Norm(s6);
        s6 = {s6[0]/tmp,s6[1]/tmp,s6[2]/tmp};
        Cross_(s5s[i], r6, s4);
        tmp = Norm(s4);
        s4 = {s4[0]/tmp,s4[1]/tmp,s4[2]/tmp};
        Cross_(s5s[i], s4, tmp_v);
        r4 = {r6[0] - d5*s5s[i][0] + a5*tmp_v[0], r6[1] - d5*s5s[i][1] + a5*tmp_v[1], r6[2] - d5*s5s[i][2] + a5*tmp_v[2]};
        rotate_by_axis_angle(s4, beta1, r4, s3);
        tmp = Norm(s3);
        s3 = {s3[0]/tmp,s3[1]/tmp,s3[2]/tmp};
        tmp = s3[1]*s3[1] + s3[0]*s3[0];
        if(tmp > SING_TOL)
            s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
        else
            s2 = {sin(q1_sing), cos(q1_sing), 0};
        J_dir(s2, s3, s4, s5s[i], s6, s7);
        sol1 = q_from_J(tmp_J);
        tmp_J.col(1) = -1*tmp_J.col(1);
        sol2 = q_from_low_J(tmp_J);
        qsols[2*i] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7s[i]};
        check_limits(qsols[2*i], 7);
        qsols[2*i+1] = {sol2[0], sol2[1], sol2[2], qsols[2*i][3], qsols[2*i][4], qsols[2*i][5], qsols[2*i][6]};
        check_limits(qsols[2*i+1], 3);
    }
    for (int i = 2*n_sols; i < 8; ++i) {
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    }
    return 2*n_sols;
}


// FUNCTIONS FOR SWIVEL ANGLE

array<double,2> theta_err_from_q7(const double q7, 
                                  const double theta, 
                                  const Eigen::Vector3d& i_E_O, 
                                  const array<double,3>& k_E_O, 
                                  Eigen::Vector3d& i_6_O, 
                                  const array<double,3>& n1_O, 
                                  const array<double,3>& r_7O_O, 
                                  const array<double,3>& u_7O_O) {
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    i_6_O = tmp_R * i_E_O;
    array<double,3> s6 = Cross(k_E_O, i_6_O);
    //r6 = r_7O_O - a7 * i_6_O
    array<double,3> r6 = {r_7O_O[0] - a7*i_6_O[0], r_7O_O[1] - a7*i_6_O[1], r_7O_O[2] - a7*i_6_O[2]};
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    if(tmp*tmp>1)
        return array<double,2>{1e10, 1e10};
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double,3> k_C_O = {-r6[0]/l, -r6[1]/l, -r6[2]/l};
    array<double,3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
    array<double,3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    double sa2, ca2;
    sa2 = sin(alpha2);
    ca2 = cos(alpha2);
    tmp = -rz * ca2 / (ry * sa2);
    if (tmp*tmp > 1)
        return array<double,2>{1e10, 1e10};
    tmp = asin(tmp);
    double v[3] = {-sa2*cos(tmp), -sa2*sin(tmp), -ca2};
    array<array<double,3>,2> s5s;
    s5s[0] = {i_C_O[0]*v[0] + j_C_O[0]*v[1] + k_C_O[0]*v[2],
              i_C_O[1]*v[0] + j_C_O[1]*v[1] + k_C_O[1]*v[2],
              i_C_O[2]*v[0] + j_C_O[2]*v[1] + k_C_O[2]*v[2]};
    tmp = 2*sa2*cos(tmp);        
    s5s[1] = {s5s[0][0] + tmp*i_C_O[0],
              s5s[0][1] + tmp*i_C_O[1],
              s5s[0][2] + tmp*i_C_O[2]};
    array<double,2> errs;
    array<double,3> s4, r4, n2_O;
    for(int i=0; i<2; i++){
        s4 = Cross(s5s[i], r6);
        tmp = Norm(s4);
        s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
        Cross_(s5s[i], s4, r4);
        r4 = {r6[0] - d5*s5s[i][0] + a5*r4[0], 
              r6[1] - d5*s5s[i][1] + a5*r4[1], 
              r6[2] - d5*s5s[i][2] + a5*r4[2]};
        Cross_(r_7O_O, r4, n2_O);
        tmp = Dot(n2_O, s4);
        if(tmp<0) 
            n2_O = {-n2_O[0], -n2_O[1], -n2_O[2]};
        errs[i] = theta - signed_angle(n1_O, n2_O, u_7O_O);
        errs[i]*=errs[i];
    }
    return errs;
}

array<array<double,7>,2> franka_ik_q7_one_sol(const double q7, 
                                             const Eigen::Vector3d& i_E_O, 
                                             const array<double,3>& k_E_O, 
                                             Eigen::Vector3d& i_6_O, 
                                             const array<double,3>& r_7O_O,
                                             const unsigned int branch,
                                             const double q1_sing) {
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    i_6_O = tmp_R * i_E_O;
    array<double,3> s6 = Cross(k_E_O, i_6_O);
    array<double,3> r6 = {r_7O_O[0] - a7*i_6_O[0], r_7O_O[1] - a7*i_6_O[1], r_7O_O[2] - a7*i_6_O[2]};
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double,3> k_C_O = {-r6[0]/l, -r6[1]/l, -r6[2]/l};
    array<double,3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
    array<double,3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double,3>,4> s5s;
    double sa2, ca2;
    sa2 = sin(alpha2);
    ca2 = cos(alpha2);
    tmp = -rz * ca2 / (ry * sa2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    tmp = asin(tmp);
    double v[3] = {-sa2*cos(tmp), -sa2*sin(tmp), -ca2};
    array<double,3> s5;
    s5 = {i_C_O[0]*v[0] + j_C_O[0]*v[1] + k_C_O[0]*v[2],
          i_C_O[1]*v[0] + j_C_O[1]*v[1] + k_C_O[1]*v[2],
          i_C_O[2]*v[0] + j_C_O[2]*v[1] + k_C_O[2]*v[2]};
    if(branch==1){
        tmp = 2*sa2*cos(tmp);
        s5 = {s5[0] + tmp*i_C_O[0],
              s5[1] + tmp*i_C_O[1],
              s5[2] + tmp*i_C_O[2]};
    }
    array<array<double,7>,2> sols;
    array<double,3> s4, r4, s3, s2;
    array<double,6> sol1;
    array<double,3> sol2;
    s4 = Cross(s5, r6);
    tmp = Norm(s4);
    s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
    r4 = Cross(s5, s4);
    r4 = {r6[0] - d5*s5[0] + a5*r4[0], 
          r6[1] - d5*s5[1] + a5*r4[1], 
          r6[2] - d5*s5[2] + a5*r4[2]};
    R_axis_angle(s4, beta1);
    s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
          tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
          tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
    tmp = Norm(s3);
    s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
    tmp = s3[1]*s3[1] + s3[0]*s3[0];
    if(tmp > SING_TOL) 
        s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
    else
        s2 = {sin(q1_sing), cos(q1_sing), 0};
    J_dir(s2, s3, s4, s5, s6, k_E_O);
    sol1 = q_from_J(tmp_J);
    tmp_J.col(1) = -1*tmp_J.col(1);
    sol2 = q_from_low_J(tmp_J);
    sols[0] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7};
    check_limits(sols[0], 7);
    sols[1] = {sol2[0], sol2[1], sol2[2], sols[0][3], sols[0][4], sols[0][5], sols[0][6]} ;
    check_limits(sols[1], 3); 
    return sols;
}


void franka_ik_q7_one_sol_arr(const double q7,
                              const Eigen::Vector3d& i_E_O,
                              const array<double, 3>& k_E_O,
                              Eigen::Vector3d& i_6_O,
                              const array<double, 3>& r_7O_O,
                              const unsigned int branch,
                              array<array<double, 7>, 8>& qsols,
                              unsigned int ind,
                              const double q1_sing) {
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    i_6_O = tmp_R * i_E_O;
    array<double, 3> s6 = Cross(k_E_O, i_6_O);
    array<double, 3> r6 = { r_7O_O[0] - a7 * i_6_O[0], r_7O_O[1] - a7 * i_6_O[1], r_7O_O[2] - a7 * i_6_O[2] };
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double, 3> k_C_O = { -r6[0] / l, -r6[1] / l, -r6[2] / l };
    array<double, 3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = { i_C_O[0] / tmp, i_C_O[1] / tmp, i_C_O[2] / tmp };
    array<double, 3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double, 3>, 4> s5s;
    double sa2, ca2;
    sa2 = sin(alpha2);
    ca2 = cos(alpha2);
    tmp = -rz * ca2 / (ry * sa2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    tmp = asin(tmp);
    double v[3] = { -sa2 * cos(tmp), -sa2 * sin(tmp), -ca2 };
    array<double, 3> s5;
    s5 = { i_C_O[0] * v[0] + j_C_O[0] * v[1] + k_C_O[0] * v[2],
          i_C_O[1] * v[0] + j_C_O[1] * v[1] + k_C_O[1] * v[2],
          i_C_O[2] * v[0] + j_C_O[2] * v[1] + k_C_O[2] * v[2] };
    if (branch == 1) {
        tmp = 2 * sa2 * cos(tmp);
        s5 = { s5[0] + tmp * i_C_O[0],
              s5[1] + tmp * i_C_O[1],
              s5[2] + tmp * i_C_O[2] };
    }
    array<double, 3> s4, r4, s3, s2;
    array<double, 6> sol1;
    array<double, 3> sol2;
    s4 = Cross(s5, r6);
    tmp = Norm(s4);
    s4 = { s4[0] / tmp, s4[1] / tmp, s4[2] / tmp };
    r4 = Cross(s5, s4);
    r4 = { r6[0] - d5 * s5[0] + a5 * r4[0],
          r6[1] - d5 * s5[1] + a5 * r4[1],
          r6[2] - d5 * s5[2] + a5 * r4[2] };
    R_axis_angle(s4, beta1);
    s3 = { tmp_R(0,0) * r4[0] + tmp_R(0,1) * r4[1] + tmp_R(0,2) * r4[2],
          tmp_R(1,0) * r4[0] + tmp_R(1,1) * r4[1] + tmp_R(1,2) * r4[2],
          tmp_R(2,0) * r4[0] + tmp_R(2,1) * r4[1] + tmp_R(2,2) * r4[2] };
    tmp = Norm(s3);
    s3 = { s3[0] / tmp, s3[1] / tmp, s3[2] / tmp };
    tmp = s3[1] * s3[1] + s3[0] * s3[0];
    if (tmp > SING_TOL)
        s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
    else
        s2 = { sin(q1_sing), cos(q1_sing), 0 };
    J_dir(s2, s3, s4, s5, s6, k_E_O);
    sol1 = q_from_J(tmp_J);
    tmp_J.col(1) = -1 * tmp_J.col(1);
    sol2 = q_from_low_J(tmp_J);
    qsols[2*ind] = { sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7 };
    check_limits(qsols[2 * ind], 7);
    qsols[2*ind + 1] = { sol2[0], sol2[1], sol2[2], qsols[2*ind][3], qsols[2*ind][4], qsols[2*ind][5], qsols[2*ind][6] };
    check_limits(qsols[2 * ind + 1], 3);
}

vector<array<double,7>> franka_ik_swivel(const array<double, 3>& r, const array<double, 9>& ROE, const double theta, const double q1_sing=0, const unsigned int n_points=600) {
    array<double,3> k_E_O = {ROE[2], ROE[5], ROE[8]};
    //array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    //r_7O_O = r_EO_O - dE * k_E_O
    array<double,3> r_7O_O = {r[0] - dE*k_E_O[0], r[1] - dE*k_E_O[1], r[2] - d1 - dE*k_E_O[2]};
    double tmp = sqrt(r_7O_O[1]*r_7O_O[1] + r_7O_O[0]*r_7O_O[0]);
    if (tmp < SING_TOL){
        cout << "ERROR: n1_O is undefined";
        return vector<array<double,7>>();
    }
    array<double,3> n1_O = {r_7O_O[1]/tmp, -r_7O_O[0]/tmp, 0};
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    Eigen::Vector3d i_6_O;
    tmp = Norm(r_7O_O);
    array<double,3> u_7O_O = {r_7O_O[0]/tmp, r_7O_O[1]/tmp, r_7O_O[2]/tmp};
    double q7_step = (q_up[6] - q_low[6])/(n_points-1);
    double q7;


    array<array<double,2>,MAX_N_POINTS> Errs;
    array<array<unsigned int, 2>,MAX_N_POINTS> close_cases;
    array<double,MAX_N_POINTS> q7s;
    unsigned int n_close_cases = 0;
    for(int i=0; i<n_points; i++){
        q7s[i] = q_low[6] + i*q7_step;
        Errs[i] = theta_err_from_q7(q7s[i], theta, i_E_O, k_E_O, i_6_O, n1_O, r_7O_O, u_7O_O);
        if(Errs[i][0] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 0;
            n_close_cases+=1;
        }
        if(Errs[i][1] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 1;
            n_close_cases+=1;
        }
    }


    /*
    array<array<double,2>,MAX_N_POINTS> Errs;
    array<double,MAX_N_POINTS> q7s;
    array<array<unsigned int, 2>,MAX_N_POINTS> close_cases;

    for(int i=0; i<n_points; i++)
        q7s[i] = q_low[6] + i*q7_step;

    const int max_num_threads = 64;
    array<array<array<unsigned int, 2>, MAX_N_POINTS>, max_num_threads> local_buffers;
    array<int, max_num_threads> local_counts = {0};

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int& local_count = local_counts[thread_id];
        auto& local_buffer = local_buffers[thread_id];
        #pragma omp for
        for(unsigned int i=0; i<n_points; i++){
            Errs[i] = theta_err_from_q7(q7s[i], theta, i_E_O, k_E_O, i_6_O, n1_O, r_7O_O, u_7O_O);
            if(Errs[i][0] < ERR_THRESH)
                local_buffer[local_count++] = {i,0};
            if(Errs[i][1] < ERR_THRESH)
                local_buffer[local_count++] = {i,1};
        }
    }

    unsigned int n_close_cases = 0;
    for(int t = 0; t < max_num_threads; t++) {
        for(int j = 0; j < local_counts[t]; j++) {
            close_cases[n_close_cases++] = local_buffers[t][j];
        }
    }
    */



    array<unsigned int, 2> min=close_cases[0];
    vector<array<unsigned int, 2>> best;
    for(int i=1; i<n_close_cases; i++){
        // identify repeated cases i.e. cases where several consecutive solutions passed the threshold
        if(close_cases[i][0]==close_cases[i-1][0]+1){
            if(Errs[close_cases[i][0]][close_cases[i][1]] < Errs[min[0]][min[1]]){
                min = close_cases[i];
                }
        }
        else{
            best.push_back(min);
            min = close_cases[i];
        }
    }
    best.push_back(min);
    size_t n_sols = best.size();
    double e0, e1, e2, e3, q71, q72, q7_opt;
    vector<array<double,7>> sols(2*n_sols);
    array<array<double,7>,2>tmp_sols;
    unsigned int ind = 0;
    for(auto m:best){
        q7_opt = q7s[m[0]];
        if(m[0]>1 &&  m[0]<n_points-2){
            if(Errs[m[0]-2][m[1]]<ERR_THRESH && Errs[m[0]+2][m[1]]<ERR_THRESH){
                if(Errs[m[0]+1][m[1]]<Errs[m[0]-1][m[1]]){
                    // 0=i-1, 1=i, 2=i+1, 3=i+2
                    e0 = Errs[m[0]-1][m[1]];
                    e1 = Errs[m[0]][m[1]];
                    e2 = Errs[m[0]+1][m[1]];
                    e3 = Errs[m[0]+2][m[1]];
                    q71 = q7s[m[0]];
                    q72 = q7s[m[0]+1];
                }
                else{
                    // 0=i-2, 1=i-1, 2=i, 3=i+1
                    e0 = Errs[m[0]-2][m[1]];
                    e1 = Errs[m[0]-1][m[1]];
                    e2 = Errs[m[0]][m[1]];
                    e3 = Errs[m[0]+1][m[1]];
                    q71 = q7s[m[0]-1];
                    q72 = q7s[m[0]];
                }
                tmp = ((e1-e0)*q71 - (e3-e2)*q72 + (e2-e1)*q7_step)/(e1-e0-e3+e2);
                if(tmp>q71 && tmp<q72)
                    q7_opt = tmp;
            }
        } 
        tmp_sols = franka_ik_q7_one_sol(q7_opt, i_E_O, k_E_O, i_6_O, r_7O_O, m[1], q1_sing);
        sols[ind] = tmp_sols[0];
        ind++;
        sols[ind] = tmp_sols[1];
        ind++;
    }
    return sols;
}


unsigned int franka_ik_swivel_arr(const array<double, 3>& r, 
                                  const array<double, 9>& ROE, 
                                  const double theta, 
                                  array<array<double, 7>, 8>& qsols,
                                  const double q1_sing = 0, 
                                  const unsigned int n_points = 600) {
    array<double, 3> k_E_O = { ROE[2], ROE[5], ROE[8] };
    //array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    //r_7O_O = r_EO_O - dE * k_E_O
    array<double, 3> r_7O_O = { r[0] - dE * k_E_O[0], r[1] - dE * k_E_O[1], r[2] - d1 - dE * k_E_O[2] };
    double tmp = sqrt(r_7O_O[1] * r_7O_O[1] + r_7O_O[0] * r_7O_O[0]);
    if (tmp < SING_TOL) {
        cout << "ERROR: n1_O is undefined";
        for (int i = 0; i < 8; i++)
            fill(qsols[i].begin(), qsols[i].end(), NAN);
        return 0;
    }
    array<double, 3> n1_O = { r_7O_O[1] / tmp, -r_7O_O[0] / tmp, 0 };
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    Eigen::Vector3d i_6_O;
    tmp = Norm(r_7O_O);
    array<double, 3> u_7O_O = { r_7O_O[0] / tmp, r_7O_O[1] / tmp, r_7O_O[2] / tmp };
    double q7_step = (q_up[6] - q_low[6]) / (n_points - 1);
    double q7;

    array<array<double, 2>, MAX_N_POINTS> Errs;
    array<array<unsigned int, 2>, MAX_N_POINTS> close_cases;
    array<double, MAX_N_POINTS> q7s;
    unsigned int n_close_cases = 0;
    for (int i = 0; i < n_points; i++) {
        q7s[i] = q_low[6] + i * q7_step;
        Errs[i] = theta_err_from_q7(q7s[i], theta, i_E_O, k_E_O, i_6_O, n1_O, r_7O_O, u_7O_O);
        if (Errs[i][0] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 0;
            n_close_cases += 1;
        }
        if (Errs[i][1] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 1;
            n_close_cases += 1;
        }
    }

    array<unsigned int, 2> min = close_cases[0];
    //array<array<unsigned int, 2>, 16> best; //setting 16 as maximum number of solutions
    //unsigned int num_best = 0;
    vector<array<unsigned int, 2>> best;
    for (int i = 1; i < n_close_cases; i++) {
        // identify repeated cases i.e. cases where several consecutive solutions passed the threshold
        if (close_cases[i][0] == close_cases[i - 1][0] + 1) {
            if (Errs[close_cases[i][0]][close_cases[i][1]] < Errs[min[0]][min[1]]) {
                min = close_cases[i];
            }
        }
        else {
            //best[num_best++] = min;
            best.push_back(min);
            min = close_cases[i];
            //if (num_best == 15) break;
        }
    }
    best.push_back(min);
    //size_t n_sols = best.size();
    unsigned int n_sols = static_cast<unsigned int>(best.size());
    //best[num_best++] = min;
    if (n_sols > 4) {
        cout << "WARNING: Number of solutions is" << 2 * n_sols << "- Only the first 8 solutions found will be returned.";
        n_sols = 4;
    } 
    double e0, e1, e2, e3, q71, q72, q7_opt;
    //vector<array<double, 7>> sols(2 * n_sols);
    //array<array<double, 7>, 2>tmp_sols;
    //unsigned int ind = 0;
    array<unsigned int, 2> m;
    for (int i=0; i < n_sols; i++) {
        m = best[i];
        q7_opt = q7s[m[0]];
        if (m[0] > 1 && m[0] < n_points - 2) {
            if (Errs[m[0] - 2][m[1]] < ERR_THRESH && Errs[m[0] + 2][m[1]] < ERR_THRESH) {
                if (Errs[m[0] + 1][m[1]] < Errs[m[0] - 1][m[1]]) {
                    // 0=i-1, 1=i, 2=i+1, 3=i+2
                    e0 = Errs[m[0] - 1][m[1]];
                    e1 = Errs[m[0]][m[1]];
                    e2 = Errs[m[0] + 1][m[1]];
                    e3 = Errs[m[0] + 2][m[1]];
                    q71 = q7s[m[0]];
                    q72 = q7s[m[0] + 1];
                }
                else {
                    // 0=i-2, 1=i-1, 2=i, 3=i+1
                    e0 = Errs[m[0] - 2][m[1]];
                    e1 = Errs[m[0] - 1][m[1]];
                    e2 = Errs[m[0]][m[1]];
                    e3 = Errs[m[0] + 1][m[1]];
                    q71 = q7s[m[0] - 1];
                    q72 = q7s[m[0]];
                }
                tmp = ((e1 - e0) * q71 - (e3 - e2) * q72 + (e2 - e1) * q7_step) / (e1 - e0 - e3 + e2);
                if (tmp > q71 && tmp < q72)
                    q7_opt = tmp;
            }
        }
        franka_ik_q7_one_sol_arr(q7_opt, i_E_O, k_E_O, i_6_O, r_7O_O, m[1], qsols, i, q1_sing);
    }
    for (int i = 2*n_sols; i < 8; ++i) {
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    }
    return 2*n_sols;
}


double franka_swivel(const array<double,7>& q){
    // swivel angle for a configuration q (not coded for speed purposes)
    array<Eigen::Matrix4d, 9> Ts = franka_fk_all_frames(q);
    array<double,3> r4 = {Ts[3](0,3), Ts[3](1,3), Ts[3](2,3)-d1};
    array<double,3> r7 = {Ts[6](0,3), Ts[6](1,3), Ts[6](2,3)-d1};
    array<double,3> s4 = {Ts[3](0,2), Ts[3](1,2), Ts[3](2,2)};
    double tmp = sqrt(r7[1]*r7[1] + r7[0]*r7[0]);
    if (tmp < SING_TOL){
        cout << "ERROR: n1_O is undefined";
        return NAN;
    }
    array<double,3> n1_O = {r7[1]/tmp, -r7[0]/tmp, 0};
    array<double,3> n2_O = Cross(r7, r4);
    tmp = Dot(n2_O, s4);
    if(tmp<0) 
        n2_O = {-n2_O[0], -n2_O[1], -n2_O[2]};
    tmp = Norm(r7);
    return signed_angle(n1_O, n2_O, array<double,3>{r7[0]/tmp,r7[1]/tmp,r7[2]/tmp});
}

















//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


vector<array<double,7>> franka_J_ik_q7(const array<double, 3>& r,
                                       const array<double, 9>& ROE,
                                       const double q7,
                                       vector<array<array<double, 6>,7>>& Jsols,
                                       const bool joint_angles = false,
                                       const bool Jacobian_ee_is_8=false,
                                       const double q1_sing = PI/2) {
    // ROE is in row-first
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    array<double,3> k_E_O = {ROE[2], ROE[5], ROE[8]};
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    Eigen::Vector3d i_6_O = tmp_R * i_E_O;
    array<double,3> s6 = Cross(k_E_O, i_6_O);
    array<double,3> r6 = {r[0] - dE*k_E_O[0] - a7*i_6_O[0], r[1] - dE*k_E_O[1] - a7*i_6_O[1], r[2] - d1 - dE*k_E_O[2] - a7*i_6_O[2]};
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    if(tmp>1){
        if((tmp - 1) * (tmp - 1) < SING_TOL){
            tmp = 1;
        }
        else{
            cout << "ERROR: unable to assembly kinematic chain";
            return vector<array<double,7>>();
        }
    }
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double,3> k_C_O = {-r6[0]/l, -r6[1]/l, -r6[2]/l};
    array<double,3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
    array<double,3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double,3>,4> s5s;
    double sa2, ca2;
    int n_alphs = 1;
    unsigned int n_sols = 0;
    if (d3+d5 < l && l < b1+b2) n_alphs = 2;
    double v[3];
    for(int i=0; i<n_alphs; i++){
        sa2 = sin(alpha2);
        ca2 = cos(alpha2);
        tmp = -rz * ca2 / (ry * sa2);
        if (tmp*tmp > 1)
            continue;
        tmp = asin(tmp);
        v[0] = -sa2*cos(tmp);
        v[1] = -sa2*sin(tmp);
        v[2] = -ca2;
        s5s[n_sols] = {i_C_O[0]*v[0] + j_C_O[0]*v[1] + k_C_O[0]*v[2],
                       i_C_O[1]*v[0] + j_C_O[1]*v[1] + k_C_O[1]*v[2],
                       i_C_O[2]*v[0] + j_C_O[2]*v[1] + k_C_O[2]*v[2]};
        tmp = 2*sa2*cos(tmp);
        //s5[n_sols+1] = s5s[n_sols] + (2*sa2*cos(tmp)*i_C_O);
        s5s[n_sols+1] = {s5s[n_sols][0] + tmp*i_C_O[0],
                        s5s[n_sols][1] + tmp*i_C_O[1],
                        s5s[n_sols][2] + tmp*i_C_O[2]};
        n_sols+=2;
        alpha2 = beta2 - actmp;
    }
    Jsols.resize(2*n_sols);
    vector<array<double,7>> sols;
    if(joint_angles)
        sols.resize(2*n_sols);
    array<double,6> sol1;
    array<double,3> sol2;
    array<double,3> s4, r4, s3, s2, s5;
    for(int i=0; i<n_sols; i++){
        s5 = s5s[i];
        s4 = Cross(s5, r6);
        tmp = Norm(s4);
        s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
        //r4 = r6 - d5 * s5 + a5 * Cross(s5, s4);
        r4 = Cross(s5, s4);
        r4 = {r6[0] - d5*s5[0] + a5*r4[0], r6[1] - d5*s5[1] + a5*r4[1], r6[2] - d5*s5[2] + a5*r4[2]};
        //s3 = R_axis_angle(s4, beta1) * r4;
        R_axis_angle(s4, beta1);
        s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
              tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
              tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
        tmp = Norm(s3);
        s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
        tmp = s3[1]*s3[1] + s3[0]*s3[0];
        if(tmp > SING_TOL) {
            s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
        }
        else{
            s2 = {sin(q1_sing), cos(q1_sing), 0};
        }
        /*
        cout<<endl<<"case"<<2*i+1<<endl;
        cout<<"s1 = (0,0,1)"<<endl;
        cout<<"r1 = ("<<-r[0]<<","<<-r[1]<<","<<d1-r[2]<<endl;
        cout<<"s2 = ("<<s2[0]<<","<<s2[1]<<","<<s2[2]<<endl;
        cout<<"r2 = ("<<-r[0]<<","<<-r[1]<<","<<d1-r[2]<<endl;
        cout<<"s3 = ("<<s3[0]<<","<<s3[1]<<","<<s3[2]<<endl;
        cout<<"r3 = ("<<-r[0]<<","<<-r[1]<<","<<d1-r[2]<<endl;
        cout<<"s4 = ("<<s4[0]<<","<<s4[1]<<","<<s4[2]<<endl;
        cout<<"r4 = ("<<r4[0]-r[0]<<","<<r4[1]-r[1]<<","<<r4[2]+d1-r[2]<<endl;
        cout<<"s5 = ("<<s5[0]<<","<<s5[1]<<","<<s5[2]<<endl;
        cout<<"r5 = ("<<r6[0]-r[0]<<","<<r6[1]-r[1]<<","<<r6[2]+d1-r[2]<<endl;
        cout<<"s6 = ("<<s6[0]<<","<<s6[1]<<","<<s6[2]<<endl;
        cout<<"r6 = ("<<r6[0]-r[0]<<","<<r6[1]-r[1]<<","<<r6[2]+d1-r[2]<<endl;
        cout<<"s7 = ("<<k_E_O[0]<<k_E_O[1]<<k_E_O[2]<<endl<<endl;
        */
        save_J_sol(s2, s3, s4, s5, s6, k_E_O, r4, r6, r, Jsols, i, Jacobian_ee_is_8);
        if(joint_angles){
            J_dir(s2, s3, s4, s5, s6, k_E_O);
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1*tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            sols[2*i] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7} ;
            check_limits(sols[2*i], 7);
            sols[2*i+1] = {sol2[0], sol2[1], sol2[2], sols[2*i][3], sols[2*i][4], sols[2*i][5], sols[2*i][6]} ;
            check_limits(sols[2*i+1], 3);
        }
    }
    return sols;
}



















unsigned int franka_J_ik_q7_arr(const array<double, 3>& r,
                                const array<double, 9>& ROE,
                                const double q7,
                                array<array<array<double, 6>, 7>, 8>& Jsols,
                                array<array<double, 7>, 8>& qsols,
                                const bool joint_angles = false,
                                const char Jacobian_ee = 'E',
                                const double q1_sing = PI / 2) {
    // ROE is in row-first
    // Jacobian_ee can only be {'E', 'F', '8', '6'}
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    array<double, 3> k_E_O = { ROE[2], ROE[5], ROE[8] };
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    Eigen::Vector3d i_6_O = tmp_R * i_E_O;
    array<double, 3> s6 = Cross(k_E_O, i_6_O);
    array<double, 3> r6 = { r[0] - dE * k_E_O[0] - a7 * i_6_O[0], r[1] - dE * k_E_O[1] - a7 * i_6_O[1], r[2] - d1 - dE * k_E_O[2] - a7 * i_6_O[2] };
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    if (tmp > 1) {
        if ((tmp - 1) * (tmp - 1) < SING_TOL) {
            tmp = 1;
        }
        else {
            cout << "ERROR: unable to assembly kinematic chain";
            for (int i = 0; i < 8; ++i) {
                fill(qsols[i].begin(), qsols[i].end(), NAN);
                for (auto& row : Jsols[i])
                    fill(row.begin(), row.end(), NAN);
            }
            return 0;
        }
    }
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double, 3> k_C_O = { -r6[0] / l, -r6[1] / l, -r6[2] / l };
    array<double, 3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = { i_C_O[0] / tmp, i_C_O[1] / tmp, i_C_O[2] / tmp };
    array<double, 3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double, 3>, 4> s5s;
    double sa2, ca2;
    int n_alphs = 1;
    unsigned int n_sols = 0;
    if (d3 + d5 < l && l < b1 + b2) n_alphs = 2;
    double v[3];
    for (int i = 0; i < n_alphs; i++) {
        sa2 = sin(alpha2);
        ca2 = cos(alpha2);
        tmp = -rz * ca2 / (ry * sa2);
        if (tmp * tmp > 1)
            continue;
        tmp = asin(tmp);
        v[0] = -sa2 * cos(tmp);
        v[1] = -sa2 * sin(tmp);
        v[2] = -ca2;
        s5s[n_sols] = { i_C_O[0] * v[0] + j_C_O[0] * v[1] + k_C_O[0] * v[2],
                       i_C_O[1] * v[0] + j_C_O[1] * v[1] + k_C_O[1] * v[2],
                       i_C_O[2] * v[0] + j_C_O[2] * v[1] + k_C_O[2] * v[2] };
        tmp = 2 * sa2 * cos(tmp);
        //s5[n_sols+1] = s5s[n_sols] + (2*sa2*cos(tmp)*i_C_O);
        s5s[n_sols + 1] = { s5s[n_sols][0] + tmp * i_C_O[0],
                        s5s[n_sols][1] + tmp * i_C_O[1],
                        s5s[n_sols][2] + tmp * i_C_O[2] };
        n_sols += 2;
        alpha2 = beta2 - actmp;
    }
    //Jsols.resize(2 * n_sols);
    //vector<array<double, 7>> sols;
    array<double, 6> sol1;
    array<double, 3> sol2;
    array<double, 3> s4, r4, s3, s2, s5;
    for (int i = 0; i < n_sols; i++) {
        s5 = s5s[i];
        s4 = Cross(s5, r6);
        tmp = Norm(s4);
        s4 = { s4[0] / tmp, s4[1] / tmp, s4[2] / tmp };
        //r4 = r6 - d5 * s5 + a5 * Cross(s5, s4);
        r4 = Cross(s5, s4);
        r4 = { r6[0] - d5 * s5[0] + a5 * r4[0], r6[1] - d5 * s5[1] + a5 * r4[1], r6[2] - d5 * s5[2] + a5 * r4[2] };
        //s3 = R_axis_angle(s4, beta1) * r4;
        R_axis_angle(s4, beta1);
        s3 = { tmp_R(0,0) * r4[0] + tmp_R(0,1) * r4[1] + tmp_R(0,2) * r4[2],
              tmp_R(1,0) * r4[0] + tmp_R(1,1) * r4[1] + tmp_R(1,2) * r4[2],
              tmp_R(2,0) * r4[0] + tmp_R(2,1) * r4[1] + tmp_R(2,2) * r4[2] };
        tmp = Norm(s3);
        s3 = { s3[0] / tmp, s3[1] / tmp, s3[2] / tmp };
        tmp = s3[1] * s3[1] + s3[0] * s3[0];
        if (tmp > SING_TOL) {
            s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
        }
        else {
            s2 = { sin(q1_sing), cos(q1_sing), 0 };
        }
        save_J_sol(s2, s3, s4, s5, s6, k_E_O, r4, r6, r, Jsols, i, Jacobian_ee);
        if (joint_angles) {
            J_dir(s2, s3, s4, s5, s6, k_E_O);
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1 * tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            qsols[2 * i] = { sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7 };
            check_limits(qsols[2 * i], 7);
            qsols[2 * i + 1] = { sol2[0], sol2[1], sol2[2], qsols[2 * i][3], qsols[2 * i][4], qsols[2 * i][5], qsols[2 * i][6] };
            check_limits(qsols[2 * i + 1], 3);
        }
    }
    for (int i = 2*n_sols; i < 8; ++i) {
        for (auto& row : Jsols[i]) 
            fill(row.begin(), row.end(), NAN);  
    }
    for (int i = joint_angles? 2*n_sols : 0; i<8; i++)
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    return 2*n_sols;
}


vector<array<double,7>> franka_J_ik_q4(const array<double, 3>& r,
                                       const array<double, 9>& ROE,
                                       const double q4,
                                       vector<array<array<double, 6>,7>>& Jsols,
                                       const bool joint_angles = false,
                                       const bool Jacobian_ee_is_8 = false,
                                       const double q1_sing = PI/2,
                                       const double q7_sing=0) {
    array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    // ROE is in row-first
    array<double,3> tmp_v = {r_EO_O[1]*ROE[8]-r_EO_O[2]*ROE[5],
                            r_EO_O[2]*ROE[2]-r_EO_O[0]*ROE[8],
                            r_EO_O[0]*ROE[5]-r_EO_O[1]*ROE[2]};
    if(tmp_v[0]*tmp_v[0] + tmp_v[1]*tmp_v[1] + tmp_v[2]*tmp_v[2] < SING_TOL)
        return franka_J_ik_q7(r, ROE, q7_sing, Jsols, joint_angles, Jacobian_ee_is_8, q1_sing);
    array<double,3> r_O7O_O = {r_EO_O[0]-dE*ROE[2], r_EO_O[1]-dE*ROE[5], r_EO_O[2]-dE*ROE[8]};
    array<double,3> r_O7O_E = {ROE[0]*r_O7O_O[0] + ROE[3]*r_O7O_O[1] + ROE[6]*r_O7O_O[2],
                               ROE[1]*r_O7O_O[0] + ROE[4]*r_O7O_O[1] + ROE[7]*r_O7O_O[2],
                               ROE[2]*r_O7O_O[0] + ROE[5]*r_O7O_O[1] + ROE[8]*r_O7O_O[2]};
    double alpha = q4 + beta1 + beta2 - PI;
    double lo2 = b1*b1 + b2*b2 - 2*b1*b2*cos(alpha);
    double lp2 = lo2 - r_O7O_E[2]*r_O7O_E[2];
    if (lp2*lp2 < SING_TOL) lp2 = 0;
    if (lp2 < 0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double gamma2 = beta2 + asin(b1*sin(alpha)/sqrt(lo2));
    double cg2 = cos(gamma2), sg2 = sin(gamma2);
    double Lp2 = r_O7O_E[0]*r_O7O_E[0] + r_O7O_E[1]*r_O7O_E[1], phi = atan2(-r_O7O_E[1], -r_O7O_E[0]);
    double tmp = (Lp2 + a7*a7 - lp2)/(2*sqrt(Lp2)*a7);
    if((tmp-1)*(tmp-1) < SING_TOL)
        tmp = 1.0;
    if(tmp > 1.0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double psi = acos(tmp), ry, rz;
    double q7s[2] = {-phi - psi - 3*PI/4, -phi + psi - 3*PI/4};
    double gammas[2] = {0,0};
    size_t ind = 0;
    array<double,3> s2, s3, s4, s5, s6, r4, r6, i_C_O, j_C_O, k_C_O;
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    array<double,6> sol1;
    array<double,3> sol2;
    vector<array<double,7>> sols(0);
    Jsols.clear();
    for (auto q7 : q7s){
        tmp_v = {cos(-q7 + 3*PI/4), sin(-q7 + 3*PI/4), 0};
        s6 = {ROE[0]*tmp_v[0] + ROE[1]*tmp_v[1], ROE[3]*tmp_v[0] + ROE[4]*tmp_v[1], ROE[6]*tmp_v[0] + ROE[7]*tmp_v[1]};
        tmp_v = {-a7*cos(-q7 + PI/4), -a7*sin(-q7 + PI/4), 0};
        r6 = {ROE[0]*tmp_v[0] + ROE[1]*tmp_v[1], ROE[3]*tmp_v[0] + ROE[4]*tmp_v[1], ROE[6]*tmp_v[0] + ROE[7]*tmp_v[1]};
        r6 = {r6[0]+r_O7O_O[0], r6[1]+r_O7O_O[1], r6[2]+r_O7O_O[2]};
        tmp = Norm(r6);
        k_C_O = {-r6[0]/tmp, -r6[1]/tmp, -r6[2]/tmp};
        Cross_(k_C_O, s6, i_C_O);
        tmp = Norm(i_C_O);
        i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
        Cross_(k_C_O, i_C_O, j_C_O);
        ry = s6[0]*j_C_O[0] + s6[1]*j_C_O[1] + s6[2]*j_C_O[2];
        rz = s6[0]*k_C_O[0] + s6[1]*k_C_O[1] + s6[2]*k_C_O[2];
        tmp = -rz*cg2/(ry*sg2);
        if (tmp*tmp > 1) continue;
        tmp = asin(tmp);
        gammas[0] = tmp;
        gammas[1] = PI - tmp;
        for (auto gamma : gammas){
            tmp_v = {-sg2*cos(gamma), -sg2*sin(gamma), -cg2};
            s5 = {i_C_O[0]*tmp_v[0] + j_C_O[0]*tmp_v[1] + k_C_O[0]*tmp_v[2],
                  i_C_O[1]*tmp_v[0] + j_C_O[1]*tmp_v[1] + k_C_O[1]*tmp_v[2],
                  i_C_O[2]*tmp_v[0] + j_C_O[2]*tmp_v[1] + k_C_O[2]*tmp_v[2]};
            Cross_(s5, r6, s4);
            tmp = Norm(s4);
            s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
            Cross_(s5, s4, r4);
            r4 = {r6[0] - d5*s5[0] + a5*r4[0], r6[1] - d5*s5[1] + a5*r4[1], r6[2] - d5*s5[2] + a5*r4[2]};
            R_axis_angle(s4, beta1);
            s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
                  tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
                  tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
            tmp = Norm(s3);
            s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
            tmp = s3[1]*s3[1] + s3[0]*s3[0];
            if(tmp > SING_TOL)
                s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
            else
                s2 = {sin(q1_sing), cos(q1_sing), 0};
            pushback_J_sol(s2, s3, s4, s5, s6, s7, r4, r6, r, Jsols, Jacobian_ee_is_8);
            if(joint_angles){
                J_dir(s2, s3, s4, s5, s6, s7);
                sol1 = q_from_J(tmp_J);
                tmp_J.col(1) = -1*tmp_J.col(1);
                sol2 = q_from_low_J(tmp_J);
                sols.push_back(array<double,7>{sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7});
                check_limits(sols.back(), 7);
                ind = sols.size() - 1;
                sols.push_back(array<double,7>{sol2[0], sol2[1], sol2[2], sols[ind][3], sols[ind][4], sols[ind][5], sols[ind][6]});
                check_limits(sols.back(), 3);
            }
        }
    }
    return sols;
}


unsigned int franka_J_ik_q4_arr(const array<double, 3>& r,
                                const array<double, 9>& ROE,
                                const double q4,
                                array<array<array<double, 6>, 7>, 8>& Jsols,
                                array<array<double, 7>, 8>& qsols,
                                const bool joint_angles = false,
                                const char Jacobian_ee = 'E',
                                const double q1_sing = PI / 2,
                                const double q7_sing = 0) {
    array<double, 3> r_EO_O = { r[0], r[1], r[2] - d1 };
    // ROE is in row-first
    array<double, 3> tmp_v = { r_EO_O[1] * ROE[8] - r_EO_O[2] * ROE[5],
                            r_EO_O[2] * ROE[2] - r_EO_O[0] * ROE[8],
                            r_EO_O[0] * ROE[5] - r_EO_O[1] * ROE[2] };
    if (tmp_v[0] * tmp_v[0] + tmp_v[1] * tmp_v[1] + tmp_v[2] * tmp_v[2] < SING_TOL)
        return franka_J_ik_q7_arr(r, ROE, q7_sing, Jsols, qsols, joint_angles, Jacobian_ee, q1_sing);
    array<double, 3> r_O7O_O = { r_EO_O[0] - dE * ROE[2], r_EO_O[1] - dE * ROE[5], r_EO_O[2] - dE * ROE[8] };
    array<double, 3> r_O7O_E = { ROE[0] * r_O7O_O[0] + ROE[3] * r_O7O_O[1] + ROE[6] * r_O7O_O[2],
                               ROE[1] * r_O7O_O[0] + ROE[4] * r_O7O_O[1] + ROE[7] * r_O7O_O[2],
                               ROE[2] * r_O7O_O[0] + ROE[5] * r_O7O_O[1] + ROE[8] * r_O7O_O[2] };
    double alpha = q4 + beta1 + beta2 - PI;
    double lo2 = b1 * b1 + b2 * b2 - 2 * b1 * b2 * cos(alpha);
    double lp2 = lo2 - r_O7O_E[2] * r_O7O_E[2];
    if (lp2 * lp2 < SING_TOL) lp2 = 0;
    if (lp2 < 0) {
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
            for (auto& row : Jsols[i])
                fill(row.begin(), row.end(), NAN);
        }
        return 0;
    }
    double gamma2 = beta2 + asin(b1 * sin(alpha) / sqrt(lo2));
    double cg2 = cos(gamma2), sg2 = sin(gamma2);
    double Lp2 = r_O7O_E[0] * r_O7O_E[0] + r_O7O_E[1] * r_O7O_E[1], phi = atan2(-r_O7O_E[1], -r_O7O_E[0]);
    double tmp = (Lp2 + a7 * a7 - lp2) / (2 * sqrt(Lp2) * a7);
    if ((tmp - 1) * (tmp - 1) < SING_TOL)
        tmp = 1.0;
    if (tmp > 1.0) {
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
            for (auto& row : Jsols[i])
                fill(row.begin(), row.end(), NAN);
        }
        return 0;
    }
    double psi = acos(tmp), ry, rz;
    double q7s[2] = { -phi - psi - 3 * PI / 4, -phi + psi - 3 * PI / 4 };
    double gammas[2] = { 0,0 };
    size_t ind = 0;
    array<double, 3> s2, s3, s4, s5, s6, r4, r6, i_C_O, j_C_O, k_C_O;
    array<double, 3> s7 = { ROE[2],ROE[5],ROE[8] };
    array<double, 6> sol1;
    array<double, 3> sol2;
    for (auto q7 : q7s) {
        tmp_v = { cos(-q7 + 3 * PI / 4), sin(-q7 + 3 * PI / 4), 0 };
        s6 = { ROE[0] * tmp_v[0] + ROE[1] * tmp_v[1], ROE[3] * tmp_v[0] + ROE[4] * tmp_v[1], ROE[6] * tmp_v[0] + ROE[7] * tmp_v[1] };
        tmp_v = { -a7 * cos(-q7 + PI / 4), -a7 * sin(-q7 + PI / 4), 0 };
        r6 = { ROE[0] * tmp_v[0] + ROE[1] * tmp_v[1], ROE[3] * tmp_v[0] + ROE[4] * tmp_v[1], ROE[6] * tmp_v[0] + ROE[7] * tmp_v[1] };
        r6 = { r6[0] + r_O7O_O[0], r6[1] + r_O7O_O[1], r6[2] + r_O7O_O[2] };
        tmp = Norm(r6);
        k_C_O = { -r6[0] / tmp, -r6[1] / tmp, -r6[2] / tmp };
        Cross_(k_C_O, s6, i_C_O);
        tmp = Norm(i_C_O);
        i_C_O = { i_C_O[0] / tmp, i_C_O[1] / tmp, i_C_O[2] / tmp };
        Cross_(k_C_O, i_C_O, j_C_O);
        ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
        rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
        tmp = -rz * cg2 / (ry * sg2);
        if (tmp * tmp > 1) continue;
        tmp = asin(tmp);
        gammas[0] = tmp;
        gammas[1] = PI - tmp;
        for (auto gamma : gammas) {
            tmp_v = { -sg2 * cos(gamma), -sg2 * sin(gamma), -cg2 };
            s5 = { i_C_O[0] * tmp_v[0] + j_C_O[0] * tmp_v[1] + k_C_O[0] * tmp_v[2],
                  i_C_O[1] * tmp_v[0] + j_C_O[1] * tmp_v[1] + k_C_O[1] * tmp_v[2],
                  i_C_O[2] * tmp_v[0] + j_C_O[2] * tmp_v[1] + k_C_O[2] * tmp_v[2] };
            Cross_(s5, r6, s4);
            tmp = Norm(s4);
            s4 = { s4[0] / tmp, s4[1] / tmp, s4[2] / tmp };
            Cross_(s5, s4, r4);
            r4 = { r6[0] - d5 * s5[0] + a5 * r4[0], r6[1] - d5 * s5[1] + a5 * r4[1], r6[2] - d5 * s5[2] + a5 * r4[2] };
            R_axis_angle(s4, beta1);
            s3 = { tmp_R(0,0) * r4[0] + tmp_R(0,1) * r4[1] + tmp_R(0,2) * r4[2],
                  tmp_R(1,0) * r4[0] + tmp_R(1,1) * r4[1] + tmp_R(1,2) * r4[2],
                  tmp_R(2,0) * r4[0] + tmp_R(2,1) * r4[1] + tmp_R(2,2) * r4[2] };
            tmp = Norm(s3);
            s3 = { s3[0] / tmp, s3[1] / tmp, s3[2] / tmp };
            tmp = s3[1] * s3[1] + s3[0] * s3[0];
            if (tmp > SING_TOL)
                s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
            else
                s2 = { sin(q1_sing), cos(q1_sing), 0 };
            save_J_sol(s2, s3, s4, s5, s6, s7, r4, r6, r, Jsols, ind, Jacobian_ee);
            if (joint_angles) {
                J_dir(s2, s3, s4, s5, s6, s7);
                sol1 = q_from_J(tmp_J);
                tmp_J.col(1) = -1 * tmp_J.col(1);
                sol2 = q_from_low_J(tmp_J);
                qsols[2*ind] = { sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7 };
                check_limits(qsols[2*ind], 7);
                qsols[2*ind + 1] = { sol2[0], sol2[1], sol2[2], qsols[2 * ind][3], qsols[2 * ind][4], qsols[2 * ind][5], qsols[2 * ind][6] };
                check_limits(qsols[2*ind + 1], 3);
            }
            ind++;
        }
    }
    for (int i = 2 * ind; i < 8; ++i) {
        for (auto& row : Jsols[i])
            fill(row.begin(), row.end(), NAN);
    }
    for (int i = joint_angles ? 2 * ind : 0; i < 8; i++)
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    return 2 * ind;
}




vector<array<double,7>> franka_J_ik_q6_parallel(const array<double, 3>& r,
                                                const array<double, 3>& r_EO_O,
                                                const array<double, 9>& ROE,
                                                const int sgn,
                                                vector<array<array<double, 6>,7>>& Jsols,
                                                const bool joint_angles,
                                                const bool Jacobian_ee_is_8,
                                                const double q1_sing) {
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    array<double,3> r_QO_O = {r_EO_O[0] + (-dE + sgn*d5)*s7[0], r_EO_O[1] + (-dE + sgn*d5)*s7[1], r_EO_O[2] + (-dE + sgn*d5)*s7[2]};
    array<double,3> r_OQ_Q = {-ROE[0]*r_QO_O[0] - ROE[3]*r_QO_O[1] - ROE[6]*r_QO_O[2],
                              -ROE[1]*r_QO_O[0] - ROE[4]*r_QO_O[1] - ROE[7]*r_QO_O[2],
                              -ROE[2]*r_QO_O[0] - ROE[5]*r_QO_O[1] - ROE[8]*r_QO_O[2]};
    double tmp = b1*b1 - r_OQ_Q[2]*r_OQ_Q[2];
    if (tmp*tmp < SING_TOL)
        tmp = 0;
    if(tmp<0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double lp = sqrt(tmp);
    array<double,3> r_OpQ_Q = {r_OQ_Q[0], r_OQ_Q[1], 0};
    double l_OpQ = sqrt(r_OQ_Q[0]*r_OQ_Q[0] + r_OQ_Q[1]*r_OQ_Q[1]);
    double alphas[2], Ls[2];
    double q7;
    Ls[0] = a5+lp,
    Ls[1] = a5-lp;
    array<double,3> tmp_v, r_6pQ_Q, i_4_Q, r_4Q_Q, s6_Q, r6_Q, s4_Q, s3_Q, s2, s3, s4, s5, s6, r4, r6;
    Eigen::Matrix<double, 3, 4> partial_J_Q, partial_J_O;
    Eigen::Matrix<double, 3, 2> rs;
    Eigen::Matrix3d ROQ;
    ROQ << ROE[0], ROE[1], ROE[2],
           ROE[3], ROE[4], ROE[5],
           ROE[6], ROE[7], ROE[8];
    const array<double,3> k{{0,0,1}};
    array<double,3> s5_Q{{0,0,-1.0*sgn}};
    Jsols.clear();
    vector<array<double,7>> sols(0);
    int tmp_sgn;
    size_t ind;
    array<double,6> sol1;
    array<double,3> sol2;
    for(auto L:Ls){
        tmp = (-L*L + a7*a7 + l_OpQ*l_OpQ) / (2*a7*l_OpQ);
        if((tmp-1)*(tmp-1) < SING_TOL)
            tmp = 1;
        else if((tmp+1)*(tmp+1) < SING_TOL)
            tmp = -1;
        if(tmp*tmp>1)
            continue;
        alphas[0] = acos(tmp);
        alphas[1] = -acos(tmp);
        for(auto alpha:alphas){
            rotate_by_axis_angle(k, alpha, r_OpQ_Q, r_6pQ_Q);
            r_6pQ_Q = {a7*r_6pQ_Q[0]/l_OpQ, a7*r_6pQ_Q[1]/l_OpQ, a7*r_6pQ_Q[2]/l_OpQ};
            i_4_Q = {r_OpQ_Q[0] - r_6pQ_Q[0], r_OpQ_Q[1] - r_6pQ_Q[1], r_OpQ_Q[2] - r_6pQ_Q[2]};
            tmp = Norm(i_4_Q);
            tmp_sgn = L<0 ? -1:1;
            i_4_Q = {tmp_sgn*i_4_Q[0]/tmp, tmp_sgn*i_4_Q[1]/tmp, tmp_sgn*i_4_Q[2]/tmp};
            r_4Q_Q = {r_6pQ_Q[0] + a5*i_4_Q[0], r_6pQ_Q[1] + a5*i_4_Q[1], r_6pQ_Q[2] + a5*i_4_Q[2]};
            Cross_(r_6pQ_Q, k, s6_Q);
            r6_Q = {r_6pQ_Q[0], r_6pQ_Q[1], r_6pQ_Q[2] -sgn*d5};
            rs << r_4Q_Q[0], r6_Q[0],
                  r_4Q_Q[1], r6_Q[1],
                  r_4Q_Q[2], r6_Q[2];
            rs = ROQ*rs; // r_4Q_O, r_6Q_O
            r4 = {rs(0,0)+r_QO_O[0], rs(1,0)+r_QO_O[1], rs(2,0)+r_QO_O[2]};
            r6 = {rs(0,1)+r_QO_O[0], rs(1,1)+r_QO_O[1], rs(2,1)+r_QO_O[2]};
            Cross_(i_4_Q, s5_Q, s4_Q);
            //s4_Q = {i_4_Q[1]*s5_Q[2], -i_4_Q[0]*s5_Q[2], 0};
            tmp_v = {r_4Q_Q[0] - r_OQ_Q[0], r_4Q_Q[1] - r_OQ_Q[1], r_4Q_Q[2] - r_OQ_Q[2]};
            rotate_by_axis_angle(s4_Q, beta1, tmp_v, s3_Q);
            tmp = Norm(s3_Q);
            //s3_Q = {s3_Q[0]/tmp,s3_Q[1]/tmp,s3_Q[2]/tmp};
            partial_J_Q << s3_Q[0]/tmp, s4_Q[0], s5_Q[0], s6_Q[0],
                           s3_Q[1]/tmp, s4_Q[1], s5_Q[1], s6_Q[1],
                           s3_Q[2]/tmp, s4_Q[2], s5_Q[2], s6_Q[2];
            partial_J_O = ROQ*partial_J_Q;
            s3 = {partial_J_O(0,0), partial_J_O(1,0), partial_J_O(2,0)};
            s4 = {partial_J_O(0,1), partial_J_O(1,1), partial_J_O(2,1)};
            s5 = {partial_J_O(0,2), partial_J_O(1,2), partial_J_O(2,2)};
            s6 = {partial_J_O(0,3), partial_J_O(1,3), partial_J_O(2,3)};
            //k_6_Q = -r_6pQ_Q / norm(r_6pQ_Q)
            //q7 = np.arctan2(cross(k_6_Q, [1, 0, 0])[2], np.dot(k_6_Q, [1, 0, 0])) + PI / 4
            q7 = atan2(r_6pQ_Q[1], -r_6pQ_Q[0]) + PI / 4;
            tmp = s3[1]*s3[1] + s3[0]*s3[0];
            if(tmp > SING_TOL)
                s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
            else
                s2 = {sin(q1_sing), cos(q1_sing), 0};
            pushback_J_sol(s2, s3, s4, s5, s6, s7, r4, r6, r, Jsols, Jacobian_ee_is_8);
            if(joint_angles){
                J_dir(s2, s3, s4, s5, s6, s7);
                sol1 = q_from_J(tmp_J);
                tmp_J.col(1) = -1*tmp_J.col(1);
                sol2 = q_from_low_J(tmp_J);
                sols.push_back(array<double,7>{sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7});
                check_limits(sols.back(), 7);
                ind = sols.size() - 1;
                sols.push_back(array<double,7>{sol2[0], sol2[1], sol2[2], sols[ind][3], sols[ind][4], sols[ind][5], sols[ind][6]});
                check_limits(sols.back(), 3);
            }
        }
    }
    return sols;
}

unsigned int franka_J_ik_q6_parallel_arr(const array<double, 3>& r,
                                         const array<double, 3>& r_EO_O,
                                         const array<double, 9>& ROE,
                                         const int sgn,
                                         array<array<array<double, 6>, 7>, 8>& Jsols,
                                         array<array<double, 7>, 8>& qsols,
                                         const bool joint_angles,
                                         const char Jacobian_ee,
                                         const double q1_sing) {
    array<double, 3> s7 = { ROE[2],ROE[5],ROE[8] };
    array<double, 3> r_QO_O = { r_EO_O[0] + (-dE + sgn * d5) * s7[0], r_EO_O[1] + (-dE + sgn * d5) * s7[1], r_EO_O[2] + (-dE + sgn * d5) * s7[2] };
    array<double, 3> r_OQ_Q = { -ROE[0] * r_QO_O[0] - ROE[3] * r_QO_O[1] - ROE[6] * r_QO_O[2],
                              -ROE[1] * r_QO_O[0] - ROE[4] * r_QO_O[1] - ROE[7] * r_QO_O[2],
                              -ROE[2] * r_QO_O[0] - ROE[5] * r_QO_O[1] - ROE[8] * r_QO_O[2] };
    double tmp = b1 * b1 - r_OQ_Q[2] * r_OQ_Q[2];
    if (tmp * tmp < SING_TOL)
        tmp = 0;
    if (tmp < 0) {
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
            for (auto& row : Jsols[i])
                fill(row.begin(), row.end(), NAN);
        }
        return 0;
    }
    double lp = sqrt(tmp);
    array<double, 3> r_OpQ_Q = { r_OQ_Q[0], r_OQ_Q[1], 0 };
    double l_OpQ = sqrt(r_OQ_Q[0] * r_OQ_Q[0] + r_OQ_Q[1] * r_OQ_Q[1]);
    double alphas[2], Ls[2];
    double q7;
    Ls[0] = a5 + lp,
        Ls[1] = a5 - lp;
    array<double, 3> tmp_v, r_6pQ_Q, i_4_Q, r_4Q_Q, s6_Q, r6_Q, s4_Q, s3_Q, s2, s3, s4, s5, s6, r4, r6;
    Eigen::Matrix<double, 3, 4> partial_J_Q, partial_J_O;
    Eigen::Matrix<double, 3, 2> rs;
    Eigen::Matrix3d ROQ;
    ROQ << ROE[0], ROE[1], ROE[2],
        ROE[3], ROE[4], ROE[5],
        ROE[6], ROE[7], ROE[8];
    const array<double, 3> k{ {0,0,1} };
    array<double, 3> s5_Q{ {0,0,-1.0 * sgn} };
    int tmp_sgn;
    unsigned int ind = 0;
    array<double, 6> sol1;
    array<double, 3> sol2;
    for (auto L : Ls) {
        tmp = (-L * L + a7 * a7 + l_OpQ * l_OpQ) / (2 * a7 * l_OpQ);
        if ((tmp - 1) * (tmp - 1) < SING_TOL)
            tmp = 1;
        else if ((tmp + 1) * (tmp + 1) < SING_TOL)
            tmp = -1;
        if (tmp * tmp > 1)
            continue;
        alphas[0] = acos(tmp);
        alphas[1] = -acos(tmp);
        for (auto alpha : alphas) {
            rotate_by_axis_angle(k, alpha, r_OpQ_Q, r_6pQ_Q);
            r_6pQ_Q = { a7 * r_6pQ_Q[0] / l_OpQ, a7 * r_6pQ_Q[1] / l_OpQ, a7 * r_6pQ_Q[2] / l_OpQ };
            i_4_Q = { r_OpQ_Q[0] - r_6pQ_Q[0], r_OpQ_Q[1] - r_6pQ_Q[1], r_OpQ_Q[2] - r_6pQ_Q[2] };
            tmp = Norm(i_4_Q);
            tmp_sgn = L < 0 ? -1 : 1;
            i_4_Q = { tmp_sgn * i_4_Q[0] / tmp, tmp_sgn * i_4_Q[1] / tmp, tmp_sgn * i_4_Q[2] / tmp };
            r_4Q_Q = { r_6pQ_Q[0] + a5 * i_4_Q[0], r_6pQ_Q[1] + a5 * i_4_Q[1], r_6pQ_Q[2] + a5 * i_4_Q[2] };
            Cross_(r_6pQ_Q, k, s6_Q);
            r6_Q = { r_6pQ_Q[0], r_6pQ_Q[1], r_6pQ_Q[2] - sgn * d5 };
            rs << r_4Q_Q[0], r6_Q[0],
                r_4Q_Q[1], r6_Q[1],
                r_4Q_Q[2], r6_Q[2];
            rs = ROQ * rs; // r_4Q_O, r_6Q_O
            r4 = { rs(0,0) + r_QO_O[0], rs(1,0) + r_QO_O[1], rs(2,0) + r_QO_O[2] };
            r6 = { rs(0,1) + r_QO_O[0], rs(1,1) + r_QO_O[1], rs(2,1) + r_QO_O[2] };
            Cross_(i_4_Q, s5_Q, s4_Q);
            //s4_Q = {i_4_Q[1]*s5_Q[2], -i_4_Q[0]*s5_Q[2], 0};
            tmp_v = { r_4Q_Q[0] - r_OQ_Q[0], r_4Q_Q[1] - r_OQ_Q[1], r_4Q_Q[2] - r_OQ_Q[2] };
            rotate_by_axis_angle(s4_Q, beta1, tmp_v, s3_Q);
            tmp = Norm(s3_Q);
            //s3_Q = {s3_Q[0]/tmp,s3_Q[1]/tmp,s3_Q[2]/tmp};
            partial_J_Q << s3_Q[0] / tmp, s4_Q[0], s5_Q[0], s6_Q[0],
                s3_Q[1] / tmp, s4_Q[1], s5_Q[1], s6_Q[1],
                s3_Q[2] / tmp, s4_Q[2], s5_Q[2], s6_Q[2];
            partial_J_O = ROQ * partial_J_Q;
            s3 = { partial_J_O(0,0), partial_J_O(1,0), partial_J_O(2,0) };
            s4 = { partial_J_O(0,1), partial_J_O(1,1), partial_J_O(2,1) };
            s5 = { partial_J_O(0,2), partial_J_O(1,2), partial_J_O(2,2) };
            s6 = { partial_J_O(0,3), partial_J_O(1,3), partial_J_O(2,3) };
            //k_6_Q = -r_6pQ_Q / norm(r_6pQ_Q)
            //q7 = np.arctan2(cross(k_6_Q, [1, 0, 0])[2], np.dot(k_6_Q, [1, 0, 0])) + PI / 4
            q7 = atan2(r_6pQ_Q[1], -r_6pQ_Q[0]) + PI / 4;
            tmp = s3[1] * s3[1] + s3[0] * s3[0];
            if (tmp > SING_TOL)
                s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
            else
                s2 = { sin(q1_sing), cos(q1_sing), 0 };
            save_J_sol(s2, s3, s4, s5, s6, s7, r4, r6, r, Jsols, ind, Jacobian_ee);
            if (joint_angles) {
                J_dir(s2, s3, s4, s5, s6, s7);
                sol1 = q_from_J(tmp_J);
                tmp_J.col(1) = -1 * tmp_J.col(1);
                sol2 = q_from_low_J(tmp_J);
                //sols.push_back(array<double, 7>{sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7});
                //check_limits(sols.back(), 7);
                //ind = sols.size() - 1;
                //sols.push_back(array<double, 7>{sol2[0], sol2[1], sol2[2], sols[ind][3], sols[ind][4], sols[ind][5], sols[ind][6]});
                //check_limits(sols.back(), 3);
                qsols[2*ind] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7 };
                check_limits(qsols[2 * ind], 7);
                qsols[2*ind + 1] = {sol2[0], sol2[1], sol2[2], qsols[2*ind][3], qsols[2*ind][4], qsols[2*ind][5], qsols[2*ind][6] };
                check_limits(qsols[2 * ind + 1], 3);
            }
            ind++;
        }
    }
    for (int i = 2 * ind; i < 8; ++i) {
        for (auto& row : Jsols[i])
            fill(row.begin(), row.end(), NAN);
    }
    for (int i = joint_angles ? 2*ind : 0; i < 8; i++)
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    return 2*ind;
}


vector<array<double,7>> franka_J_ik_q6(const array<double, 3>& r,
                                       const array<double, 9>& ROE,
                                       const double q6,
                                       vector<array<array<double, 6>,7>>& Jsols,
                                       const bool joint_angles = false,
                                       const bool Jacobian_ee_is_8 = false,
                                       const double q1_sing=PI/2,
                                       const double q7_sing=0) {
    array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    // ROE is in row-first
    array<double,3> tmp_v = {r_EO_O[1]*ROE[8]-r_EO_O[2]*ROE[5],
                            r_EO_O[2]*ROE[2]-r_EO_O[0]*ROE[8],
                            r_EO_O[0]*ROE[5]-r_EO_O[1]*ROE[2]};
    if(tmp_v[0]*tmp_v[0] + tmp_v[1]*tmp_v[1] + tmp_v[2]*tmp_v[2] < SING_TOL)
        return franka_J_ik_q7(r, ROE, q7_sing, Jsols, joint_angles, Jacobian_ee_is_8, q1_sing);
    if(sin(q6)*sin(q6) < SING_TOL)
        return franka_J_ik_q6_parallel(r, r_EO_O, ROE, cos(q6)>=0 ? 1:-1, Jsols, joint_angles, Jacobian_ee_is_8, q1_sing);
    // NON-PARALLEL CASE:
    array<double,3> s7 = {ROE[2],ROE[5],ROE[8]};
    double gamma1 = PI-q6;
    double cg1 = cos(gamma1);
    double sg1 = sin(gamma1);
    array<double,3> r_O7O_O = {r_EO_O[0]-dE*ROE[2], r_EO_O[1]-dE*ROE[5], r_EO_O[2]-dE*ROE[8]};
    array<double,3> r_PO_O = {r_O7O_O[0] + (a7/tan(gamma1))*s7[0], r_O7O_O[1] + (a7/tan(gamma1))*s7[1], r_O7O_O[2] + (a7/tan(gamma1))*s7[2]};
    double lP = Norm(r_PO_O);
    double lC = a7/sg1;
    double Cx = -(ROE[0]*r_PO_O[0] + ROE[3]*r_PO_O[1] + ROE[6]*r_PO_O[2]);
    double Cy = -(ROE[1]*r_PO_O[0] + ROE[4]*r_PO_O[1] + ROE[7]*r_PO_O[2]);
    double Cz = -(ROE[2]*r_PO_O[0] + ROE[5]*r_PO_O[1] + ROE[8]*r_PO_O[2]);
    double c = sqrt(a5*a5 + (lC+d5)*(lC+d5));
    double tmp = (-b1*b1 + lP*lP + c*c)/(2*lP*c);
    if((tmp-1)*(tmp-1) < SING_TOL)
        tmp = 1.0;
    if(tmp > 1.0){
        cout << "ERROR: unable to assembly kinematic chain";
        return vector<array<double,7>>();
    }
    double tau = acos(tmp);
    unsigned int n_gamma_sols = 1;
    if((d3 + d5 + lC < lP) && (lP < b1 + c)) n_gamma_sols = 2;
    double gamma2s[2];
    if(d5 < -lC)
        gamma2s[0] = tau + atan(a5/(d5+lC))+PI;
    else
        gamma2s[0] = tau + atan(a5/(d5+lC));
    if(n_gamma_sols>1)
        gamma2s[1] = gamma2s[0] - 2*tau;
    array<array<double,3>,4> s5s;
    double q7s[4];
    double d, u1, u2;
    unsigned int n_sols = 0;
    for(int i=0; i<n_gamma_sols; i++){
        d = lP*cos(gamma2s[i]);
        tmp = (d + Cz*cg1) / (sqrt(Cx*Cx*sg1*sg1 + Cy*Cy*sg1*sg1));
        if((tmp-1)*(tmp-1) < SING_TOL)
            tmp = 1;
        else if((tmp+1)*(tmp+1) < SING_TOL)
            tmp = -1;
        if(tmp*tmp>1)
            continue;
        u1 = asin(tmp);
        u2 = atan2(Cx*sg1, Cy*sg1);
        q7s[n_sols] = 5*PI/4 - u1 + u2;
        tmp_v = {-sg1*cos(u1 - u2), -sg1*sin(u1 - u2), cg1};
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
        q7s[n_sols] = PI/4 + u1 + u2;
        tmp_v = {-sg1*cos(PI - u1 - u2), -sg1*sin(PI - u1 - u2), cg1};
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
    }
    array<double,3> s2, s3, s4, s6, r4, r6;
    array<double,6> sol1;
    array<double,3> sol2;
    Jsols.clear();
    Jsols.resize(2*n_sols);
    vector<array<double,7>> sols;
    if(joint_angles)
        sols.resize(2*n_sols);
    for(int i; i<n_sols; i++){
        r6 = {r_PO_O[0] - lC*s5s[i][0], r_PO_O[1] - lC*s5s[i][1], r_PO_O[2] - lC*s5s[i][2]};
        tmp_v = {r_O7O_O[0]-r6[0], r_O7O_O[1]-r6[1], r_O7O_O[2]-r6[2]};
        Cross_(s7, tmp_v, s6);
        tmp = Norm(s6);
        s6 = {s6[0]/tmp,s6[1]/tmp,s6[2]/tmp};
        Cross_(s5s[i], r6, s4);
        tmp = Norm(s4);
        s4 = {s4[0]/tmp,s4[1]/tmp,s4[2]/tmp};
        Cross_(s5s[i], s4, tmp_v);
        r4 = {r6[0] - d5*s5s[i][0] + a5*tmp_v[0], r6[1] - d5*s5s[i][1] + a5*tmp_v[1], r6[2] - d5*s5s[i][2] + a5*tmp_v[2]};
        rotate_by_axis_angle(s4, beta1, r4, s3);
        tmp = Norm(s3);
        s3 = {s3[0]/tmp,s3[1]/tmp,s3[2]/tmp};
        tmp = s3[1]*s3[1] + s3[0]*s3[0];
        if(tmp > SING_TOL)
            s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
        else
            s2 = {sin(q1_sing), cos(q1_sing), 0};
        save_J_sol(s2, s3, s4, s5s[i], s6, s7, r4, r6, r, Jsols, i, Jacobian_ee_is_8);
        if(joint_angles){
            J_dir(s2, s3, s4, s5s[i], s6, s7);
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1*tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            sols[2*i] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7s[i]};
            check_limits(sols[2*i], 7);
            sols[2*i+1] = {sol2[0], sol2[1], sol2[2], sols[2*i][3], sols[2*i][4], sols[2*i][5], sols[2*i][6]};
            check_limits(sols[2*i+1], 3);
        }
    }
    return sols;
}




unsigned int franka_J_ik_q6_arr(const array<double, 3>& r,
                                const array<double, 9>& ROE,
                                const double q6,
                                array<array<array<double, 6>, 7>, 8>& Jsols,
                                array<array<double, 7>, 8>& qsols,
                                const bool joint_angles = false,
                                const char Jacobian_ee = 'E',
                                const double q1_sing = PI / 2,
                                const double q7_sing = 0) {
    array<double, 3> r_EO_O = { r[0], r[1], r[2] - d1 };
    // ROE is in row-first
    array<double, 3> tmp_v = { r_EO_O[1] * ROE[8] - r_EO_O[2] * ROE[5],
                            r_EO_O[2] * ROE[2] - r_EO_O[0] * ROE[8],
                            r_EO_O[0] * ROE[5] - r_EO_O[1] * ROE[2] };
    if (tmp_v[0] * tmp_v[0] + tmp_v[1] * tmp_v[1] + tmp_v[2] * tmp_v[2] < SING_TOL)
        return franka_J_ik_q7_arr(r, ROE, q7_sing, Jsols, qsols, joint_angles, Jacobian_ee, q1_sing);
    if (sin(q6) * sin(q6) < SING_TOL)
        return franka_J_ik_q6_parallel_arr(r, r_EO_O, ROE, cos(q6) >= 0 ? 1 : -1, Jsols, qsols, joint_angles, Jacobian_ee, q1_sing);
    // NON-PARALLEL CASE:
    array<double, 3> s7 = { ROE[2],ROE[5],ROE[8] };
    double gamma1 = PI - q6;
    double cg1 = cos(gamma1);
    double sg1 = sin(gamma1);
    array<double, 3> r_O7O_O = { r_EO_O[0] - dE * ROE[2], r_EO_O[1] - dE * ROE[5], r_EO_O[2] - dE * ROE[8] };
    array<double, 3> r_PO_O = { r_O7O_O[0] + (a7 / tan(gamma1)) * s7[0], r_O7O_O[1] + (a7 / tan(gamma1)) * s7[1], r_O7O_O[2] + (a7 / tan(gamma1)) * s7[2] };
    double lP = Norm(r_PO_O);
    double lC = a7 / sg1;
    double Cx = -(ROE[0] * r_PO_O[0] + ROE[3] * r_PO_O[1] + ROE[6] * r_PO_O[2]);
    double Cy = -(ROE[1] * r_PO_O[0] + ROE[4] * r_PO_O[1] + ROE[7] * r_PO_O[2]);
    double Cz = -(ROE[2] * r_PO_O[0] + ROE[5] * r_PO_O[1] + ROE[8] * r_PO_O[2]);
    double c = sqrt(a5 * a5 + (lC + d5) * (lC + d5));
    double tmp = (-b1 * b1 + lP * lP + c * c) / (2 * lP * c);
    if ((tmp - 1) * (tmp - 1) < SING_TOL)
        tmp = 1.0;
    if (tmp > 1.0) {
        cout << "ERROR: unable to assembly kinematic chain";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
            for (auto& row : Jsols[i])
                fill(row.begin(), row.end(), NAN);
        }
        return 0;
    }
    double tau = acos(tmp);
    unsigned int n_gamma_sols = 1;
    if ((d3 + d5 + lC < lP) && (lP < b1 + c)) n_gamma_sols = 2;
    double gamma2s[2];
    if (d5 < -lC)
        gamma2s[0] = tau + atan(a5 / (d5 + lC)) + PI;
    else
        gamma2s[0] = tau + atan(a5 / (d5 + lC));
    if (n_gamma_sols > 1)
        gamma2s[1] = gamma2s[0] - 2 * tau;
    array<array<double, 3>, 4> s5s;
    double q7s[4];
    double d, u1, u2;
    unsigned int n_sols = 0;
    for (int i = 0; i < n_gamma_sols; i++) {
        d = lP * cos(gamma2s[i]);
        tmp = (d + Cz * cg1) / (sqrt(Cx * Cx * sg1 * sg1 + Cy * Cy * sg1 * sg1));
        if ((tmp - 1) * (tmp - 1) < SING_TOL)
            tmp = 1;
        else if ((tmp + 1) * (tmp + 1) < SING_TOL)
            tmp = -1;
        if (tmp * tmp > 1)
            continue;
        u1 = asin(tmp);
        u2 = atan2(Cx * sg1, Cy * sg1);
        q7s[n_sols] = 5 * PI / 4 - u1 + u2;
        tmp_v = { -sg1 * cos(u1 - u2), -sg1 * sin(u1 - u2), cg1 };
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
        q7s[n_sols] = PI / 4 + u1 + u2;
        tmp_v = { -sg1 * cos(PI - u1 - u2), -sg1 * sin(PI - u1 - u2), cg1 };
        column_1s_times_vec(ROE, tmp_v, s5s[n_sols]);
        n_sols++;
    }
    array<double, 3> s2, s3, s4, s6, r4, r6;
    array<double, 6> sol1;
    array<double, 3> sol2;
    for (int i; i < n_sols; i++) {
        r6 = { r_PO_O[0] - lC * s5s[i][0], r_PO_O[1] - lC * s5s[i][1], r_PO_O[2] - lC * s5s[i][2] };
        tmp_v = { r_O7O_O[0] - r6[0], r_O7O_O[1] - r6[1], r_O7O_O[2] - r6[2] };
        Cross_(s7, tmp_v, s6);
        tmp = Norm(s6);
        s6 = { s6[0] / tmp,s6[1] / tmp,s6[2] / tmp };
        Cross_(s5s[i], r6, s4);
        tmp = Norm(s4);
        s4 = { s4[0] / tmp,s4[1] / tmp,s4[2] / tmp };
        Cross_(s5s[i], s4, tmp_v);
        r4 = { r6[0] - d5 * s5s[i][0] + a5 * tmp_v[0], r6[1] - d5 * s5s[i][1] + a5 * tmp_v[1], r6[2] - d5 * s5s[i][2] + a5 * tmp_v[2] };
        rotate_by_axis_angle(s4, beta1, r4, s3);
        tmp = Norm(s3);
        s3 = { s3[0] / tmp,s3[1] / tmp,s3[2] / tmp };
        tmp = s3[1] * s3[1] + s3[0] * s3[0];
        if (tmp > SING_TOL)
            s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
        else
            s2 = { sin(q1_sing), cos(q1_sing), 0 };
        save_J_sol(s2, s3, s4, s5s[i], s6, s7, r4, r6, r, Jsols, i, Jacobian_ee);
        if (joint_angles) {
            J_dir(s2, s3, s4, s5s[i], s6, s7);
            sol1 = q_from_J(tmp_J);
            tmp_J.col(1) = -1 * tmp_J.col(1);
            sol2 = q_from_low_J(tmp_J);
            qsols[2 * i] = { sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7s[i] };
            check_limits(qsols[2 * i], 7);
            qsols[2 * i + 1] = { sol2[0], sol2[1], sol2[2], qsols[2 * i][3], qsols[2 * i][4], qsols[2 * i][5], qsols[2 * i][6] };
            check_limits(qsols[2 * i + 1], 3);
        }
    }
    for (int i = 2 * n_sols; i < 8; ++i) {
        for (auto& row : Jsols[i])
            fill(row.begin(), row.end(), NAN);
    }
    for (int i = joint_angles ? 2 * n_sols : 0; i < 8; i++)
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    return 2*n_sols;
}


void franka_J_ik_q7_one_sol(const double q7,
                            const Eigen::Vector3d& i_E_O,
                            const array<double,3>& k_E_O,
                            Eigen::Vector3d& i_6_O,
                            const array<double,3>& r_7O_O,
                            const array<double,3>& r,
                            vector<array<array<double, 6>,7>>& Jsols,
                            vector<array<double, 7>>& sols,
                            const bool joint_angles,
                            const bool Jacobian_ee_is_8,
                            const unsigned int branch,
                            const double q1_sing) {
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    i_6_O = tmp_R * i_E_O;
    array<double,3> s6 = Cross(k_E_O, i_6_O);
    array<double,3> r6 = {r_7O_O[0] - a7*i_6_O[0], r_7O_O[1] - a7*i_6_O[1], r_7O_O[2] - a7*i_6_O[2]};
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double,3> k_C_O = {-r6[0]/l, -r6[1]/l, -r6[2]/l};
    array<double,3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = {i_C_O[0]/tmp, i_C_O[1]/tmp, i_C_O[2]/tmp};
    array<double,3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double,3>,4> s5s;
    double sa2, ca2;
    sa2 = sin(alpha2);
    ca2 = cos(alpha2);
    tmp = -rz * ca2 / (ry * sa2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    tmp = asin(tmp);
    double v[3] = {-sa2*cos(tmp), -sa2*sin(tmp), -ca2};
    array<double,3> s5;
    s5 = {i_C_O[0]*v[0] + j_C_O[0]*v[1] + k_C_O[0]*v[2],
          i_C_O[1]*v[0] + j_C_O[1]*v[1] + k_C_O[1]*v[2],
          i_C_O[2]*v[0] + j_C_O[2]*v[1] + k_C_O[2]*v[2]};
    if(branch==1){
        tmp = 2*sa2*cos(tmp);
        s5 = {s5[0] + tmp*i_C_O[0],
              s5[1] + tmp*i_C_O[1],
              s5[2] + tmp*i_C_O[2]};
    }
    array<double,3> s4, r4, s3, s2;
    array<double,6> sol1;
    array<double,3> sol2;
    array<array<double,7>,2> tmp_sols;
    s4 = Cross(s5, r6);
    tmp = Norm(s4);
    s4 = {s4[0]/tmp, s4[1]/tmp, s4[2]/tmp};
    r4 = Cross(s5, s4);
    r4 = {r6[0] - d5*s5[0] + a5*r4[0],
          r6[1] - d5*s5[1] + a5*r4[1],
          r6[2] - d5*s5[2] + a5*r4[2]};
    R_axis_angle(s4, beta1);
    s3 = {tmp_R(0,0)*r4[0] + tmp_R(0,1)*r4[1] + tmp_R(0,2)*r4[2],
          tmp_R(1,0)*r4[0] + tmp_R(1,1)*r4[1] + tmp_R(1,2)*r4[2],
          tmp_R(2,0)*r4[0] + tmp_R(2,1)*r4[1] + tmp_R(2,2)*r4[2]};
    tmp = Norm(s3);
    s3 = {s3[0]/tmp, s3[1]/tmp, s3[2]/tmp};
    tmp = s3[1]*s3[1] + s3[0]*s3[0];
    if(tmp > SING_TOL)
        s2 = {-s3[1]/sqrt(tmp), s3[0]/sqrt(tmp), 0};
    else
        s2 = {sin(q1_sing), cos(q1_sing), 0};
    pushback_J_sol(s2, s3, s4, s5, s6, k_E_O, r4, r6, r, Jsols, Jacobian_ee_is_8);
    if(joint_angles){
        J_dir(s2, s3, s4, s5, s6, k_E_O);
        sol1 = q_from_J(tmp_J);
        tmp_J.col(1) = -1*tmp_J.col(1);
        sol2 = q_from_low_J(tmp_J);
        tmp_sols[0] = {sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7};
        check_limits(tmp_sols[0], 7);
        tmp_sols[1] = {sol2[0], sol2[1], sol2[2], tmp_sols[0][3], tmp_sols[0][4], tmp_sols[0][5], tmp_sols[0][6]} ;
        check_limits(tmp_sols[1], 3);
        sols.push_back(tmp_sols[0]);
        sols.push_back(tmp_sols[1]);
    }
}



void franka_J_ik_q7_one_sol_arr(const double q7,
                                const Eigen::Vector3d& i_E_O,
                                const array<double, 3>& k_E_O,
                                Eigen::Vector3d& i_6_O,
                                const array<double, 3>& r_7O_O,
                                const array<double, 3>& r,
                                array<array<array<double, 6>, 7>, 8>& Jsols,
                                array<array<double, 7>, 8>& qsols,
                                unsigned int ind,
                                const bool joint_angles,
                                const char Jacobian_ee,
                                const unsigned int branch,
                                const double q1_sing) {
    R_axis_angle(k_E_O, -(q7 - PI / 4));
    i_6_O = tmp_R * i_E_O;
    array<double, 3> s6 = Cross(k_E_O, i_6_O);
    array<double, 3> r6 = { r_7O_O[0] - a7 * i_6_O[0], r_7O_O[1] - a7 * i_6_O[1], r_7O_O[2] - a7 * i_6_O[2] };
    double l = Norm(r6);
    double tmp = (b1 * b1 - l * l - b2 * b2) / (-2 * l * b2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    double actmp = acos(tmp);
    double alpha2 = beta2 + actmp;
    array<double, 3> k_C_O = { -r6[0] / l, -r6[1] / l, -r6[2] / l };
    array<double, 3> i_C_O = Cross(k_C_O, s6);
    tmp = Norm(i_C_O);
    i_C_O = { i_C_O[0] / tmp, i_C_O[1] / tmp, i_C_O[2] / tmp };
    array<double, 3> j_C_O = Cross(k_C_O, i_C_O);
    double ry = s6[0] * j_C_O[0] + s6[1] * j_C_O[1] + s6[2] * j_C_O[2];
    double rz = s6[0] * k_C_O[0] + s6[1] * k_C_O[1] + s6[2] * k_C_O[2];
    array<array<double, 3>, 4> s5s;
    double sa2, ca2;
    sa2 = sin(alpha2);
    ca2 = cos(alpha2);
    tmp = -rz * ca2 / (ry * sa2);
    // The exception tmp*tmp>1 was already handled when Errs was generated
    tmp = asin(tmp);
    double v[3] = { -sa2 * cos(tmp), -sa2 * sin(tmp), -ca2 };
    array<double, 3> s5;
    s5 = { i_C_O[0] * v[0] + j_C_O[0] * v[1] + k_C_O[0] * v[2],
          i_C_O[1] * v[0] + j_C_O[1] * v[1] + k_C_O[1] * v[2],
          i_C_O[2] * v[0] + j_C_O[2] * v[1] + k_C_O[2] * v[2] };
    if (branch == 1) {
        tmp = 2 * sa2 * cos(tmp);
        s5 = { s5[0] + tmp * i_C_O[0],
              s5[1] + tmp * i_C_O[1],
              s5[2] + tmp * i_C_O[2] };
    }
    array<double, 3> s4, r4, s3, s2;
    array<double, 6> sol1;
    array<double, 3> sol2;
    s4 = Cross(s5, r6);
    tmp = Norm(s4);
    s4 = { s4[0] / tmp, s4[1] / tmp, s4[2] / tmp };
    r4 = Cross(s5, s4);
    r4 = { r6[0] - d5 * s5[0] + a5 * r4[0],
          r6[1] - d5 * s5[1] + a5 * r4[1],
          r6[2] - d5 * s5[2] + a5 * r4[2] };
    R_axis_angle(s4, beta1);
    s3 = { tmp_R(0,0) * r4[0] + tmp_R(0,1) * r4[1] + tmp_R(0,2) * r4[2],
          tmp_R(1,0) * r4[0] + tmp_R(1,1) * r4[1] + tmp_R(1,2) * r4[2],
          tmp_R(2,0) * r4[0] + tmp_R(2,1) * r4[1] + tmp_R(2,2) * r4[2] };
    tmp = Norm(s3);
    s3 = { s3[0] / tmp, s3[1] / tmp, s3[2] / tmp };
    tmp = s3[1] * s3[1] + s3[0] * s3[0];
    if (tmp > SING_TOL)
        s2 = { -s3[1] / sqrt(tmp), s3[0] / sqrt(tmp), 0 };
    else
        s2 = { sin(q1_sing), cos(q1_sing), 0 };
    save_J_sol(s2, s3, s4, s5, s6, k_E_O, r4, r6, r, Jsols, ind, Jacobian_ee);
    if (joint_angles) {
        J_dir(s2, s3, s4, s5, s6, k_E_O);
        sol1 = q_from_J(tmp_J);
        tmp_J.col(1) = -1 * tmp_J.col(1);
        sol2 = q_from_low_J(tmp_J);
        qsols[2*ind] = { sol1[0], sol1[1], sol1[2], sol1[3], sol1[4], sol1[5], q7 };
        check_limits(qsols[2*ind], 7);
        qsols[2*ind + 1] = { sol2[0], sol2[1], sol2[2], qsols[2*ind][3], qsols[2*ind][4], qsols[2*ind][5], qsols[2*ind][6] };
        check_limits(qsols[2*ind + 1], 3);
    }
}

vector<array<double,7>> franka_J_ik_swivel(const array<double, 3>& r,
                                           const array<double, 9>& ROE,
                                           const double theta,
                                           vector<array<array<double, 6>,7>>& Jsols,
                                           const bool joint_angles = false,
                                           const bool Jacobian_ee_is_8 = false,
                                           const double q1_sing=PI/2,
                                           const unsigned int n_points=600) {
    array<double,3> k_E_O = {ROE[2], ROE[5], ROE[8]};
    //array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    //r_7O_O = r_EO_O - dE * k_E_O
    array<double,3> r_7O_O = {r[0] - dE*k_E_O[0], r[1] - dE*k_E_O[1], r[2] - d1 - dE*k_E_O[2]};
    double tmp = sqrt(r_7O_O[1]*r_7O_O[1] + r_7O_O[0]*r_7O_O[0]);
    if (tmp < SING_TOL){
        cout << "ERROR: n1_O is undefined";
        return vector<array<double,7>>();
    }
    array<double,3> n1_O = {r_7O_O[1]/tmp, -r_7O_O[0]/tmp, 0};
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    Eigen::Vector3d i_6_O;
    tmp = Norm(r_7O_O);
    array<double,3> u_7O_O = {r_7O_O[0]/tmp, r_7O_O[1]/tmp, r_7O_O[2]/tmp};
    double q7_step = (q_up[6] - q_low[6])/(n_points-1);
    double q7;
    array<array<double,2>,MAX_N_POINTS> Errs;
    array<array<unsigned int, 2>,MAX_N_POINTS> close_cases;
    array<double,MAX_N_POINTS> q7s;
    unsigned int n_close_cases = 0;
    for(int i=0; i<n_points; i++){
        q7s[i] = q_low[6] + i*q7_step;
        Errs[i] = theta_err_from_q7(q7s[i], theta, i_E_O, k_E_O, i_6_O, n1_O, r_7O_O, u_7O_O);
        if(Errs[i][0] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 0;
            n_close_cases+=1;
        }
        if(Errs[i][1] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 1;
            n_close_cases+=1;
        }
        //cout<<"q7 = "<<q7s[i]*180/PI<<", Errs = "<<Errs[i][0]<<", "<<Errs[i][1]<<endl;
    }
    //cout<<endl;
    array<unsigned int, 2> min=close_cases[0];
    vector<array<unsigned int, 2>> best;
    for(int i=1; i<n_close_cases; i++){
        // identify repeated cases i.e. cases where several consecutive solutions passed the threshold
        if(close_cases[i][0]==close_cases[i-1][0]+1){
            if(Errs[close_cases[i][0]][close_cases[i][1]] < Errs[min[0]][min[1]]){
                min = close_cases[i];
                }
        }
        else{
            best.push_back(min);
            min = close_cases[i];
        }
    }
    best.push_back(min);
    int n_sols = best.size();
    double e0, e1, e2, e3, q71, q72, q7_opt;
    vector<array<double,7>> sols(0);
    Jsols.clear();
    array<array<double,7>,2>tmp_sols;
    unsigned int ind = 0;
    for(auto m:best){
        q7_opt = q7s[m[0]];
        if(m[0]>1 &&  m[0]<n_points-2){
            if(Errs[m[0]-2][m[1]]<ERR_THRESH && Errs[m[0]+2][m[1]]<ERR_THRESH){
                if(Errs[m[0]+1][m[1]]<Errs[m[0]-1][m[1]]){
                    // 0=i-1, 1=i, 2=i+1, 3=i+2
                    e0 = Errs[m[0]-1][m[1]];
                    e1 = Errs[m[0]][m[1]];
                    e2 = Errs[m[0]+1][m[1]];
                    e3 = Errs[m[0]+2][m[1]];
                    q71 = q7s[m[0]];
                    q72 = q7s[m[0]+1];
                }
                else{
                    // 0=i-2, 1=i-1, 2=i, 3=i+1
                    e0 = Errs[m[0]-2][m[1]];
                    e1 = Errs[m[0]-1][m[1]];
                    e2 = Errs[m[0]][m[1]];
                    e3 = Errs[m[0]+1][m[1]];
                    q71 = q7s[m[0]-1];
                    q72 = q7s[m[0]];
                }
                tmp = ((e1-e0)*q71 - (e3-e2)*q72 + (e2-e1)*q7_step)/(e1-e0-e3+e2);
                if(tmp>q71 && tmp<q72)
                    q7_opt = tmp;
            }
        }
        franka_J_ik_q7_one_sol(q7_opt, i_E_O, k_E_O, i_6_O, r_7O_O, r, Jsols, sols, joint_angles, Jacobian_ee_is_8, m[1], q1_sing);
    }
    return sols;
}




unsigned int franka_J_ik_swivel_arr(const array<double, 3>& r,
                                    const array<double, 9>& ROE,
                                    const double theta,
                                    array<array<array<double, 6>, 7>, 8>& Jsols,
                                    array<array<double, 7>, 8>& qsols,
                                    const bool joint_angles = false,
                                    const char Jacobian_ee = 'E',
                                    const double q1_sing = PI / 2,
                                    const unsigned int n_points = 600) {
    array<double, 3> k_E_O = { ROE[2], ROE[5], ROE[8] };
    //array<double,3> r_EO_O = {r[0], r[1], r[2] - d1};
    //r_7O_O = r_EO_O - dE * k_E_O
    array<double, 3> r_7O_O = { r[0] - dE * k_E_O[0], r[1] - dE * k_E_O[1], r[2] - d1 - dE * k_E_O[2] };
    double tmp = sqrt(r_7O_O[1] * r_7O_O[1] + r_7O_O[0] * r_7O_O[0]);
    if (tmp < SING_TOL) {
        cout << "ERROR: n1_O is undefined";
        for (int i = 0; i < 8; ++i) {
            fill(qsols[i].begin(), qsols[i].end(), NAN);
            for (auto& row : Jsols[i])
                fill(row.begin(), row.end(), NAN);
        }
        return 0;
    }
    array<double, 3> n1_O = { r_7O_O[1] / tmp, -r_7O_O[0] / tmp, 0 };
    Eigen::Vector3d i_E_O(ROE[0], ROE[3], ROE[6]);
    Eigen::Vector3d i_6_O;
    tmp = Norm(r_7O_O);
    array<double, 3> u_7O_O = { r_7O_O[0] / tmp, r_7O_O[1] / tmp, r_7O_O[2] / tmp };
    double q7_step = (q_up[6] - q_low[6]) / (n_points - 1);
    double q7;
    array<array<double, 2>, MAX_N_POINTS> Errs;
    array<array<unsigned int, 2>, MAX_N_POINTS> close_cases;
    array<double, MAX_N_POINTS> q7s;
    unsigned int n_close_cases = 0;
    for (int i = 0; i < n_points; i++) {
        q7s[i] = q_low[6] + i * q7_step;
        Errs[i] = theta_err_from_q7(q7s[i], theta, i_E_O, k_E_O, i_6_O, n1_O, r_7O_O, u_7O_O);
        if (Errs[i][0] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 0;
            n_close_cases += 1;
        }
        if (Errs[i][1] < ERR_THRESH)
        {
            close_cases[n_close_cases][0] = i;
            close_cases[n_close_cases][1] = 1;
            n_close_cases += 1;
        }
    }
    array<unsigned int, 2> min = close_cases[0];
    //array<array<unsigned int, 2>, 16> best; //setting 16 as maximum number of solutions
    vector<array<unsigned int, 2>> best;
    //unsigned int num_best = 0;
    for (int i = 1; i < n_close_cases; i++) {
        // identify repeated cases i.e. cases where several consecutive solutions passed the threshold
        if (close_cases[i][0] == close_cases[i - 1][0] + 1) {
            if (Errs[close_cases[i][0]][close_cases[i][1]] < Errs[min[0]][min[1]]) {
                min = close_cases[i];
            }
        }
        else {
            //best[num_best++] = min;
            best.push_back(min);
            min = close_cases[i];
            //if (num_best == 16) break;
        }
    }
    //if (num_best < 16) best[num_best++] = min;
    best.push_back(min);
    unsigned int n_sols = static_cast<unsigned int>(best.size());
    if (n_sols > 4) {
        cout << "WARNING: Number of solutions is" << 2 * n_sols << "- Only the first 8 solutions found will be returned.";
        n_sols = 4;
    }
    double e0, e1, e2, e3, q71, q72, q7_opt;
    array<unsigned int, 2> m;
    for (int i = 0; i < n_sols; i++) {
        m = best[i];
        q7_opt = q7s[m[0]];
        if (m[0] > 1 && m[0] < n_points - 2) {
            if (Errs[m[0] - 2][m[1]] < ERR_THRESH && Errs[m[0] + 2][m[1]] < ERR_THRESH) {
                if (Errs[m[0] + 1][m[1]] < Errs[m[0] - 1][m[1]]) {
                    // 0=i-1, 1=i, 2=i+1, 3=i+2
                    e0 = Errs[m[0] - 1][m[1]];
                    e1 = Errs[m[0]][m[1]];
                    e2 = Errs[m[0] + 1][m[1]];
                    e3 = Errs[m[0] + 2][m[1]];
                    q71 = q7s[m[0]];
                    q72 = q7s[m[0] + 1];
                }
                else {
                    // 0=i-2, 1=i-1, 2=i, 3=i+1
                    e0 = Errs[m[0] - 2][m[1]];
                    e1 = Errs[m[0] - 1][m[1]];
                    e2 = Errs[m[0]][m[1]];
                    e3 = Errs[m[0] + 1][m[1]];
                    q71 = q7s[m[0] - 1];
                    q72 = q7s[m[0]];
                }
                tmp = ((e1 - e0) * q71 - (e3 - e2) * q72 + (e2 - e1) * q7_step) / (e1 - e0 - e3 + e2);
                if (tmp > q71 && tmp < q72)
                    q7_opt = tmp;
            }
        }
        franka_J_ik_q7_one_sol_arr(q7_opt, i_E_O, k_E_O, i_6_O, r_7O_O, r, Jsols, qsols, i, joint_angles, Jacobian_ee, m[1], q1_sing);
    }
    for (int i = 2*n_sols; i < 8; ++i) {
        for (auto& row : Jsols[i])
            fill(row.begin(), row.end(), NAN);
    }
    for (int i = joint_angles ? 2*n_sols : 0; i < 8; i++)
        fill(qsols[i].begin(), qsols[i].end(), NAN);
    return 2*n_sols;
}













//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

/*
int main(){
    array<double,9> ROE = {-0.2619425,-0.05341576432114013,0.9636041,0.6560851000000001,0.72239677919841,0.21839260000000027,-0.7077701000000001,0.6894125959344101,-0.15418120000000007};
    array<double,3> r = {0.04447315,0.23491803,0.52787603};
    double q7 = -0.7052075689235409;
    vector<array<double,7>> sols;
    auto start = high_resolution_clock::now();
    sols = franka_ik_q7(r, ROE, q7);
    auto end = high_resolution_clock::now();
    cout << "SOLUTIONS ...................."<<endl;
    for(int i=0; i<sols.size(); i++){
        cout<<"solution "<< i+1 << endl;
        for(int j=0; j<7; j++){
            cout << "q_" << j+1 << " = " << sols[i][j]*180/PI << endl;
        }
        cout << endl;
    }
    auto duration = duration_cast<microseconds>(end - start);
    cout << endl<<"duration:" << duration.count() << endl;
}
*/