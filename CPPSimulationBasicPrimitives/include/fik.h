#ifndef FIK_H
#define FIK_H
#include <iostream>
#include <array> 
#include <vector>
#include <chrono>
#include "Eigen/Dense"
using namespace std;
using namespace std::chrono;


std::vector<std::array<double,7>> franka_J_ik_swivel(const std::array<double, 3>& r,
                                           const std::array<double, 9>& ROE,
                                           const double theta,
                                           std::vector<std::array<std::array<double, 6>,7>>& Jsols,
                                           const bool joint_angles = false,
                                           const bool Jacobian_ee_is_8 = false,
                                           const double q1_sing=1.707,
                                           const unsigned int n_points=600);   
double franka_swivel(const std::array<double,7>& q);  
Eigen::Matrix4d franka_fk(const std::array<double,7>& q, const unsigned int ee=9);


vector<std::array<double,7>> franka_J_ik_q7(
    const std::array<double, 3>& r, 
    const std::array<double, 9>& ROE, 
    const double q7, 
    vector<std::array<std::array<double, 6>,7>>& Jsols, 
    const bool joint_angles = false, 
    const bool Jacobian_ee_is_8=false,
    const double q1_sing = 1.707);

vector<std::array<double,7>> franka_ik_q7(
    const std::array<double, 3>& r, 
    const std::array<double, 9>& ROE, 
    const double q7, 
    const double q1_sing = 1.707);

std::array<std::array<double, 6>,7> partial_Jacobian(
    const std::array<double, 7>& q, 
    const unsigned int ee = 9);
  
unsigned int franka_J_ik_q7_arr(const std::array<double, 3>& r,
                                const std::array<double, 9>& ROE,
                                const double q7,
                                std::array<std::array<std::array<double, 6>, 7>, 8>& Jsols,
                                std::array<std::array<double, 7>, 8>& qsols,
                                const bool joint_angles = false,
                                const char Jacobian_ee = 'E',
                                const double q1_sing = 1.707);

unsigned int franka_ik_q7_arr(const std::array<double, 3>& r,
    const std::array<double, 9>& ROE,
    const double q7,
    std::array<std::array<double, 7>, 8>& qsols,
    const double q1_sing = 1.707);

void check_limits(std::array<double,7>& q, int n);


std::array<double, 7> J_to_q_6joints(
    const std::array<std::array<double, 6>, 7>& J, 
    const std::array<std::array<double, 3>, 3>& R, 
    const char ee = '8');

double signed_angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& s);

double signed_angle(const std::array<double,3>& v1, const std::array<double,3>& v2, const std::array<double,3>& s);

void Cross_(const std::array<double,3>& u, const Eigen::Vector3d& v, std::array<double,3>& w);
void Cross_(const std::array<double, 3>& u, const std::array<double, 3>& v, std::array<double, 3>& w);
void R_axis_angle(const Eigen::Vector3d& s, double theta);
void R_axis_angle(const std::array<double,3>& s, double theta);

vector<std::array<double,7>> franka_ik_swivel(
    const std::array<double, 3>& r, 
    const std::array<double, 9>& ROE, 
    const double theta, 
    const double q1_sing=0, 
    const unsigned int n_points=600);

std::array<std::array<double,6>,7> Jacobian(
    const std::array<double,7>& q, 
    const bool Jacobian_ee_is_8 = false);


#endif