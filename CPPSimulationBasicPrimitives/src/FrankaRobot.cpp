/**
(C) Copyright 2019 DQ Robotics Developers

This file is part of DQ Robotics.

    DQ Robotics is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DQ Robotics is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with DQ Robotics.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
- Murilo M. Marinho (murilo@nml.t.u-tokyo.ac.jp)
*/

#ifndef DQ_ROBOTICS_FRANKA_DH_H
#define DQ_ROBOTICS_FRANKA_DH_H

#include "../include/FrankaRobot.h"
#include<dqrobotics/utils/DQ_Constants.h>

namespace DQ_robotics
{

    // DQ_SerialManipulatorDH FrankaRobot::kinematics()
    // {
    //     const double pi2 = pi/2.0;

    //     MatrixXd franka_dh(5,7);
    //     franka_dh <<    0,      0,          0,          0,      0,      0,          0,
    //                 0.333,      0,      0.316,          0,  0.384,      0,      0.107,
    //                     0,      0,     0.0825,    -0.0825,      0,  0.088,     0.0003,
    //                 -pi2,    pi2,        pi2,       -pi2,    pi2,    pi2,          0,
    //                     0,      0,          0,          0,      0,      0,          0;
    //     DQ_SerialManipulatorDH franka(franka_dh);
    //     std::cout<<"kinematics running..."<<std::endl;
    //     return franka;
    // }

    MatrixXd FrankaRobot::_get_mdh_matrix()
    {
        const double pi2 = pi/2.0;
        Matrix<double,5,7> raw_franka_mdh(5,7);
        raw_franka_mdh <<  0,    0,       0,         0,         0,      0,      -pi/2,
                        0.333,  0, 3.16e-1,         0,   3.84e-1,      0,      0,
                        0,     0,       0,   8.25e-2,  -8.25e-2,      0, 8.8e-2,
                        0,  -pi2,     pi2,       pi2,      -pi2,    pi2,    pi2,
                        0,     0,       0,         0,         0,      0,      0;

        return raw_franka_mdh;
    }

    DQ FrankaRobot::_get_offset_base()
    {
        // return 1+ E_ * DQ(0, 0.596499 - 0.534905, 0, 0);
        // return 1 + E_ * 0.5 * DQ(0, 0.0413,0,0);
        return 1 + E_ * 0.5 * DQ(0, 0,0,0);
        // double z_delta = 45/180 * pi;
        // DQ z_rotation = DQ(cos(z_delta/2), 0, 0, sin(z_delta/2));
        // DQ translation = DQ(0, 0, 0, 1.07e-1);
        // return z_rotation + E_ * 0.5 * z_rotation * translation;
        

    }

    

    DQ FrankaRobot::_get_offset_flange()
    {
        // return 1+E_*0.5*k_*1.07e-1;  // 0.5??????
        // double z_delta = 45/180 * pi;
        // DQ z_rotation = DQ(cos(z_delta/2), 0, 0, sin(z_delta/2));
        // DQ translation = DQ(0, 0, 0, 1.07e-1);
        // return z_rotation + E_ * 0.5 * z_rotation * translation;
        return 1 + E_ * 0.5 * DQ(0, 0, 0, 0);
    
        
        // DQ frame = 1 + E_*0.5*k_*(1.07e-1 + 0.1034);
        // DQ transform_45deg = DQ(0.9238795, 0, 0, -0.3826834).normalize();
        // DQ frame_new = transform_45deg + 0.5 * E_ *k_*(1.07e-1 + 0.1034) * transform_45deg;

        // return frame_new;
    }


    std::tuple<const VectorXd, const VectorXd> FrankaRobot::_get_q_limits()
    {
        const VectorXd q_max_ = ((VectorXd(7) <<  2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525,  2.8973).finished());
        const VectorXd q_min_ = ((VectorXd(7) << -2.8973,-1.7628,-2.8973, -3.0718,-2.8973, -0.0175, -2.8973).finished());
        return std::make_tuple(q_min_, q_max_);
    }

    std::tuple<const VectorXd, const VectorXd> FrankaRobot::_get_q_dot_limits()
    {
        const VectorXd q_min_dot_ = ((VectorXd(7) << -2.1750, -2.1750, -2.1750, -2.1750, -2.61, -2.61, -2.61).finished());
        const VectorXd q_max_dot_ = ((VectorXd(7) <<  2.1750,  2.1750,  2.1750,  2.1750,  2.61,  2.61,  2.61).finished());
        return std::make_tuple(q_min_dot_, q_max_dot_);
    }

    DQ_SerialManipulatorMDH FrankaRobot::kinematics()
    {
        // std::cout<<"frankaemikarobot kinematics running..."<<std::endl;
        // _get_offset_base();
        DQ_SerialManipulatorMDH franka(_get_mdh_matrix());
        // std::cout << "get_offset_flange()1" << _get_offset_flange() << std::endl;
        
        franka.set_base_frame(_get_offset_base());

        DQ xe = 1 + E_*0.5*k_*1.07e-1;
        franka.set_effector(xe);

        
        franka.set_reference_frame(_get_offset_base());
        // std::cout << "get_offset_flange()2" << _get_offset_flange() << std::endl;
        franka.set_effector(_get_offset_flange());
        VectorXd q_min;
        VectorXd q_max;
        VectorXd q_dot_min;
        VectorXd q_dot_max;
        std::tie(q_min, q_max) = _get_q_limits();
        std::tie(q_dot_min, q_dot_max) = _get_q_dot_limits();
        franka.set_lower_q_limit(q_min);
        franka.set_upper_q_limit(q_max);
        franka.set_lower_q_dot_limit(q_dot_min);
        franka.set_upper_q_dot_limit(q_dot_max);
        return franka;
    }

}

#endif
