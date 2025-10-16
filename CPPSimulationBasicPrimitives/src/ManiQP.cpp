#include "../include/ManiQPHuman.h"
#include "../include/franka_ik_He.hpp"
#include "../include/fik.h"

#include <iostream>
#include <dqrobotics/interfaces/vrep/DQ_VrepInterface.h>
#include <dqrobotics/utils/DQ_Geometry.h>

#include <dqrobotics/robot_control/DQ_PseudoinverseController.h>
#include <thread>
#include "yaml-cpp/yaml.h"
#include<dqrobotics/utils/DQ_Geometry.h>
#include<dqrobotics/utils/DQ_Constants.h>

Vector3d calculate_error(MatrixXd original_transformation, VectorXd q){
    std::array<double, 7> q_arr;
    for (size_t i = 0; i < 7; i++) {
        q_arr[i] = q(i);
    }
    MatrixXd real_transform = franka_fk(q_arr);
    Matrix3d rotation_matrix = real_transform.block(0,0,3,3);
    Vector3d translation_vector = real_transform.block(0,3,3,1);
	// std::cout << "translation_vector: " << translation_vector << std::endl;
    Eigen::AngleAxisd angle_axis(rotation_matrix);
    // Get rotation angle (in radians)
    double angle = angle_axis.angle();
    // Get rotation axis (unit vector)
    Eigen::Vector3d axis = angle_axis.axis();

    
    Matrix3d original_rotation =  original_transformation.block(0,0,3,3);
    Vector3d original_translation = original_transformation.block(0,3,3,1);
    Eigen::AngleAxisd original_angle_axis(original_rotation);
    // Get rotation angle (in radians)
    double original_angle = original_angle_axis.angle();
    Eigen::Vector3d original_axis = original_angle_axis.axis();

    double translation_error = (translation_vector - original_translation).norm();
    Matrix3d R_error = (original_rotation.transpose())* rotation_matrix;
    Eigen::AngleAxisd R_error_quat(R_error);

    double error_angle = R_error_quat.angle();

    Eigen::Vector3d error_axis = R_error_quat.axis();
    double value = (R_error.trace() - 1) / 2;
    value = std::min(1.0, std::max(-1.0, value)); // Clamp to [-1,1]
    double angle_error = acos(value);

    double axis_error = (original_axis - axis).norm();

    Vector3d error_all;
    error_all << translation_error, error_angle, error_axis.norm();
    return error_all;
}

int main()
{    
	auto vi = DQ_VrepInterface();
	try {
    	if (!vi.connect("127.0.0.1", 19997,100,10))
        	throw std::runtime_error("Unable to connect to CoppeliaSim.");
    	vi.set_synchronous(true);
    	vi.start_simulation();
    	std::this_thread::sleep_for(std::chrono::milliseconds(100));
    	//----------------------------------------------------------
    	std::vector<std::string> jointnames = {"Franka_joint1", "Franka_joint2",
                                           	"Franka_joint3", "Franka_joint4",
                                           	"Franka_joint5", "Franka_joint6",
                                           	"Franka_joint7"};



		std::string file_name = "./config/parameter.yml";
		YAML::Node config = YAML::LoadFile(file_name);
		std::array<double, 3> desired_pos_array = config["data"]["desired_position"].as<std::array<double, 3>>();
  		Eigen::Map<Eigen::Vector3d> desired_pos(desired_pos_array.data());
  		//std::string path = config["data"]["joint_data_folder"].as<std::string>()+ config["data"]["task"].as<std::string>();
		//std::cout << "0000000000000" << std::endl;

		//std::string human_angle_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_pose.csv"; 
		std::string human_angle_path = "/media/panda/data/yaqi/CPPSimulationBasicPrimitives/Exp_2/simulation/push_side/push_side_2/modified_pose.csv"; 
		std::cout << "55555555555" << std::endl;
		std::cout << "human_angle_path: " << human_angle_path << std::endl;

  		
		MatrixXd human_joint_angle = load_csv(human_angle_path)/180 * pi;
		//path = config["data"]["wrist_data_folder"].as<std::string>() + config["data"]["task"].as<std::string>();
		//std::cout << "99999999999999" << std::endl;

		//std::string human_wrist_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_wrist.csv"; 
		std::string human_wrist_path = "/media/panda/data/yaqi/CPPSimulationBasicPrimitives/Exp_2/simulation/push_side/push_side_2/modified_wrist.csv";
		std::cout << "------------" << std::endl;

		MatrixXd wrist_traj = load_csv3(human_wrist_path);
		std::cout << "11111111111111" << std::endl;

		DQ_SerialManipulatorMDH robot = FrankaRobot::kinematics();		
		VectorXd q_initial  = vi.get_joint_positions(jointnames);
		DQ initial_pos_dq = robot.fkm(q_initial,5);

		DQ_SerialManipulatorMDH human = RealSenseRobot::kinematics();		
		VectorXd qh_initial  = human_joint_angle.col(0);
		DQ initial_pos_human_dq = human.fkm(qh_initial);
		std::cout << "222222222222" << std::endl;

		
		int size_all = wrist_traj.cols();

		double gain = config["data"]["gain"].as<double>();
		
		for(int i=0; i<wrist_traj.cols(); i++){
			wrist_traj.col(i) = (wrist_traj.col(i) - wrist_traj.col(0)); 
		}  

		wrist_traj.row(2) = -wrist_traj.row(2);
  		wrist_traj.row(1) = -wrist_traj.row(1);

		for(int i=0; i<wrist_traj.cols(); i++){
			wrist_traj.col(i) = (wrist_traj.col(i)  + initial_pos_dq.translation().vec3()); // robot desired cartesian translation 
		}  

		double step = 0;
		double error_angle = 9999999999999999999;
		//Vector4d error_full;
		Eigen::MatrixXd Jaco_ik(4, 7);
		Eigen::VectorXd error_full(Jaco_ik.rows());

		error_full.setConstant(1000);
		//VectorXd q_dot(7);
		Eigen::VectorXd q_dot(Jaco_ik.cols());
		q_dot.setZero();
		int n_human = 8;


		SparseMatrix<c_float> H_sparse;
    	SparseMatrix<c_float> A_sparse;
		VectorXd q_current  = vi.get_joint_positions(jointnames);

		while(true)
		{
			if (step>= human_joint_angle.cols())
			{
				step = human_joint_angle.cols() -1;
			}					
			// VectorXd q_current  = vi.get_joint_positions(jointnames);
			VectorXd q_robot_current = q_current;
			
			//manipulability tracking
			VectorXd q_human_current(8);
			q_human_current << 0, 0.0601534 ,  1.47194 ,0.0211931, 0,  -2.48352,  0 , 0;
		
			MatrixXd J_anal_human_current = human.pose_jacobian(q_human_current);
			MatrixXd J_geom_human_current = geomJac8(human, J_anal_human_current, q_human_current, n_human).block(3,0,3,8);
			MatrixXd ME_human_current = J_geom_human_current * J_geom_human_current.transpose();
			EigenSolver<MatrixXd> eigensolver_feasible_human(ME_human_current); // WHy different from Matlab eig()????
			auto EigenVec_human_feasible = eigensolver_feasible_human.eigenvalues();
			VectorXd EV_count_human = EigenVec_human_feasible.array().real();


			double human_ev = EV_count_human(0) * EV_count_human(1) * EV_count_human(2);
			
			// next timestep human manipulability ellipsoid
			VectorXd q_human_next =  human_joint_angle.col(step);
			MatrixXd J_anal_human_next = human.pose_jacobian(q_human_next);
			MatrixXd J_geom_human_next = geomJac8(human, J_anal_human_next, q_human_next, n_human).block(3,0,3,8);
			MatrixXd ME_human_next = J_geom_human_next * J_geom_human_next.transpose();

			// current timestep robot manipulability ellipsoid
			VectorXd q_robot_initial(7);
			q_robot_initial << -2.82151 ,-1.46356,  0.046521  , -3.06096 ,  2.65921 ,-0.0137487, 1.86731;
			DQ x_robot_current_6joint_dq = robot.fkm(q_robot_initial, 5);

			MatrixXd J_anal_robot_current = robot.pose_jacobian(q_robot_initial, 5);
			MatrixXd J_geom_robot_current= (geomJac6(robot, J_anal_robot_current, q_robot_initial, 6)).block(3, 0, 3, 6); 
			MatrixXd ME_robot_current = J_geom_robot_current * J_geom_robot_current.transpose();
			MatrixXd ME_robot_initial = ME_robot_current;
			
			EigenSolver<MatrixXd> eigensolver_feasible(ME_robot_current); // WHy different from Matlab eig()????
			auto EigenVec_robot_feasible = eigensolver_feasible.eigenvalues();
			VectorXd EV_count = EigenVec_robot_feasible.array().real();
		
			double robot_ev = EV_count(0) * EV_count(1) * EV_count(2);
			// Manipulability distance in manifold   
			// MatrixXd  ME_human_distance_next = logmap(ME_human_next, ME_human_current, 1);        
			// MatrixXd ME_human_distance_next = logmap(ME_human_next/0.0012 * robot_ev, ME_human_current/0.0012 * robot_ev, 1);
			MatrixXd ME_human_distance_next = logmap(ME_human_next/0.0012 * robot_ev, ME_human_current/0.0012 * robot_ev, 1);

			MatrixXd ME_robot_desire_distance_next = para_trans(ME_human_current, ME_robot_current, ME_human_distance_next);
			MatrixXd ME_robot_desire_next = expmap(ME_robot_desire_distance_next, ME_robot_current);
        
			// line to line jacobian

			DQ curr_pos_dq = robot.fkm(q_current,5);
			MatrixXd anal_J = robot.pose_jacobian(q_current,5);
			MatrixXd geom_J = geomJac6(robot, anal_J, q_current, 6).block(3,0,3,6);
			Vector3d curr_pos = curr_pos_dq.translation().vec3();
			
			// we should not use 7 joint. The 7th joint will always rotate itself to satisfy the desired direction, instead of rtate the full arm
			//Matrix<double, 4, 6> Jaco_ik;
			Eigen::MatrixXd Jaco_ik(4, 7);
			Jaco_ik.setZero();
			

			// rotate 90 degree
			DQ rot = cos(pi/4) + sin(pi/4) * (k_);

			curr_pos_dq = robot.fkm(q_current, 5);
			DQ current_ee_translation = curr_pos_dq.translation();


			DQ current_ee_direction = ((curr_pos_dq.rotation() * i_ * curr_pos_dq.rotation().conj()).normalize()).normalize(); 
			DQ current_line =  current_ee_direction + E_  * cross(current_ee_translation, current_ee_direction);


			DQ initial_ee_direction = ((initial_pos_dq.rotation() * i_ * initial_pos_dq.rotation().conj()).normalize()).normalize();

			DQ desired_human_rotation = initial_pos_human_dq.rotation().conj() * human.fkm(human_joint_angle.col(step)).rotation();
			DQ desired_human_direction = desired_human_rotation * i_ * desired_human_rotation.conj();
			DQ desired_line =  desired_human_direction + E_  * cross(DQ(current_ee_translation), desired_human_direction);

			// std::cout << "desired_human_direction" << desired_human_direction << std::endl;
			// std::cout << "error_angle" << error_angle << std::endl;
			// std::cout << "current_ee_direction: " << current_ee_direction << std::endl;

			
			MatrixXd line_jacobian = robot.line_jacobian(robot.pose_jacobian(q_current, 5), robot.fkm(q_current, 5), DQ(i_));
			
			error_angle = DQ_robotics::DQ_Geometry::line_to_line_angle(current_line, desired_line);
			error_angle = 1-1 * cos(error_angle);
			MatrixXd line_to_line_angle_jaco = robot.line_to_line_angle_jacobian(line_jacobian, current_line, desired_line);	
			if (line_to_line_angle_jaco.array().isNaN().any())
			{
				line_to_line_angle_jaco.setZero();
			}
			Jaco_ik.block(0,0,3,6) = geom_J;
			error_full.block(0,0,3,1) = (wrist_traj.col(step)  - (curr_pos_dq.translation()).vec3()) * gain;
			Jaco_ik.block(3,0,1,6)  = line_to_line_angle_jaco;
			error_full(3) = -( error_angle)  * gain;
			//q_dot.block(0,0,6,1) =  (Jaco_ik).completeOrthogonalDecomposition().pseudoInverse() * (error_full);
			Eigen::MatrixXd J_pinv = Jaco_ik.completeOrthogonalDecomposition().pseudoInverse(); // 7×4
			q_dot = J_pinv * error_full;   

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double error_now = -( error_angle) * 1; // (current_ee_rotation.vec3().transpose() * human_x_rotation.vec3()).value();     
        VectorXd q_velocity_desired = Jaco_ik.completeOrthogonalDecomposition().pseudoInverse() * error_full;
		std::cout << "q_velocity_desired: " << q_velocity_desired.transpose() << std::endl;
        // VectorXd q_velocity_desired = error_vel * line_to_line_dist_jac.completeOrthogonalDecomposition().pseudoInverse();
		VectorXd q_vel_desired7(7);
		q_vel_desired7.setZero();
		q_vel_desired7.block(0,0,6,1) = q_velocity_desired;

		double min = pi/180 * (-30);
		double max = pi/180 * 30;
        VectorXd q_desired = q_current + q_vel_desired7 * 0.01;

		std::cout << "desired -------------------------------" << std::endl;
		std::cout << "q_desired: " << q_desired.transpose() << std::endl;
		// vi.set_joint_positions(jointnames, q_desired);
		q_current = q_desired;


        std::random_device random_generator;
        std::mt19937 gen(random_generator());

        std::uniform_real_distribution<double> distribution(min, max);
        
        int jj = 0;
        double distance_to_mani = 999999;
		int SIZE = 100;
		VectorXd q_work;
		q_work.setZero();

        for (int i = 0; i < SIZE; ++i) {

			Eigen::Matrix4d TFE_for_all;
			TFE_for_all <<  0.70710678,  0.70710678,  0.,         0.,
						-0.70710678,  0.70710678,  0.,         0.,
						0.,          0.,          1.,         0.1034,
                		0.,          0.,          0.,         1.;


            double x_delta = distribution(gen); 
            double y_delta = distribution(gen);  
            double z_delta = distribution(gen);  


            // std::cout << "x_delta: " << x_delta << std::endl;
            // std::cout << "y_delta: " << y_delta << std::endl;
            // std::cout << "z_delta: " << z_delta << std::endl;
            
            DQ x_rotation = DQ(cos(x_delta/2), sin(x_delta/2), 0, 0);
            DQ y_rotation = DQ(cos(y_delta/2), 0, sin(y_delta/2), 0);
            DQ z_rotation = DQ(cos(z_delta/2), 0, 0, sin(z_delta/2));


            DQ rotation_1;
            rotation_1 = x_rotation.conj() * y_rotation * x_rotation;
            DQ rotation_2;
            rotation_2 = rotation_1 * z_rotation;

            DQ desired_rotation_robot_current = robot.fkm(q_desired).normalize().rotation();
        //    VectorXd quat_vec = (desired_rotation_robot_current * rotation_2).vec4();  
            VectorXd quat_vec = (desired_rotation_robot_current).vec4();  

            Eigen::Quaterniond quat(quat_vec(0), quat_vec(1), quat_vec(2), quat_vec(3));
            Eigen::Matrix3d rotationMatrix = quat.toRotationMatrix();

            // // Vector3d angles_ShoulderBase = rotationMatrix.eulerAngles(0, 1, 2);
            // // Eigen::AngleAxisd rollAngle(angles_ShoulderBase.row(2).value(), Eigen::Vector3d::UnitZ());
            // // Eigen::AngleAxisd yawAngle(0, Eigen::Vector3d::UnitY());
            // // Eigen::AngleAxisd pitchAngle(0, Eigen::Vector3d::UnitX());


            // // Eigen::Quaternion<double> rotm = rollAngle * yawAngle * pitchAngle;

			Eigen::Matrix<double, 4, 4> O_T_EE_matrix;

            O_T_EE_matrix.block(0,0,3,3) = rotationMatrix;
            O_T_EE_matrix.block(0,3,3,1) = robot.fkm(q_desired).normalize().translation().vec3();
			std::array<double, 7> q_des;
			for (size_t i = 0; i < 7; i++) {
				q_des[i] = q_desired(i);
			}
			O_T_EE_matrix = franka_fk(q_des);
			VectorXd rotation2_vec = (rotation_2).vec4(); 
			Eigen::Quaterniond rotation2(rotation2_vec(0), rotation2_vec(1), rotation2_vec(2), rotation2_vec(3));
            rotationMatrix = rotation2.toRotationMatrix();
			O_T_EE_matrix.block(0,0,3,3) = O_T_EE_matrix.block(0,0,3,3) * rotationMatrix;
            O_T_EE_matrix.block(0,3,3,1) = wrist_traj.col(step);

			// std::cout << "O_T_EE_matrix before : " << O_T_EE_matrix  << std::endl;

            // std::cout << "q_desired: " <<robot.fkm(q_desired).translation().vec3().transpose() << std::endl;

			// std::cout << "O_T_EE_matrix 5 before : " << franka_fk(q_des, 5)  << std::endl;

            // std::cout << "q_desired 5: " <<robot.fkm(q_desired, 5).translation().vec3().transpose() << std::endl;
			
			
			// std::cout << "O_T_EE_matrix: " << O_T_EE_matrix  * TFE_for_all << std::endl;
			Matrix<double, 4, 4> O_T_EE_matrix_cali_for_all = O_T_EE_matrix * TFE_for_all;
            // O_T_EE_matrix.block(3,3,1,1) << 1;
			// Matrix<double, 4, 4> O_T_EE_matrix_cali_for_all = O_T_EE_matrix;// * TFE_for_all;
			// std::cout << "desired translation after moving: " << O_T_EE_matrix_cali_for_all.block(0,3,3,1).transpose() << std::endl;
			// Matrix3d rotation_cali = O_T_EE_matrix_cali_for_all.block(0,0,3,3);
			// Eigen::Quaterniond quaternion(rotation_cali);
			// std::cout << "desired rotation after moving: " << quaternion << std::endl;

            // double* O_T_EE_vec = O_T_EE_matrix_cali_for_all.data();
            // std::array<double, 16> O_T_EE_ik_array;
            // std::copy(O_T_EE_vec, O_T_EE_vec + 16, O_T_EE_ik_array.begin()); 

            // double* q_robot_current_vec = q_robot_current.data();
            // // double* q_robot_current_vec =  q_robot_IK_solver.col(j-1).data();

            // std::array<double, 7> q_robot_current_array;
            // std::copy(q_robot_current_vec, q_robot_current_vec + 7, q_robot_current_array.begin());
            
            double q7 = q_robot_current[6];
            // std::array< std::array<double, 7>, 4> q_ik_solution;
            // // std::array<double, 7> q_ik_solution;
       

			Eigen::Matrix4d TFE;
			TFE <<  0.70710678,  0.70710678,  0.,         0.,
				-0.70710678,  0.70710678,  0.,         0.,
					0.,          0.,          1.,         0.1034,
					0.,          0.,          0.,         1.;
			Matrix<double, 4, 4> O_T_EE_matrix_cali = O_T_EE_matrix * TFE;
			std::array<double, 9> ROE;
			MatrixXd init_rot_trans = O_T_EE_matrix_cali.block(0,0,3,3).transpose();
			std::copy(init_rot_trans.data(), init_rot_trans.data() + 9, ROE.begin());
			// std::array<double,9> ROE = {-0.5878745000000002,-0.7990239702980599,-0.12635000000000007,-0.16720689999999983,0.2728366411870001,-0.9474186000000001,0.7914831,-0.5358366439507003,-0.293996};
			// std::array<double,9> ROE = {-0.5878745000000002,-0.7990239702980599,-0.12635000000000007,-0.16720689999999983,0.2728366411870001,-0.9474186000000001,0.7914831,-0.5358366439507003,-0.293996};
			std::array<double,3> r = {O_T_EE_matrix_cali.block(0,3,1,1).value(), O_T_EE_matrix_cali.block(1,3,1,1).value(), O_T_EE_matrix_cali.block(2,3,1,1).value()};       
			// std::vector<std::array<std::array<double, 6>,7>> Jsols;
			std::array<std::array<double,7>,8> qsols;

			// std::cout << "desired translation after moving: " << O_T_EE_matrix_cali.block(0,3,3,1).transpose() << std::endl;
			Matrix3d rotation_cali = O_T_EE_matrix_cali.block(0,0,3,3);
			Eigen::Quaterniond quaternion(rotation_cali);
			// std::cout << "desired rotation after moving: " << quaternion << std::endl;

			unsigned int nsols;

			double* initial_q_vec = q_desired.data();
			std::array<double, 7> initial_q_vec_array;
			std::copy(initial_q_vec, initial_q_vec + 7, initial_q_vec_array.begin());
			double swivel_theta = franka_swivel(initial_q_vec_array);
			// q_ik_solution = franka_ik_swivel(r, ROE, swivel_theta);


			nsols = franka_ik_q7_arr(r, ROE, q7, qsols);
            std::vector<std::array<double,7>> q_ik_solution;
            q_ik_solution = std::vector<std::array<double, 7>>(qsols.begin(), qsols.end());
			
			// for (const auto& arr : q_ik_solution) {
			// 	std::cout << "[ ";
			// 	for (double val : arr) {
			// 		std::cout << val << " ";
			// 	}
			// 	std::cout << "]" << std::endl;
			// }
            // std::cout << "time: " << duration << " nanoseconds" << std::endl;
              
            // }
			for (size_t i = 0; i < q_ik_solution.size(); ++i) {
                // std::cout << "Solution " << i + 1 << ": ";
                for (double val : q_ik_solution[i]) {
                    // std::cout << val << " ";
                }
                // std::cout << std::endl;
            }
            
            for (int jjj = 0; jjj<q_ik_solution.size(); jjj++){


                if  (!std::isnan(q_ik_solution[jjj][0]) & !std::isnan(q_ik_solution[jjj][1]) & !std::isnan(q_ik_solution[jjj][3])
                & !std::isnan(q_ik_solution[jjj][2]) & !std::isnan(q_ik_solution[jjj][4]) & !std::isnan(q_ik_solution[jjj][5])
                & !std::isnan(q_ik_solution[jjj][6])){
                    // if  (!std::isnan(q_ik_solution[0]) & !std::isnan(q_ik_solution[1]) & !std::isnan(q_ik_solution[3])
                    //     & !std::isnan(q_ik_solution[2]) & !std::isnan(q_ik_solution[4]) & !std::isnan(q_ik_solution[5])
                    //     & !std::isnan(q_ik_solution[6])){
					 Eigen::VectorXd vec_i = Eigen::Map<const Eigen::VectorXd>(q_ik_solution[jjj].data(), 7);
					Vector3d error_all = calculate_error(O_T_EE_matrix_cali_for_all, vec_i);
					
					double translation_error = error_all(0);
					double angle_error = error_all(1);
					double axis_error = error_all(2); 
					// std::cout << "translation_error: " << translation_error << std::endl;
					// std::cout << "angle_error: " << angle_error << std::endl;
					// std::cout << "axis_error: " << axis_error << std::endl;
					//////////////////////////////////////////////////////////////////
                    
                    VectorXd q_feasible(7);
                        
                    Eigen::Map< Eigen::Matrix<double, 7, 1> > q_ik_one_solution(q_ik_solution[jjj].data());
                    // Eigen::Map< Eigen::Matrix<double, 7, 1>> q_ik_one_solution(q_ik_solution.data());

                    q_feasible = q_ik_one_solution;

                    // abort();


                    MatrixXd J_anal = robot.pose_jacobian(q_feasible, 5);
                    MatrixXd J_geom = (geomJac6(robot, J_anal, q_feasible, 6)).block(3, 0, 3, 6); 
                    MatrixXd ME_robot_feasible = J_geom*J_geom.transpose();
                    MatrixXd ME_feasible_to_initial = logmap(ME_robot_feasible, ME_robot_initial,1);

                    MatrixXd ME_ik_dist;
					ME_ik_dist = logmap(ME_robot_feasible, ME_robot_initial, 1);

                     
                    MatrixXd ME_compare = logmap(ME_robot_desire_distance_next, ME_ik_dist, 1);
                    double cost =abs((ME_compare.norm())) +  0.5 * abs((q_robot_current.block(0,0,3,1) - q_feasible.block(0,0,3,1)).norm());//abs((ME_compare.norm())) + 0.5 * abs((q_robot_current - q_feasible).norm());
					
					// std::cout << "feasible: " << q_feasible.transpose() << std::endl;
					// std::cout << "feasible translation: " << robot.fkm(q_feasible).translation() << std::endl;
					// std::cout << "feasible rotation: " << robot.fkm(q_feasible).rotation() << std::endl;

					std::array<double, 7> q_arr;
					for (size_t i = 0; i < 7; i++) {
						q_arr[i] = q_feasible(i);
					}

				
					MatrixXd real_transform = franka_fk(q_arr);
					// std::cout << "real_transform: " << real_transform << std::endl;
					// abort();

					real_transform = real_transform * TFE_for_all.inverse();
					Matrix3d rotation_matrix = real_transform.block(0,0,3,3);
					Vector3d translation_vector = real_transform.block(0,3,3,1);
					// std::cout << "feasible translation_vector: " << translation_vector << std::endl;
					MatrixXd J_anal7 = robot.pose_jacobian(q_feasible);
                    MatrixXd J_geom7 = (geomJac7(robot, J_anal7, q_feasible, 7)); 
					// std::cout << "J_geom7: " <<  J_geom7 << std::endl;
					
					std::array<std::array<double,6>,7> Jsols;
					
					Jsols = Jacobian(q_arr, true);
					std::array<double, 42> Jsols_flat;
					for (size_t i = 0; i < 7; ++i)
						for (size_t j = 0; j < 6; ++j)
							Jsols_flat[i * 6 + j] = Jsols[i][j];

					// Now map it using Eigen
					Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> J(Jsols_flat.data());

					// std::cout << "J_geom7: " <<  J.transpose() << std::endl;
					// std::cout << "-----------------------------------------------------------" << std::endl;

					// std::cout << "cost: " << cost << std::endl;
					// std::cout << "q_feasible: " << q_feasible.transpose() << std::endl;
					// std::cout << "-----------------------------------------------------------" << std::endl;
					// vi.set_joint_positions(jointnames, q_feasible);
    				// std::this_thread::sleep_for(std::chrono::milliseconds(1000));




					


                    if (distance_to_mani > cost){
                        q_work = q_feasible;
                        distance_to_mani = cost;  
                    }
                    jj = jj + 1;
                }
            }
        }
		std::cout << "Step: " << step << std::endl;
		if (q_work.norm() < 0.1) {
			q_work = q_robot_current;
		}
		q_work = q_desired;
		std::cout << "work ----------------------------- " << std::endl;
 		std::cout << "q_desired: " << q_desired.transpose() << std::endl;
 		std::cout << "q_work: " << q_work.transpose() << std::endl;
		std::cout << "cost final: " << (q_work.block(0,0,3,1) - q_desired.block(0,0,3,1)).norm() << std::endl;
		// vi.set_joint_positions(jointnames, q_work);
    	// std::this_thread::sleep_for(std::chrono::milliseconds(1000));
		

		// abort();
		// q_current = q_work;
		// abort();






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			Matrix<double, 6, 6> H_zero;
			H_zero.setZero();
			Matrix<c_float, 4, 1> lb;
			Matrix<c_float, 4, 1> ub;
			Matrix<c_float, 4, 6> A;
			Matrix<c_float, 4, 1> lb_lower;
			Matrix<c_float, 4, 1> ub_lower;
			lb_lower << -0.03, -0.03, -0.03, 0;
			ub_lower << 0.03, 0.03, 0.03, 0;

			lb = error_full.block(0,0,4,1);// + lb_lower;
			ub = error_full.block(0,0,4,1);// + ub_lower;
			A = Jaco_ik;
			// std::cout << "step: " << step << std::endl;

			// tracking human manip
			Matrix<double, 6, 6> weight_identity;
			weight_identity.setIdentity();

			// std::cout << "11111111111111111111111111111" << std::endl;

			MatrixXd J_geom_robot = (geomJac6(robot, J_anal_robot_current, q_robot_current, 6)).block(0, 0, 6, 6);
			// std::cout << "---------------------------" << std::endl;

			MatrixXd J_anal_desired = robot.pose_jacobian(q_work, 5);

			// std::cout << "99999999999999999999999999999999999999999999999999999999" << std::endl;

			MatrixXd J_geom_desired = (geomJac6(robot, J_anal_desired, q_work, 6)).block(0, 0, 6, 6);

			// std::cout << "88888888888888888888888888888888888" << std::endl;


			MatrixXd ME_robot = J_geom_robot * J_geom_robot.transpose();
			// std::cout << "1111111111111" << std::endl;

			MatrixXd ME_desired = J_geom_desired * J_geom_desired.transpose();
			// std::cout << "2222222222222222" << std::endl;

			MatrixXd ME_dist_human = logmap(ME_desired, ME_robot, 1);
			// std::cout << "3333333333333333333" << std::endl;

			Tensor<double, 3> J_gradient_robot3D(3,6,6);
			// std::cout << "3444444444444444" << std::endl;

			J_gradient_robot3D = jacobianEst6D(q_robot_current, 6, robot);
			// std::cout << "55555555555555" << std::endl;

			MatrixXd Jacobian_Mani3D_pt = redManipulabilityJacobian6D(J_geom_robot, J_gradient_robot3D);

			// std::cout << "222222222222222222222222222222222222222222" << std::endl;



			VectorXd vectorization_ME_robot_desire_distance_pt = spd2vec_vec(ME_dist_human);

			c_float K_qp = 2; 
			Matrix<double, 1, 6> f_zero;
			f_zero.setZero();
			
			Matrix<c_float, 6, 6> H =  (Jacobian_Mani3D_pt.transpose() * Jacobian_Mani3D_pt) +  0.1 * weight_identity;//Jaco_ik.transpose() * Jaco_ik ;
			Matrix<c_float, 1, 6> f = - K_qp * (vectorization_ME_robot_desire_distance_pt.transpose() * Jacobian_Mani3D_pt);//  - K_qp * error_full.transpose() * Jaco_ik;

			H_sparse = H.sparseView();
			H_sparse.pruned(1e-9);

			// std::cout << "333333333333333333333333333333333333333333333333" << std::endl;


			lb = error_full.block(0,0,4,1) + lb_lower;
			ub = error_full.block(0,0,4,1) + ub_lower;
			A = Jaco_ik;

			A_sparse = A.sparseView();			

			OsqpEigen::Solver solver;
			solver.settings()->setVerbosity(false); // print outptu or not
			solver.settings()->setAlpha(1.5); // ADMM relaxation parameter/step size/penalty parameter
			solver.data()->setNumberOfVariables(6);
			//eigenvalue (6) + robot_x_diffing(3) + aixs tracking(1) + limits(7)
			solver.data()->setNumberOfConstraints(4); 
			solver.data()->setHessianMatrix(H_sparse);
			solver.data()->setGradient(f.transpose());
			solver.data()->setLinearConstraintsMatrix(A_sparse);
			solver.data()->setLowerBound(lb);
			solver.data()->setUpperBound(ub);
			solver.initSolver();
			solver.solveProblem();
			q_dot.block(0,0,6,1) = solver.getSolution();
			VectorXd q_pos = q_dot * 0.01 + q_current;
			q_current = q_pos;
			q_current = q_desired;
		
			// vi.set_joint_target_velocities(jointnames, q_dot);
			// vi.trigger_next_simulation_step();
			// vi.wait_for_simulation_step_to_end();
			std::cout << "final---------------------" << std::endl;
			vi.set_joint_positions(jointnames, q_desired);
	    	// std::this_thread::sleep_for(std::chrono::milliseconds(1000));

			if(190 <= step) {
				abort();
			}
			step = step + 1;
		


		}
		std::cout << "finished" << std::endl;
    	//---------------------------------------------------------
    	vi.stop_simulation();
    	vi.disconnect();
		


	} catch (std::exception& e) {
    	std::cout<<e.what()<<std::endl;
    	vi.stop_simulation();
    	vi.disconnect();
	}
	return 0;
}
