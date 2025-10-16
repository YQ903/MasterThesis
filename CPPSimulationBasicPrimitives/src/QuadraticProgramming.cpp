#include "../include/QuadraticProgramming.h"
#include <iostream>
#include <dqrobotics/interfaces/vrep/DQ_VrepInterface.h>
#include <dqrobotics/utils/DQ_Geometry.h>

#include <dqrobotics/robot_control/DQ_PseudoinverseController.h>
#include <thread>
#include "yaml-cpp/yaml.h"
#include<dqrobotics/utils/DQ_Geometry.h>
#include<dqrobotics/utils/DQ_Constants.h>



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
  		std::string path = config["data"]["joint_data_folder"].as<std::string>()+ config["data"]["task"].as<std::string>();
		std::string human_angle_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_pose.csv"; 
  		
		MatrixXd human_joint_angle = load_csv(human_angle_path)/180 * pi;
		path = config["data"]["wrist_data_folder"].as<std::string>() + config["data"]["task"].as<std::string>();
		std::string human_wrist_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_wrist.csv"; 
		MatrixXd wrist_traj = load_csv3(human_wrist_path);

		DQ_SerialManipulatorMDH robot = FrankaRobot::kinematics();		
		VectorXd q_initial  = vi.get_joint_positions(jointnames);
		DQ initial_pos_dq = robot.fkm(q_initial,5);

		DQ_SerialManipulatorMDH human = RealSenseRobot::kinematics();		
		VectorXd qh_initial  = human_joint_angle.col(0);
		DQ initial_pos_human_dq = human.fkm(qh_initial);

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
		Vector4d error_full;
		error_full.setConstant(1000);
		VectorXd q_dot(7);
		q_dot.setZero();

		SparseMatrix<c_float> H_sparse;
    	SparseMatrix<c_float> A_sparse;


		while(true)
		{
			if (step>= human_joint_angle.cols())
			{
				step = human_joint_angle.cols() -1;
			}					
			VectorXd q_current  = vi.get_joint_positions(jointnames);

			DQ curr_pos_dq = robot.fkm(q_current,5);
			MatrixXd anal_J = robot.pose_jacobian(q_current,5);
			MatrixXd geom_J = geomJac6(robot, anal_J, q_current, 6).block(3,0,3,6);
			Vector3d curr_pos = curr_pos_dq.translation().vec3();

			// we should not use 7 joint. The 7th joint will always rotate itself to satisfy the desired direction, instead of rtate the full arm
			Matrix<double, 4, 6> Jaco_ik;
			Jaco_ik.setZero();

			// rotate 90 degree
			DQ rot = cos(pi/4) + sin(pi/4) * (k_);

			curr_pos_dq = robot.fkm(q_current, 5);
			DQ current_ee_translation = curr_pos_dq.translation();

			DQ current_ee_direction = ((curr_pos_dq.rotation() * k_ * curr_pos_dq.rotation().conj()).normalize()).normalize(); 
			DQ current_line =  current_ee_direction + E_  * cross(current_ee_translation, current_ee_direction);


			DQ initial_ee_direction = ((initial_pos_dq.rotation() * k_ * initial_pos_dq.rotation().conj()).normalize()).normalize();

			DQ desired_human_rotation = initial_pos_human_dq.rotation().conj() * human.fkm(human_joint_angle.col(step)).rotation();
			DQ desired_human_direction = desired_human_rotation * k_ * desired_human_rotation.conj();
			DQ desired_line =  desired_human_direction + E_  * cross(DQ(current_ee_translation), desired_human_direction);

			std::cout << "desired_human_direction" << desired_human_direction << std::endl;
			std::cout << "error_angle" << error_angle << std::endl;
			std::cout << "current_ee_direction: " << current_ee_direction << std::endl;

			
			MatrixXd line_jacobian = robot.line_jacobian(robot.pose_jacobian(q_current, 5), robot.fkm(q_current, 5), DQ(k_));
			
			error_angle = DQ_robotics::DQ_Geometry::line_to_line_angle(current_line, desired_line);
			error_angle = 1-1 * cos(error_angle);
			MatrixXd line_to_line_angle_jaco = robot.line_to_line_angle_jacobian(line_jacobian, current_line, desired_line);	
			if (line_to_line_angle_jaco.array().isNaN().any())
			{
				line_to_line_angle_jaco.setZero();
			}

			Matrix<double, 6, 6> H_zero;
			H_zero.setIdentity();
			Matrix<c_float, 6, 6> H = H_zero;// + numbla * k_for_slack;// + numbla* Identity_6_dim;// + numbla * (J_geom_robot_current6D7Joint.block(3,0,3,6)).transpose() * J_geom_robot_current6D7Joint.block(3,0,3,6);// + numbla * Identity_6_dim;// + 0.0001 * Identity_9_dim.block(0,0,6,7)/ (dt * dt);// + 0.000005 * Identity_9_dim.block(0,0,6,7) / (dt * dt);
			H_sparse = H.sparseView();
			H_sparse.pruned(1e-9); // set those who smaller than 0.01 as zero

			Matrix<double, 1, 6> f_zero;
			f_zero.setZero();
			Matrix<c_float, 1, 6> f  = f_zero;
			// Constraints:
			Matrix<c_float, 4, 1> lb;
			Matrix<c_float, 4, 1> ub;
			Matrix<c_float, 4, 6> A;

			Jaco_ik.block(0,0,3,6) = geom_J;
			error_full.block(0,0,3,1) = (wrist_traj.col(step)  - (curr_pos_dq.translation()).vec3()) * gain;
			Jaco_ik.block(3,0,1,6)  = line_to_line_angle_jaco;
			error_full(3) = -(error_angle) * 1;
			std::cout << "error_angle: " << error_angle << std::endl;


			lb = error_full.block(0,0,4,1);
			ub = error_full.block(0,0,4,1);
			A = Jaco_ik;
			std::cout << "step: " << step << std::endl;

			// let 7th joint be consistent with 6th joint
			// DQ current_ee_direction_k = ((curr_pos_dq.rotation() * -j_ * curr_pos_dq.rotation().conj()).normalize()).normalize(); 
			// DQ current_line_k =  current_ee_direction_k + E_  * cross(current_ee_translation, current_ee_direction_k);

			// DQ desired_5th_k_rotation = robot.fkm(q_current, 4).rotation();
			// DQ desired_5th_k_direction = desired_5th_k_rotation * k_ * desired_5th_k_rotation.conj();
			// DQ desired_line_k =  desired_5th_k_direction + E_  * cross(DQ(current_ee_translation), desired_5th_k_direction);

			// std::cout << "k_current_ee_direction" << current_ee_direction_k << std::endl;
			// std::cout << "K_desired_7th_k_direction: " << desired_5th_k_direction << std::endl;
			// abort();
			
			// MatrixXd line_jacobian_k = robot.line_jacobian(robot.pose_jacobian(q_current, 5), robot.fkm(q_current, 5), DQ(-j_));
			
			// double error_angle_k = DQ_robotics::DQ_Geometry::line_to_line_angle(current_line_k, desired_line_k);
			// error_angle_k = 1-1 * cos(error_angle_k);
			// MatrixXd line_to_line_angle_jaco_k = robot.line_to_line_angle_jacobian(line_jacobian_k, current_line_k, desired_line_k);	
			// if (line_to_line_angle_jaco_k.array().isNaN().any())
			// {
			// 	line_to_line_angle_jaco_k.setZero();
			// }



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
		

			vi.set_joint_target_velocities(jointnames, q_dot);
			vi.trigger_next_simulation_step();
			vi.wait_for_simulation_step_to_end();
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
