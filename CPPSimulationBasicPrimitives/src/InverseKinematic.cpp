#include "../include/InverseKinematic.h"
#include <iostream>
#include <dqrobotics/interfaces/vrep/DQ_VrepInterface.h>
#include <dqrobotics/utils/DQ_Geometry.h>

#include <dqrobotics/robot_control/DQ_PseudoinverseController.h>
#include <thread>
#include "yaml-cpp/yaml.h"
#include<dqrobotics/utils/DQ_Geometry.h>
#include<dqrobotics/utils/DQ_Constants.h>
#include <fstream> 
#include <iomanip> 



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
		
		// 
		std::string file_name = "../config/parameter.yml";
		// std::cout<< "before yaml\n";
		YAML::Node config = YAML::LoadFile(file_name);
		// std::cout<<"after yaml\n";
		std::array<double, 3> desired_pos_array = config["data"]["desired_position"].as<std::array<double, 3>>();
  		Eigen::Map<Eigen::Vector3d> desired_pos(desired_pos_array.data());
  		//std::string path = config["data"]["joint_data_folder"].as<std::string>()+ config["data"]["task"].as<std::string>();
		//std::string human_angle_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_pose.csv"; 
		std::string human_angle_path = "/media/panda/data/yaqi/CPPSimulationBasicPrimitives/Exp_2/simulation/push_side/push_side_2/modified_pose.csv"; 

		MatrixXd human_joint_angle = load_csv(human_angle_path)/180 * pi;
		//path = config["data"]["wrist_data_folder"].as<std::string>() + config["data"]["task"].as<std::string>();
		//std::string human_wrist_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_wrist.csv"; 
		std::string human_wrist_path = "/media/panda/data/yaqi/CPPSimulationBasicPrimitives/Exp_2/simulation/push_side/push_side_2/modified_wrist.csv";
		MatrixXd wrist_traj = load_csv3(human_wrist_path);

		int n_human = 8;  // Huma
		size_t col = human_joint_angle.cols();
		MatrixXd human_joint_angle_next(n_human, col); // Joint next step angles
		human_joint_angle_next.block(0, 0, n_human, col-1) = human_joint_angle.block(0, 1, n_human, col-1);
		human_joint_angle_next.col(col-1) = human_joint_angle.col(col-1);

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
		//Vector4d error_full;
		Eigen::MatrixXd Jaco_ik(4,7);
		Eigen::VectorXd error_full(Jaco_ik.rows());
		error_full.setConstant(1000);
		//VectorXd q_dot(7);
		Eigen::VectorXd q_dot(Jaco_ik.cols());
		q_dot.setZero();

		// [ADD] Open CSV for robot translational manipulability (upper triangle of 3x3 SPD)
		std::ofstream robot_me_csv("robot_manipulability.csv");  
		robot_me_csv.setf(std::ios::fixed);                      
		robot_me_csv << std::setprecision(8);                    
		//robot_me_csv << "step,"
        //     << "m11,m12,m13,"
        //     << "m21,m22,m23,"
        //     << "m31,m32,m33\n";       

		std::ofstream human_me_csv("human_manipulability.csv");   // [ADD]
		human_me_csv.setf(std::ios::fixed);                       // [ADD]
		human_me_csv << std::setprecision(8);                     // [ADD]

		int a = 0;




		
		while(a == 0)
		{
			//if (step>= human_joint_angle.cols())
			//{
			//	step = human_joint_angle.cols() -1;
			//}
			// ADD
			if (step >= human_joint_angle.cols()) { a = 1; continue; }
					
			VectorXd q_current  = vi.get_joint_positions(jointnames);

			DQ curr_pos_dq = robot.fkm(q_current,5);
			MatrixXd anal_J = robot.pose_jacobian(q_current,5);
			MatrixXd geom_J = geomJac6(robot, anal_J, q_current, 6).block(3,0,3,6);

            // [ADD] Robot translational manipulability: M_r = J_v J_v^T (3x3, SPD)
			Eigen::Matrix3d M_r = geom_J * geom_J.transpose(); 
			robot_me_csv << M_r(0,0) << "," << M_r(0,1) << "," << M_r(0,2) << ","
             << M_r(1,0) << "," << M_r(1,1) << "," << M_r(1,2) << ","
             << M_r(2,0) << "," << M_r(2,1) << "," << M_r(2,2) << "\n";

			const int step_i = static_cast<int>(step);                      // [ADD] (if not already present)

			// --- HUMAN manipulability  ---
			MatrixXd J_an_h  = human.pose_jacobian(human_joint_angle_next.col(step_i));          // [ADD]
			MatrixXd J_lin_h = geomJac6(human, J_an_h, human_joint_angle_next.col(step_i), 6)    // [ADD]
								  .block(3,0,3,J_an_h.cols());                                    // [ADD]
			Eigen::Matrix3d M_h = J_lin_h * J_lin_h.transpose();                                  // [ADD]
			 
			// write 9 entries (row-major), same format as robot
			human_me_csv << M_h(0,0) << "," << M_h(0,1) << "," << M_h(0,2) << ","
						 << M_h(1,0) << "," << M_h(1,1) << "," << M_h(1,2) << ","
						 << M_h(2,0) << "," << M_h(2,1) << "," << M_h(2,2) << "\n";              // [ADD]




			Vector3d curr_pos = curr_pos_dq.translation().vec3();

			// we should not use 7 joint. The 7th joint will always rotate itself to satisfy the desired direction, instead of rtate the full arm
			//Matrix<double, 4, 6> Jaco_ik;
			//Eigen::MatrixXd Jaco_ik(4,7);
			Jaco_ik.setZero();

			// rotate 90 degree
			DQ rot = cos(pi/4) + sin(pi/4) * (k_);

			curr_pos_dq = robot.fkm(q_current, 5);
			DQ current_ee_translation = curr_pos_dq.translation();

			DQ current_ee_direction = ((curr_pos_dq.rotation() * k_ * curr_pos_dq.rotation().conj()).normalize()).normalize(); 
			DQ current_line =  current_ee_direction + E_  * cross(current_ee_translation, current_ee_direction);
			DQ initial_ee_direction = ((initial_pos_dq.rotation() * k_ * initial_pos_dq.rotation().conj()).normalize()).normalize(); 
			DQ desired_human_rotation = initial_pos_human_dq.rotation().conj() * human.fkm(human_joint_angle_next.col(step_i)).rotation();
			DQ desired_human_direction = desired_human_rotation * k_ * desired_human_rotation.conj();
			
			DQ human_x_dq = human.fkm(human_joint_angle_next.col(step_i));
			DQ desired_line =  desired_human_direction + E_  * cross(human_x_dq.translation(), desired_human_direction);
			MatrixXd line_jacobian = robot.line_jacobian(robot.pose_jacobian(q_current, 5), robot.fkm(q_current, 5), DQ(k_));
			error_angle = DQ_robotics::DQ_Geometry::line_to_line_angle(current_line, desired_line);
			error_angle = 1-1 * cos(error_angle);
			MatrixXd line_to_line_angle_jaco = robot.line_to_line_angle_jacobian(line_jacobian, current_line, desired_line);	
			if (line_to_line_angle_jaco.array().isNaN().any())
			{
				line_to_line_angle_jaco.setZero();
			}

			// Task:: check human direction!!!

			std::cout << "desired_human_direction" << desired_human_direction << std::endl;
			std::cout << "error_angle" << error_angle << std::endl;
			std::cout << "current_ee_direction: " << current_ee_direction << std::endl;
			// abort();
			
			Jaco_ik.block(0,0,3,6) = geom_J;
			error_full.block(0,0,3,1) = (wrist_traj.col(step_i)  - (curr_pos_dq.translation()).vec3()) * gain;
			Jaco_ik.block(3,0,1,6)  = line_to_line_angle_jaco;
			error_full(3) = -(error_angle);
			//q_dot.block(0,0,6,1) =  (Jaco_ik).completeOrthogonalDecomposition().pseudoInverse() * (error_full);
			// compute pseudo‐inverse into its own matrix
            Eigen::MatrixXd J_pinv = Jaco_ik.completeOrthogonalDecomposition().pseudoInverse();
			
			std::cerr
			  << "Jaco_ik dims: "    << Jaco_ik.rows()    << "×" << Jaco_ik.cols()    << ", "
			  << "J_pinv dims: "    << J_pinv.rows()    << "×" << J_pinv.cols()    << ", "
			  << "error_full dims: " << error_full.rows() << "×" << error_full.cols()
			  << std::endl;

			Eigen::VectorXd q_dot(J_pinv.rows());  // size = number of rows of J_pinv (i.e. #columns of Jaco_ik)
			q_dot = J_pinv * error_full;


			vi.set_joint_target_velocities(jointnames, q_dot);
			vi.trigger_next_simulation_step();
			vi.wait_for_simulation_step_to_end();
			step = step + 1;
		}
		std::cout << "finished" << std::endl;
		//ADD
		if (robot_me_csv.is_open()) robot_me_csv.close();  // [ADD]
		if (human_me_csv.is_open()) human_me_csv.close();  //ADD
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
