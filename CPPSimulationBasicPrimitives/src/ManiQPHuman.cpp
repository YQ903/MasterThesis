#include "../include/ManiQPHuman.h"
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
  		//std::string path = config["data"]["joint_data_folder"].as<std::string>()+ config["data"]["task"].as<std::string>();
		//std::string human_angle_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_pose.csv"; 
		std::string human_angle_path = "/media/panda/data/yaqi/CPPSimulationBasicPrimitives/Exp_2/simulation/push_side/push_side_2/modified_pose.csv"; 

		MatrixXd human_joint_angle = load_csv(human_angle_path)/180 * pi;
		//path = config["data"]["wrist_data_folder"].as<std::string>() + config["data"]["task"].as<std::string>();
		//std::string human_wrist_path = "/home/yuhe/Desktop/project/1-GettingReadyToBeTaught/HumanDataProcessingOffline/data/3-applied_data/" + path + "/modified_wrist.csv"; 
		std::string human_wrist_path = "/media/panda/data/yaqi/CPPSimulationBasicPrimitives/Exp_2/simulation/push_side/push_side_2/modified_wrist.csv";
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
		int n_human = 8;


		SparseMatrix<c_float> H_sparse;
    	SparseMatrix<c_float> A_sparse;


		while(true)
		{
			if (step>= human_joint_angle.cols())
			{
				step = human_joint_angle.cols() -1;
			}					
			VectorXd q_current  = vi.get_joint_positions(jointnames);
			
			//manipulability tracking
			VectorXd q_human_current(8);
			q_human_current <<0, 0.0601534 ,  1.47194 ,0.0211931, 0,  -2.48352,  0 , 0;
		
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
			VectorXd q_robot_current(7);
			q_robot_current << -2.82151 ,-1.46356,  0.046521  , -3.06096 ,  2.65921 ,-0.0137487, 1.86731;
			DQ x_robot_current_6joint_dq = robot.fkm(q_robot_current, 5);

			MatrixXd J_anal_robot_current = robot.pose_jacobian(q_robot_current, 5);
			MatrixXd J_geom_robot_current= (geomJac6(robot, J_anal_robot_current, q_robot_current, 6)).block(3, 0, 3, 6); 
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
			Matrix<double, 4, 6> Jaco_ik;
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

			std::cout << "desired_human_direction" << desired_human_direction << std::endl;
			std::cout << "error_angle" << error_angle << std::endl;
			std::cout << "current_ee_direction: " << current_ee_direction << std::endl;

			
			MatrixXd line_jacobian = robot.line_jacobian(robot.pose_jacobian(q_current, 5), robot.fkm(q_current, 5), DQ(i_));
			
			error_angle = DQ_robotics::DQ_Geometry::line_to_line_angle(current_line, desired_line);
			error_angle = 1-1 * cos(error_angle);
			MatrixXd line_to_line_angle_jaco = robot.line_to_line_angle_jacobian(line_jacobian, current_line, desired_line);	
			if (line_to_line_angle_jaco.array().isNaN().any())
			{
				line_to_line_angle_jaco.setZero();
			}



			Matrix<double, 6, 6> H_zero;
			H_zero.setZero();
			Matrix<c_float, 4, 1> lb;
			Matrix<c_float, 4, 1> ub;
			Matrix<c_float, 4, 6> A;
			Matrix<c_float, 4, 1> lb_lower;
			Matrix<c_float, 4, 1> ub_lower;
			lb_lower << -0.0, -0.0, -0.0, 0;
			ub_lower << 0.0, 0.0, 0.0, 0;

			Jaco_ik.block(0,0,3,6) = geom_J;
			error_full.block(0,0,3,1) = (wrist_traj.col(step)  - (curr_pos_dq.translation()).vec3()) * gain;
			Jaco_ik.block(3,0,1,6)  = line_to_line_angle_jaco;
			error_full(3) = -( error_angle) * 1;
			std::cout << "error_angle: " << error_angle << std::endl;
	
			lb = error_full.block(0,0,4,1);// + lb_lower;
			ub = error_full.block(0,0,4,1);// + ub_lower;
			A = Jaco_ik;
			std::cout << "step: " << step << std::endl;

			// tracking human manip
			Matrix<double, 6, 6> weight_identity;
			weight_identity.setIdentity();

			MatrixXd ME_robot_begin_human = J_geom_robot_current * J_geom_robot_current.transpose();
			MatrixXd ME_dist_human = logmap(ME_human_next/0.0012 * robot_ev, ME_robot_begin_human, 1);
			Tensor<double, 3> J_gradient_robot3D(3,6,6);
			J_gradient_robot3D = jacobianEst3D(q_robot_current, 6, robot);
			MatrixXd Jacobian_Mani3D_pt = redManipulabilityJacobian3D(J_geom_robot_current, J_gradient_robot3D);


			VectorXd vectorization_ME_robot_desire_distance_pt = spd2vec_vec(ME_dist_human);

			c_float K_qp = 2; 
			Matrix<double, 1, 6> f_zero;
			f_zero.setZero();
			
			Matrix<c_float, 6, 6> H =  (Jacobian_Mani3D_pt.transpose() * Jacobian_Mani3D_pt) +  0.01 * weight_identity;//Jaco_ik.transpose() * Jaco_ik ;
			Matrix<c_float, 1, 6> f = - K_qp * (vectorization_ME_robot_desire_distance_pt.transpose() * Jacobian_Mani3D_pt);//  - K_qp * error_full.transpose() * Jaco_ik;

			H_sparse = H.sparseView();
			H_sparse.pruned(1e-9);

			lb = error_full.block(0,0,4,1);// + lb_lower;
			ub = error_full.block(0,0,4,1);// + ub_lower;
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
