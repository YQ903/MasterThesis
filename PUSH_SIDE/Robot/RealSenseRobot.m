%% Contributors to this file:
%     Yuhe Gong - yuhe.gong@nottingham.ac.uk

%% File:
%     Input: the real sense shoulder joints q1, q2, q3
%            elbow flexion joint: q4
%            zero joint: q5, q6, the wrist joint always set 0
%     Base frame/ Reference frame: same as Franka Emika Robots

%% Code: 
classdef RealSenseRobot
    methods (Static)
        function real_sense_kinematics = kinematics()

          
            % D-H of real_sense: alpha before theta
            % shoulder xyz at 234, elbow at 6
            %real_sense_DH_theta   = [pi/2, 0,      pi/2,   pi/2,-pi/2,  0,      -pi/2,      0, 0,    0, 0];
            %real_sense_DH_d       = [0,    0,      0,      0,   0,      0,      0,          0, 0,    0, 0];
            %real_sense_DH_a       = [0,    0,      0,      0,   0,      0.2443, 0.2566,     0, 0,    0, 0];
            %real_sense_DH_alpha   = [-pi/2,pi/2,   pi/2,   pi/2,pi/2,   pi/2,      pi,      0, 0, pi/2, 0];

            real_sense_DH_theta   = [pi/2, 0,      pi/2,   pi/2,-pi/2,  0,      -pi/2,      0];
            real_sense_DH_d       = [0,    0,      0,      0,   0,      0,      0,          0];
            real_sense_DH_a       = [0,    0,      0,      0,   0,      0.2443, 0.2566,     0];
            real_sense_DH_alpha   = [-pi/2,pi/2,   pi/2,   pi/2,pi/2,   pi/2,      pi,      0];



            
            real_sense_DH_type    = repmat(DQ_SerialManipulatorMDH.JOINT_ROTATIONAL,1,8);
            
            real_sense_DH_matrix = [ real_sense_DH_theta;
                                     real_sense_DH_d;
                                     real_sense_DH_a;
                                     real_sense_DH_alpha;
                                     real_sense_DH_type;]; % D&H parameters matrix for the arm model
            
            real_sense_kinematics = DQ_SerialManipulatorMDH(real_sense_DH_matrix); % Defines robot model using dual quaternions
%           
            xb = 1 + DQ.E * 0.5 * DQ([0, 0, 0, 0]);
            xc = 1 + DQ.E * 0.5 * DQ([0, 0, 0, 0.0107]); 
            real_sense_kinematics.set_reference_frame(xb);
            real_sense_kinematics.set_base_frame(xb);

        end
    end    
end