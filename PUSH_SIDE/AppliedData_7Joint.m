close all
clear

addpath('Functions/DQ_robotics_matlab/');
scriptFullPath = mfilename('fullpath');
[projectRoot,~,~] = fileparts(scriptFullPath);
addpath( genpath( fullfile(projectRoot,'Functions','DQ_robotics_matlab') ) );
rehash toolboxcache

addpath('Functions/DQ_robotics_matlab/utils/');
addpath('Functions/mani_qp_plot_matlab/');
addpath('Functions/mani_qp_plot_matlab/function');
addpath('Robot/');
include_namespace_dq
% data_name = "push_side_1";
data_name = "Exp_2";
M_all = csvread("data_2/2-processed_data/" + data_name + "/all_angles.csv");
position_all = csvread("data_2/2-processed_data/" + data_name + "/all_cartesian_translation.csv");
%BaseWrist = csvread("data_2/2-processed_data/" + data_name + "/Base_to_Wrist.csv");
rsr = RealSenseRobotTrans7DoF;
fk = rsr.kinematics();
position = position_all;
M = M_all(:,1:11);

ee_position = position(:,end-2:end);


%initial vector wrist to hand
v_start = position(1, 10:12) - position(1, 7:9);
v_start = v_start/norm(v_start);

%vector wrist to hand from frame 2
v_end = position(2:end, 10:12) - position(2:end, 7:9);
rowendNorms = vecnorm(v_end,2,2);
v_end = v_end ./ rowendNorms;

len = 1;

elbow_shouder = position(:, 1:3) - position(:, 4:6);
L1 = sqrt(sum(elbow_shouder.^2,2));


wrist_elbow = position(:,7:9) - position(:,4:6);
L2 =  sqrt(sum(wrist_elbow.^2,2));
s = size(M);
s1 = s(1);

x = zeros(s1, 3);
mani = zeros(s1, 3);

offset = [0,180,180,180,0,90, 0,0,0,0,0] - M(1,:);

M(1, 4) = M(1, 4) + 90;
c_initial = M(1, :)/180 * pi;


M_7joints = zeros(s1, 7); 

M_7joints(:,1:3) = M(:,2:4);
M_7joints(:,4) = M(:,6);
%M_7joints(:,5:7) = 0;
M_7joints(:,5:7) = M(:,9:11);
M_7joints(:, :)= M_7joints(:, :)/180 * pi;
  
%transfer the initial wrist-to-hand vector back to base
R = eul2rotm([M_7joints(1,1), M_7joints(1,2), M_7joints(1,3)+ 90/180 * pi], 'XYZ');   
angles_yzx = rotm2eul(R, 'YZX');
M_7joints(1,1:3) = angles_yzx(1:3);

%FK to get the ee_pose   
c = M_7joints(1, :); 
ee = fk.fkm(c);
a_dq         = ee; 
%extract rotation and reverse the initial wrist-to-hand vector to base
q8           = vec8(a_dq); 
qvec = q8(1:4); 
q = quaternion(qvec(1), qvec(2), qvec(3), qvec(4));  
q = normalize(q);
base_vec = rotatepoint(conj(q), v_start);  


for j = 2:s1
    R = eul2rotm([M_7joints(j,1), M_7joints(j,2), M_7joints(j,3)+ 90/180 * pi], 'XYZ');
    %elbowplot(j,:) = R *[0.0, 0.0, 0.24]';
    angles_yzx = rotm2eul(R, 'YZX');
    M_7joints(j,1:2) = angles_yzx(1:2);
    M_7joints(j,3) = angles_yzx(3);
    % M_7joints(j, 2) = M_7joints(j, 2);
 

    c = M_7joints(j, :); 
    ee_pose = fk.fkm(c);
    x(j, :) = vec3(ee_pose.translation);
    %extract rotation and rotate the base wrist-to-hand vector to the end
    q_8 = vec8(ee_pose); 
    q_vec = q_8(1:4); 
    q_ = quaternion(q_vec(1), q_vec(2), q_vec(3), q_vec(4));  
    q_ = normalize(q_);
    v_comend(j,:) = rotatepoint(q_, base_vec);

    
end

figure
plot(ee_position(1:199,:) - ee_position(1,:))
title("ee transition human data")

figure
plot(x(4:end,:) - x(4,:))
title("ee transition calculated")

figure
plot(v_end(2:198,:))
title("collect unit vector end")

figure
plot(v_comend(2:end,:))
title("calculated unit vector end")

%figure;
%plot(elbowplot);




