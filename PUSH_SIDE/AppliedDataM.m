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
M_all = csvread("data_3/2-processed_data/" + data_name + "/all_angles.csv");
position_all = csvread("data_3/2-processed_data/" + data_name + "/all_cartesian_translation.csv");
%BaseWrist = csvread("data_2/2-processed_data/" + data_name + "/Base_to_Wrist.csv");
rsr = RealSenseRobotTrans7DoF;
fk = rsr.kinematics();
position = position_all;
M = M_all(:,1:11);

% data: "push_front";
% position = position_all(1:210, :);
% M = M_all(1:210, :);

% data: "move_up";
% position = position_all(50:130,:);
% M = M_all(50:130,:);

wrist_position = position(:,end-2:end);
%plot(wrist_position - wrist_position(1,:))
plot(position(:, 7:9)-position(1, 7:9))
len = 1;

elbow_shouder = position(:, 1:3) - position(:, 4:6);
L1 = sqrt(sum(elbow_shouder.^2,2));

 

wrist_elbow = position(:,7:9) - position(:,4:6);
L2 =  sqrt(sum(wrist_elbow.^2,2));
%fk = rsr.kinematics();

%a = wrist_elbow(1,:)
%a = a/norm(a);
%b = [1,0,0]/norm([1,0,0])
%v = cross(a',b')
 %s = norm(v)
 %c = a * b'
 %vx = [0, - v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0]
 %R = eye(3,3) + vx + vx^2 * (1-c)/(s^2)
 %R * a'
 %sss = size(wrist_elbow)
 %sss111 = sss(1)
 %for p = 1:1:sss111
 %    a150 = wrist_elbow(p,:);
 %    a150 = a150/norm(a150);
 %    wrist_elbow(p,:) = wrist_elbow(p,:) / norm(wrist_elbow(p,:));
 %    wrist_elbow(p,:) = a150 * R ;

 %end
 %plot(wrist_elbow)
s = size(M);
s1 = s(1);

x = zeros(s1, 3);
mani = zeros(s1, 3);

offset = [0,180,180,180,0,90, 0,0,0,0,0] - M(1,:);

M_7joints = zeros(s1, 7); 

M_7joints(:,1:3) = M(:,2:4);
M_7joints(:,4) = M(:,6);
%M_7joints(:,5:7) = 0;
M_7joints(:,5:7) = M(:,9:11);
M_7joints(:, :)= M_7joints(:, :)/180 * pi;

%c_initial = M(1, :)/180 * pi;
%ee_pose_initail = fk.fkm(c_initial);
%r_initail = ee_pose_initail.rotation();
%f = r_initail.conj() + 0.5 * E_ * r_initail.conj() * DQ([0,0,0,0]);
%fk.set_base_frame(f);
%fk.set_reference_frame(f);
mani = zeros(s1, 9);

for j = 1:s1
    R = eul2rotm([M_7joints(j,1), M_7joints(j,2), M_7joints(j,3)+ 90/180 * pi], 'XYZ');
    %elbowplot(j,:) = R *[0.0, 0.0, 0.24]';
    angles_yzx = rotm2eul(R, 'YZX');
    M_7joints(j,1:2) = angles_yzx(1:2);
    M_7joints(j,3) = angles_yzx(3);
    c = M_7joints(j, :); 

    ee_pose = fk.fkm(c);
    ee_pose.translation
    ee_rotation = ee_pose.rotation();
    l_direction = ee_rotation * i_ * ee_rotation.conj();
    a = ee_pose.translation().vec3();
    a(2) = -a(2);
    a(3) = -a(3);
    b = fk.fkm(c,6).translation().vec3();
    b(2) = -b(2);
    b(3) = -b(3);
    l_direction = a  - b;
    %l_back = r_initail.conj() * l_direction * r_initail;
    wrist_elbow(j,:) = l_direction/norm(l_direction);
    x(j, :) = vec3(ee_pose.translation);
    
    jaco = geomJ(fk,c,7);
    mani_elli = jaco(4:6, :) * jaco(4:6, :)';
    mani(j,:) = reshape(mani_elli.', 1, 9);
    %mani_elli = jaco(4:6, :) * transpose(jaco(4:6, :));
    %mani(j, :) = transpose(eig(mani_elli));
    
end



n  = size(M_7joints,2);
W  = 5; 
Npos = size(position,1);
Nq   = size(M_7joints,1);
N    = min(Npos, Nq);
mani_h = nan(N,9);                   % 1x9 flattened 3x3 manipulability per frame

for t = 1:N
    i0  = max(1, t-W);
    i1  = min(N, t+W);                     % clamp to N
    idx = i0:i1;              % neighborhood around frame t
    dQ = M_7joints(idx,:) - M_7joints(t,:);   % (2W+1) x n
    dX = position(idx,10:12) - position(t,10:12);% (2W+1) x 3

    Jv = (dQ \ dX)';                  % 3 x n translational Jacobian

    S  = Jv * Jv.';                   % 3x3 manipulability ellipsoid
    S  = (S+S.')/2;                   % symmetrize (numerical hygiene)
    mani_h(t,:) = reshape(S.',1,9);   % row-major: m11 m12 m13 m21 ...
end



wrist_elbow(:,2) = wrist_elbow(:,2)
wrist_elbow(:,3) = wrist_elbow(:,3)
%plot(wrist_elbow)
s = size(M);
s1 = s(1);
s2 = s(2);
M_extend = zeros(s1*len, s2);

wrist_extend = zeros(s1*len, 3);
initial_wrist = x(1, :);


csvwrite('data_3/3-applied_data/' + data_name + '/manipulabilities.csv', mani)
csvwrite('data_3/3-applied_data/' + data_name + '/manipulabilities_human.csv', mani_h)

% plot(x)
% legend("cartesian x", "cartesian y", "cartesian z")


% plot(M(:,2:4)/180 * 2 * pi)
% legend("shoulder x", "shoulder y", "shouder z")
% title("shoulder angles")
% 
% figure
% plot(M(:,6)/180 * 2 * pi)
% legend("elbow")
% title("elbow angles")



% %%%%%%%%%%%%%%%%%%%%

% s = size(position);
% s1 = s(1);
% 
% wrist = M(:,7:7);


%%
% ind = 1;
% M_ind = ind+6;
% offset = position(1:1,M_ind:M_ind) - x(1:1,ind:ind);
% 
% subplot(2, 3, ind)
% plot(position(:,M_ind:M_ind))
% title("wrist x - calculated")
% subplot(2, 3, ind + 3)
% plot(x(:,ind:ind))
% title("wrist x - recorded")
%legend("recorded", "calculated by fkm")

%%

% for i = 1:3
%      subplot(1, 3, i)
%      plot(mani(:,i:i))
%  end
%  sgtitle("eigen values of kinematic mani - move up")
