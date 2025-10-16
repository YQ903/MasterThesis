robot = FrankaEmikaPandaRobot.kinematics();


dq_jaco = readmatrix('human_data/human_Jacobian/geometrical_Jacobian_Drill_1.csv');
s = size(dq_jaco);
human_jaco_vol = zeros(s(1),6);

robot_pos = readmatrix("human_data/robot/q_position_mean_traj.csv");
s2 = size(robot_pos);
robot_jaco_vol = zeros(s2(1),6);

C8 = diag([-1 ones(1,3) -1 ones(1,3)]');
C4m = -C8(1:4,1:4);  
CJ4_2_J3= [zeros(3,1) eye(3)];

for i = 1:s
    i_jaco_col = dq_jaco(i,:);
    i_jaco_8D = reshape(i_jaco_col,8,8);
    i_jaco = i_jaco_8D([2:4, 6:8],:)/ 1.6;
    human_mani = i_jaco * transpose(i_jaco);
    eigen_human = transpose(svd(i_jaco));
    human_jaco_vol(i,:) = eigen_human;



     
    J = zeros(6,7);
    Jacob = robot.pose_jacobian(robot_pos(i, :));
    xm = robot.fkm(robot_pos(i, :));
    J(1:3,:) = CJ4_2_J3*2*haminus4(xm.P')*Jacob(1:4,:); 
    J(4:6,:) = CJ4_2_J3*2*( hamiplus4(xm.D)*C4m*Jacob(1:4,:) + haminus4(xm.P')*Jacob(5:8,:)  );
    
    robot_jaco = J;
    %robot.pose_jacobian();

    robot_mani = robot_jaco * transpose(robot_jaco);
    eigen_robot = transpose(svd(J));
    robot_jaco_vol(i,:) = eigen_robot;

    % geo_jaco = geomJ(robot,position_real(i,:));
    % geo_jaco_vol = z
end

subplot(2,3,1)
plot(human_jaco_vol(:,1))
hold on
plot(robot_jaco_vol(:,1))
legend('human','robot')
title('eigenvalue 1')

subplot(2,3,2)
plot(human_jaco_vol(:,2))
hold on
plot(robot_jaco_vol(:,2))
legend('human','robot')
title('eigenvalue 2')

subplot(2,3,3)
plot(human_jaco_vol(:,3))
hold on
plot(robot_jaco_vol(:,3))
legend('human','robot')
title('eigenvalue 3')

subplot(2,3,4)
plot(human_jaco_vol(:,4))
hold on
plot(robot_jaco_vol(:,4))
legend('human','robot')
title('eigenvalue 4')

subplot(2,3,5)
plot(human_jaco_vol(:,5))
hold on
plot(robot_jaco_vol(:,5))
legend('human','robot')
title('eigenvalue 5')

subplot(2,3,6)
plot(human_jaco_vol(:,6))
hold on
plot(robot_jaco_vol(:,6))
legend('human','robot')
title('eigenvalue 6')
