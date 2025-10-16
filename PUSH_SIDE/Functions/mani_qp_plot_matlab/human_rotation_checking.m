robot = FrankaEmikaPandaRobot.kinematics();


wrist_pose = readmatrix('human_data/human_Jacobian/wrist_pose_shoulder_Drill_1.csv');
s = size(wrist_pose);
human_x_position = zeros(s(1),4);

robot_pos = readmatrix("human_data/robot/q_position_mean_traj.csv");
s2 = size(robot_pos);
robot_x_position = zeros(s2(1),4);

for i = 1:s
    human_x_position(i,:) = transpose(vec4(DQ(wrist_pose(i,:)).rotation)) * 0.01;
    Cartesian_pose = robot.fkm(robot_pos(i,:));
    robot_x_position(i,:) = transpose(vec4(Cartesian_pose.rotation));

    % geo_jaco = geomJ(robot,position_real(i,:));
    % geo_jaco_vol = z
end

subplot(1,4,1)
plot(human_x_position(:,1))
hold on
plot(robot_x_position(:,1))
legend('human','robot')
title('x axis (human cartesian/100)')

subplot(1,4,2)
plot(human_x_position(:,2))
hold on
plot(robot_x_position(:,2))
legend('human','robot')
title('y axis(human cartesian/100)')

subplot(1,4,3)
plot(human_x_position(:,3))
hold on
plot(robot_x_position(:,3))
legend('human','robot')
title('z axis(human cartesian/100)')

subplot(1,4,4)
plot(human_x_position(:,4))
hold on
plot(robot_x_position(:,4))
legend('human','robot')
title('z axis(human cartesian/100)')



