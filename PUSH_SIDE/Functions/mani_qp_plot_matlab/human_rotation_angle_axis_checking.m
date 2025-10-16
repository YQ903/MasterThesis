robot = FrankaEmikaPandaRobot.kinematics();


wrist_pose = readmatrix('human_data/human_Jacobian/wrist_pose_shoulder_Drill_1.csv');
s = size(wrist_pose);
human_x_rotation_angle = zeros(s(1),1);
human_x_rotation_axis = zeros(s(1),3);

robot_pos = readmatrix("human_data/robot/q_position_mean_traj.csv");
s2 = size(robot_pos);
robot_x_rotation_angle = zeros(s2(1),1);
robot_x_rotation_axis = zeros(s2(1),3);

for i = 1:s
    human_x_rotation_angle(i,:) = (DQ(wrist_pose(i,:)).rotation_angle);
    human_x_rotation_axis(i,:) = transpose(vec3((DQ(wrist_pose(i,:)).rotation_axis)));


    Cartesian_pose = robot.fkm(robot_pos(i,:));
    robot_x_rotation_angle(i,:) = (DQ(Cartesian_pose.rotation_angle));
    robot_x_rotation_axis(i,:) = transpose(vec3((DQ(Cartesian_pose.rotation_axis))));

    % geo_jaco = geomJ(robot,position_real(i,:));
    % geo_jaco_vol = z
end

subplot(2,2,1)
plot(human_x_rotation_angle(:,1))
hold on
plot(robot_x_rotation_angle(:,1))
legend('human','robot')
title('rotation angle')

subplot(2,2,2)
plot(human_x_rotation_axis(:,1))
hold on
plot(robot_x_rotation_axis(:,1))
legend('human','robot')
title('rotation x axis')

subplot(2,2,3)
plot(human_x_rotation_axis(:,2))
hold on
plot(robot_x_rotation_axis(:,2))
legend('human','robot')
title('rotation y axis')

subplot(2,2,4)
plot(human_x_rotation_axis(:,3))
hold on
plot(robot_x_rotation_axis(:,3))
legend('human','robot')
title('rotation z axis')
