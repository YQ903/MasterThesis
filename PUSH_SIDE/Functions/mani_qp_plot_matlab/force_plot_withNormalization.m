% Guo Yu
% Plot results from csv file
% Please change the base's reference frame in FrankaEmikaPandaRobot.m!
% (try to combine several loops together
% addpath(genpath('..\..\..\dqrobotics-toolbox-matlab'));
% position_guid = readmatrix('joint_position_exam_force_traj.csv'); 
% position_real = readmatrix('joint_position_real_joint_limit.csv'); 
% velocity_guid = readmatrix('joint_velocity_exam_force.csv'); 
% velocity_real = readmatrix('joint_velocity_real_joint_limit.csv'); 
force_data = readmatrix('/home/yre/Desktop/franka_ros/admittance/data/force_recording/force_recording.csv');
F_ext_threshold_array = [5, 10, 15, 20];
s7 = size(F_ext_threshold_array);
s8 = s7(2);
for k = 1:s8
    subplot(2, 2, k)
    force = vertcat(force_data, zeros(1000, 6));
    F_ext_threshold = F_ext_threshold_array(1,k);
    s2 = size(force);
    s = s2(1);
    plot(force(3000:s,1))
    legend("original foce")
    fc_array = [10, 1, 0.8, 0.5];
    s3 = size(fc_array);
    for j = 1:s3(2)   
        Ts = 0.001;
        fc = fc_array(1,j);
        force = vertcat(force_data, zeros(1000, 6));
        hold on
        a = 2 * 3.14 * fc * Ts / (2 * 3.14 * fc * Ts + 1);
        % Low Pass Filter: y_n = a x_n + (1 - a) * y_{n-1} = y_{n-1} + a * (x_n - y_{n-1})
        
        for i = 1:s
            if (i == 1)
                force(i,:) = force(1,:);
                if (norm(force(i,:)) > F_ext_threshold)
                   force(i,:) = (norm(force(i,:)) - F_ext_threshold)/ norm(force(i,:))*force(i,:);
                   F_ext_fil_last = force(i,:);
                else
                   force(i,:) =  zeros(size(force(i,:)));
                   F_ext_fil_last = force(i,:);
                end
            else
                if (norm(force(i,:)) > F_ext_threshold)
                   force(i,:) = (norm(force(i,:)) - F_ext_threshold)/ norm(force(i,:))*force(i,:);
                else
                   force(i,:) =  zeros(size(force(i,:)));
                end
                force(i,:) = F_ext_fil_last + a * (force(i,:) - F_ext_fil_last);
                F_ext_fil_last = force(i,:);
            end
        end 
        plot(force(3000:s,1),'DisplayName','fc='+string(fc))
    end
    
    title('force with deadzone'+string(F_ext_threshold))
end
