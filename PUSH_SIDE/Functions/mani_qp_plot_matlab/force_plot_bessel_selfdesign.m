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

force = vertcat(force_data, zeros(1000, 6));

F_ext_threshold = 5;

s2 = size(force);
s = s2(1);
fc_array = [1];
s3 = size(fc_array);

for i = 1:s
    if (i == 1)
        force(i,:) = force(1,:);
        if (norm(force(i,:)) > F_ext_threshold)
           force(i,:) = (norm(force(i,:)) - F_ext_threshold)/ norm(force(i,:))*force(i,:);
           F_ext_fil_last = force(i,:);
        else
           force(i,:) =  zeros(size(force(i,:)));
        end 
    end 
end 


for j = 1:1
     %   Design a 5-th order Bessel analog filter (fc = 100 Hz) and convert it to a
    %   digital filter.
 
    Fs =1000;                            % Sampling frequency
    N = 2; 
    fc = 20;
    [z,p,k] = besself(N,fc);          % Bessel analog filter design
    [num,den]=zp2tf(z,p,k);             % Convert to transfer function form
    [numd,dend]=bilinear(num,den,Fs);   % Analog to Digital conversion
%     [sos] = zp2sos(zd,pd,kd);           % Convert to SOS form
%     fvtool(sos)                         % Visualize the digital filter
    samples = 6002;
    signal = force(:, 1);
    dt = 1/Fs;
    time = (0:samples-1)*dt;
    signal_filtered_bessel = filtfilt(numd,dend,signal);
    plot(signal)
    hold on
    plot(signal_filtered_bessel)
    hold on


    Ts = 0.001;
    fc = fc_array(1,j);
    force = vertcat(force_data, zeros(1000, 6));
    hold on
    a = 2 * 3.14 * fc * Ts / (2 * 3.14 * fc * Ts + 1);
    % Low Pass Filter: y_n = a x_n + (1 - a) * y_{n-1} = y_{n-1} + a * (x_n - y_{n-1})
    
    for i = 1:s
        if (i == 1)
            force(i,:) = force(1,:);
            F_ext_fil_last = force(i,:);
        else
            force(i,:) = F_ext_fil_last + a * (force(i,:) - F_ext_fil_last);
            F_ext_fil_last = force(i,:);
        end
    end 
    plot(force(1:s,1),'DisplayName','fc='+string(fc))
end

title("force without deadzone")