%% INTEGRATION OF ACCELERATION

clc; %reset command window
clearvars; %reset variables
results=importdata('/Users/robintherond/Documents/Polytech/PIFE/instantWave/Tests/test_robot_cut.csv', ' ');

%Global variables
samplePeriod = 0.02; %s
order_butter = 2; %order of the HP filter
filtCutOff = 1; %cut off frequency of HP filter

size2=size(results);
size_file=size2(1);

acc_x=results(:,26);
acc_y=results(:,27);
acc_z=results(:,28);

gyr_x=results(:,29);
gyr_y=results(:,30);
gyr_z=results(:,31);


acc=[acc_x,acc_y,acc_z];
gyr=[gyr_x,gyr_y,gyr_z];

x=1:1:size_file;


R = zeros(3,3,length(gyr));     % rotation matrix describing sensor relative to Earth

ahrs = MahonyAHRS('SamplePeriod', samplePeriod, 'Kp', 1);

for i = 1:length(gyr)
    ahrs.UpdateIMU(gyr(i,:) * (pi/180), acc(i,:));	% gyroscope units must be radians
    R(:,:,i) = quatern2rotMat(ahrs.Quaternion)';    % transpose because ahrs provides Earth relative to sensor
end

%% Calculate 'tilt-compensated' accelerometer

% tcAcc = zeros(size(acc));  % accelerometer in Earth frame
% 
% for i = 1:length(acc)
%     tcAcc(i,:) = R(:,:,i) * acc(i,:)';
% end

figure('Name','Acceleration unfiltered')
hold on;
plot(x,acc_x);
plot(x,acc_y);
plot(x,acc_z);
xlabel('sample');
ylabel('g');
legend('X', 'Y', 'Z');




%% Calculate linear acceleration in Earth frame (subtracting gravity)

linAcc = acc - [zeros(length(acc), 1), zeros(length(acc), 1), ones(length(acc), 1)];
linAcc = linAcc * 9.81;     % convert from 'g' to m/s/s

%% Plot
figure('NumberTitle', 'off', 'Name', 'Linear Acceleration');
hold on;
plot(linAcc(:,1), 'r');
plot(linAcc(:,2), 'g');
plot(linAcc(:,3), 'b');
xlabel('sample');
ylabel('m/s^2');
title('Linear acceleration');
legend('X', 'Y', 'Z');

%% Calculate linear velocity (integrate acceleartion)

linVel = zeros(size(linAcc));

for i = 2:length(linAcc)
    linVel(i,:) = linVel(i-1,:) + linAcc(i,:) * samplePeriod;
end

%% Integration of acceleration using polynomial methods

linVel2 = zeros(size(linAcc));
linVel_X = zeros(size(linAcc(:,1)));
linVel_Y = zeros(size(linAcc(:,2)));
linVel_Z = zeros(size(linAcc(:,3)));

% for i=2:1:length(linAcc)-1
% 
% interv=(i):1:(i+1);
% 
% %Integrate x axis
% [p,~,mu] = polyfit(interv,[linAcc(i,1),linAcc(1+i,1)],1);
% f=polyval(p,interv,[],mu); %return coeff of the polynomial p in desending powers.
% integral_X = polyint(f); %return coeff of the primitive of P
% integral_X = polyval(integral_X,i+1);
% linVel_X(i+1) = linVel_X(i) + integral_X;
% 
% 
% %Integrate Y axis
% [p,~,mu] = polyfit(interv,[linAcc(i,2),linAcc(1+i,2)],1);
% f=polyval(p,interv,[],mu);
% integral_Y = polyint(p);
% integral_Y = polyval(integral_Y,i+1);
% linVel_Y(i+1) = linVel_Y(i) + integral_Y;
% 
% %Integrate Z axis
% [p,~,mu] = polyfit(interv,[linAcc(i,3),linAcc(1+i,3)],1);
% f=polyval(p,interv,[],mu);
% integral_Z = polyint(p);
% integral_Z = polyval(integral_Z,i+1);
% linVel_Z(i+1) = linVel_Z(i) + integral_Z;
% 
% linVel2=[linVel_X, linVel_Y, linVel_Z];
% end;

x=1:1:size(linAcc(:,1));
x=x'
order_poly_integ = 100 %length(linAcc)

%Polynomial extrapolation
[p,~,mu] = polyfit(x,linAcc(:,1),order_poly_integ);
f=polyval(p,x,[],mu); %return coeff of the polynomial p in desending powers.
integral_X = polyint(p); %return coeff of the primitive of P.

%Actual integration
for i = 2:length(linAcc)
    linVel(i,:) = linVel(i-1,:) + linAcc(i,:) * samplePeriod;
end 

figure('Name','Velocity after integration')
plot(x,f,'o');
hold on ;
plot(x,linAcc(:,1));
hold off;


%% Plot classical integration
figure('NumberTitle', 'off', 'Name', 'Linear Velocity');
hold on;
plot(linVel(:,1), 'r');
plot(linVel(:,2), 'g');
plot(linVel(:,3), 'b');
xlabel('sample');
ylabel('m/s');
title('Linear velocity');
legend('X', 'Y', 'Z');

%% Plot polynomial integration
figure('NumberTitle', 'off', 'Name', 'Linear Velocity after Poly integr');
hold on;
plot(linVel2(:,1), 'r');
plot(linVel2(:,2), 'g');
plot(linVel2(:,3), 'b');
xlabel('sample');
ylabel('m/s');
title('Linear velocity 2');
legend('X', 'Y', 'Z');

%% High-pass filter linear velocity to remove drift

[b, a] = butter(order_butter, (2*filtCutOff)/(1/samplePeriod), 'high');
linVelHP = filtfilt(b, a, linVel);

%% Plot
% figure('NumberTitle', 'off', 'Name', 'High-pass filtered Linear Velocity');
% hold on;
% plot(linVelHP(:,1), 'r');
% plot(linVelHP(:,2), 'g');
% plot(linVelHP(:,3), 'b');
% xlabel('sample');
% ylabel('m/s');
% title('High-pass filtered linear velocity');
% legend('X', 'Y', 'Z');

%% Calculate linear position (integrate velocity)

linPos = zeros(size(linVelHP));

for i = 2:length(linVelHP)
    linPos(i,:) = linPos(i-1,:) + linVelHP(i,:) * samplePeriod;
end

% Plot
figure('NumberTitle', 'off', 'Name', 'Linear Position');
hold on;
plot(linPos(:,1), 'r');
plot(linPos(:,2), 'g');
plot(linPos(:,3), 'b');
xlabel('sample');
ylabel('m');
title('Linear position');
legend('X', 'Y', 'Z');

%% High-pass filter linear position to remove drift

[b, a] = butter(order_butter, (2*filtCutOff)/(1/samplePeriod), 'high');
linPosHP = filtfilt(b, a, linPos);

% Plot
figure('NumberTitle', 'off', 'Name', 'High-pass filtered Linear Position');
hold on;
plot(linPosHP(:,1), 'r');
plot(linPosHP(:,2), 'g');
plot(linPosHP(:,3), 'b');
xlabel('sample');
ylabel('m');
title('High-pass filtered linear position');
legend('X', 'Y', 'Z');

%figure('Name','3D trajectory')
%scatter3(linPosHP(:,1),linPosHP(:,2),linPosHP(:,3))


x=1:1:size(results(:,36));
figure('Name','Height');
hold on;
plot(x,results(:,36));
[p,~,mu] = polyfit(x,results(:,36)',200);
f=polyval(p,x,[],mu);
plot(x,f);

% f=f';
% 
% for i=1:1:size_file
% linPosHP(i,3)=linPosHP(i,3)*0.5+f(i,1)*0.5;
% end
% 
% %% Plot After Complementary filter
% figure('NumberTitle', 'off', 'Name', 'Complementary filter');
% hold on;
% plot(linPosHP(:,1), 'r');
% plot(linPosHP(:,2), 'g');
% plot(linPosHP(:,3), 'b');
% xlabel('sample');
% ylabel('m');
% title('High-pass filtered linear position + Complementary filter');
% legend('X', 'Y', 'Z');

SamplePlotFreq = 2;

% SixDOFanimation(linPosHP, R, ...
%                  'SamplePlotFreq', SamplePlotFreq, 'Trail', 'Off', ...
%                  'Position', [9 39 1280 720], ...
%                  'AxisLength', 0.1, 'ShowArrowHead', false, ...
%                  'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, 'Title', 'Unfiltered',...
%                  'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq)); 