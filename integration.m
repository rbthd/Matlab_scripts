%% Global variables

clc; %reset command window
clearvars; %reset variables
results=importdata('/Users/robintherond/Documents/Polytech/PIFE/instantWave/Tests/test_barometer_slow_rob.csv', ' ');
scrsz = get(groot,'ScreenSize');

samplePer = 0.02; %s

orderButter = 2; %order of the HP filter
orderPolyBaro = 10;
orderButterBaro = 2;

filtCutOff = 0.3; %cut off frequency of HP filter
filtCutOffLP = 1;
SamplePlotFreq = 2; %For the Animation

length_data = length(results);

acc=[results(:,26),results(:,27),results(:,28)];
gyr=[results(:,29),results(:,30),results(:,31)];
baro=results(:,36);

R = zeros(3,3,length(gyr));     % rotation matrix describing sensor relative to Earth

interv=0.02:0.02:0.02*length(results(:,36));

%reference function of motion
h0=baro(1); %First barometer output -> altitude reference
height_ref=h0+0.08*sin(2*pi*(1/7.75)*interv);


ahrs = MahonyAHRS('SamplePeriod', samplePer, 'Kp', 1);

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


linAcc = ComputelinAcc(acc);
linVel = BasicIntAcc(linAcc, samplePer);
linVelHP = HPFilter(linVel,orderButter,filtCutOff,samplePer);
linPos = BasicIntVel(linVel, samplePer);
linPosHP = HPFilter(linPos,orderButter,filtCutOff,samplePer);

%plotAcc(interv,acc);
plotLinAcc(interv,linAcc);
%plotLinVel(interv,linVel);
%plotLinVelHP(interv,linVelHP);
%plotLinPos(interv,linPos);
plotLinPosHP(interv,linPosHP,height_ref);

%plotBarometer(baro,interv, orderPolyBaro, orderButterBaro, filtCutOffLP, samplePer)




% SixDOFanimation(linPosHP, R, ...
%                  'SamplePlotFreq', SamplePlotFreq, 'Trail', 'Off', ...
%                  'Position', [9 39 1280 720], ...
%                  'AxisLength', 0.1, 'ShowArrowHead', false, ...
%                  'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, 'Title', 'Unfiltered',...
%                  'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePer) / SamplePlotFreq)); 



%% Digital Butterworth High-pass filter
function [output] = HPFilter(input, orderButter, filtCutOff,samplePer)
[b, a] = butter(orderButter, (2*filtCutOff)/(1/samplePer), 'high');
output = filtfilt(b, a, input);
end

%% Digital butterworth Low-pass filter
function [output] = LPFilter(input, orderButter, filtCutOff, samplePer)
[b, a] = butter(orderButter, (2*filtCutOff)/(1/samplePer), 'low');
output = filtfilt(b, a, input);
end

%% Calculate linear acceleration in Earth frame (substracting gravity)
function[linAcc] = ComputelinAcc(acc)
linAcc = acc - [zeros(length(acc), 1), zeros(length(acc), 1), ones(length(acc), 1)]; %z = z-1;
linAcc = linAcc * 9.81;     % convert from 'g' to m/s/s
end

%% Calculate linear velocity (integrate acceleartion)
function[linVel] = BasicIntAcc(linAcc, samplePer)
linVel=zeros(size(linAcc));
for i = 2:length(linAcc)
    linVel(i,:) = linVel(i-1,:) + linAcc(i,:) * samplePer;
end
end

%% Calculate linear position (integrate velocity)
function[linPos] = BasicIntVel(linVel, samplePer)
linPos=zeros(size(linVel));
for i = 2:length(linVel)
    linPos(i,:) = linPos(i-1,:) + linVel(i,:) * samplePer;
end
end

%% Plot unfiltered acceleration
function plotAcc(interv,acc)
figure('Name','Acceleration unfiltered')
hold on;
plot(interv,acc(:,1));
plot(interv,acc(:,2));
plot(interv,acc(:,3));
xlabel('time (s)');
ylabel('g');
legend('X', 'Y', 'Z');
end

%% Plot Linear Acceleration
function plotLinAcc(interv,linAcc)
figure('NumberTitle', 'off', 'Name', 'Linear Acceleration');
hold on;
plot(interv,linAcc(:,1), 'r');
plot(interv,linAcc(:,2), 'g');
plot(interv,linAcc(:,3), 'b');
xlabel('time(s)');
ylabel('m/s^2');
title('Linear acceleration');
legend('X', 'Y', 'Z');
end

%% Plot Linear Volecity
function plotLinVel(interv,linVel)
figure('NumberTitle', 'off', 'Name', 'Linear Velocity');
hold on;
plot(interv,linVel(:,1), 'r');
plot(interv,linVel(:,2), 'g');
plot(interv,linVel(:,3), 'b');
xlabel('sample');
ylabel('m/s');
title('Linear velocity');
legend('X', 'Y', 'Z');
end

%% Plot high pass filtered linear velocity
function plotLinVelHP(interv,linVelHP)
figure('NumberTitle', 'off', 'Name', 'High-pass filtered Linear Velocity');
hold on;
plot(interv,linVelHP(:,1), 'r');
plot(interv,linVelHP(:,2), 'g');
plot(interv,linVelHP(:,3), 'b');
xlabel('sample');
ylabel('m/s');
title('High-pass filtered linear velocity');
legend('X', 'Y', 'Z');
end

%% Plot unfiltered linear position
function plotLinPos(interv,linPos)
figure('NumberTitle', 'off', 'Name', 'Linear Position');
hold on;
plot(interv,linPos(:,1), 'r');
plot(interv,linPos(:,2), 'g');
plot(interv,linPos(:,3), 'b');
xlabel('sample');
ylabel('m');
title('Linear position');
legend('X', 'Y', 'Z');
end

%% Plot High-pass filtered Linear Position
function plotLinPosHP(interv, linPosHP, f)
% Plot
figure('NumberTitle', 'off', 'Name', 'High-pass filtered Linear Position');
hold on;
plot(interv,linPosHP(:,1), 'r', 'LineWidth',2);
plot(interv,linPosHP(:,2), 'g', 'LineWidth',2);
plot(interv,linPosHP(:,3), 'b', 'LineWidth',2);
plot(interv,f, 'm','LineWidth', 2); %Reference function for comparison
f=f+0.6;
plot(interv,f, 'c','LineWidth', 2); %Reference offsted function for comparison

xlabel('time (s)');
ylabel('m');
title('High-pass filtered linear position');
legend('X', 'Y', 'Z','Reference function','Offsted Reference function');
axis([0 23 -0.8 0.4])
hold off;

figure('Name','Resampling');
plot(interv,f,'c','LineWidth',3); %Ref func
x=0.05:0.05:0.05*length(linPosHP);
hold on;
plot(interv,linPosHP(:,3), 'b', 'LineWidth',1);
plot(x-(400*0.04),linPosHP(:,3),'r','LineWidth',3);
hold off;
legend('Reference function','Real outputs','Resampled & time-delayed outputs');
axis([0 23 -0.2 0.2])
title('Z position after resampling and time-delay');
xlabel('time (s)');
ylabel('altitude (m)');

end


%% Plot barometer outputs

function plotBarometer(baro,interv,orderPolyBaro, orderButterBaro, filtCutOffLP, samplePer)
scrsz = get(groot,'ScreenSize');
figure('Name','Altitude from barometer','Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);
hold on;

%REFERENCE Function 
h0=baro(1); %First barometer output -> altitude reference
height_ref=h0+0.1*sin(2*pi*(1/7.75)*interv);
plot(interv,height_ref,'m*'); %Reference function

%Unfiltered data
plot(interv,baro,'Color','k'); %Unfiltered output

%LOW PASS Filter
altitudeLP = LPFilter(baro, orderButterBaro, filtCutOffLP, samplePer);
plot(interv,altitudeLP,'Color','b','LineWidth',3); %Low pass filtered output

%Sliding average
for i = 1:length(altitudeLP)
baro_avg(i) = mean(altitudeLP(max([1 i-50]):i));
end
plot(interv,baro_avg,'Color','y','LineWidth',4); %Averaged output


%Polynomial extrapolation
[p,~,mu] = polyfit(interv,baro_avg,orderPolyBaro); %scaling and centering
f=polyval(p,interv,[],mu);
plot(interv,f,'LineWidth',4,'Color','r');
xlabel('Time (s)');
ylabel('Height (m)');


legend({'Unfiltered','Low pass filtered','Sliding average','Polynomial extrapolation'},'FontSize',15);
hold off;

end