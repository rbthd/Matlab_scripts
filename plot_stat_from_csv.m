clc; %reset command window
clear all; %reset variables
%results=importdata('test_position_settings2.csv', ' ');
%results=importdata('/Volumes/Home Directory/instantWave/minimu9/tests/test3.csv', ' ');
results=importdata('/Users/robintherond/Documents/Polytech/PIFE/instantWave/Tests/test_WAVE_robot.csv', ' ');

%% Plotting stat from Euler angles 
scrsz = get(groot,'ScreenSize');
figure('Name' ,'Normal distribution of Euler Angles','Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);

subplot(321);
[mu,s,muci,sci] = normfit(results(:,23));
norm_yaw = normpdf(results(:,23),mu, s);
scatter(results(:,23),norm_yaw,'.');
dim = [.15 .8 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('YAW');
xlabel('degrees');

subplot(323);
[mu,s,muci,sci] = normfit(results(:,24));
norm_pitch = normpdf(results(:,24),mu, s);
scatter(results(:,24),norm_pitch,'.');
dim = [.15 .5 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('PITCH');
xlabel('degrees');

subplot(325);
[mu,s,muci,sci] = normfit(results(:,25));
norm_roll = normpdf(results(:,25),mu, s);
scatter(results(:,25),norm_roll,'.');
dim = [.15 .2 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('ROLL');
xlabel('degrees');

%% Plotting stat from Euler angles 

subplot(322);
[mu,s,muci,sci] = normfit(results(:,26));
norm_acc_x = normpdf(results(:,26),mu, s);
scatter(results(:,26),norm_acc_x,'.');
dim = [.8 .8 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('ACCELERO X');
xlabel('g');

subplot(324);
[mu,s,muci,sci] = normfit(results(:,27));
norm_acc_y = normpdf(results(:,27),mu, s);
scatter(results(:,27),norm_acc_y,'.');
dim = [.8 .5 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('ACCELERO Y');
xlabel('g');

subplot(326);
[mu,s,muci,sci] = normfit(results(:,28));
norm_acc_z = normpdf(results(:,28),mu, s);
scatter(results(:,28),norm_acc_z,'.');
dim = [.8 .2 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('ACCELERO Z');
xlabel('g');

%% Plotting magnitude of accelero vector

size2=size(results)
size=size2(1);
mag_acc=zeros(size,1);
x=(1:1:size);

for i=1:1:size
    mag_acc(i,1)=sqrt((results(i,26)*results(i,26))+(results(i,27)*results(i,27))+(results(i,28)*results(i,28)));
end

figure('Name' ,'Magnitude of acceleration vector','Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);
subplot(221);
plot(x,mag_acc);
xlabel('Samples');
ylabel('Amplitude of Acceleration vector');
subplot(222);
[mu,s,muci,sci] = normfit(mag_acc);
normal_acc = normpdf(mag_acc,mu, s);
max_normal_acc = max(normal_acc)
for i=1:1:size
    normal_acc(i,1) = normal_acc(i,1)/max_normal_acc;
end
scatter(mag_acc,normal_acc,'.');
xlabel('Amplitude of Acceleration vector');
ylabel('Normalized Density of probability');
norm_acc_avg=mean(mag_acc);
std_deviation=std(mag_acc);
dim = [.4 .82 .1 .1];
str={'Average = ' num2str(norm_acc_avg),'Standard deviation = ' num2str(std_deviation)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');





%% Plotting magnitude of magneto vector


mag_magneto=zeros(size,1);
x=(1:1:size);

for i=1:1:size
    mag_magneto(i,1)=sqrt((results(i,29)*results(i,29))+(results(i,30)*results(i,30))+(results(i,31)*results(i,31)));
end

figure('Name' ,'Magnitude of magnetometer vector','Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);
subplot(221);
plot(x,mag_magneto);
subplot(222);
[mu,s,muci,sci] = normfit(mag_magneto);
norm_magneto = normpdf(mag_magneto,mu, s);
scatter(mag_magneto,norm_magneto,'.');
norm_acc_avg=mean(mag_magneto);
std_deviation=std(mag_magneto);
str={'Average = ' num2str(norm_acc_avg),'Standard deviation = ' num2str(std_deviation)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');



