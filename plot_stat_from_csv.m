clc; %reset command window
clear all; %reset variables
results=importdata('/Users/robintherond/Documents/Polytech/PIFE/instantWave/Tests/test_nonMoving2.csv', ' ');

%The column number might vary depending on the csv
%file you use. Verify the order if you use data from
%another program than minimu9-ahrs in 'all' mode.

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

%% Plotting stat from Accelerometers  

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
normal_mag = normpdf(mag_acc,mu, s);
max_normal_acc = max(normal_mag)
for i=1:1:size
    normal_mag(i,1) = normal_mag(i,1)/max_normal_acc;
end
scatter(mag_acc,normal_mag,'.');
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


%% Plotting stat from Gyrometers  
x=1:1:length(results(:,29));
figure('Name','Stats from gyrometers','Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)])
subplot(321);
[mu,s,muci,sci] = normfit(results(:,29));
norm_gyr_x = normpdf(results(:,29),mu, s);
plot(results(:,29),norm_gyr_x,'.');
dim = [.37 .8 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('GYROMETER X');
xlabel('°/s');
ylabel('Density of probability');
subplot(322)
scatter(x,results(:,29))
title('GYROMETER X');
xlabel('samples');
ylabel('°/s');


subplot(323);
[mu,s,muci,sci] = normfit(results(:,30));
norm_gyr_y = normpdf(results(:,30),mu, s);
plot(results(:,30),norm_gyr_y,'.');
dim = [.37 .5 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
xlabel('°/s');
ylabel('Density of probability');
title('GYROMETER Y');

subplot(324)
scatter(x,results(:,30))
xlabel('samples');
ylabel('°/s');
title('GYROMETER Y');

subplot(325);
[mu,s,muci,sci] = normfit(results(:,31));
norm_gyr_z = normpdf(results(:,31),mu, s);
plot(results(:,31),norm_gyr_z,'.');
dim = [.37 .2 .1 .1];
str={'Average = ' num2str(mu),'Standard deviation = ' num2str(s)};
annotation('textbox',dim,'String',str , 'FitBoxToText','on');
title('GYROMETER Z');

xlabel('°/s');
ylabel('Density of probability');
subplot(326)
scatter(x,results(:,31))
title('GYROMETER Z');
xlabel('samples');
ylabel('°/s');

title('GYROMETER Z');
xlabel('°/s');
