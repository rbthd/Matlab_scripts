clc; %reset command window
clf; %reset figures
clear all; %reset variables
period=20; %milliseconds
results=importdata('test_MAGNETO.csv', ' ');

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);
subplot(221);
scatter3(results(:,7),results(:,8),results(:,9), '.');
title('Magnetometers measurement en 3D view - Uncalibrated')
hold on
plot3([-5000 5000],[0 0], [0 0],'r',[0 0], [-5000 5000], [0 0], 'r', [0 0],[0 0], [-5000 15000], 'r');
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')

subplot(222);
scatter3(results(:,26),results(:,27),results(:,28), '.');
title('Magnetometers measurement en 3D view Calibrated')
hold on
plot3([0 0],[-10 10], [0 0],'r',[-10 10], [0 0], [0 0], 'r', [0 0],[0 0], [-10 10], 'r');
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')

subplot(223);
grid on
scatter3(results(:,7),results(:,8),results(:,9), '.');
hold on
plot3([-2000 6000],[0 0], [0 0],'r',[0 0], [-4000 4000], [0 0], 'r', [0 0],[0 0], [0 9000], 'r');

axis([-2000 6000 -4000 4000 0 9000]);

title('Magnetometers measurement en 3D view - Uncalibrated')
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')

subplot(224);
scatter3(results(:,26),results(:,27),results(:,28), '.');
hold on 
plot3([0 0],[-10 10], [0 0],'r',[-10 10], [0 0], [0 0], 'r', [0 0],[0 0], [-10 10], 'r');
axis([-1.5 1.5 -1.5 1.5]);
title('Magnetometers measurement en 3D view - Calibrated')
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')



scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);
subplot(334);
scatter(results(:,26),results(:,27));
hold on
plot([0 0],[-10 10],'r',[-10 10], [0 0],'r');
title('Magneto readings - XY view - Calibrated');
xlabel('X axis');
ylabel('Y axis');



subplot(337)
scatter(results(:,26),results(:,27));
hold on
plot([0 0],[-10 10],'r',[-10 10], [0 0],'r');axis([-1.5 1.5 -1.5 1.5]);
title('Magneto readings - XY view - Calibrated');
xlabel('X axis');
ylabel('Y axis');

subplot(335)
scatter(results(:,26),results(:,28));
hold on
plot([0 0],[-10 10],'r',[-10 10], [0 0],'r');
title('Magneto readings - XZ view - Calibrated');
xlabel('X axis');
ylabel('Z axis');

subplot(338)
scatter(results(:,26),results(:,28));
axis([-1.5 1.5 -1.5 1.5]);
hold on
plot([0 0],[-10 10],'r',[-10 10], [0 0],'r');
title('Magneto readings - XZ view - Calibrated');
xlabel('X axis');
ylabel('Z axis');

subplot(336)
scatter(results(:,27),results(:,28));
hold on
plot([0 0],[-10 10],'r',[-10 10], [0 0],'r');
title('Magneto readings - YZ view - Calibrated');
xlabel('Y axis');
ylabel('Z axis');

subplot(339)
scatter(results(:,27),results(:,28));
axis([-1.5 1.5 -1.5 1.5]);
hold on
plot([0 0],[-10 10],'r',[-10 10], [0 0],'r');
title('Magneto readings - YZ view - Calibrated');
xlabel('Y axis');
ylabel('Z axis');

subplot(331)
scatter(results(:,7),results(:,8));
hold on
plot([0 0],[-6000 4000],'r',[-5000 10000], [0 0],'r');
title({'Magneto readings - XY view ';'Not calibrated'});
xlabel('X axis');
ylabel('Y axis');

subplot(332)
scatter(results(:,7),results(:,9));
hold on
plot([0 0],[-5000 15000],'r',[-5000 10000], [0 0],'r');
title({'Magneto readings - XZ view ';'Not calibrated'});
xlabel('X axis');
ylabel('Z axis');

subplot(333)
scatter(results(:,8),results(:,9));
hold on
plot([0 0],[-5000 15000],'r',[-5000 5000], [0 0],'r');
title({'Magneto readings - YZ view ';'Not calibrated'});
xlabel('Y axis');
ylabel('Z axis');


size_tmp = size(E)
size_tmp2 = size(F)

size_E = size_tmp(1,1)
size_F = size_tmp2(1,1)

yaw_E = E(:,14);
pitch_E = E(:,15);
roll_E = E(:,16);

yaw_F = F(:,14);
pitch_F = F(:,15);
roll_F = F(:,16);

for k=1:1:size_E
   if (yaw_E(k,1) > 100)
      yaw_E(k,1) = 180 - yaw_E(k,1);
   end
end
for m=1:1:size_F
   if (yaw_F(k,1) > 100)
      yaw_F(k,1) = 180 - yaw_F(k,1);
   end
end

x=1:1:size_E;
figure (10);
subplot(221);
plot(x,yaw_E);
subplot(222);
plot(x,pitch_E);
subplot(223);
plot(x,roll_E);

x=1:1:size_F;
figure (11);
subplot(221);
plot(x,yaw_F);
subplot(222);
plot(x,pitch_F);
subplot(223);
plot(x,roll_F);

subplot(221);
scatter3(C(:,1),C(:,2),C(:,3));
subplot(222);
scatter3(D(:,1),D(:,2),D(:,3));
subplot(223);
scatter3(D(:,29),D(:,30),D(:,31));

acc_xA = A(:,26);
acc_yA = A(:,27);
acc_zA = A(:,28);
acc_xB = B(:,26);
acc_yB = B(:,27);
acc_zB = B(:,28);

mag_xA = A(:,29);
mag_yA = A(:,30);
mag_zA = A(:,31);
mag_xB = B(:,29);
mag_yB = B(:,30);
mag_zB = B(:,31);
mag_xD = D(:,29);
mag_yD = D(:,30);
mag_zD = D(:,31);

quat1B = B(:,19);
quat2B = B(:,20);
quat3B = B(:,21);
quat4B = B(:,22);

yawA = A(:,23);
pitchA = A(:,24);
rollA = A(:,25);
yawB = B(:,23);
pitchB = B(:,24);
rollB = B(:,25);

%raw_mag_x = A(:,1);
%raw_mag_y = A(:,2);
%raw_mag_z = A(:,3);


for i=1:1:size(A)
mag2XA=mag_xA(i,1)*mag_xA(i,1);
mag2YA=mag_yA(i,1)*mag_yA(i,1);
mag2ZA=mag_zA(i,1)*mag_zA(i,1);

mag2A=mag2XA+mag2YA+mag2ZA;
magA(i,1)=sqrt(mag2A);

mag2XB=mag_xB(i,1)*mag_xB(i,1);
mag2YB=mag_yB(i,1)*mag_yB(i,1);
mag2ZB=mag_zB(i,1)*mag_zB(i,1);

mag2B=mag2XB+mag2YB+mag2ZB;
magB(i,1)=sqrt(mag2B);

mag2XD=mag_xD(i,1)*mag_xD(i,1);
mag2YD=mag_yD(i,1)*mag_yD(i,1);
mag2ZD=mag_zD(i,1)*mag_zD(i,1);

mag2D=mag2XD+mag2YD+mag2ZD;
magD(i,1)=sqrt(mag2D);

end
x=1:1:size(magD);
%plot3(raw_mag_x,raw_mag_y,raw_mag_z);
%plot(x,magnitude);
figure(6);
plot(x,magD);
% x=1:1:size(A)
% figure(3)
% subplot(221);
% plot(x,yawA);
% subplot(222);
% plot(x,pitchA);
% subplot(223);
% plot(x,rollA);
% subplot(224);
% scatter(x,magA,'.');

x=1:1:size(B);

% figure(4);
% subplot(221);
% plot(x,yawB);
% subplot(222);
% plot(x,pitchB);
% subplot(223);
% plot(x,rollB);
% x=1:1:size(magB);
% subplot(224);
% scatter(x,magB,'.');

x=1:1:size(B);
figure(5);
subplot(221);
plot(x,quat1B);
subplot(222);
plot(x,quat2B);
subplot(223);
plot(x,quat3B);
subplot(224);
plot(x,quat4B);

size_vect=size(x);
size=size_vect(1,2);

vit_x=zeros(size,1);
vit_y=zeros(size,1);
vit_z=zeros(size,1);

pos_x=zeros(size,1);
pos_y=zeros(size,1);
pos_z=zeros(size,1);

%%%% Integration %%%%%
for i=1:1:size-1
    
    vit_x(i,1)=(acc_x(i+1,1)-acc_x(i,1))/period;
    vit_y(i,1)=(acc_y(i+1,1)-acc_y(i,1))/period;
    vit_z(i,1)=(acc_z(i+1,1)-acc_z(i,1))/period;
    pos_x(i,1)=(vit_x(i+1,1)-vit_x(i,1))/period;
    pos_y(i,1)=(vit_y(i+1,1)-vit_y(i,1))/period;
    pos_z(i,1)=(vit_z(i+1,1)-vit_z(i,1))/period;
end
%plot3(pos_x,pos_y,pos_z);
   
         
