clc; %reset command window
clear all; %reset variables
%results=importdata('test_position_settings2.csv', ' ');
results=importdata('test_crot2.csv', ' ');

size2=size(results)
size_file=size2(1)

scrsz = get(groot,'ScreenSize');
figure('Name' ,'YAW ','Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)]);



x=1:1:size_file;
subplot(311)
plot(x,results(:,25));
subplot(312);
plot(x,results(:,26));
subplot(313);
plot(x,results(:,27));
