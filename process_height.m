clc;
clear all;
results=importdata('/Users/robintherond/Documents/Polytech/PIFE/instantWave/Tests/height/test_height.csv', ' ');
size2=size(results);
size_file=size2(1);

height=results(:,1);
x=1:1:size_file;
figure(1)
subplot(211)
plot(x,height)

min_h=min(height);
max_h=max(height);

slope=max_h-min_h;



for i=1:5:size_file-9  
    height(i)
    height(i+1)
    
    i=i+10;
end
subplot(212)
plot(x,height)