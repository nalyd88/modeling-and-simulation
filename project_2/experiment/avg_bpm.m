close all
clear all
clc


load('data_hb.mat')

[b a] = butter(2, 0.06); 

time    = VarName2;
sensor1 = VarName3*5/1000;
sensor2 = VarName4*5/1000; 
sensor3 = VarName5*5/1000;

time = time/1e6;


sensor1_filt = filtfilt(b,a,sensor1);
sensor2_filt = filtfilt(b,a,sensor2);
sensor3_filt = filtfilt(b,a,sensor3);

mean1 = mean(sensor1_filt);
mean2 = mean(sensor2_filt);
mean3 = mean(sensor3_filt);
var1 = var(sensor1_filt);
var2 = var(sensor2_filt);
var3 = var(sensor3_filt);

sensor1_filt = sensor1_filt + 0.5*ones(size(sensor1_filt));
sensor2_filt = sensor1_filt + 0.5*ones(size(sensor1_filt));
sensor3_filt = sensor1_filt + 0.5*ones(size(sensor1_filt));


[pks1,locs1] = findpeaks(sensor1_filt,'MinPeakHeight',5,'MinPeakDistance',1.5);
[pks2,locs2] = findpeaks(sensor2_filt,'MinPeakHeight',5.5,'MinPeakDistance',1.5);
[pks3,locs3] = findpeaks(sensor3_filt,'MinPeakHeight',5.5,'MinPeakDistance',1.5);

clc
BPM1 = length(pks1)/((time(end)-time(1))/60);
BPM2 = length(pks2)/((time(end)-time(1))/60);
BPM3 = length(pks3)/((time(end)-time(1))/60);

disp('\n \n BPM from sensor 1: \n')
disp(BPM1)

disp('\n \n BPM from sensor 2: \n')
disp(BPM2)

disp('\n \n BPM from sensor 3: \n')
disp(BPM3)

    
figure(1)
plot(time,sensor1_filt,'Linewidth',2)
hold on
plot(time(locs1),pks1,'ro','Linewidth',2)
xlabel('time (sec)')
ylabel('Voltage (V)')
hold on

figure(2)
plot(time,sensor2_filt,'Linewidth',2)
hold on
plot(time(locs2),pks2,'ro','Linewidth',2)
xlabel('time (sec)')
ylabel('Voltage (V)')

figure(3)
plot(time,sensor3_filt,'Linewidth',2)
hold on
plot(time(locs3),pks3,'ro','Linewidth',2)
xlabel('time (sec)')
ylabel('Voltage (V)')
