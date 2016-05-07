close all
clear all
clc


load('data_hb.mat')

[b a] = butter(2, 0.03); 

time    = VarName2;
sensor1 = VarName3*5/1000;
sensor2 = VarName4*5/1000; 
sensor3 = VarName5*5/1000;

time = time/1e6;


N = length(time);
mean_dt = mean(diff(time));
local_t = 10;
local_N = round(local_t/mean_dt);

BPM1 = 0;
BPM2 = 0;
BPM3 = 0;
jCount = 1;

amp_gain = 2;

plotCount =1;
for i=1:1:N
  sensor1_store(jCount) = sensor1(i);
  sensor2_store(jCount) = sensor2(i);
  sensor3_store(jCount) = sensor3(i);

  if mod(i,local_N)==0
    % one second has passed
  
    sensor1_store_filt = amp_gain*filtfilt(b,a,sensor1_store) + 0.5*ones(size(sensor1_store));
    sensor2_store_filt = amp_gain*filtfilt(b,a,sensor2_store) + 0.5*ones(size(sensor1_store));
    sensor3_store_filt = amp_gain*filtfilt(b,a,sensor3_store) + 0.5*ones(size(sensor1_store));
  
    [pks1,locs1] = findpeaks(sensor1_store_filt,'MinPeakHeight',8,'MinPeakDistance',1.5);
    [pks2,locs2] = findpeaks(sensor2_store_filt,'MinPeakHeight',6,'MinPeakDistance',1.5);
    [pks3,locs3] = findpeaks(sensor3_store_filt,'MinPeakHeight',8,'MinPeakDistance',1.5);
    
    BPM1 = length(pks1)/(local_t/60); 
    BPM2 = length(pks2)/(local_t/60); 
    BPM3 = length(pks3)/(local_t/60); 
    
    
    plotCount = plotCount + 1;
    
    delete sensor1_store sensor2_store sensor3_store
    jCount = 1;
  end
  
  BPM1_store(i) = BPM1; 
  BPM2_store(i) = BPM2; 
  BPM3_store(i) = BPM3; 
  
  jCount = jCount + 1;
  
end

sensor1_filt = filtfilt(b,a,sensor1);
sensor2_filt = filtfilt(b,a,sensor2);
sensor3_filt = filtfilt(b,a,sensor3);

sensor1_filt = amp_gain*(sensor1_filt) + 0.5*ones(size(sensor1_filt));
sensor2_filt = amp_gain*(sensor2_filt) + 0.5*ones(size(sensor1_filt));
sensor3_filt = amp_gain*(sensor3_filt) + 0.5*ones(size(sensor1_filt));


[pks1,locs1] = findpeaks(sensor1_filt,'MinPeakHeight',8,'MinPeakDistance',1.5);
[pks2,locs2] = findpeaks(sensor2_filt,'MinPeakHeight',6,'MinPeakDistance',1.5);
[pks3,locs3] = findpeaks(sensor3_filt,'MinPeakHeight',8,'MinPeakDistance',1.5);

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

figure(40)
plot(time,BPM1_store,'Linewidth',2)
hold on
xlabel('time (sec)')
ylabel('BPM 1')

figure(50)
plot(time,BPM2_store,'Linewidth',2)
hold on
xlabel('time (sec)')
ylabel('BPM 2')

figure(60)
plot(time,BPM3_store,'Linewidth',2)
hold on
xlabel('time (sec)')
ylabel('BPM 3')
