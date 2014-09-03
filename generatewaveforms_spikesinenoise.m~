% Generate inputs for checking distributions
% Indices are as follows
% (1) sine wave
% (2) sin+noise
% (3) sin+brownnoise
% (4) noise
% (5) brownnoise
% (6) upward spikes
% (7) downward spikes

clear all; close all;

Fs = 1e3;
duration=10;
t = linspace(0,duration,duration*Fs);
noise = randn(1,length(t))*0.1;
noisebrown = cumsum(noise);

% Sine wave with and without noise
y(:,1,1) = 5*sin(10*pi*t);
y(:,2,2) = 5*sin(10*pi*t)+noise;
y(:,3,3) = 5*sin(10*pi*t)+noisebrown;

% Noise only
y(:,4,4) = noise;
y(:,5,5) = noisebrown;

% Spiking
timeStepS = 1/Fs;                  
spikesPerS = 50;                  
durationS = duration;                  
times = [0:timeStepS:durationS-timeStepS];
spikes = zeros(1, length(times));
for train = 1:1
vt = rand(size(times));
y(:,6,6) = (spikesPerS*timeStepS) > vt;
end
y(:,7,7) = -y(:,6);

Xd = y;
clear y;
F_rand = [1 2 3 4 5 6 7];
k_rand = [1 1 1 1 1 1 1];

save('/Users/joshsalvi/Downloads/output/artificialwaveforms/waves.mat');
