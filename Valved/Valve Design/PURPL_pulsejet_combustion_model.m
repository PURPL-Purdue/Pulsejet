% Define parameters
close all;
freq = 60; % 60 Hz
period = 1/freq;
amp = 10; %amplitude

% Generate pulse train
sim_t = 0.1; %simulation time (s)
fs = 100000; % Sampling frequency
t = 0:1/fs:sim_t; % \time period for simulation
d = 0:period:sim_t; % Pulse locations

width = 0.1; %width of base of guasian pulse in %
pulse1 = pulstran(t,d,'rectpuls',period/2)*amp; % Creates the pulse train (redundant)
pulse2 = pulstran(t,d,'gauspuls',freq*(1/width),1)*amp;% Creates the pulse train (useful)
pulse2(pulse2 < mean(pulse2)*(1/width)) = 0; %set everyting below threshold to zero
figure()
hold on
plot(t,pulse1,"b") %redundant
plot(t,linspace(mean(pulse2),mean(pulse2),numel(t)),"--c")
plot(t,pulse2,"r")