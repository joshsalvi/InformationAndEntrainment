function [Sig Time]=HilbertAnalysis()
clear all
%close all

TestT = 200;
w = 2*pi/TestT;
%100 periods
t_range = 100*TestT;
%100 sample points per period
sample_no_per_period = 100;
Time = 0:TestT/sample_no_per_period:t_range;

allplots = 1;
if allplots == 1
plotnum = 8;

for nplot = 1:plotnum
    
if nplot == 1
%Noise
DetSig = zeros(1,numel(Time));
elseif nplot == 2
%Sinusoidal function
DetSig = cos(w.*Time);
elseif nplot == 3
%Symmetric square wave
duty = 50; %Percentange of period signal is positive
DetSig = square(w.*Time,duty);
elseif nplot == 4
%Spiky square wave
duty = 30; %Percentange of period signal is positive
%duty = 70;
DetSig = square(w.*Time,duty);
elseif nplot == 5
%Symmetric sawtooth wave
DetSig = sawtooth(w.*Time,0.5);
elseif nplot == 6
%Asymmetric sawtooth wave
DetSig = sawtooth(w.*Time,0);
%DetSig = sawtooth(w.*Time,1);
elseif nplot == 7
%Asymmetric sawtooth plus square wave signal
DetSig = (sawtooth(w.*Time,0) + square(w.*Time))/2;
%DetSig = (sawtooth(w.*Time,1) - square(w.*Time))/2;
elseif nplot == 8
%Big jump relaxation oscillator shaped signal
SigP = sawtooth(w.*Time,1) - square(w.*Time);
SigP = SigP.*(SigP>=0);
SigN = sawtooth(w.*Time + pi,0) + square(w.*Time + pi);
SigN = SigN.*(SigN<=0);
%SigP = sawtooth(w.*Time,0) + square(w.*Time);
%SigP = SigP.*(SigP>=0);
%SigN = sawtooth(w.*Time + pi,1) - square(w.*Time + pi);
%SigN = SigN.*(SigN<=0);
DetSig = (SigP + SigN)/2;
elseif nplot == 9
%FHN relaxation oscillator shaped signal
SigP = sawtooth(w.*Time,0);
SigP = SigP.*(SigP>=0);
SigN = sawtooth(w.*Time + pi,1);
SigN = SigN.*(SigN<=0);
DetSig = (SigP + SigN)/2;
elseif nplot == 10
%HB relaxation oscillator shaped signal
SigP = sawtooth(w.*Time,1);
SigP = SigP.*(SigP>=0);
SigN = sawtooth(w.*Time + pi,0);
SigN = SigN.*(SigN<=0);
DetSig = (SigP + SigN)/2;
end

Sig = DetSig + 0.1.*randn(length(Time),1)';
figure;
plot(Time,Sig,'r')
hold on
plot(Time,DetSig,'k');
h = gca;
axis([0 10*TestT 1.1*min(Sig) 1.1*max(Sig)])

xlabel('Time','FontSize',24); ylabel('x','FontSize',24);

DetSig = hilbert(DetSig);
Sig = hilbert(Sig);

minR = min(real(Sig));
maxR = max(real(Sig));
minI = min(imag(Sig));
maxI = max(imag(Sig));

%Bin width using the Freedman/Diaconis rule
BinSizeR = 2*iqr(real(Sig))/length(Sig)^(1/3);
NRbins = floor((maxR - minR)/BinSizeR) + 2;
BinSizeI = 2*iqr(imag(Sig))/length(Sig)^(1/3);
NIbins = floor((maxI - minI)/BinSizeI) + 2;

AHist = zeros(NRbins,NIbins);
Rpos = (minR-BinSizeR/2):BinSizeR:(maxR+BinSizeR/2);
Ipos = (minI-BinSizeI/2):BinSizeI:(maxI+BinSizeI/2);

for R = 1:NRbins
    for I = 1:NIbins
    AHist(R,I) = sum(real(Sig)<(Rpos(R)+BinSizeR/2) & real(Sig)>=(Rpos(R)-BinSizeR/2) & imag(Sig)<(Ipos(I)+BinSizeI/2) & imag(Sig)>=(Ipos(I)-BinSizeI/2));
    end
end
figure
surf(Rpos,Ipos,AHist');
view([0 90])
hs = gca;
%colormap(flipud(jet))
shading interp
hold on

period_pts = (numel(Time)-2*sample_no_per_period):(numel(Time)-sample_no_per_period);
%plot 2d curve in plane at height equal to max of Hilbert distribution
plot3(real(DetSig(period_pts)),imag(DetSig(period_pts)),max(max(AHist)).*ones(1,numel(period_pts)),'k','LineWidth',2);
axis square
hold on

%Create one period's worth of data for the vector plot
DetSigR = real(DetSig(period_pts));
DetSigI = imag(DetSig(period_pts));
DetSigRg = gradient(DetSigR);
DetSigIg = gradient(DetSigI);
%Remove the points with very large gradients
VecPlot_pts = find(abs(DetSigRg) < 5*std(DetSigRg) & abs(DetSigIg) < 5*std(DetSigIg));
%Reduce the density of vectors
for r = 1:5
VecPlot_pts = VecPlot_pts(find(mod(find(VecPlot_pts<numel(DetSigRg)),2)));
end
vec_scale = 1;
%plot 2d vector plot in plane at height equal to max of Hilbert distribution
quiver3(DetSigR(VecPlot_pts),DetSigI(VecPlot_pts),max(max(AHist)).*ones(1,numel(VecPlot_pts)),DetSigRg(VecPlot_pts),DetSigIg(VecPlot_pts),zeros(1,numel(VecPlot_pts)),vec_scale,'k','LineWidth',2);
xlabel('x','FontSize',24); ylabel('H(x)','FontSize',24);

if nplot == 1
title(h,'White Noise Signal','FontSize',12); 
title(hs,'White Noise Distribution','FontSize',12);
elseif nplot == 2
title(h,'Sinusoidal Signal','FontSize',12); 
title(hs,'Sinusoidal Distribution','FontSize',12);
elseif nplot == 3
title(h,'Square Wave Signal','FontSize',12); 
title(hs,'Square Wave Distribution','FontSize',12);
elseif nplot == 4
title(h,'Spiky Square Wave Signal','FontSize',12); 
title(hs,'Spiky Square Wave Distribution','FontSize',12);
elseif nplot == 5
title(h,'Symmetric Sawtooth Signal','FontSize',12); 
title(hs,'Symmetric Sawtooth Distribution','FontSize',12);
elseif nplot == 6
title(h,'Asymmetric Sawtooth Signal','FontSize',12); 
title(hs,'Asymmetric Sawtooth Distribution','FontSize',12);
elseif nplot == 7
title(h,'Asymmetric Sawtooth Plus Square Wave Signal','FontSize',12); 
title(hs,'Asymmetric Sawtooth Plus Square Wave Distribution','FontSize',12);
elseif nplot == 8
title(h,'Big Jump Relaxation Oscillator Shaped Signal','FontSize',12); 
title(hs,'Big Jump Relaxation Oscillator Shaped Distribution','FontSize',12);
elseif nplot == 9
title(h,'FHN Relaxation Oscillator Shaped Signal','FontSize',12); 
title(hs,'FHN Relaxation Oscillator Shaped Distribution','FontSize',12);
elseif nplot == 10
title(h,'HB Relaxation Oscillator Shaped Signal','FontSize',12); 
title(hs,'HB Relaxation Oscillator Shaped Distribution','FontSize',12);
end

clear Sig DetSig
end
end

%%%%%%%Simulations%%%%%%%%%

TimeS = 0:TestT/sample_no_per_period:2*t_range;
Time = TimeS(find(TimeS == t_range):find(TimeS == 2*t_range));

plotnum = 5;
for nplot = 1:5
    
noise_type = 2;
if noise_type == 1
%Add noise to deterministic trajectory
    
options = odeset('MStateDependence','none','MaxStep', 1);

if nplot == 1
    sol = ode45(@HBmodel,TimeS,[1 1],options);
elseif nplot == 2
    sol = ode45(@HBmodel2,TimeS,[1 1],options);
elseif nplot == 3
    sol = ode45(@HBmodel3,TimeS,[1 1],options);
elseif nplot == 4
    sol = ode45(@FHN,TimeS,[1 1],options);
elseif nplot == 5
    sol = ode45(@VdP,TimeS,[1 1],options);
end

[Y,YP] = deval(sol,TimeS);

%Find the Hilbert transform of the steady state signal
DetSig = Y(1,find(TimeS == t_range):find(TimeS == 2*t_range));

%Noise amplitude is 0.1 times the oscillation amplitude
%Normally distributed noise mean = 0, std = 0.1*max(DetSig).
Sig = DetSig + 0.1.*randn(length(Time),1)';

elseif noise_type == 2
%Stochastic Ito integration

if nplot == 1
    [Xdet Xem] = HBmodelEM(TimeS);
elseif nplot == 2
    [Xdet Xem] = HBmodel2EM(TimeS);
elseif nplot == 3
    [Xdet Xem] = HBmodel3EM(TimeS);
elseif nplot == 4
    [Xdet Xem] = FHNEM(TimeS);
elseif nplot == 5
    [Xdet Xem] = VdPEM(TimeS);
end

%Find the Hilbert transform of the steady state signal
DetSig = Xdet(1,find(TimeS == t_range):find(TimeS == 2*t_range));
Sig = Xem(1,find(TimeS == t_range):find(TimeS == 2*t_range));

maxDetSig = max(DetSig);
DetSig = DetSig/maxDetSig(1,1);
Sig = Sig/maxDetSig(1,1);
end

figure;
plot(Time,Sig,'r')
hold on
plot(Time,DetSig,'k');
h = gca;
axis([Time(end)-10*TestT Time(end) 1.1*min(Sig) 1.1*max(Sig)])

xlabel('Time','FontSize',24); ylabel('x','FontSize',24);

DetSig = hilbert(DetSig);
Sig = hilbert(Sig);

minR = min(real(Sig));
maxR = max(real(Sig));
minI = min(imag(Sig));
maxI = max(imag(Sig));

%Bin width using the Freedman/Diaconis rule
BinSizeR = 2*iqr(real(Sig))/length(Sig)^(1/3);
NRbins = floor((maxR - minR)/BinSizeR) + 2;
BinSizeI = 2*iqr(imag(Sig))/length(Sig)^(1/3);
NIbins = floor((maxI - minI)/BinSizeI) + 2;

AHist = zeros(NRbins,NIbins);
Rpos = (minR-BinSizeR/2):BinSizeR:(maxR+BinSizeR/2);
Ipos = (minI-BinSizeI/2):BinSizeI:(maxI+BinSizeI/2);

for R = 1:NRbins
    for I = 1:NIbins
    AHist(R,I) = sum(real(Sig)<(Rpos(R)+BinSizeR/2) & real(Sig)>=(Rpos(R)-BinSizeR/2) & imag(Sig)<(Ipos(I)+BinSizeI/2) & imag(Sig)>=(Ipos(I)-BinSizeI/2));
    end
end
figure
surf(Rpos,Ipos,AHist');
view([0 90])
hs = gca;
%colormap(flipud(jet))
shading interp
hold on

period_pts = (numel(Time)-2*sample_no_per_period):(numel(Time)-sample_no_per_period);
%plot 2d curve in plane at height equal to max of Hilbert distribution
plot3(real(DetSig(period_pts)),imag(DetSig(period_pts)),max(max(AHist)).*ones(1,numel(period_pts)),'k','LineWidth',2);
axis square
hold on

%Create one period's worth of data for the vector plot
DetSigR = real(DetSig(period_pts));
DetSigI = imag(DetSig(period_pts));
DetSigRg = gradient(DetSigR);
DetSigIg = gradient(DetSigI);
%Remove the points with very large gradients
VecPlot_pts = find(abs(DetSigRg) < 5*std(DetSigRg) & abs(DetSigIg) < 5*std(DetSigIg));
%Reduce the density of vectors
for r = 1:5
VecPlot_pts = VecPlot_pts(find(mod(find(VecPlot_pts<numel(DetSigRg)),2)));
end
vec_scale = 1;
%plot 2d vector plot in plane at height equal to max of Hilbert distribution
quiver3(DetSigR(VecPlot_pts),DetSigI(VecPlot_pts),max(max(AHist)).*ones(1,numel(VecPlot_pts)),DetSigRg(VecPlot_pts),DetSigIg(VecPlot_pts),zeros(1,numel(VecPlot_pts)),vec_scale,'k','LineWidth',2);
xlabel('x','FontSize',24); ylabel('H(x)','FontSize',24);

if nplot == 1
title(h,'HB Model Signal','FontSize',12); 
title(hs,'HB Model Distribution','FontSize',12);
elseif nplot == 2
title(h,'HB Model Different Parameters Signal','FontSize',12); 
title(hs,'HB Model Different Parameters Distribution','FontSize',12);
elseif nplot == 3
title(h,'Alternate HB Model Signal','FontSize',12); 
title(hs,'Alternate HB Model Distribution','FontSize',12);
elseif nplot == 4
title(h,'FHN Signal','FontSize',12); 
title(hs,'FHN Distribution','FontSize',12);
elseif nplot == 5
title(h,'VdP Signal','FontSize',12); 
title(hs,'VdP Distribution','FontSize',12);
end

clear Sig DetSig

end
end

%Negative Stiffness
function dx = HBmodel(t,x)
%Parameter set one

%Increasing a increases the fraction of a period during which
%adaptation occurs
a = 5;
%b > 1 has unbounded solutions
b = 1;
tau = 25;
Fc = 0;
gamma = 1;
k = 1.5;

dx = zeros(2,1); 
dx(1) = (-k*x(1) + a*(x(1)-x(2)) - (x(1)-x(2))^3 + Fc)/gamma;
dx(2) = (b*x(1) - x(2))/tau;
 
end

function dx = HBmodel2(t,x)
%Parameter set two

%Increasing a increases the fraction of a period during which
%adaptation occurs
a = 3.5;
%b > 1 has unbounded solutions
b = 0.5;
tau = 60;
Fc = 0;
gamma = 0.5;
k = 2.5;

dx = zeros(2,1); 
dx(1) = (-k*x(1) + a*(x(1)-x(2)) - (x(1)-x(2))^3 + Fc)/gamma;
dx(2) = (b*x(1) - x(2))/tau;
 
end

function dx = HBmodel3(t,x)
%Alternate HB model

%Increasing a increases the fraction of a period during which
%adaptation occurs
a = 2;
b = -2;
tau = 30;
Fc = 0;
gamma = 1;
k = 1;

dx = zeros(2,1); 
dx(1) = (-k*x(1) + a*x(1) - (x(1)-x(2))^3 + Fc)/gamma;
dx(2) = (b*x(1) - x(2))/tau;
 
end

%Negative Stiffness
function dx = FHN(t,x)

a = 4;
b = 0.1;
tau = 20;
Fc = 0;
gamma = 1;
k = 3.5;

dx = zeros(2,1); 
dx(1) = (-k*x(1) + a*(x(1)-x(2)) - x(1)^3 + Fc)/gamma;
dx(2) = (b*x(1) - x(2))/tau;
 
end

%Negative Damping
function dx = VdP(t,x)

%Period is greater than 2*pi/sqrt(k/m)
m = 2;
%Increasing gamma increases nonlinearity
gamma = -5;
k = 0.1;
Fc = 0;

dx = zeros(2,1); 
dx(1) = x(2);
dx(2) = (-gamma*(1-x(1)^2)*x(2) - k*x(1) + Fc)/m;
 
end
