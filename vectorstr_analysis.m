% This code analyzes hair-bundle data in order to calculate the Mutual
% Information, Vector Strength, and the comparisons between them. It also
% performs an optional check to see if each parameter converges.

%% Center the input (optional)
for j = 1:a
for i = 1:(logdata.data(j,iter_loc))
for p = 1:numfreq
    input_z(:,p,i,j)=smooth(input(:,p,i,j),length(input(:,p,i,j))/N);
    input_center(:,p,i,j) = input(:,p,i,j) - input_z(:,p,i,j);
end
end
end
  
%% Generate a sine wave of the same frequency and lag as the stimulus
clear y;
freq(1:numfreq) =freq_stim;
i = 1;
t = linspace(0,3,length(input(:,i,1,1)));
for i = 1:length(freq)
    y(:,i) = sin(freq(i)*2*pi*t);
end

%% Calculate the vector strength, associated statistics, and mutual information

% Initialize variables
VS(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
rayleightest_p(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kde(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kde_p_chi2(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kde_p(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
selfMI(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_norm(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kdedispl_p(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kdedispl(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kdephase(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;
MI_kdephase_p(numfreq,(logdata.data(j+ramp_loc1-1,iter_loc)),a) = 0;

dwnspl = 100;            % downsampling factor: ensure it is still greater than the Nyquist
%{
for i = 1:length(freq)
    [~,selfMI_kdedispl_p(i),~,selfMI_kdedispl(i)] = mutualinfostatkde(y(1:stim_time(i),i),y(1:stim_time(i),i),[],10,[],2^10,1,2,dwnspl);
    [~,selfMI_kdephase_p(i),~,selfMI_kdephase(i)] = mutualinfostatkdephase(y(1:stim_time(i),i),y(1:stim_time(i),i),[],10,[],2^10,1,2,dwnspl);
end
%}

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    clear x11 y11
    x11 = stim_center(1:stim_time(p),p,i,j);
    y11 = y(1:length(x11),p);
    
    % Vector strength and Rayleigh test for nonuniformity
    [VS(p,i,j), rayleightest_p(p,i,j)] = vscalc2(x11,y11,1);        
    
    % Mutual information from displacement histograms
    %[MI_displ(p,i,j) MI_displ_p(p,i,j)] = mutualinformation2(x11,y11,0);
    
    % Mutual information from phase
    %MI_phase(p,i,j) = mutualinformation3(x11,y11,0);  
    
    % Mutual information from kernel density estimates
    %[MI_kde(p,i,j), MI_kde_p_chi2(p,i,j)] = mutualinformation4(x11,y11,0);
    %[~, MI_kde_p(p,i,j)] = mutualinfostat(x11,y11,[],10,[]);
    
    % Normalized mutual information
    %selfMI(p,i,j) = mutualinformation4(y11,y11,0);
    %MI_norm(p,i,j) = MI_kde(p,i,j)/selfMI(p,i,j);
    %[~, MI_norm_p(p,i,j)] = mutualinfostatnorm(x11,y11,[],10,[]);
    
    % Use another method for mutual information
    %MI_mut(p,i,j) = mutualinfo(x11,y11);
    %[~,MI_mut_p(p,i,j)] = mutualinfostat(x11,y11,[],100,[]);
    
    % Alternatively, one can solve for the MI using only one function
    [~,MI_kdedispl_p(p,i,j),~,MI_kdedispl(p,i,j)] = mutualinfostatkde(x11,y11,[],10,[],2^10,1,2,dwnspl);
    [~,MI_kdephase_p(p,i,j),~,MI_kdephase(p,i,j)] = mutualinfostatkdephase(x11,y11,[],10,[],2^10,1,2,dwnspl);
    %MI_kderapid_norm(p,i,j) = MI_kderapid(p,i,j)/selfMI_kderapid(p);
end
end
end

% Find the averages
averageyn = 1;
% average all?
averageallyn = 0;
if averageallyn == 1;
    avgind = 1:a;
else
    avgind = [1 2 3];
end

if averageyn ==1
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    %VS_mean(p,i) = mean(VS(p,i,avgind));
    %VS_sem(p,i) = std(VS(p,i,avgind))/sqrt(length(avgind));
    %MI_displ_mean(p,i) = mean(MI_displ(p,i,avgind));
    %MI_displ_sem(p,i) = std(MI_displ(p,i,avgind))/sqrt(length(avgind));   
    %MI_phase_mean(p,i) = mean(MI_phase(p,i,avgind));
    %MI_phase_sem(p,i) = std(MI_phase(p,i,avgind))/sqrt(length(avgind));
    %MI_norm_mean(p,i) = mean(MI_norm(p,i,avgind));
    %MI_norm_sem(p,i) = std(MI_norm(p,i,avgind))/sqrt(length(avgind));
    %MI_mut_mean(p,i) = mean(MI_mut(p,i,avgind));
    %MI_mut_sem(p,i) = std(MI_mut(p,i,avgind))/sqrt(length(avgind));
    %MI_kde_mean(p,i) = mean(MI_kde(p,i,avgind));
    %MI_kde_sem(p,i) = std(MI_kde(p,i,avgind))/sqrt(length(avgind));
    MI_kdedispl_mean(p,i) = mean(MI_kdedispl(p,i,avgind));
    MI_kdedispl_sem(p,i) = std(MI_kdedispl(p,i,avgind))/sqrt(length(avgind));
    MI_kdephase_mean(p,i) = mean(MI_kdephase(p,i,avgind));
    MI_kdephase_sem(p,i) = std(MI_kdephase(p,i,avgind))/sqrt(length(avgind));
end
end
end
end

%% Does VS converge?

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
for k = 1:100           % 1%-100% of signal
    clear xhilb xhilb_eiphi yhilb yhilb_eiphi x11 y11
    
    % Step through window lengths
    x11 = stim_center(1:stim_time(p)*k/100,p,i,j);
    y11 = y(1:length(x11),p);
    
    % Vector strength calculation
    [VS_conv(k,p,i,j), rayleightest_p(k,p,i,j)] = vscalc2(x11,y11,1);
    
end
end
end
end

for p = 1:numfreq
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for k = 1:100
VS_conv_mean(k,p,i) = mean(VS_conv(k,p,i,:));
VS_conv_sem(k,p,i) = std(VS_convVSc(k,p,i,:))/sqrt(a);
end
end 
end

figure;
set(0,'DefaultAxesColorOrder',autumn(numfreq));
for p = 1:numfreq
    %errorbar(300*(1:100),meanMIconv_all(:,p,i),semMIconv_all(:,p,i));
    plot(300*(1:100)/Fs,VS_conv_mean(:,p,i));
    xlabel('Window Size (sec)');ylabel('Vector Strength');
    hold all;
end

%% Does MI converge?

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
for k = 1:100       % 1%-100% of signal
    clear x11 y11

    % Step through window lengths
    x11 = stim_center(1:stim_time(p)*k/100,p,i,j);
    y11 = y(1:length(x11),p);    
    
    MI_displ_conv(k,p,i,j) = mutualinformation2(x11,y11,0);
    MI_phase_conv(k,p,i,j) = mutualinformation3(x11,y11,0);  
    MI_kde_conv(k,p,i,j)   = mutualinformation4(x11,y11,0);  
    
end
end
end
end

% Find averages
for p = 1:numfreq
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for k = 1:100
MI_displ_conv_mean(k,p,i) = mean(MI_displ_conv(k,p,i,:));
MI_displ_conv_sem(k,p,i) = std(MI_displ_conv(k,p,i,:))/sqrt(a);
MI_phase_conv_mean(k,p,i) = mean(MI_phase_conv(k,p,i,:));
MI_phase_conv_sem(k,p,i) = std(MI_phase_conv(k,p,i,:))/sqrt(a);
end
end 
end

figure;
set(0,'DefaultAxesColorOrder',autumn(numfreq));
for p = 1:numfreq
    plot(300*(1:100)/Fs,MI_displ_conv_mean(:,p,i));
    xlabel('Window Size (sec)');ylabel('Mutual Information (bits) from displacement');
    hold all;
end
figure;
set(0,'DefaultAxesColorOrder',autumn(numfreq));
for p = 1:numfreq
    plot(300*(1:100)/Fs,MI_phase_conv_mean(:,p,i));
    xlabel('Window Size (sec)');ylabel('Mutual Information (bits) from phase');
    hold all;
end

%% Plot mutual information and vector strength together
% Note: you will want to calculate the mean and sem of the cases listed
% first.

% SemiLog? (1=yes)
logyn = 1;

figure;
[hAx,hLine1,hLine2]=plotyy(Fcs(:,1),VS_mean,Fcs(:,1),MI_displ_mean);
hold(hAx(1),'on');errorbar(hAx(1),Fcs(:,1),VS_mean,VS_sem);
hold(hAx(2),'on');errorbar(hAx(2),Fcs(:,1),MI_displ_mean,MI_displ_sem);
ylabel(hAx(1),'Vector Strength') % left y-axis
ylabel(hAx(2),'Mutual Information (bits) from displacement') % right y-axis
if logyn == 1
    set(hAx(1),'xscale','log')
    set(hAx(2),'xscale','log')
end
%{
figure;
[hAx,hLine1,hLine2]=plotyy(amp_stim,VS_mean,amp_stim,MI_phase_mean);
hold(hAx(1),'on');errorbar(hAx(1),amp_stim,VS_mean,VS_sem);
hold(hAx(2),'on');errorbar(hAx(2),amp_stim,MI_phase_mean,MI_phase_sem);
ylabel(hAx(1),'Vector Strength') % left y-axis
ylabel(hAx(2),'Mutual Information (bits) from phase') % right y-axis
%}


figure;
[hAx,hLine1,hLine2]=plotyy(Fcs(:,1),VS_mean,Fcs(:,1),MI_norm_mean);
hold(hAx(1),'on');errorbar(hAx(1),Fcs(:,1),VS_mean,VS_sem);
hold(hAx(2),'on');errorbar(hAx(2),Fcs(:,1),MI_norm_mean,MI_norm_sem);
ylabel(hAx(1),'Vector Strength') % left y-axis
ylabel(hAx(2),'Normalized Mutual Information') % right y-axis
if logyn == 1
    set(hAx(1),'xscale','log')
    set(hAx(2),'xscale','log')
end
