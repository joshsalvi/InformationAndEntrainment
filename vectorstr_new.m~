%%
  for j = 1:a
  for i = 1:(logdata.data(j,iter_loc))
 for p = 1:numfreq
input_z(:,p,i,j)=smooth(input(:,p,i,j),length(input(:,p,i,j))/N);
input_center(:,p,i,j) = input(:,p,i,j) - input_z(:,p,i,j);
 end
  end
  end
  
  %%
  clear y;
  freq(1:numfreq) =freq_stim;
  i = 1;
  t = linspace(0,3,length(input(:,i,1,1)));
  for i = 1:length(freq)
      y(:,i) = sin(freq(i)*2*pi*t);
  end
  
  
  %%
  
  warning off;
  
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    xhilb = hilbert(stim_center(1:stim_time(p),p,i,j));
    %xdelt = angle(xhilb);
    xhilb_eiphi = xhilb./abs(xhilb);                      % vector strength in new way, normalize all lengths to 1
    
    
    yhilb = hilbert(y(1:length(xhilb),p));
    %cdelt2 = angle(chilb2);
    yhilb_eiphi = yhilb./abs(yhilb);                      % vector strength in new way, normalize all lengths to 1
   
    
    VS2(p,i,j) = abs(sum((xhilb_eiphi./yhilb_eiphi))/length(xhilb));
    
    
    % RAYLEIGH TEST
    %[VS_rayleigh_p(p,i,j) VS_rayleigh_z(p,i,j)] = circ_rtest(angle(xhilb_eiphi./yhilb_eiphi));
    N = length(xhilb)*10/4;
    VS2_n(p,i,j) = VS2(p,i,j)*N;
    VS_rayleigh_p(p,i,j) = exp(sqrt(1+4*N+4*(N^2-VS2_n(p,i,j)^2))-(1+2*N));
    
    % V TEST
    N = length(xhilb)*10/4;
    VS2_n(p,i,j) = VS2(p,i,j)*N;
    VS_vtest_V(p,i,j) = VS2_n(p,i,j)*cos(mean(angle(xhilb_eiphi./yhilb_eiphi))-circ_mean(angle(xhilb_eiphi./yhilb_eiphi)));
    VS_vtest_mu(p,i,j) = VS_vtest_V(p,i,j)*sqrt(2/N);
    VS_vtest_p(p,i,j) = 1 - normcdf(VS_vtest_mu(p,i,j));
    
    %O OMNIBUS TEST
    [VS_omnibus_p(p,i,j)] = circ_otest(angle(xhilb_eiphi./yhilb_eiphi),1);
    
    % Calculate mutual information
    MI(p,i,j) = mutualinformation2(stim_center(1:stim_time(p),p,i,j),y(1:length(xhilb),p));
    
end
end
end

%% Does VS converge?

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
for k = 1:100
    clear xhilb xhilb_eiphi yhilb yhilb_eiphi
    xhilb = hilbert(stim_center(1:stim_time(p)*k/100,p,i,j));
    xhilb_eiphi = xhilb./abs(xhilb);
    yhilb = hilbert(y(1:length(xhilb),p));
    yhilb_eiphi = yhilb./abs(yhilb);
    VSc(k,p,i,j) =  abs(sum((xhilb_eiphi./yhilb_eiphi))/length(xhilb));
end
end
end
end

for p = 1:numfreq
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for k = 1:100
meanVSconv_all(k,p,i) = mean(VSc(k,p,i,:));semVSconv_all(k,p,i) = std(VSc(k,p,i,:))/sqrt(4);
end
end 
end

figure;
set(0,'DefaultAxesColorOrder',autumn(numfreq));
for p = 1:numfreq
    %errorbar(300*(1:100),meanMIconv_all(:,p,i),semMIconv_all(:,p,i));
    plot(300*(1:100)/Fs,meanVSconv_all(:,p,i));
    xlabel('Window Size (sec)');ylabel('Vector Strength');
    hold all;
end
%% Does MI converge?

for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
for k = 1:100
MI2(k,p,i,j) = mutualinformation2(stim_center(1:stim_time(p)*k/100,p,i,j),y(1:stim_time(p)*k/100,p));
end
end
end
end

for p = 1:numfreq
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for k = 1:100
meanMIconv_all(k,p,i) = mean(MI2(k,p,i,:));semMIconv_all(k,p,i) = std(MI2(k,p,i,:))/sqrt(4);
end
end 
end

figure;
set(0,'DefaultAxesColorOrder',autumn(numfreq));
for p = 1:numfreq
    %errorbar(300*(1:100),meanMIconv_all(:,p,i),semMIconv_all(:,p,i));
    plot(300*(1:100)/Fs,meanMIconv_all(:,p,i));
    xlabel('Window Size (sec)');ylabel('Mutual Information');
    hold all;
end


%% PLOT mutual information and vector strength

[hAx,hLine1,hLine2]=plotyy(amp_stim,meanVS2_all,amp_stim,meanMI_all)
ylabel(hAx(1),'Vector Strength') % left y-axis
ylabel(hAx(2),'Mutual Information (bits)') % right y-axis


%%

for i = 1:10
    for j = 1:5
        MEAN_VS2_278(i,j) = mean(VS2(i,j,[2 7 8]));SEM_VS2_278(i,j) = std(VS2(i,j,[2 7 8]))/sqrt(3);
    end
end

%% Calculate Entrainment


for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
    
    entrainment(p,i,j) = ((abs(stim_fft1(p,i,j)))/sqrt(2))/std(stim_center(1:round(cycles/freq_stim(p)*1e4),p,i,j));;
end
end
end
