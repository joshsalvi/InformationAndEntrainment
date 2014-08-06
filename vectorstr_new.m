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
    
end
end
end


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
