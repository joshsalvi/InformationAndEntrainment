p=6;i=1;j=1;

clear rb2 spcount rb sphase sphasehilb xhilb dhilb rdelt rb Xdcrs xdelt 
close all
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
  
  freq =freq_stim;
  i = 1;
  t = linspace(0,3,length(input(:,i,1,1)));
  for i = 1:length(freq)
      y(:,i) = sin(freq(i)*2*pi*t);
  end
  

  %%
for j = 1:a
for i = 1:(logdata.data(j+ramp_loc1-1,iter_loc))
for p = 1:numfreq
% Divide Phase into bins

N = 20;                                                 % Number of bins


dhilb = hilbert(Delta_in_center(1:stim_time(p),p,i,j));                        % Find Hilbert transform of input for phase information
rdelt = angle(dhilb);
%{
rb = min(rdelt):(abs(min(rdelt))+abs(max(rdelt)))/N:max(rdelt);
rb2 = min(rdelt):(abs(min(rdelt))+abs(max(rdelt)))/(N-1):max(rdelt);

rb = rb/max(rb)*pi;
rb2 = rb2/max(rb2)*pi;

rdbin = zeros(length(rdelt),N);

for l = 1:N
    rdbin(:,l) = (rdelt >= rb(l) & rdelt <= rb(l+1));   % Divide phase into bins in time
end

% Find zero crossings of signal
warning off
Xdcrs = find(diff(stim_center(1:stim_time(p),p,i,j)>0)~=0)+1;    % all zero crossings
Xdcrs = find(diff(stim_center(1:stim_time(p),p,i,j)>0)==1)+1;    % positive zero crossings
%}
xhilb = hilbert(stim_center(1:stim_time(p),p,i,j));
xdelt = angle(xhilb);

chilb = hilbert(input_center(1:stim_time(p),p,i,j));
cdelt = angle(chilb);

chilb2 = hilbert(y(1:length(cdelt),p));
cdelt2 = angle(chilb2);


sphasehilb2 = -(xdelt - rdelt);
sphasehilb4 = xdelt;
sphasehilb = -(cdelt - xdelt);
sphasehilb3 = -(cdelt2(1:length(xdelt)/2) - xdelt(1:length(xdelt)/2));

%{
for w = 1:N
    clear spind
    spind(:,w) = find(rdbin(Xdcrs,w)==1);
    spcount(w) = length(spind(:,w));
end

sphase=zeros(max(spcount),N);

for w = 1:N
    sphase(1:spcount(w),w) = rb2(w);
end

sphase(sphase==0)=[];


% Rayleigh's Test

[pval(p,i,j) zval(p,i,j)] = circ_rtest(sphase);
vec_length(p,i,j) = circ_r(sphase);
vec_dir(p,i,j) = circ_mean(sphase);
vec_std(p,i,j) = circ_std(sphase);
vec_sem = vec_std/sqrt(N);
%}
[pvalhilb(p,i,j) zvalhilb(p,i,j)] = circ_rtest(sphasehilb);
vec_lengthhilb(p,i,j) = circ_r(sphasehilb);
vec_dirhilb(p,i,j) = circ_mean(sphasehilb);
vec_stdhilb(p,i,j) = circ_std(sphasehilb);

[pvalhilb2(p,i,j) zvalhilb2(p,i,j)] = circ_rtest(sphasehilb2);
vec_lengthhilb2(p,i,j) = circ_r(sphasehilb2);
vec_dirhilb2(p,i,j) = circ_mean(sphasehilb2);
vec_stdhilb2(p,i,j) = circ_std(sphasehilb2);

[pvalhilb3(p,i,j) zvalhilb3(p,i,j)] = circ_rtest(sphasehilb3);
vec_lengthhilb3(p,i,j) = circ_r(sphasehilb3);
vec_dirhilb3(p,i,j) = circ_mean(sphasehilb3);
vec_stdhilb3(p,i,j) = circ_std(sphasehilb3);

[pvalhilb4(p,i,j) zvalhilb4(p,i,j)] = circ_rtest(sphasehilb4);
vec_lengthhilb4(p,i,j) = circ_r(sphasehilb4);
vec_dirhilb4(p,i,j) = circ_mean(sphasehilb4);
vec_stdhilb4(p,i,j) = circ_std(sphasehilb4);


% Plot polar coordinates
%{
figure(j)
subplot(logdata.data(j,iter_loc),numfreq,p+(i-1)*numfreq)
g=polar(rb2,spcount);   
x = get(g,'Xdata');
y = get(g,'Ydata');
%g=patch(x,y,'r');

figure(a+j)
subaxis(logdata.data(j,iter_loc),numfreq,p+(i-1)*numfreq,'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
set(gca, 'LooseInset', get(gca,'TightInset'))
polar(0,0.2); hold on;
circ_plot(sphase,'hist',[],N,true,true,'linewidth',2,'color','r');  % histogram with normalized count
%}
%{
figure(2*a+j)
subaxis(logdata.data(j,iter_loc),numfreq,p+(i-1)*numfreq,'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
set(gca, 'LooseInset', get(gca,'TightInset'))
polar(0, 0.2); hold on;
circ_plot(sphasehilb,'hist',[],N,true,true,'linewidth',2,'color','r');  % histogram with normalized count
%}
%{
figure(2*a+j)
subaxis(logdata.data(j,iter_loc),numfreq,p+(i-1)*numfreq,'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
set(gca, 'LooseInset', get(gca,'TightInset'))
polar(0, 0.2); hold on;
circ_plot(sphasehilb3,'hist',[],N,true,true,'linewidth',2,'color','r');  % histogram with normalized count
%}
end
end
end

