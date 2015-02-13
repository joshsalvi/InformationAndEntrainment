function [numpeak, ptau, tauescape, taum, anclusteredsum, at, bn] = kldivisions(y,t,dt,taumin,taumax,Ntau,plotyn)
%
% This function calculates the number of low-energy states (LES) in a time series
% by employing the Kullback-Leibler divergence and k-means clustering. Each
% state is determined by finding a PDF with a width that has been optimized
% to minimize KL divergence between other window sizes. The PDFs of the
% signal are then clustered by KL divergence between them. To determine if
% a cluster denotes an LES, the escape times (the amount of time the signal
% spends in that state) is compared with the window used for each pdf in
% that state. If the windows are significantly smaller than the escape
% times, then the state is LES. 
%
% [numpeak, ptau, tauescape, taum, anclusteredsum, at bn] = kldivisions(y,t,dt,taumin,taumax,Ntau,plotyn)
%
%
% Inputs:
% y : 1xN signal
% t : 1xN time vector
% dt : step size for the PDFs. For a PDF at every point (excluding the edges), dt=1.
% taumin : minimum window size to be surveyed (in points)
% taumax : maximum window size to be surveyed (in points)
% Ntau : number of windows to survey between taumin and taumax
% plotyn : plot the data? (1=yes; 0=no)
%
% Outputs:
% numpeak : number of peaks in your signal for a p-value threshold of 1E-4
% ptau : p values for each of the clusters (the windows versus the escape
% times)
% tauescape : mean escape time for each cluster
% taum : mean window size for each cluster
% anclusteredsum : PDF of each cluster
% at : PDF of entire signal
% bn : x-axis for PDF
%
% Please contact jsalvi@rockefeller.edu with questions.
%


%t=linspace(0,30,5e2);
%y=sin(2*pi*t);
%noise = randn(1,length(t));
%y=y+noise;
%y=noise;
M=dt;
tau=linspace(taumin,taumax,Ntau);
tn=1:M:length(t);tn=round(tn);
fs2=0.01;
Fs = 1/(t(2)-t(1));

for j = 1:length(tn)
clear yn tau2
tau2=tau;ind=1;ind2=1:1:length(tau2);
iter=1;
    while length(tau2) > 1
        clear a b c d yn kldiv2
    for k = 1:length(tau2)
        if tn(j) > tau2(k)/2 && length(y)-tn(j) > tau2(k)/2
        if floor(tau2(k)/2) > 1
            yn{k} = y(tn(j)-floor(tau2(k)/2):tn(j)+floor(tau2(k)/2));
        else 
            yn{k} = 0;
        end
        [a{k}, b{k}] = ksdensity(yn{k},1.5*min(y):fs2:1.5*max(y));
        a{k}=a{k}./sum(a{k});
        a{k}(a{k}==0)=1e-323;
        elseif tn(j) <= tau2(k)/2
        if floor(tau2(k)/2) > 1
            yn{k} = y(1:round(tn(j)+floor(tau2(k)/2)));
        else 
            yn{k} = 0;
        end
        [a{k}, b{k}] = ksdensity(yn{k},1.5*min(y):fs2:1.5*max(y));
        a{k}=a{k}./sum(a{k});
        a{k}(a{k}==0)=1e-323;
        else
        if floor(tau2(k)/2) > 1
            yn{k} = y(tn(j)-floor(tau2(k)/2):end);
        else 
            yn{k} = 0;
        end
        [a{k}, b{k}] = ksdensity(yn{k},1.5*min(y):fs2:1.5*max(y));
        a{k}=a{k}./sum(a{k});
        a{k}(a{k}==0)=1e-323;  
        end
    end
    
            
    for k = 1:length(tau2)
        kldiv2(k) = sum(a{ind}.*(log(a{ind})-log(a{ind2(k)})));
    end
    [c,d]=ksdensity(kldiv2);
    c=c./sum(c);
    [pkc,trc]=PTDetect(c,1e-10);
    if isempty(trc) == 1
        clear yn
        tau2 = tau2(findnearest(kldiv2,0));tau2=tau2(1);
        if tn(j) > tau2/2 && length(y)-tn(j) > tau2/2
            yn = y(tn(j)-floor(tau2/2):tn(j)+floor(tau2/2));
        elseif tn(j) <= tau2/2 && length(y)-tn(j) > tau2/2
            yn = y(1:round(tn(j)+floor(tau2/2)));
        else
            yn = y(tn(j)-floor(tau2/2):end);
        end
        break;
    end
    pkmax=pkc(c(pkc) == max(c(pkc)));
    if isempty(pkmax) ==0
        pkmax=pkmax(1);
        klmax=d(pkmax);
        pkmin = trc(trc>=pkmax);
    else
        pkmax=0;
        klmax=max(d);
    end
    if isempty(pkmin) ==0
        pkmin=pkmin(1);
        klmin=d(pkmin);
        indt=find(kldiv2<=klmin);
    else
        pkmin=min(trc);
        klmin=d(pkmin);
        indt= find(kldiv2>=klmin);   
    end
    taumax=tau2(findnearest(kldiv2,klmax));
    taumax=taumax(1);
    tau2=tau2(indt);
    ind=find(tau2==taumax);
    if isempty(ind) == 1
        ind=1;
    end
    ind2=1:1:length(tau2);
    iter=iter+1;
    end
    ydiv{j}=yn;
    taudiv(j)=tau2/Fs;
    iterdiv(j)=iter;
    [adiv{j}, bdiv{j}] = ksdensity(ydiv{j},1.5*min(y):fs2:1.5*max(y));
end



for k = 1:length(tn)
    if isempty(ydiv{k}) ==0
    [an{k}, bn{k}] = ksdensity(ydiv{k},1.5*min(y):fs2:1.5*max(y));
    an{k}=an{k}./sum(an{k});
    else
        [an{k}, bn{k}] = ksdensity(zeros(1,100),1.5*min(y):fs2:1.5*max(y));
    end
    an{k}(an{k}==0)=1e-323;
end

for j = 1:length(tn)
    for k = 1:length(tn)
        klcluster(j,k) = sum(an{j}.*(log(an{j})-log(an{k})));
    end
end

if plotyn ==1
    figure;
    scatter(tn,taudiv);
end
if plotyn ==1
    set(0,'DefaultAxesColorOrder',cool(length(tn)));
    figure;
    for k = 1:length(tn)
    scatter(tn,klcluster(:,k));hold all
    end
end

N=3;
idxcluster = kmeans(klcluster,N);
for n = 1:N
idxclustermeans(n) = min(mean(klcluster(idxcluster==n)));
end
idxclusterlowmean = find(idxclustermeans == min(idxclustermeans));
for n = 1:N
    indc = find(idxcluster==n);
    for h = 1:length(indc)
        anclustered{n}(:,h) = an{indc(h)};
    end
    anclusteredsum{n} = sum(anclustered{n},2);
    anclusteredsum{n} = anclusteredsum{n}./sum(anclusteredsum{n});
end

if plotyn ==1
figure;
set(0,'DefaultAxesColorOrder',cool(N));
for n = 1:N
    plot(bn{1},anclusteredsum{n}); hold all;
end
[at bt] = ksdensity(y,1.5*min(y):fs2:1.5*max(y));
at=at./sum(at);
hold on;
plot(bt,at,'k--');

figure;
silhouette(klcluster,idxcluster)
end
% Calculate escape time for each time point (how long does it take to
% escape the cluster? Then determine whether that point is an LES by
% determining whether tau < tau_escape ... This will define the correct
% number of states.
%taudiv=taudiv/Fs;
taux = cell(N);
for n = 1:N
    divdiff{n} = diff(idxcluster==n);
    divdiffenter{n} = find(divdiff{n}==1);
    divdiffexit{n} = find(divdiff{n}==-1);
    if isempty(divdiffenter{n}) == 0
    if length(divdiffenter{n}) >= length(divdiffexit{n})
        if divdiffenter{n}(1) < divdiffexit{n}(1)
            for m = 1:length(divdiffexit{n})
                taux{n}(m) = tn(divdiffexit{n}(m)) - tn(divdiffenter{n}(m));
            end
        elseif length(divdiffexit{n}) > 1
            for m = 1:length(divdiffexit{n})-1
                taux{n}(m) = tn(divdiffexit{n}(m+1)) - tn(divdiffenter{n}(m));
            end
        else
            for m = 1:1
                taux{n}(m) = tn(divdiffexit{n}(m+1)) - tn(divdiffenter{n}(m));
            end
        end
    else
            if divdiffenter{n}(1) < divdiffexit{n}(1)
            for m = 1:length(divdiffenter{n})
                taux{n}(m) = tn(divdiffexit{n}(m)) - tn(divdiffenter{n}(m));
            end
            elseif length(divdiffenter{n}) > 1
            for m = 1:length(divdiffenter{n})-1
                taux{n}(m) = tn(divdiffexit{n}(m+1)) - tn(divdiffenter{n}(m));
            end
            else
                for m = 1:1
                taux{n}(m) = tn(divdiffexit{n}(m+1)) - tn(divdiffenter{n}(m));
            end
            end
    end
        taux{n} = taux{n}./Fs;
        tauescape(n) = mean(taux{n});
    else
        tauescape(n) = 0;
    end
end

for n = 1:N
    clear tind
    tind = find(idxcluster==n);
    taum(n) = mean(taudiv(tind));
    [htau{n},ptau(n),citau{n},statstau{n}]=ttest(taudiv(tind),taux{n},'Tail','Left');
end
numpeak=sum(ptau<1e-4);
disp(['Escape times = ' num2str(tauescape)]);
disp(['Mean window size = ' num2str(taum)]);
disp(['p values = ' num2str(ptau)]);
disp(['Number of peaks = ' num2str(numpeak)]);

for n = 1:N
    for k = 1:N
    [hks(n,k), pkst(n,k), kstat(n,k)] = kstest2(anclusteredsum{n},anclusteredsum{k});
    end
end
pkst(pkst==1)=0;
pkstind=pkst<1e-10; % threshold for p-value in kstest

