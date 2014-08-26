function I = miupperbound(y,Fs)
% This function calculates the upper bound of mutual information for two
% signals, a stimulus X and a response Y. Multiple presentations are
% required for this to work, where Y is an m x n matrix, m the length
% of the signal in the time domain and n the number of presentations.
%
%  I = mutual information rate
%  y = m x n matrix of responses
%  Fs = scan rate (Hz)
%  
%
% Joshua D Salvi
% jsalvi@rockefeller.edu
%

sizeY = size(y);

% Find the individual noise in each response in the frequency domain
yavg = mean(y,2);       % find average
NFFT = (2^4)*2^nextpow2(sizeY(1));
[yavgpsd favg] = pwelch(yavg,[],[],NFFT,Fs);    % average psd
for i = 1:sizeY(2)
    y(:,i) = y(:,i) - yavg;         % find individual noise
    [ypsd(:,i),f(:,i)] = pwelch(y(:,i),[],[],NFFT,Fs);
end
ypsdavg = mean(ypsd,2);

% Find the SNR
SNR = ypsdavg./yavgpsd;

% Find the mutual information rate
I = sum(log2(1+SNR)*(f(2,1)-f(1,1)));

end



