function [MIfreq MIfreq_p] = mutualinfospec(y,Fs);
% This function calculates the mutual information from a dataset, y, that
% is organized into an MxN matrix, where M is the length of time and N is
% the number of trials
%
% [MIfreq MIfreq_p] = mutualinfospec(y,Fs);
%
%

sizeY = size(y);

Ymean = mean(y,2);

for i = 1:sizeY(2)
    yn(:,i) = y(:,i) - Ymean;
    [yp(:,i) yf(:,i)] = pwelch(yn(:,i),[],[],[],Fs);
end
[ymp ymf] = pwelch(Ymean,[],[],[],Fs);

ypmean = mean(yp,2);
ysnr = ymp./ypmean;
ind = ysnr~=0;
MIfreq = sum(log2(ysnr(ind)));
end

