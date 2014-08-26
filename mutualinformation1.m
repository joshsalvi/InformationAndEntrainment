function I = mutualinformation1(x,y)
% This function calculates the mutual information for two signals.
%
%  I = mutualinformation1(x,y)
%
%  I : mutual information (bits)
%  x,y : input signals
%  
%
% I recommend using the Freedman-Diaconis rule for calculating the
% appropriate bin size.
%
% Joshua D Salvi, jsalvi@rockefeller.edu
%

% Freedman-Diaconis rule to calculate bin size
bwx = 2*iqr(x)/length(x)^(1/3);
nbx = ceil((max(x) - min(x))/bwx);
bwy = 2*iqr(y)/length(y)^(1/3);
nby = ceil((max(y) - min(y))/bwy);

% Calculate mutual information
I = entropy1(x,nbx) + entropy1(y,nby) - jointentropy1(x,y,[nbx nby]);
