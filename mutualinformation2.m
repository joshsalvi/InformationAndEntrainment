function I = mutualinformation2(x,y)
% This function calculates the mutual information for two signals using
% a 2D histogram.
%
%  I = mutualinformation2(x,y)
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

% 2D histogram
[nh,co] = hist3([x y],[nbx nby]);
nh = nh./sum(nh);

% Marginal probabilities
Inh1 = log2(sum(nh,1));
Inh2 = log2(sum(nh,2));

% Mutual information from the joint and marginal probabilities
I = sum(sum(nh.*bsxfun(@minus,bsxfun(@minus,log2(nh),Inh1),Inh2)));

end
