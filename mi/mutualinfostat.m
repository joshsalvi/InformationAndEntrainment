function [h p] = mutualinfostat(x,y,varargin)
% This function calculates a p-value for the mutual information value from
% two vectors. The function takes the vectors, reshuffles them "iter"
% times to calculate a probability distribution, and it then calculates a
% p-value based upon the input mutual information.
%
% [h p] = mutualinfostat(x,y,I,iter,alpha)
%
%       I : mutual information input (if empty, will calculate it for you)
%       x : data set 1 (Mx1)
%       y : data set 2 (Mx1)
%       iter : number of iterations (default = 1E4)
%       alpha : alpha level for the statistical test (default=0.05)
%       h : reject or accept the null hypothesis that the mutual
%       information comes from a random distribution as calculated here (1
%       = reject at alpha level, 0 = do not reject)
%       p : p-value associated with the statistical test
%
% [h p] = mutualinfostat(x,y,[],[],[]);
%
%
%   Joshua D. Salvi
%   jsalvi@rockefeller.edu

if isempty(varargin{2})
    iter = 1E4;
else
    iter = varargin{2};
end
if isempty(varargin{3})
    alpha = 0.05;
else
    alpha = varargin{3};
end

% Randomly shuffle the X and Y data and calculate MI
for i = 1:iter
    %MIS(i) = mutualinformation4(Shuffle(x),Shuffle(y),0);
    clear z x1 y1
    z = Shuffle(vertcat(x,y));
    x1 = z(1:length(z)/2);
    y1 = z(length(z)/2+1:length(z));
    MIS(i) = mutualinfo(x1,y1);
end
% Create a kernel density estimate 
[a b] = ksdensity(MIS,0:1e-4:10);
a=a./sum(a);
% Calculate mutual information if not an input
if isempty(varargin{1})
    I = mutualinfo(x,y);
else
    I = varargin{1};
end

% Find p-value for one-tailed test
p = sum(a(findnearest(b,I):end));
% Significant
if p < alpha
    h = 1;
else
    h = 0;
end
end
