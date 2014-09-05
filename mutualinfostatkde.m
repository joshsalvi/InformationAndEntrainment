function [h p MIS I] = mutualinfostatkde(x,y,varargin)
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

if iscolumn(x) == 0
    x = x';
end
if iscolumn(y) == 0
    y = y';
end

% Randomly shuffle the X and Y data and calculate MI - create surrogates
% Pethel et al, Entropy (2014) 16:2839-2849
for i = 1:iter
    
    [fx wx ux vx] = trans_count(x,0);       % zeroth order markov model
    [fy wy uy vy] = trans_count(y,0);
    xn = whittle_surrogate(fx,wx,ux,vx);    % generate surrogates
    yn = whittle_surrogate(fy,wy,uy,vy);
    MIS(i) = mutualinformation4(xn,yn,0);
    clear z fx fy wx wy ux uy vx vy xn yn
    %{
    z = Shuffle(vertcat(x,y));
    x1 = z(1:length(z)/2);
    y1 = z(length(z)/2+1:length(z));
    MIS(i) = rapidmi(x1,y1);
    %}
end
% Create a kernel density estimate 
[a b] = ksdensity(MIS,0:max(MIS)/10000:max(MIS)*2);
a=a./sum(a);
% Calculate mutual information if not an input
if isempty(varargin{1})
    I = rapidmi(x,y);
else
    I = varargin{1};
end

% Find p-value for single-tailed test
p = sum(a(findnearest(b,I):end));
% Significant
if p < alpha
    h = 1;
else
    h = 0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ir = rapidmi(m,n)

if iscolumn(m) == 0
    m = m';
end
if iscolumn(n) == 0
    n = n';
end

[bwmn,dmn,meshmnm,meshmnn]=kde2d([m n],2^5);
dmn=abs(dmn);dmn = dmn./sum(sum(dmn));dmnlog = log2(dmn);
dmnlog(dmnlog==inf | dmnlog==-inf)=0;

Inh1 = log2(sum(dmn,1));Inh2 = log2(sum(dmn,2));
Inh1(Inh1==inf | Inh1==-inf)=0;Inh2(Inh2==inf | Inh2==-inf)=0;

Ir = sum(sum(dmn.*bsxfun(@minus,bsxfun(@minus,dmnlog,Inh1),Inh2)));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
