function [h p MIS I] = mutualinfostatkdephase(x,y,varargin)
% This function calculates a p-value for the mutual information value from
% two vectors. The function takes the vectors, reshuffles them "iter"
% times to calculate a probability distribution, and it then calculates a
% p-value based upon the input mutual information.
%
% [h p] = mutualinfostatkdephase(x,y,I,iter,alpha,bins,markovreps,maxorder)
%
%       I : mutual information input (if empty, will calculate it for you)
%       x : data set 1 (Mx1)
%       y : data set 2 (Mx1)
%       iter : number of iterations (default = 1E4)
%       alpha : alpha level for the statistical test (default=0.05)
%       bins : number of bins for meshgrid
%       markovreps : number of times to repeat the Markov process (default
%       = 100)
%       maxorder : maximum Markov order to be used (default = 2)
%       h : reject or accept the null hypothesis that the mutual
%       information comes from a random distribution as calculated here (1
%       = reject at alpha level, 0 = do not reject)
%       p : p-value associated with the statistical test
%
% [h p] = mutualinfostat(x,y,[],[],[],[],[]);
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
if isempty(varargin{4})
    bins = 2^4;
else
    bins = varargin{4};
end
if isempty(varargin{5})
    markovreps = 100;
else
    markovreps = varargin{5};
end
if isempty(varargin{6})
    maxorder = 2;
else
    maxorder = varargin{6};
end
if isempty(varargin{7})
    downsample = 5;
else
    downsample = varargin{7};
end


if iscolumn(x) == 0
    x = x';
end
if iscolumn(y) == 0
    y = y';
end


% Randomly shuffle the X and Y data and calculate MI - create surrogates
% Pethel et al, Entropy (2014) 16:2839-2849
% Determine the order of your Markov model, up to order 2
x=x(round(1:downsample:end));y=y(round(1:downsample:end));              % Downsample by a factor of scan rate / filter freq (10 kHz / 2 kHz)

xh = hilbert(x);                                     % find analytic signal
xph = atan2(imag(xh),real(xh));                      % instantaneous phase
yh = hilbert(y);
yph = atan2(imag(yh),real(yh));                      


[px] = MarkovOrderTests(xph,markovreps,maxorder);     % determine markov order
[py] = MarkovOrderTests(yph,markovreps,maxorder);
rx=find(px>0.05);ry=find(py>0.05);
if isempty(rx) == 1         % if maximum markov order is not large enough, set to the next highest order
    rx = maxorder+2;
end
if isempty(ry) == 1
    ry = maxorder+2;
end
rx=rx(1);ry=ry(1);    % select the smallest Markov order that was found from the above algorithm

[fx, wx, ux, vx] = trans_count(xph,rx-1);       % Nth order markov model
[fy, wy, uy, vy] = trans_count(yph,ry-1);

for i = 1:iter
    clear z xn yn
    xnph = whittle_surrogate(fx,wx,ux,vx);                  % generate surrogates
    ynph = whittle_surrogate(fy,wy,uy,vy);
    %xnh = hilbert(xn);                                     % find analytic signal
    %xnph = atan2(imag(xnh),real(xnh));                     % instantaneous phase
    %ynh = hilbert(yn);
    %ynph = atan2(imag(ynh),real(ynh)); 
    
    if length(xnph) > length(ynph)
        MIS(i) = rapidmi(xnph(1:length(ynph)),ynph,bins);
    elseif length(xnph) <= length(ynph)
        MIS(i) = rapidmi(xnph,ynph(1:length(xnph)),bins);
    end
end

minlength=min([length(xnph),length(ynph)]);
% Create a kernel density estimate 
[a, b] = ksdensity(MIS,round(0:max(MIS)/10000:max(MIS)*2));
a=a./sum(a);
% Calculate mutual information if not an input
if isempty(varargin{1})
    I = rapidmi(xph,yph,bins);
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
function Ir = rapidmi(m,n,bins)

if iscolumn(m) == 0
    m = m';
end
if iscolumn(n) == 0
    n = n';
end

[bwmn,dmn,meshmnm,meshmnn]=kde2d([m n],bins);
dmn=abs(dmn);dmn = dmn./sum(sum(dmn));dmnlog = log2(dmn);
dmnlog(dmnlog==inf | dmnlog==-inf)=0;

Inh1 = log2(sum(dmn,1));Inh2 = log2(sum(dmn,2));
Inh1(Inh1==inf | Inh1==-inf)=0;Inh2(Inh2==inf | Inh2==-inf)=0;

Ir = sum(sum(dmn.*bsxfun(@minus,bsxfun(@minus,dmnlog,Inh1),Inh2)));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
