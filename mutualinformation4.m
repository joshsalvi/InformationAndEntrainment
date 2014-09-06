function [I p] = mutualinformation4(x,y,bins)
% This function calculates the mutual information for two signals using
% a kernel density estimator.
%
%  I = mutualinformation2(x,y)
%
%  I : mutual information (bits)
%  x,y : input signals
%  bins : how many bins would you like to use? (default = 2^4)
%  p : p-value for the mutual information value calculated here
%  
%
% I recommend using the Freedman-Diaconis rule for calculating the
% appropriate bin size.
%
% Joshua D Salvi, jsalvi@rockefeller.edu
%

if iscolumn(x) == 0
    x = x';
end
if iscolumn(y) == 0
    y = y';
end

if isempty(bins) == 1
    bins = 2^4;
end

% Kernel density estimate
%{
[bwx,dx,meshx,cdfx]=kde1d(x,2^10);
[bwy,dy,meshy,cdfy]=kde1d(y,2^10);
dx = abs(dx);
dx = dx./sum(sum(dx));
dxlog = log2(dx);
dxlog(dxlog==inf | dxlog==-inf)=0;
dy = abs(dy);
dy = dy./sum(sum(dy));
dylog = log2(dy);
dylog(dylog==inf | dylog==-inf)=0;
%}

% 2D kernel density estimate
%[xy]=gkde2([x y]);
%dxy = xy.pdf;
[bwxy,dxy,meshxyx,meshxyy]=kde2d([x y],bins);
dxy=abs(dxy);
dxy = dxy./sum(sum(dxy));
dxylog = log2(dxy);
dxylog(dxylog==inf | dxylog==-inf)=0;

% Marginal probabilities
Inh1 = log2(sum(dxy,1));
Inh2 = log2(sum(dxy,2));
Inh1(Inh1==inf | Inh1==-inf)=0;
Inh2(Inh2==inf | Inh2==-inf)=0;

% Mutual information from the joint and marginal probabilities
I = sum(sum(dxy.*bsxfun(@minus,bsxfun(@minus,dxylog,Inh1),Inh2)));

% Calculate the p value
% Mutual information is equal to the G-test statistic divided by 2N where
% N is the sample size. The G-test is also roughly equal to a chi-squared
% distribution. 
% http://en.wikipedia.org/wiki/Mutual_information
% http://www.biostathandbook.com/chigof.html#chivsg
df = (length(meshxyx)-1)*(length(meshxyy)-1);
Gstat = I*2*length(meshxyx)*length(meshxyy);
p = gammainc(Gstat/2,df/2,'upper');     % p-value


end
