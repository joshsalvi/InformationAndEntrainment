function I = mutualinformation3(x,y,plotyn)
% This function calculates the mutual information for two signals using
% a 2D histogram of phase.
%
%  I = mutualinformation2(x,y)
%
%  I : mutual information (bits)
%  x,y : input signals
%  plotyn : plot the probability distributions? (1 = yes)
%  
%
% I recommend using the Freedman-Diaconis rule for calculating the
% appropriate bin size, as implemented in this script.
%
% Joshua D Salvi, jsalvi@rockefeller.edu
%

if iscolumn(x) == 0
    x = x';
end
if iscolumn(y) == 0
    y = y';
end

xh = hilbert(x);
xph = angle(xh);                      % normalize all lengths to 1
yh = hilbert(y);
yph = angle(yh);                      % normalize all lengths to 1

% Freedman-Diaconis rule to calculate bin size
bwx = 2*iqr(xph)/length(xph)^(1/3);
nbx = ceil((max(xph) - min(xph))/bwx);
bwy = 2*iqr(yph)/length(yph)^(1/3);
nby = ceil((max(yph) - min(yph))/bwy);

% 2D histogram
[nh,co] = hist3([xph yph],[nbx nby]);
nh = nh./sum(sum(nh));
nhlog = log2(nh);
nhlog(nhlog==inf | nhlog==-inf)=0;

% Marginal probabilities
Inh1 = log2(sum(nh,1));
Inh2 = log2(sum(nh,2));
Inh1(Inh1==inf | Inh1==-inf)=0;
Inh2(Inh2==inf | Inh2==-inf)=0;

% Mutual information from the joint and marginal probabilities
I = sum(sum(nh.*bsxfun(@minus,bsxfun(@minus,nhlog,Inh1),Inh2)));

if plotyn == 1
    figure;
    imagesc(co{2},co{1},nh); title('Joint Probability');
    xlabel('Y');ylabel('X');colorbar
    
    figure;
    subplot(1,2,1);plot(co{2},sum(nh,1)); title('Marginal Probability (Y)');
    subplot(1,2,2);plot(co{1},sum(nh,2)); title('Marginal Probability (X)');
end

end
