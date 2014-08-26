function H = entropy1(y,bins)
% This function calculates the marginal entropy of a signal.
%
%  H = entropy1(y,bins)
%
%  H : marginal entropy
%  y : input signal
%  bins : bin size for histograms
%
% I recommend using the Freedman-Diaconis rule for calculating the
% appropriate bin size.
%
% Joshua D Salvi, jsalvi@rockefeller.edu
%


% Calculate a histogram
[ nh xo ] = hist(y,bins);

% Normalize the histogram to find a pdf
nh = nh / sum(nh);
dxo = xo(2) - xo(1);

% Calculate the entropy
ind = nh ~= 0;
H = -sum(nh(ind).*log2(nh(ind)).*dxo);

end
