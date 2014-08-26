function H = jointentropy1(x,y,bins)
% This function calculates the joint entropy for two signals.
%
%  H = jointentropy1(x,y,bins)
%
%  H : joint entropy
%  x,y : input signals
%  bins : bin size for histograms (two-element array)
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

% Calculate the 2D histogram
[ nh , co ] = hist3([x y],bins);
hx = co{1,1};
hy = co{1,2};

% Normalize the histogram to make a pdf
nh = nh./sum(sum(nh));
dxo = hx(2) - hx(1);
dyo = hy(2) - hy(1);

% Calculate the joint entropy
ind = nh ~= 0;
H = -dxo*dyo*sum(nh(ind).*log2(nh(ind)));

end
