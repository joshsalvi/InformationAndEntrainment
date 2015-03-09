function [Ij, Ir, Ii, Dkl] = rapidentropyRI(x,bins)
%
% [Ij, Ir, Ii, Dkl] = rapidentropyRI(x,bins)
%
% Ij = joint entropy for real/imaginary 2D phase portrait
% Ir = entropy for 1D histogram (real)
% Ii = entropy for 1D histogram (imag)
% Dkl = KL divergence between real and imaginary distributions
% x = data (input), 1D vector
% bins = # of bins, set to 0 for freedman-diaconis rule 
%
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%

if bins == 0
    bins = freedmandiaconis(x);
end

xr = real(hilbert(x));
xi = imag(hilbert(x));
xp = angle(hilbert(x));
if iscolumn(xr) == 0
    xr = xr';
end

if iscolumn(xi) == 0
    xi = xi';
end
    
if iscolumn(xp) == 0
    xp = xp';
end

[~,dmn,~,~]=kde2d([xr xi],bins);
[~,dmp]=kde1d(xp,2^8);dmp=abs(dmp);dmp=dmp./sum(dmp);
dmn(dmn<0)=0;
dmn = dmn./sum(sum(dmn));dmnlog = log2(dmn);
dmnlog(dmnlog==inf | dmnlog==-inf)=0;

dmnr = sum(dmn,1);dmni = sum(dmn,2);
Inh1 = log2(dmnr);Inh2 = log2(dmni);
Inh3 = log2(dmnr./dmni');
Inh1(Inh1==inf | Inh1==-inf)=0;Inh2(Inh2==inf | Inh2==-inf)=0;
Inh3(Inh3==inf | Inh3==-inf)=0;

Ir = -sum(bsxfun(@times,dmnr,Inh1));
Ii = -sum(bsxfun(@times,dmni,Inh2));
Ij = -sum(sum(bsxfun(@times,dmn,dmnlog)));

% KL divergence between real and imaginary distributions (Dkl(R||I))
Dkl = sum(bsxfun(@times,dmnr,Inh3));

end
