function [VS, rayleigh_p] = vscalc2(x,y,rayleighyn)
% This function calculates the vector strength for two signals using
% each function's analytic signal.
%
%  [VS, rayleigh_p, rayleigh_stat] = vscalc(x,y,1)
%
%  VS : vector strength
%  rayleigh_stat : rayleigh statistic
%  x,y : input signals
%  rayleighyn : run rayleigh test? (1=yes)
%  
% The function does not call angle() so that it may avoid wrapping
% artifacts when calculating the instantaneous phase.
%
% Joshua D Salvi, jsalvi@rockefeller.edu
%

% Calculate the analytic signal using the Hilbert transform.
xhilb = hilbert(x);
xhilb_eiphi = xhilb./abs(xhilb);    % normalize all lengths to 1
yhilb = hilbert(y);
yhilb_eiphi = yhilb./abs(yhilb);    % normalize all lengths to 1
   
% Calculate vector strength.
VS = abs(sum((xhilb_eiphi./yhilb_eiphi))/length(xhilb));

if rayleighyn == 1
    N = length(xhilb)*10/4;
    VS2_n = VS*N;
    rayleigh_p = exp(sqrt(1+4*N+4*(N^2-VS2_n^2))-(1+2*N));
end

end
