function [VS, rayleigh_z, rayleigh_p, rayleigh_h] = vscalcMN(varargin)
% This function calculates the vector strength for two signals using
% each function's analytic signal for different m:n ratios of phase locking.
%
%  [VS, rayleigh_z, rayleigh_p, rayleigh_h] = vscalcMN(x,y,rayleighyn,alphalevel,m,n)
%
%  VS : vector strength
%  rayleigh_stat : rayleigh statistic
%  rayleigh_h : accept/reject null (1=accept,0=reject)
%  x,y : input signals
%  rayleighyn : run rayleigh test? (1=yes)
%  alphalevel : alphalevel for rayleigh test (default=0.01)
%  m,n : calculate m:n phase locking (default=1:1)
%  
% The function does not call angle() so that it may avoid wrapping
% artifacts when calculating the instantaneous phase.
%
% Joshua D Salvi, jsalvi@rockefeller.edu
%

narginchk(2, 6)

x=varargin{1};
y=varargin{2};

if nargin < 3
    rayleighyn=0;
    rayleigh_z='NaN';
    rayleigh_p='NaN';
    rayleigh_h='NaN';
else
    if isempty(varargin{3})==0
        rayleighyn=varargin{3};
    else
        rayleighyn=0;
    end
end

if nargin < 4
    alphalevel=0.01;
else
    if isempty(varargin{4})==0
        alphalevel=varargin{4};
    else
        alphalevel=0.01;
    end
end

if nargin < 5
    m=1;n=1;
else
    if isempty(varargin{5})==0
        m=varargin{5};
    else
        m=1;
    end
end

if nargin < 6
    n=1;
else
    if isempty(varargin{6})==0
        n=varargin{6};
    else
        n=1;
    end
end

if iscolumn(x)==0
    x=x';
end
if iscolumn(y)==0
    y=y';
end
if isempty(alphalevel)==1
    alphalevel=0.01;
end

if alphalevel > 0.2
    disp('Alpha level too large. Restrict to 0.2.');
end

% Calculate the analytic signal using the Hilbert transform.
xhilb = hilbert(x);
xp = atan2(imag(xhilb),real(xhilb));
yhilb = hilbert(y);
yp = atan2(imag(yhilb),real(yhilb));
N = length(xhilb);
   
% Calculate vector strength.
clear i
VS = abs(1/N*sum(exp(i*(m*xp-n*yp))));

if rayleighyn == 1
    N = length(xhilb)*10/4;
    VS2_n = VS*N;
    rayleigh_z = N*VS^2;           % http://ethologie.unige.ch/pdf-libres/test.de.rayleigh.et.v-test.pdf
    rayleigh_h = exp(sqrt(1+4*N+4*(N^2-VS2_n^2))-(1+2*N));
    rayleigh_p = exp(-rayleigh_z);      % Grun and Rotter, Analysis of Parallel Spike Trains; AND Fisher (1993) Statistical Analysis of Circular Data

    if alphalevel == 0.1
        if rayleigh_z >= 2.30
            rayleigh_h = 0;
        else
            rayleigh_h = 1;
        end
    elseif alphalevel == 0.05
        if rayleigh_z >= 3
            rayleigh_h = 0;
        else
            rayleigh_h = 1;
        end
    elseif alphalevel == 0.02
        if rayleigh_z >= 3.91
            rayleigh_h = 0;
        else
            rayleigh_h = 1;
        end
    elseif alphalevel == 0.01
        if rayleigh_z >= 4.61
            rayleigh_h = 0;
        else
            rayleigh_h = 1;
        end
    elseif alphalevel == 0.001
        if rayleigh_z >= 6.91
            rayleigh_h = 0;
        else
            rayleigh_h = 1;
        end
    else
        if rayleigh_p < alphalevel
            rayleigh_h = 0;
        else
            rayleigh_h = 1;
        end
    end
    %}

end

end
