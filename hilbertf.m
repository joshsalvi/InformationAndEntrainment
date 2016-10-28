function x = hilbertf(xr,Fs,f0,ninterp)
% Discrete-time analytic signal via Hilbert transform
%
%   X = hilbertf(Xr) computes the so-called discrete-time analytic signal
%   X = Xr + i*Xi such that Xi is the Hilbert transform of real vector Xr.
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%   If Xr is a matrix, then HILBERT operates along the columns of Xr.
%
%   hilbertf(Xr,Fs,F0) computes the Hilbert transform for the frequency range 
%   of interest F0 = [F1 f2], given a sample rate Fs.
%   
%
%   For a discrete-time analytic signal X, the last half of fft(X) is zero, 
%   and the first (DC) and center (Nyquist) elements of fft(X) are purely real.
%
%   EXAMPLE:
%          Xr = [1 2 3 4];
%          X = hilbert(Xr)
%          % produces X=[1+1i 2-1i 3-1i 4+1i] such that Xi=imag(X)=[1 -1 -1 1] is the
%          % Hilbert transform of Xr, and Xr=real(X)=[1 2 3 4].  Note that the last half
%          % of fft(X)=[10 -4+4i -2 0] is zero (in this example, the last half is just
%          % the last element).  Also note that the DC and Nyquist elements of fft(X)
%          % (10 and -2) are purely real.
%
%   See also FFT, IFFT.

%   Copyright 1988-2008 The MathWorks, Inc.

%   References:
%     [1] Alan V. Oppenheim and Ronald W. Schafer, Discrete-Time
%     Signal Processing, 2nd ed., Prentice-Hall, Upper Saddle River, 
%     New Jersey, 1998.
%
%     [2] S. Lawrence Marple, Jr., Computing the discrete-time analytic 
%     signal via FFT, IEEE Transactions on Signal Processing, Vol. 47, 
%     No. 9, September 1999, pp.2600--2603.
if exist('ninterp')
    t1=1:length(xr);t1q=linspace(t1(1),t1(end),ninterp*length(t1));
    xr = interp1(t1,xr,t1q,'spline');
end
if nargin<2, Fs=[]; f0=[]; end
if nargin==2
    disp('Need a frequency of interest as third argument.');
    f0=[];
end

if ~isreal(xr)
  warning(message('signal:hilbert:Ignore'))
  xr = real(xr);
end
% Work along the first nonsingleton dimension
[xr,nshifts] = shiftdim(xr);

n=(size(xr,1));
x = fft(xr,n,1); % n-point FFT over columns.

if nargin >= 3 && isempty(Fs)==0
    f=Fs/2*linspace(0,1,length(x)/2+1);f=[f f];
    if length(f0) == 2
        f01=findnearest(f,f0(1));f01(f01<1)=1;
        f02=findnearest(f,f0(2));f02(f02>n)=n;
    elseif length(f0) == 1
        f0 = [f0 f0];
        f01=findnearest(f,f0(1));
        f02=findnearest(f,f0(2));
    else
        disp('Need a frequency range. Exiting.')
        return;
    end
    if isempty(f01)==0 && isempty(f02)==0 && length(f01)==2 && length(f02)==2
        f01=f01([1 end]);f02=f02([1 end]);
        x=[x(f01(1):f02(1)); x(f01(2):f02(2))];
    else
        disp('Could not find frequency of interest.');
    end
end
n=length(x);
h  = zeros(n,~isempty(x)); % nx1 for nonempty. 0x0 for empty.


if n > 0 && 2*fix(n/2) == n
  % even and nonempty
  h([1 n/2+1]) = 1;
  h(2:n/2) = 2;
elseif n>0
  % odd and nonempty
  h(1) = 1;
  h(2:(n+1)/2) = 2;
end
x = ifft(x.*h(:,ones(1,size(x,2))));

% Convert back to the original shape.
x = shiftdim(x,-nshifts);

if exist('ninterp')
    x = x(1:ninterp:end);
end
end


function [r,c,V] = findnearest(srchvalue,srcharray,bias)

if nargin<2
    error('Need two inputs: Search value and search array')
elseif nargin<3
    bias = 0;
end
% find the differences
srcharray = srcharray-srchvalue;
if bias == -1   % only choose values <= to the search value    
    srcharray(srcharray>0) =inf;        
elseif bias == 1  % only choose values >= to the search value    
    srcharray(srcharray<0) =inf;        
end

% give the correct output
if nargout==1 | nargout==0
    
    if all(isinf(srcharray(:)))
        r = [];
    else
        r = find(abs(srcharray)==min(abs(srcharray(:))));
    end 
        
elseif nargout>1
    if all(isinf(srcharray(:)))
        r = [];c=[];
    else
        [r,c] = find(abs(srcharray)==min(abs(srcharray(:))));
    end
    
    if nargout==3
        V = srcharray(r,c)+srchvalue;
    end
end
end
