function [varargout] = analyticinformation(varargin)
%
% Analytic Information
% --------------------
%
% (1) AI = analyticinformation(x)
%
% Computes the analytic information (AI) for vector, x. The analytic
% information is a normalized value using the mutualinfo() function.
%
%
% (2) [AI, VS, rayleigh_z, rayleigh_pval] = analyticinformation(x,y)
%
% Computes the analytic information for x and the analytic information for
% y, in which x and y are two vectors of the same length.
% Outputs the analytic information for each (AI = [AI_x , AI_y]), the
% vector strength (synchronization index) of x and y, and the Rayleigh
% statistic (rayleigh_z) and associated p-value (rayleigh_pval) for the
% vector strength.
%
% Coder: Joshua D. Salvi
% Email: jsalvi@rockefeller.edu
% Year: 2016
%

% Initialization
for j = 1:nargout
    varargout{j} = NaN;
end

if nargin < 1
    disp('Requires at least one input argument. Exiting.');
    return    
elseif nargin == 1
    x = varargin{1};
    % Calculate AI for x
    varargout{1} = mutualinfo(x,hilbert(x));
elseif nargin == 2
    x = varargin{1};
    y = varargin{2};
    if length(x) == length(y)
        N = length(x);
        % Calculate AI for each variable (x and y)
        varargout{1} = [mutualinfo(x,hilbert(x)), mutualinfo(y,hilbert(y))];
        % Calculate VS for x and y
        varargout{2} = VScalc(x,y);
        % Calculate Rayleigh R and z values
        R = N * varargout{2};
        varargout{3} = R^2 / N;
        % Estimate the p-value from Rayleigh z
        varargout{4} = exp(-varargout{3});
    else
        disp('X and Y must be of the same length. Exiting.');
        return
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% VScalc(x,y) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VS] = VScalc(x,y)

if iscolumn(x)==0
    x=x';
end
if iscolumn(y)==0
    y=y';
end
xhilb = hilbert(x);
xhilb_eiphi = xhilb./abs(xhilb);    
yhilb = hilbert(y);
yhilb_eiphi = yhilb./abs(yhilb);
VS = abs(sum((xhilb_eiphi./yhilb_eiphi))/length(xhilb));

end

    
