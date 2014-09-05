
function z = embed2(x,d)
% rearranges a one dimensional array x into
% overlapping rows of length d
x = x(:);
N = length(x);
i = (1:d);
j = (0:N-d)';
a = i(ones(N-d+1,1),:);
b = j(:,ones(1,d));
z = x(a+b);
end
