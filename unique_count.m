
function [y c] = unique_count(x)
% Unique_count sorts and counts the rows of x
% y lists the unique rows and c their counts
if size(x,1) == 1;
    x = x(:);
end
y = sortrows(x);
d = any(y(1:end-1,:)-y(2:end,:),2);
d = [d; true];
ind = (1:size(y,1));
g=ind(d);
n = length(g);
c(1) = g(1);
c(2:n) = g(2:n) - g(1:n-1);
c = c';
y = y(g,:);

end
