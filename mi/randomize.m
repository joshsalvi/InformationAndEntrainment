function W = randomize(A)
[m,n] = size(A);
E = A(:);
W(1) = E(1);
E(1) =[];
N = m*n;
while length(E) > 0
    K = length(W);
    RandInd = randi(length(E),1);
    for j = 1: K 
       P(j) =  E(RandInd) ~= W(j); 
    end
    if all(P) 
    W =[W,E(RandInd)];
    E(RandInd) =[];
    end
end  
W = reshape(W,m,n);
