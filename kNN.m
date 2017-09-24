function [k, acc ] = kNN(sig, type, max_k)
%KNN Summary of this function goes here
%   Detailed explanation goes here

pd = squareform(pdist(sig));

k= 1:max_k;
acc=zeros(1,max_k);

for i=1:max_k
    z = zeros(length(type),1);
    for j = 1:length(type)
        [~, nn] = sort(pd(:,j),1,'ascend');
        z(j) = mode(type(nn(2:(i+1),:)),1);
    end
acc(i) = sum(z == type) / length(type);
end

end

