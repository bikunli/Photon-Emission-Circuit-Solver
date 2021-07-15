function [Bigram] = Tableau2Bigram(Tableau)
% Any tableau has a Bigram, which is the output of this function
% Bigram: there are n rows, each row is [l(g),r(g)]
n = numel(find(~all(Tableau == 0,2))); % non trivial row of tableau
Bigram = zeros(n,2); %
for i_r = 1:n
    Bigram(i_r,1) = find(Tableau(i_r,1:n)+Tableau(i_r,n+(1:n))~=0, 1, 'first');
    Bigram(i_r,2) = find(Tableau(i_r,1:n)+Tableau(i_r,n+(1:n))~=0, 1, 'last');
end
end