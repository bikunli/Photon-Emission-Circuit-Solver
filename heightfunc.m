function [h] = heightfunc(Tableau,option)
% input an arbitary tableau
% output a vector h, representing h(x)
% h(1) = h(x=1)
% This function uses the method introduced  in paper: New J. Phys. 7, 170 (2005)
% Version 4.0, Date: 06/09/2021
[n,n2] = size(Tableau);
if n2 ~= 2*n
    error('The size of Tableau_photon is improper!');
end


h = zeros(1,n);

if nargin == 1 % general case
    rv_ind = n:-1:1;
    Tableau_reversed = Tableau(:,[rv_ind,n+rv_ind]);
    [T_echelon,~] = ECHELON_TABLEAU(Tableau_reversed);
    B = Tableau2Bigram(T_echelon);
    
    for i = 1:n
        h(i) = i - numel( find(B(:,1)>= (n-i+1)) );
    end
elseif nargin == 2
    switch option
        case 'pure' % if pure state
            [T_echelon,~] = ECHELON_TABLEAU(Tableau);
            B = Tableau2Bigram(T_echelon);
            for i = 1:n
                h(i) = n - i - numel( find(B(:,1)> i) );
            end
        otherwise
            h = heightfunc(Tableau);
    end
end


% % % % % input an arbitary tableau
% % % % % output a vector h, representing h(x)
% % % % % h(1) = h(x=1)
% % % % % Version 2.0, Date: 01/27/2021
% % % % [n,n2] = size(Tableau);
% % % % if n2 ~= 2*n
% % % %     error('The size of Tableau_photon is improper!');
% % % % end
% % % % h = zeros(1,n);
% % % % [T_clipped,~] = CLIPPING_TABLEAU(Tableau);
% % % % B = Tableau2Bigram(T_clipped);
% % % % for i = 1:n-1
% % % %     h(i) = 1/2*numel(intersect(find(B(:,1)<=i),find(B(:,2)>i)));
% % % % end
% % % % end