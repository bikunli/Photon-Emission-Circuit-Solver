function [Tableau_echelon,R] = ECHELON_TABLEAU(Tableau)
% The echelon form Algorithm introduced in paper: New J. Phys. 7, 170 (2005)
% Input:
% Tableau -- The Tabluea matrix (Lx2L) to be par-trace
% Output:
% Tableau_echelon -- The Tableau under echelon gauge
% R: a binary matrix that transforms rows of tableau: Tabluea_echelon = R * Tableau
% -- Version: v2.1  Date: 04/11/2021 --

[L,Lcol] = size(Tableau);
R = logical(eye(L));
if 2*L ~= Lcol
    warning('The input Tableau does not have proper size!');
end

% pick all non-zero rows, and put them to the top of the tableau
nz_row_resort_ind = [find(~all(Tableau == 0,2));find(all(Tableau == 0,2))];
Tableau = Tableau(nz_row_resort_ind,:);
R = R(nz_row_resort_ind,:);

Tableau_alt = Tableau(:,reshape([1:L;(1:L)+L],2*L,1)); % switch to the form of alternating indices
nz_row_num = numel(find(~all(Tableau==0,2))); % non-zero rows number


if nz_row_num > 1
    for i = 1:(nz_row_num-1)
        % i is the current working row/column
        % find the position of pivot elements:
        [rpiv, cpiv] = find(Tableau_alt(i:nz_row_num,:) == 1,1);
        rpiv = rpiv + i - 1; % if rpiv = 1, then output rpiv = i;
        if isempty(cpiv) == 1
            warning('Encounter empty case! \n');% which shouldn't happen
            break;   % if cpiv column is empty, doing nothing
        elseif isempty(cpiv) == 0
            % Here goes if cpiv column is NOT empty:
            if rpiv ~= i % swap rows: moves the pivot row to i-th row
                Tableau_alt([i,rpiv],:) = Tableau_alt([rpiv,i],:);
                R([i,rpiv],:) = R([rpiv,i],:);
            end
            % index vector of rows to-be-eliminated
            ind_2belim = find(Tableau_alt((i+1):nz_row_num,cpiv)==1) + i;
            if isempty(ind_2belim) ~= 1 % doing rows elimination
                Tableau_alt(ind_2belim,:) = mod(Tableau_alt(ind_2belim,:) + Tableau_alt(i,:),2);
                R(ind_2belim,:) = mod(R(ind_2belim,:) + R(i,:),2);
            end
        end
    end
elseif nz_row_num == 0
    warning('The tableau is empty!');
elseif nz_row_num == 1
    warning('The tableau has only 1 non-trivial row!');
end

% Return to conventional Tableau colomn order
Tableau_echelon = Tableau_alt(:,[2*(1:L)-1,2*(1:L)]);
R = logical(R);

end