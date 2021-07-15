function [Generators] = GaugeTransformation(Generators,R)
% This function do gauge transformation (row transformation).
% Generators is a structure variable, which is updated in the output:
% Generator.Tableau is the tableau representation of generators
% Generator.SignVector is a binary vector records the sign: (-1)^SignVector
% R is a binary n x n matrix that complete the row transformation
% Version: v1.0, Date: 01/17/2021


SignVector_temp = mod( R*Generators.SignVector,2 );
n = size(Generators.Tableau,1);
for k = 1:n
    g_mat = Generators.Tableau(logical(R(k,:)),:);
    ExtraSign = 0; % phase factor: (-1)^ExtraSign
    m = size(g_mat,1);
    % The follow steps mutiplies the i-th row with the reference row,
    % It will find if a minus sign should be obtained.
    rr = g_mat(1,:); % Reference row: the first row
    for i = 2:m
        r4ins = g_mat(i,:); % row for inspection: will be multipied with the rr row.
        for j = 1:n
            % update the phase factor: (-1)^ph
            ExtraSign = ExtraSign + ((rr(j)==1)*(rr(j+n)==1)*(r4ins(j+n) - r4ins(j)) ...
                + (rr(j)==1)*(rr(j+n)==0)*(r4ins(j+n)*(2*r4ins(j)-1)) ...
                + (rr(j)==0)*(rr(j+n)==1)*(r4ins(j)*(1-2*r4ins(j+n))))/2;
        end
        rr = mod(rr + r4ins,2);
    end
    if ExtraSign ~= round(ExtraSign) % ph must be an interger !
        warning('ExtraSign is not an interger!');    % which shouldn't happen if g_mat is legal
    else
        ExtraSign = round(ExtraSign);
    end
    Generators.SignVector(k) = mod(SignVector_temp(k) + ExtraSign,2);
end
Generators.Tableau = mod( R*Generators.Tableau,2 );
end