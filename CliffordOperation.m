function [Generators] = CliffordOperation(Generators,gate_str,site)
% This function operates Generators by specific clifford gate
% Generators is a structure variable:
% Generator.Tableau is the tableau representation
% Generator.SignVector is a binary vector records the sign: (-1)^SignVector
% gate_str can be 'I','H','P','PH','HP','HPH','X','Y','Z','CNOT'.......
% ** but 'SWAP' is not included (too trivial)
% (more gate_str see the dictionary: OtherLC)
% Other operation should be realized repeatedly using this function
% site is a scalar, if the operation acts on single qubit gate
% site is a vector [i,j], if the operation is a control gate CNOT_{i,j}
% This function needs to call function: Qmatrix() (included at the end)
% Version v1.0,  Date: 01/17/2021

n = size(Generators.Tableau,1);
OtherLC = {'HP','PH','HPH', ...
    'HX','PX','HPX','PHX','HPHX', ...
    'HY','PY','HPY','PHY','HPHY', ...
    'HZ','PZ','HPZ','PHZ','HPHZ'};
switch gate_str
    case {'I',''}
    case 'H'
        Generators.Tableau = mod(Generators.Tableau * Qmatrix(n, gate_str ,site),2);
        Generators.SignVector = mod(Generators.SignVector + Generators.Tableau(:,site).*Generators.Tableau(:,n+site),2);
    case 'P'
        Generators.Tableau = mod(Generators.Tableau * Qmatrix(n, gate_str ,site),2);
        Generators.SignVector = mod(Generators.SignVector + Generators.Tableau(:,site).*Generators.Tableau(:,n+site),2);
    case 'X'
        Generators.SignVector = mod(Generators.SignVector + Generators.Tableau(:,n+site),2);
    case 'Y'
        Generators.SignVector = mod(Generators.SignVector + Generators.Tableau(:,site) + Generators.Tableau(:,n+site),2);
    case 'Z'
        Generators.SignVector = mod(Generators.SignVector + Generators.Tableau(:,site),2);
    case OtherLC
        for i_str = 1:length(gate_str)
            Generators = CliffordOperation(Generators,gate_str(i_str),site);
        end
    case {'CNOT','CX'}
        Generators.Tableau = mod(Generators.Tableau * Qmatrix(n,'CNOT',site(1), site(2)),2);
        Generators.SignVector = mod(Generators.SignVector ...
            + Generators.Tableau(:,site(1)).*Generators.Tableau(:,n+site(2)) ...
            .*(Generators.Tableau(:,site(2)) + Generators.Tableau(:,n+site(1)) + 1),2);
    otherwise
        fprintf(1,['gate_str = %s .\n'],gate_str);
        error('Unexpected string of variable ''gate_str''  ');
end
end

function [Q]=Qmatrix(N,type,i,j)
% Q matrix is a matrix that right-multiply on the tableau matrix
% N is the qubit number
% type is the operation type
% i,j are sites to be acted on.
if strcmp(type,'CZ') == 1 || strcmp(type,'CNOT') == 1 || strcmp(type,'CX') == 1
    if nargin < 4
        error('Insufficient input arguments!');
    end
end
Q = eye(2*N);

switch type
    case 'I'
        Q([i,N+i],[i,N+i]) = eye(2);
    case 'H'
        Q([i,N+i],[i,N+i]) = [0,1;1,0];
    case 'P'
        Q([i,N+i],[i,N+i]) = [1,1;0,1];
    case 'HP'
        Q([i,N+i],[i,N+i]) = [0,1;1,1];
    case 'PH'
        Q([i,N+i],[i,N+i]) = [1,1;1,0];
    case 'HPH'
        Q([i,N+i],[i,N+i]) = [1,0;1,1];
    case {'CNOT','CX'}
        Q(i,j) = 1;
        Q(N+j,N+i) = 1;
    case 'CZ'
        Q(j,[N+i]) = 1;
        Q(i,[N+j]) = 1;
    otherwise
        error('Unexpected string of variable ''type''! ');
end

end
