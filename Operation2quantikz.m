function [] = Operation2quantikz(Operation,n_p,n_e,filename)
% This function will print the quantikz code for quantum circuit, 
% based on the recipe Operation
% n_p is the photons number
% n_e is the emitter qubits number
% The final output will be a .txt file: filename.txt, in root directory

n = n_e + n_p;

CC = cell(n,1); % The cell variable that prepares everything to be printed
IniState0_p = '\\gate[style={fill=red!20}]{\\ket{0}}';
IniState0_e = '\\gate[style={fill=green!20}]{\\ket{0}}'; % defines the style of 'initialization box'


c = 1; % initialize the column
for i = 1:n_p
    %     CC{i,c} = ['\\lstick{$p_{',num2str(i),'}:\\ket{0}$}'];
    CC{i,c} = ['\\lstick{$p_{',num2str(i),'}$}'];
end
for i = (1:n_e)+n_p
    %     CC{i,c} = ['\\lstick{$e_{',num2str(i),'}:\\ket{0}$}'];
    CC{i,c} = ['\\lstick{$e_{',num2str(i-n_p),'}$}'];
end

c = c + 1;
for i = 1:n_p
    CC{i,c} = IniState0_p;
end
for i = (1:n_e) + n_p
    CC{i,c} = IniState0_e;
end


% prepare W0
Op_type = Operation.emitters.W0.type;
Op_site = Operation.emitters.W0.site;
c = c + 1;
for k = size(Op_site,2):-1:1
    if isempty(Op_type{1,k}) == 0
        switch Op_type{1,k}
            case {'I',''}
            case {'CNOT','CX'}
                if size(CC,2) >= c
                    if all(cellfun('isempty',CC(:,c)) == 1) == 0
                        c = c + 1;
                    end
                end
                for i = (1:n_e)+n_p
                    if i == Op_site{1,k}(1)
                        d = Op_site{1,k}(2) - Op_site{1,k}(1);
                        CC{i,c} = ['\\ctrl{',num2str(d),'}'];
                    elseif i == Op_site{1,k}(2)
                        CC{i,c} = '\\targ{}';
                    else
                        CC{i,c} = '\\qw';
                    end
                end
                c = c + 1;
            otherwise  % nontrivial  single qubit gate ( local clifford ) 
                if size(CC,2) >= c
                    if isempty( CC{Op_site{1,k},c} ) == 0
                        c = c + 1;
                    end
                end
                for s = 1:length(Op_type{1,k})
                    CC{Op_site{1,k},c-1 + s} = ['\\gate{',Op_type{1,k}(length(Op_type{1,k})+1-s),'}'];
                end
                c = c - 1 + length(Op_type{1,k});
        end
    end
end


for j = 1:n_p
    % Emission CNOT
    if size(CC,2) >= c
        if all(cellfun('isempty',CC(:,c)) == 1) == 0
            c = c + 1;
        end
    end
    for i = 1:n
        if i == Operation.emitters.EmissionSite{j,1}
            d = j - i;
            CC{i,c} = ['\\ctrl{',num2str(d),'}'];
        elseif i == j
            CC{i,c} = '\\targ{}';
        else
            CC{i,c} = '\\qw';
        end
    end
    c = c + 1;
    
    % Ue_j & Up_j
    Op_type = [Operation.emitters.Ue.type(j,:),Operation.photons.Up.type(j,1)];
    Op_site = [Operation.emitters.Ue.site(j,:),j];
    for k = size(Op_site,2):-1:1
        if isempty(Op_type{1,k}) == 0
            switch Op_type{1,k}
                case {'I',''} % trivial
                case 'CNOT'
                    if size(CC,2) >= c
                        if all(cellfun('isempty',CC(:,c)) == 1) == 0
                            c = c + 1;
                        end
                    end
                    for i = (1:n_e)+n_p
                        if i == Op_site{1,k}(1)
                            d = Op_site{1,k}(2) - Op_site{1,k}(1);
                            CC{i,c} = ['\\ctrl{',num2str(d),'}'];
                        elseif i == Op_site{1,k}(2)
                            CC{i,c} = '\\targ{}';
                        else
                            CC{i,c} = '\\qw';
                        end
                    end
                    c = c + 1;
                otherwise
                    if size(CC,2) >= c
                        if isempty( CC{Op_site{1,k},c} ) == 0
                            c = c + 1;
                        end
                    end
                    for s = 1:length(Op_type{1,k})
                        CC{Op_site{1,k},c-1 + s} = ['\\gate{',Op_type{1,k}(length(Op_type{1,k})+1-s),'}'];
                    end
                    c = c - 1 + length(Op_type{1,k});
            end
        end
    end
    
    
    
    % measurement
    if isempty(Operation.emitters.MeasurementSite{j,1}) == 0
        if size(CC,2) >= c
            if all(cellfun('isempty',CC(:,c)) == 1) == 0
                c = c + 1;
            end
        end
        for i = 1:n
            if i == Operation.emitters.MeasurementSite{j,1}
                d = j - i;
                CC{i,c} = ['\\meter{} \\vcw{',num2str(d),'}'];
            elseif i == j
                CC{i,c} = '\\gate{X}';
            else
                CC{i,c} = '\\qw';
            end
        end
        c = c + 1;
        
        %-----re-initialize emitter qubit-----
        for i = 1:n
            if i == Operation.emitters.MeasurementSite{j,1}
                CC{i,c} = IniState0_e;
            else
                CC{i,c} = '\\qw';
            end
        end
        c = c + 1;
        for i = 1:n
            if i == Operation.emitters.MeasurementSite{j,1}
                CC{i,c} = ['\\gate{H}'];
            else
                CC{i,c} = '\\qw';
            end
        end
        c = c + 1;
        %----------
        
        
        
        % W
        Op_type = Operation.emitters.W.type(j,:);
        Op_site = Operation.emitters.W.site(j,:);
        for k = size(Op_site,2):-1:1
            if isempty(Op_type{1,k}) == 0
                switch Op_type{1,k}
                    case {'I',''}
                    case 'CNOT'
                        if size(CC,2) >= c
                            if all(cellfun('isempty',CC(:,c)) == 1) == 0
                                c = c + 1;
                            end
                        end
                        for i = (1:n_e)+n_p
                            if i == Op_site{1,k}(1)
                                d = Op_site{1,k}(2) - Op_site{1,k}(1);
                                CC{i,c} = ['\\ctrl{',num2str(d),'}'];
                            elseif i == Op_site{1,k}(2)
                                CC{i,c} = '\\targ{}';
                            else
                                CC{i,c} = '\\qw';
                            end
                        end
                        c = c + 1;
                    otherwise
                        if size(CC,2) >= c
                            if isempty( CC{Op_site{1,k},c} ) == 0
                                c = c + 1;
                            end
                        end
                        CC{Op_site{1,k},c} = ['\\gate{',Op_type{1,k},'}'];
                end
            end
        end
    end
end


% complete \qw
[ind1,ind2] = find(cellfun('isempty',CC));
for i = 1:length(ind1)
    CC{ind1(i),ind2(i)} = '\\qw';
end

% trim around the re-initialize state
[reini_ind1,reini_ind2] = find(strcmp(CC,IniState0_e));
for i = 1:length(reini_ind1)
    %     CC{reini_ind1(i),reini_ind2(i)-1} = '';  % delete the \\qw before |0>
    if strcmp(CC{reini_ind1(i),reini_ind2(i)+1},'\\gate{H}') && strcmp(CC{reini_ind1(i),reini_ind2(i)+2},'\\gate{H}')
        CC{reini_ind1(i),reini_ind2(i)+1} = '\\qw'; % delete the repetive H gate
        CC{reini_ind1(i),reini_ind2(i)+2} = '\\qw';
    end
end

% trim long tail of \qw
for i = 1:size(CC,1)
    ind = find(~strcmp(CC(i,:),'\\qw'),1,'last');
    for j = (ind+2):size(CC,2)
        CC{i,j} = '' ;
    end
end



% ------ print ----
fileID = fopen(['quantikz_',filename,'.txt'],'w');
fprintf(fileID, '\\begin{tikzpicture} \n');
fprintf(fileID, '\t\\node[scale = 0.5]{	\n');
fprintf(fileID, '\t\t\\begin{quantikz} [row sep={0.7cm,between origins},column sep=0.12cm]\n');
for i_r = 1:size(CC,1)
    fprintf(fileID, '\t\t\t');
    for i_c = 1:size(CC,2)
        fprintf(fileID, [' & ',CC{i_r,i_c}]);
    end
    fprintf(fileID, ' \\\\   \n');
end
fprintf(fileID, '\t\t\\end{quantikz} \n');
fprintf(fileID, '\t}; \n');
fprintf(fileID, '\\end{tikzpicture}');
fclose(fileID);

fprintf(1,['\t** The quantikz code for latex has been saved as: quantikz_',filename,'.txt ** \n']);
end

