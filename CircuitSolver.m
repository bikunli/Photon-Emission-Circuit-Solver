function [Operation,InvOperation,Stat] = CircuitSolver(Generators)
% version: 02/22/2021, v2.1

Tableau_Phi = Generators.Tableau;  % load any Tableau here
SignVector = Generators.SignVector; % Signs for generators in Tableau

n_p = size(Generators.Tableau,1);
h0 = heightfunc(Generators.Tableau); % the original height function
n_e = max(h0);
n = n_p + n_e;

%-------structure variables that record protocol------------
% Operation is the inverted operation on stabilizer operators,
% but Operation is a normal operation on state. (Finally written on report)
Operation.photons.Up.type = cell(n_p,1);
Operation.emitters.Ue.type = cell(n_p,1); % n_p rounds
Operation.emitters.Ue.site = cell(n_p,1);
Operation.emitters.W.type = cell(n_p,1);
Operation.emitters.W.site = cell(n_p,1);
Operation.emitters.W0.type = cell(1,1);
Operation.emitters.W0.site = cell(1,1);
Operation.emitters.EmissionSite = cell(n_p,1);
Operation.emitters.MeasurementSite = cell(n_p,1); % intermidiate measurement, finall measurement not included
% Operation.photons.Nozzlesite = cell(n_e,1); % nozzle photons for PMS

% InvOperation is the inverted operation on state, (call this on simulator)
% but InvOperation is a normal operation on stablizer operators.
InvOperation.photons.Up.type = cell(n_p,1);
InvOperation.emitters.Ue.type = cell(n_p,1); % n_p rounds
InvOperation.emitters.Ue.site = cell(n_p,1);
InvOperation.emitters.W.type = cell(n_p,1);
InvOperation.emitters.W.site = cell(n_p,1);
InvOperation.emitters.W0.type = cell(1,1);
InvOperation.emitters.W0.site = cell(1,1);
InvOperation.emitters.EmissionSite = cell(n_p,1);
InvOperation.emitters.MeasurementSite = cell(n_p,1); % intermidiate measurement, finall measurement not included
%---------------------------------------

Tableau_temp = zeros(n,2*n); % the final Tableau, which is the first tableau in photon absorption process
Tableau_temp(1:n_p, [1:n_p, n+(1:n_p)]) = Tableau_Phi; % evaluate the subblock
Tableau_temp(n_p + (1:n_e),n + n_p + (1:n_e)) = eye(n_e); % g = +Z for all emitters
G_temp.Tableau = Tableau_temp; % target final Generator
G_temp.SignVector = [SignVector ; zeros(n_e,1)];


for k = 0:(n_p-1)
    j = n_p - k;
    [~,R] = ECHELON_TABLEAU(G_temp.Tableau); % echelon gauge
    G_temp = GaugeTransformation(G_temp,R);
    bigram = Tableau2Bigram(G_temp.Tableau);
    h = heightfunc(G_temp.Tableau); % PAI (1)
    
    if j > 1
        dh = h(j) - h(j-1);
    else
        dh = h(j);
    end
    switch dh
        case {0,1}      % PAI (2a)
            % skip
        case -1         % PAI (2b)
            % %             fprintf(1,'Try to fix the descending point ... j = %2d\n', j);
            b = find(bigram(:,1)>n_p, 1, 'first');
            g_b = G_temp.Tableau(b,:);
            % find the nontrivial emitter position on g_b
            e_v_b = find(g_b(n_p+(1:n_e)) + g_b(n+n_p+(1:n_e)) ~= 0);
            mu = n_p + e_v_b(1); % can be other choices
            
            gate_ct = 1;
            switch numel(e_v_b)
                case 1 % turn every Pauli matrix as X
                    s = mu; % site of operation
                    if g_b([s, n + s]) == [1,0]  % if it is X
                        % do nothing
                        gate_str = 'I';
                        gate_str_inv = 'I';
                    elseif g_b([s, n + s]) == [1,1] % if it is Y
                        gate_str = 'P';
                        gate_str_inv = 'PZ';
                    elseif g_b([s, n + s]) == [0,1] % if it is Z
                        gate_str = 'H';
                        gate_str_inv = 'H';
                    end
                    G_temp = CliffordOperation(G_temp,gate_str,s);
                    Operation.emitters.W.type{j,gate_ct} = gate_str;
                    Operation.emitters.W.site{j,gate_ct} = s;
                    InvOperation.emitters.W.type{j,gate_ct} = gate_str_inv;
                    InvOperation.emitters.W.site{j,gate_ct} = s;
                    gate_ct = gate_ct + 1;
                otherwise % numel(e_v_b) > 1, % first make them all X
                    % turn every Pauli matrix as X
                    for i_e = 1:numel(e_v_b)
                        s = n_p + e_v_b(i_e);
                        if g_b([s, n + s]) == [1,0]  % if it is X
                            % do nothing
                            gate_str = 'I';
                            gate_str_inv = 'I';
                        elseif g_b([s, n + s]) == [1,1] % if it is Y
                            gate_str = 'P';
                            gate_str_inv = 'PZ';
                        elseif g_b([s, n + s]) == [0,1] % if it is Z
                            gate_str = 'H';
                            gate_str_inv = 'H';
                        end
                        G_temp = CliffordOperation(G_temp,gate_str,s);
                        Operation.emitters.W.type{j,gate_ct} = gate_str;
                        Operation.emitters.W.site{j,gate_ct} = s;
                        InvOperation.emitters.W.type{j,gate_ct} = gate_str_inv;
                        InvOperation.emitters.W.site{j,gate_ct} = s;
                        gate_ct = gate_ct + 1;
                    end
                    % Control gates on emitters: CNOT X-X pairs becomes X-I:
                    e_v_b_remain = setdiff(e_v_b,mu - n_p);
                    for i_e = 1:numel(e_v_b_remain)
                        s = n_p + e_v_b_remain(i_e);
                        gate_str = 'CNOT';
                        gate_str_inv = 'CNOT';
                        G_temp = CliffordOperation(G_temp,gate_str,[mu,s]);
                        Operation.emitters.W.type{j,gate_ct} = gate_str;
                        Operation.emitters.W.site{j,gate_ct} = [mu,s];
                        InvOperation.emitters.W.type{j,gate_ct} = gate_str_inv;
                        InvOperation.emitters.W.site{j,gate_ct} = [mu,s];
                        gate_ct = gate_ct + 1;
                    end
            end
            % turn (-)X as (+)X, if the sign is incorrect
            if G_temp.SignVector(b) == 1
                gate_str = 'Z';
                gate_str_inv = 'Z';
                G_temp = CliffordOperation(G_temp,gate_str,mu);
                Operation.emitters.W.type{j,gate_ct} = gate_str;
                Operation.emitters.W.site{j,gate_ct} = mu;
                InvOperation.emitters.W.type{j,gate_ct} = gate_str_inv;
                InvOperation.emitters.W.site{j,gate_ct} = mu;
                gate_ct = gate_ct + 1;
            end
            % turn (+)I-X as (+)X-X by CNOT, which will 'boost' h(x)
            InvOperation.emitters.MeasurementSite{j} = mu;
            Operation.emitters.MeasurementSite{j} = mu;
            G_temp = CliffordOperation(G_temp,'CNOT',[mu,j]);
% %             % refresh echelon gauge so we can find g_a afterwards
            [~,R] = ECHELON_TABLEAU(G_temp.Tableau); % echelon gauge
            G_temp = GaugeTransformation(G_temp,R);
            bigram = Tableau2Bigram(G_temp.Tableau);
    end
    
    
    a = find(bigram(:,1) == j, 1, 'last');  % find the index of g_a
    g_a = G_temp.Tableau(a,:); %
    % LC on photon,
    % % turn the Pauli matrix to Z:
    if g_a([j, n + j]) == [1,0] % if it is X
        gate_str = 'H';
        gate_str_inv = 'H';
    elseif g_a([j, n + j]) == [1,1] % if it is Y
        gate_str = 'PH';
        gate_str_inv = 'HPZ';
    elseif g_a([j, n + j]) == [0,1] % if it is Z
        gate_str = 'I';
        gate_str_inv = 'I';
    end
    G_temp = CliffordOperation(G_temp,gate_str,j);
    Operation.photons.Up.type{j} = gate_str;
    InvOperation.photons.Up.type{j} = gate_str_inv;
    
    e_v_a = find(g_a(n_p+(1:n_e)) + g_a(n+n_p+(1:n_e)) ~= 0); % find the nontrivial emitter position on g_a
    eta = n_p + e_v_a(1); % can be other choices
    
    Operation.emitters.EmissionSite{j} = eta;
    InvOperation.emitters.EmissionSite{j} = eta;
    
    % Clifford operations on emitters: find U_e
    gate_ct = 1; % counter for history recording
    switch numel(e_v_a)
        case 1
            % turn the Pauli matrix to Z
            s = eta;
            if g_a([s, n + s]) == [1,0] % if it is X
                gate_str = 'H';
                gate_str_inv = 'H';
            elseif g_a([s, n + s]) == [1,1] % if it is Y
                gate_str = 'PH';
                gate_str_inv = 'HPZ';
            elseif g_a([s, n + s]) == [0,1] % if it is Z
                % do nothing
                gate_str = 'I';
                gate_str_inv = 'I';
            end
            G_temp = CliffordOperation(G_temp,gate_str,s);
            Operation.emitters.Ue.type{j,gate_ct} = gate_str;
            Operation.emitters.Ue.site{j,gate_ct} = s;
            InvOperation.emitters.Ue.type{j,gate_ct} = gate_str_inv;
            InvOperation.emitters.Ue.site{j,gate_ct} = s;
            gate_ct = gate_ct + 1;
        otherwise % if there are more than one non-trivial Paulimatrices
            % first make them all Z
            for i_e = 1:numel(e_v_a)  % LC on emitters
                s = n_p + e_v_a(i_e);
                g_a_emitterpauli = g_a([s, n + s]);
                if g_a_emitterpauli == [1,0] % if it is X
                    gate_str = 'H';
                    gate_str_inv = 'H';
                elseif g_a_emitterpauli == [1,1] % if it is Y
                    gate_str = 'PH';
                    gate_str_inv = 'HPZ';
                elseif g_a_emitterpauli == [0,1] % if it is Z
                    % do nothing
                    gate_str = 'I';
                    gate_str_inv = 'I';
                end
                G_temp = CliffordOperation(G_temp,gate_str,s);
                Operation.emitters.Ue.type{j,gate_ct} = gate_str;
                Operation.emitters.Ue.site{j,gate_ct} = s;
                InvOperation.emitters.Ue.type{j,gate_ct} = gate_str_inv;
                InvOperation.emitters.Ue.site{j,gate_ct} = s;
                gate_ct = gate_ct + 1;
            end
            % Control gates within emitters: CNOT Z-Z pairs becomes Z-I
            e_v_a_remain = setdiff(e_v_a, eta - n_p); % target emitters qubit
            for i_e = 1:numel(e_v_a_remain)
                s = n_p + e_v_a_remain(i_e);
                gate_str = 'CNOT';
                gate_str_inv = 'CNOT';
                G_temp = CliffordOperation(G_temp,gate_str,[s,eta]);
                Operation.emitters.Ue.type{j,gate_ct} = gate_str;
                Operation.emitters.Ue.site{j,gate_ct} = [s,eta];
                InvOperation.emitters.Ue.type{j,gate_ct} = gate_str_inv;
                InvOperation.emitters.Ue.site{j,gate_ct} = [s,eta];
                gate_ct = gate_ct + 1;
            end
    end
    
    % Correct the minus sign of g_a = (-)Z-Z to (+)Z-Z, by flipping the emitter:
    % (flipping the photon qubit seems available)
    if G_temp.SignVector(a) == 1
        gate_str = 'X';
        gate_str_inv = 'X';
        G_temp = CliffordOperation(G_temp,gate_str,eta);
        Operation.emitters.Ue.type{j,gate_ct} = gate_str;
        Operation.emitters.Ue.site{j,gate_ct} = eta;
        InvOperation.emitters.Ue.type{j,gate_ct} = gate_str_inv;
        InvOperation.emitters.Ue.site{j,gate_ct} = eta;
        gate_ct = gate_ct + 1;
    end
    
    %%%%%
    % Absorption (inverted emission):
    G_temp = CliffordOperation(G_temp,'CNOT',[eta,j]);
    %%%%
    
    
    R = eye(n);
    R(:,a) =  (G_temp.Tableau(:,j) + G_temp.Tableau(:,n+j))~=0 ;
    G_temp = GaugeTransformation(G_temp,R);
    
end

%-----------------------------------
% Find W0:

[~,R] = ECHELON_TABLEAU(G_temp.Tableau); % echelon gauge
G_temp = GaugeTransformation(G_temp,R);
% bigram = Tableau2Bigram(G_temp.Tableau);
gate_ct = 1;
for i_r = (1:n_e) + n_p
    g_c = G_temp.Tableau(i_r,:);
    e_v_c = find(g_c(n_p+(1:n_e)) + g_c(n+n_p+(1:n_e)) ~= 0); % find the nontrivial position on g_c
    kappa = n_p + e_v_c(1);
    switch numel(e_v_c)
        case 1
            % turn Pauli matrix as Z
            if g_c([kappa, n + kappa]) == [1,0] % if it is X
                gate_str = 'H';
                gate_str_inv = 'H';
            elseif g_c([kappa, n + kappa]) == [1,1] % if it is Y
                gate_str = 'PH';
                gate_str_inv = 'HPZ';
            elseif g_c([kappa, n + kappa]) == [0,1] % if it is Z
                % do nothing
                gate_str = 'I';
                gate_str_inv = 'I';
            else
                error('Unexpected value!')
            end
            G_temp = CliffordOperation(G_temp,gate_str,kappa);
            Operation.emitters.W0.type{1,gate_ct} = gate_str;
            Operation.emitters.W0.site{1,gate_ct} = kappa;
            InvOperation.emitters.W0.type{1,gate_ct} = gate_str_inv;
            InvOperation.emitters.W0.site{1,gate_ct} = kappa;
            gate_ct = gate_ct + 1;
        otherwise  % numel(e_v_c) > 1
            for i_e = 1:numel(e_v_c)  % LC on emitters
                % turn every Pauli matrix as X
                s = n_p + e_v_c(i_e);
                g_c_emitterpauli = g_c([s, n + s]);
                % turn Pauli matrix as X
                if g_c_emitterpauli == [1,0] % if it is X
                    % do nothing
                    gate_str = 'I';
                    gate_str_inv = 'I';
                elseif g_c_emitterpauli == [1,1] % if it is Y
                    gate_str = 'P';
                    gate_str_inv = 'PZ';
                elseif g_c_emitterpauli == [0,1] % if it is Z
                    gate_str = 'H';
                    gate_str_inv = 'H';
                end
                G_temp = CliffordOperation(G_temp,gate_str,s);
                Operation.emitters.W0.type{1,gate_ct} = gate_str;
                Operation.emitters.W0.site{1,gate_ct} = s;
                InvOperation.emitters.W0.type{1,gate_ct} = gate_str_inv;
                InvOperation.emitters.W0.site{1,gate_ct} = s;
                gate_ct = gate_ct + 1;
            end
            % Control gates on emitters: CNOT X-X pairs becomes X-I
            e_v_c_remain = setdiff(e_v_c,kappa - n_p); % target emitters qubit
            for i_e = 1:numel(e_v_c_remain)
                s = n_p + e_v_c_remain(i_e);
                gate_str = 'CNOT';
                gate_str_inv = 'CNOT';
                G_temp = CliffordOperation(G_temp,gate_str,[kappa,s]);
                Operation.emitters.W0.type{1,gate_ct} = gate_str;
                Operation.emitters.W0.site{1,gate_ct} = [kappa,s];
                InvOperation.emitters.W0.type{1,gate_ct} = gate_str_inv;
                InvOperation.emitters.W0.site{1,gate_ct} = [kappa,s];
                gate_ct = gate_ct + 1;
            end
            
            % Recover X as Z
            gate_str = 'H';
            gate_str_inv = 'H';
            G_temp = CliffordOperation(G_temp,gate_str,kappa);
            Operation.emitters.W0.type{1,gate_ct} = gate_str;
            Operation.emitters.W0.site{1,gate_ct} = kappa;
            InvOperation.emitters.W0.type{1,gate_ct} = gate_str_inv;
            InvOperation.emitters.W0.site{1,gate_ct} = kappa;
            gate_ct = gate_ct + 1;
    end
    
    % doing row transformation, which eliminate all other Zs in the same column
    R = logical(eye(n));
    R(:,i_r) = all(G_temp.Tableau(:,[kappa,n+kappa]) == [0,1],2);
    G_temp = GaugeTransformation(G_temp,R);
    
end

[~,R] = ECHELON_TABLEAU(G_temp.Tableau); % echelon gauge
G_temp = GaugeTransformation(G_temp,R);
% so far the Tableau is all Z, only the sign of emitters need to be fixed:
for i_e = (1:n_e) + n_p
    if G_temp.SignVector(i_e) == 1
        gate_str = 'X';
        gate_str_inv = 'X';
    else
        gate_str = 'I';
        gate_str_inv = 'I';
    end
    G_temp = CliffordOperation(G_temp,gate_str,i_e);
    Operation.emitters.W0.type{1,gate_ct} = gate_str;
    Operation.emitters.W0.site{1,gate_ct} = i_e;
    InvOperation.emitters.W0.type{1,gate_ct} = gate_str_inv;
    InvOperation.emitters.W0.site{1,gate_ct} = i_e;
    gate_ct = gate_ct + 1;
end




if all(G_temp.Tableau == [zeros(n),eye(n)],'all') && all(G_temp.SignVector == [zeros(n,1)],'all')
    % fprintf('The Generators set of ''photon + emitter'' has recovered as ''standard generator set''! \n');
    fprintf(1,'\t** The protocol is SOLVED correctly! **\n');
else
    % fprintf('The Generators set of ''photon + emitter'' is NOT yet a ''standard generator set''! \n');
    error('\t** The protocol is NOT solved correctly! **\n');
end

%----
Stat.PhotonsNumber = n_p;
Stat.EmittersNumber = n_e;
Stat.HeightFunc = h0;
Stat.OperationNumber.Ue = numel( find( ...
    ~cellfun('isempty', Operation.emitters.Ue.type) .* ~strcmp(Operation.emitters.Ue.type,'I')));
Stat.OperationNumber.W = numel( find( ...
    ~cellfun('isempty', Operation.emitters.W.type) .* ~strcmp(Operation.emitters.W.type,'I')));
Stat.OperationNumber.W0 = numel( find( ...
    ~cellfun('isempty', Operation.emitters.W0.type) .* ~strcmp(Operation.emitters.W0.type,'I')));
Stat.OperationNumber.AllEmitterUnitaries = Stat.OperationNumber.Ue + Stat.OperationNumber.W + Stat.OperationNumber.W0;
Stat.OperationNumber.Measurement = numel(find(~cellfun('isempty',Operation.emitters.MeasurementSite)));
Stat.OperationNumber.CNOT =  numel(find(strcmp(Operation.emitters.Ue.type,'CNOT'))) ...
    + numel(find(strcmp(Operation.emitters.W.type,'CNOT'))) + numel(find(strcmp(Operation.emitters.W0.type,'CNOT')));
%----

end
