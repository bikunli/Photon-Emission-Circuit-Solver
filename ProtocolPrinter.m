function [] = ProtocolPrinter(Operation,Stat,filename)

n_e = Stat.EmittersNumber;
n_p = Stat.PhotonsNumber;
n_mu = Stat.OperationNumber.Measurement;
n_CNOT = Stat.OperationNumber.CNOT;
n_LC = Stat.OperationNumber.AllEmitterUnitaries - n_CNOT;

fileID = fopen(['ProtocolReport_',filename,'.txt'],'w');
fprintf(fileID,'************************************************\n');
fprintf(fileID,'This report shows the solution of protocol.\n');
fprintf(fileID,[datestr([datetime('now')]),'\n']);
fprintf(fileID,'************************************************\n');
fprintf(fileID,'\n');
fprintf(fileID, '(The following protocol is executed in time order :) \n');
fprintf(fileID, '1 <= j <= n_p labels photons. (n_p = %1d) \n',n_p);
fprintf(fileID, '1 <= i <= n_e labels emitters. (n_e = %1d) \n',n_e);
fprintf(fileID, 'In this protocol solution, there are: \n');
fprintf(fileID, '\t%1d local measurements on emitters.\n', n_mu);
fprintf(fileID, '\t%1d CNOT gates on emitters.\n', n_CNOT);
fprintf(fileID, '\t%1d single qubit gates on emitters. (Caution: redundant repetive gates include (~ total measurement number).) \n', n_LC);
fprintf(fileID,'\n');
fprintf(fileID,'========================================================\n');
fprintf(fileID, 'Perform following operations in time order: \n');
fprintf(fileID, 'Prepare all emitters in |0> state, perform following operations on state vector: \n');
fprintf(fileID,'\n');

n_g = size(Operation.emitters.W0.type,2);
for q = 1:n_g
    gate_str = Operation.emitters.W0.type{n_g - q + 1};
    gate_site = Operation.emitters.W0.site{n_g - q + 1};
    if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
        switch  gate_str
            case 'CNOT' % if control gate
                fprintf(fileID, ['\t',gate_str,' on emitters (i_c,i_t) = (%1d,%1d). \n'], ...
                    gate_site(1)-n_p,gate_site(2)-n_p);
            otherwise % if local clifford
                fprintf(fileID, ['\t',gate_str,' on emitter i = %1d. \n'],gate_site-n_p);
        end
    end
end


for j = 1:n_p
    fprintf(fileID,'--------------\n');
    fprintf(fileID,'* Emit the photon j = %1d, with emitter i = %1d . \n', j, Operation.emitters.EmissionSite{j}-n_p );
    %%%%%%%photon
    gate_str = Operation.photons.Up.type{j};
    if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
        fprintf(fileID, ['\t',Operation.photons.Up.type{j},' on photon j = %1d \n'],j);
    end
    %%%%%%%
    n_g = size(Operation.emitters.Ue.type,2);
    for q = 1:n_g
        gate_str = Operation.emitters.Ue.type{j,n_g - q + 1};
        gate_site = Operation.emitters.Ue.site{j,n_g - q + 1};
        if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
            switch  gate_str
                case 'CNOT' % if control gate
                    fprintf(fileID, ['\t',gate_str,' on emitters (i_c,i_t) = (%1d,%1d). \n'], ...
                        gate_site(1)-n_p,gate_site(2)-n_p);
                otherwise % if local clifford
                    fprintf(fileID, ['\t',gate_str,' on emitter i = %1d. \n'],gate_site-n_p);
            end
        end
    end
    
    % if measurement happens
    if isempty(Operation.emitters.MeasurementSite{j}) == 0 % If there is measurement after emission
        fprintf(fileID, ['\t@ Measure the emitter i = %1d on z-axis, with outcome s_j = 0 or 1. \n'], ...
            Operation.emitters.MeasurementSite{j}-n_p);
        fprintf(fileID, ['\t\t (if s_j = 1, then X on photon j = %1d and emitter i = %1d) \n'], ...
            j,Operation.emitters.MeasurementSite{j}-n_p);
        

       
        % find the last non-empty string, which is the list length of W.
        n_g = find(~cellfun(@isempty,Operation.emitters.W.type(j,:)),1,'last');
        gate_str_temp = Operation.emitters.W.type{j,n_g};
        gate_site_temp = Operation.emitters.W.site{j,n_g};
        ImpendingH = (strcmp(gate_str_temp,'H'))*isequal(gate_site_temp, Operation.emitters.MeasurementSite{j});
        
        switch ImpendingH
            case 0
                % This is an EXTRA Hadamard gate that is needed but not recorded :
                fprintf(fileID, ['\tH on emitter i = %1d. \n'],Operation.emitters.MeasurementSite{j}-n_p);
            case 1
                % The above mentioned EXTRA Hadamard gate is not displayed, if there is an
                % impending Hadamard gate cancels it.
                n_g = n_g - 1; % skip the printing of this impending Hadamard gate.
        end

        for q = 1:n_g
            gate_str = Operation.emitters.W.type{j,n_g - q + 1};
            gate_site = Operation.emitters.W.site{j,n_g - q + 1};
            if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
                switch  gate_str
                    case 'CNOT' % if control gate
                        fprintf(fileID, ['\t',gate_str,' on emitters (i_c,i_t) = (%1d,%1d). \n'], ...
                            gate_site(1)-n_p,gate_site(2)-n_p);
                    otherwise % if local clifford
                        fprintf(fileID, ['\t',gate_str,' on emitter i = %1d. \n'],gate_site-n_p);
                end
            end
        end
    end  
end

fprintf(fileID,'========================================================\n');
fclose(fileID);

fprintf(1,['\n\n \t *** The protocol has been saved as ','ProtocolReport_',filename,'.txt  ***.\n \n']');
end