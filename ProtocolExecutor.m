function [G_f_recover, G_Phi_recover] = ProtocolExecutor(InvOperation,n_e,n_p)
% feed this function with InvOperation, it gives a final tableau
% all measurement outcomes in this function are random
% G_f_recover is the generators set outcome of 'photons + emitters'
% G_Phi_recover is the generators set outcome of photons only

rng shuffle
n = n_e + n_p;
G_s.Tableau = [zeros(n),eye(n)]; % standard Genrators set: all +Z (See Gottesmann's paper)
G_s.SignVector = zeros(n,1);
G_f = G_s;

n_g = size(InvOperation.emitters.W0.type,2); % number of gates in W0
for q = 1:n_g
    gate_str = InvOperation.emitters.W0.type{n_g + 1 - q};
    if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
        gate_site = InvOperation.emitters.W0.site{n_g + 1 - q};
        G_f = CliffordOperation(G_f,gate_str,gate_site);
    end
end

for j = 1:n_p
    G_f = CliffordOperation(G_f,'CNOT',[InvOperation.emitters.EmissionSite{j}, j]);
    %----- photons LC ----- (non trivial LC is before pi-pulse from LOCC)
    gate_str = InvOperation.photons.Up.type{j};
    G_f = CliffordOperation(G_f, gate_str,j);
    %-----------------
    n_g = size(InvOperation.emitters.Ue.type,2);
    for q = 1:n_g  % Operation: Ue
        gate_str = InvOperation.emitters.Ue.type{j,n_g + 1 - q};
        if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
            gate_site = InvOperation.emitters.Ue.site{j,n_g + 1 - q};
            G_f = CliffordOperation(G_f,gate_str,gate_site);
        end
    end
    
    if isempty(InvOperation.emitters.MeasurementSite{j}) == 0
        s = randi([0,1]); % a random outcome
        mu = InvOperation.emitters.MeasurementSite{j}; % Operation: measurement
        G_f = LocalMeasurement(G_f, s, 'z', mu); % measure on z-axis, with outcome s
        if s == 1 % classical communication
            G_f = CliffordOperation(G_f,'X',mu);%%%%  refresh the emitter as |0> before further emission
            G_f = CliffordOperation(G_f,'X',j); % flip the photon depends on outcome s
        end
        % Operation to recover stabilizer operator +X for the emitter just measured
        G_f = CliffordOperation(G_f,'H',mu); 
        n_g = size(InvOperation.emitters.W.type,2); % W_j can also be H, which just cancel each other 
        for q = 1:n_g % Operation: W_j
            gate_str = InvOperation.emitters.W.type{j,n_g + 1 - q};
            if isempty(gate_str) == 0 && strcmp(gate_str,'I') ~= 1
                gate_site = InvOperation.emitters.W.site{j,n_g + 1 - q};
                G_f = CliffordOperation(G_f,gate_str,gate_site);
            end
        end 
    end
end

G_f_recover = G_f;

[~,R] = ECHELON_TABLEAU(G_f.Tableau);
G_f = GaugeTransformation(G_f,R);
G_Phi_recover.Tableau = G_f.Tableau(1:n_p,[1:n_p,n+(1:n_p)]);
G_Phi_recover.SignVector = G_f.SignVector(1:n_p,1);

end