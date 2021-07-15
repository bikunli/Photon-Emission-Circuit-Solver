function [Generators] = LocalMeasurement(Generators, Sign, direction, site)
% This function make measurement on the Generators
% direction can be 'x','y','z', or a row of vector (tableau representation)
% If direction is one of 'x','y' and 'z', site should be specified.
% Sign is either 0 or 1, which represents the measurement outcome '+/-'
% Version: v1.0, Date: 01/17/2021
n = size(Generators.Tableau,1);
g_new = zeros(1,2*n);
switch direction
    case 'x'
        g_new(1,[site,n+site]) = [1,0];
    case 'y'
        g_new(1,[site,n+site]) = [1,1];
    case 'z'
        g_new(1,[site,n+site]) = [0,1];
    otherwise
        g_new = direction;
        if numel(direction) ~= 2*n
            error('Input direction error');
        end
end
P = [zeros(n),eye(n);eye(n),zeros(n)];
AnticommuteTorF = logical(mod(Generators.Tableau * P * g_new',2));
n_A = numel(find(AnticommuteTorF == 1));
switch n_A
    case 0
        % if the measurement is trivial, that doesn't change state
        % then ignore the paradox of contradictory Sign, maintain the
        % previous Generator.
    case 1
        Generators.Tableau(AnticommuteTorF,:) = g_new;
        Generators.SignVector(AnticommuteTorF) = Sign;
    otherwise % if more than one anti-commute g, gauge transform it and go to case 1
        R = logical(eye(n));
        ind1st = find(AnticommuteTorF == 1,1,'first');
        R(:,ind1st) = AnticommuteTorF;
        Generators = GaugeTransformation(Generators,R);
        Generators = LocalMeasurement(Generators, Sign, direction, site);
end
end