function [y,z,dy,dz] = ssfp_simulator(alpha, TR, T1, T2, B0, B1)

Npulse = length(alpha);

%% B0 and B1 correction + some precalculations
% B1:
alpha = B1 * alpha;
% B0:
alpha = 2 * asin(sqrt(sin(alpha/2).^2 ./ (cos(alpha/2).^2 .* cos(B0*TR/2).^2 + sin(alpha/2).^2)));

%% 
if nargout > 2
    dy = zeros(Npulse,Npulse,2);
end

% signal itself
y(:,1) = sin(alpha) ./ (T1/T2 + 1 - cos(alpha) * (T1/T2 - 1));
% derivative wrt. T1
y(:,2) = (T2*sin(alpha).*(cos(alpha) - 1))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^2;
% derivative wrt. T2
y(:,3) = -(T1*sin(alpha).*(cos(alpha) - 1))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^2;

% Dummy
z = y(:,1).' ./ tan(alpha/2);
    
% Calculate the derivatives wrt. alpha only if we need it
if nargout > 2
%     dyda = cos(alpha) ./ (T1/T2 + 1 - cos(alpha) * (T1/T2 - 1)) - ...
%            sin(alpha) ./ (T1/T2 + 1 - cos(alpha) * (T1/T2 - 1)).^2 .* sin(alpha) .* (T1/T2 - 1);
       
    dyda = (T2*(T2 - T1 + T1*cos(alpha) + T2*cos(alpha)))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^2;
    dydT1da = (T2*(cos(alpha) - 1).*(3*T2 - T1 + T1*cos(alpha) + 3*T2*cos(alpha)))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^3;
    dydT2da = -(T1*(cos(alpha) - 1).*(3*T2 - T1 + T1*cos(alpha) + 3*T2*cos(alpha)))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^3;
       
    
    %% covert from spherical to Cartesian coordinates
    dy(:,:,1) = diag(dyda);
    dy(:,:,2) = diag(dydT1da);
    dy(:,:,3) = diag(dydT2da);
    dz(:,:) = 0;
end
end