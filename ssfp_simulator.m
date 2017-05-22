function [y,z,dy,dz] = ssfp_simulator(alpha, T1, T2)


%% Dimensions
Npulse = length(alpha);

%% Loop over different T1, T2 values
dy = zeros(Npulse,Npulse,2);

% signal itself
y(:,1) = sin(alpha) ./ (T1/T2 + 1 - cos(alpha) * (T1/T2 - 1));
% derivative wrt. T1
y(:,2) = (T2*sin(alpha).*(cos(alpha) - 1))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^2;
% derivative wrt. T2
y(:,3) = -(T1*sin(alpha).*(cos(alpha) - 1))./(T1 + T2 - T1*cos(alpha) + T2*cos(alpha)).^2;

% Dummy
z = 0;
    
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