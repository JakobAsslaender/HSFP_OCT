function [y,z,dy,dz] = spgr_simulator(alpha, TR, T1)


%% Dimensions
Npulse = length(alpha);

%% Loop over different T1, T2 values
dy = zeros(Npulse,Npulse,2);

% signal itself
y(:,1) = sin(alpha) ./ (1 - T1/TR .* log(cos(alpha)));
% derivative wrt. T1
% y(:,2) = sin(alpha) .* (- log(cos(alpha))/TR) ./ (1 - T1 * log(cos(alpha))/TR).^2;
y(:,2) = (TR*log(cos(alpha)).*sin(alpha))./(TR - T1*log(cos(alpha))).^2;

% Dummy
z = 0;
    
% Calculate the derivatives wrt. alpha only if we need it
if nargout > 2
    dyda = cos(alpha) ./ (1 - T1 .* log(cos(alpha))/TR) - ...
           sin(alpha) .* (T1 .* tan(alpha)/TR) ./ (1 - T1 .* log(cos(alpha))/TR).^2;
       
%     dydT1da = cos(alpha) .* (- log(cos(alpha))/TR) ./ (1 - T1 * log(cos(alpha))/TR).^2 + ...
%               sin(alpha) .* (      tan(alpha) /TR) ./ (1 - T1 * log(cos(alpha))/TR).^2 - 2 * ...
%               sin(alpha) .* (- log(cos(alpha))/TR) ./ (1 - T1 * log(cos(alpha))/TR).^3 .* (T1 * tan(alpha)/TR);
          
    dydT1da = -(TR.*(TR.*sin(alpha).^2 - TR.*log(cos(alpha)).*cos(alpha).^2 + T1.*log(cos(alpha)).*sin(alpha).^2 + ...
                T1.*log(cos(alpha)).^2.*cos(alpha).^2))./(cos(alpha).*(TR - T1.*log(cos(alpha))).^3);
    
    %% covert from spherical to Cartesian coordinates
    dy(:,:,1) = diag(dyda);
    dy(:,:,2) = diag(dydT1da);
    dz(:,:) = 0;
end
end