function [y,z] = hybrid_state_simulator(theta, TR, T1, T2, r0, B0, B1)


%% Dimensions
Nspin  = length(T1);

%% B0 and B1 correction + some precalculations
% B1:
theta = B1 * theta;
% B0:
theta = asin(sqrt(sin(theta).^2 ./ (cos(theta).^2 .* cos(B0*TR/2).^2 + sin(theta).^2)));

stheta = sin(theta);
ctheta = cos(theta);

if ischar(r0)
    periodic = true;
    if strcmp(r0, 'per')
        beta = 1;
    elseif strcmp(r0, 'anti')
        % B1 correction is only implemented in the calculation of r and in none of the
        % derivatives.
        beta = - sqrt((stheta(end) * sin(B1*pi/2)^2)^2 + (ctheta(end) * cos(B1*pi))^2);
    else
        error('r0 must either be a number, ''per'' or ''anti''');
    end
else
    periodic = false;
end

%% Loop over different T1, T2 values
for n = Nspin:-1:1
    a = exp(- TR * cumsum(stheta.^2/T2(n) + ctheta.^2/T1(n)));
    int_cta  = cumsum(ctheta./a);
    
    if periodic
        r0 =  TR/T1(n) * beta * a(end)/(1 - beta * a(end)) * int_cta(end);
    end
    
    r = a .* (r0 + TR/T1(n) * int_cta);
    y(:,n) = r .* stheta;
    z(:,n) = r .* ctheta;   
    
end
end