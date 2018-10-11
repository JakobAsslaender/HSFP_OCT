function [y,z,dy,dz] = radial_relaxation_simulator_anti_periodic(theta, TR, T1, T2, B0, B1)


%% Dimensions
Nspin  = length(T1);
Npulse = length(theta);

%% B0 and B1 correction + some precalculations
% B1:
theta = B1 * theta;
% B0:
theta = asin(sqrt(sin(theta).^2 ./ (cos(theta).^2 .* cos(B0*TR/2).^2 + sin(theta).^2)));

stheta = sin(theta);
ctheta = cos(theta);

% B1 correction is only implemented in the calculation of r and in none of the
% derivatives.
B1cor = sqrt((stheta(end) * sin(B1*pi/2)^2)^2 + (ctheta(end) * cos(B1*pi))^2);

if nargout > 2
    dy = zeros(Npulse,Npulse,3,Nspin);
end


%% Loop over different T1, T2 values
for n = Nspin:-1:1
    %% Calculate the spin evolution
    a = exp(- TR * cumsum(stheta.^2/T2(n) + ctheta.^2/T1(n)));   
    int_cta  = cumsum(ctheta./a);
    r = TR/T1(n) .* a .* (int_cta - a(end)/(a(end)+1/B1cor) * int_cta(end));
    
    %% Calculate the derivative of the spin evolution wrt. T1
    dadT1 =  a .* TR .* cumsum(ctheta.^2/T1(n)^2);
    dint_dta_dT1 = -cumsum(ctheta .* dadT1./a./a);
    drdT1 = TR/T1(n) .* ((dadT1 - a/T1(n)) .* (int_cta - a(end)/(a(end)+1) .* int_cta(end)) + ...
                         a .* (dint_dta_dT1 - ...
                               dadT1(end)/(a(end)+1) .* int_cta(end) .* (1 + a(end)/(a(end)+1)) - ...
                                a   (end)/(a(end)+1) .* sum(- ctheta .* dadT1./a./a)));
                              
        
    %% Calculate the derivative of the spin evolution wrt. T2
    dadT2 = TR .* a  .* cumsum(stheta.^2/T2(n)^2);
    dint_dta_dT2 = -cumsum(ctheta .* dadT2./a./a);
    drdT2 = TR/T1(n) .* (dadT2 .* (int_cta - a(end)/(a(end)+1) .* int_cta(end)) + ...
                          a    .* (dint_dta_dT2 - ...
                                   dadT2(end)/(a(end)+1) .* int_cta(end) .* (1 + a(end)/(a(end)+1)) - ...
                                    a   (end)/(a(end)+1) .* sum(-dadT2./a./a .* ctheta)));

    
    %% covert from spherical to Cartesian coordinates
    y(:,3,n) = drdT2 .* stheta;
    y(:,2,n) = drdT1 .* stheta;
    y(:,1,n) =  r    .* stheta;
    z = (r .* ctheta).';
    
    % Calculate the derivatives wrt. theta only if we need it
    if nargout > 2        
        %% Calculate the derivative of the spin evolution wrt. theta
        dadtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * a;
        dadtheta = triu(dadtheta);
        
        drdtheta = bsxfun(@times, dadtheta, int_cta.' - a(end)/(a(end)+1) * int_cta(end));
        
        drdtheta = drdtheta - (stheta./a + diag(dadtheta).'./a./a .* ctheta).' * a;
        drdtheta = triu(drdtheta);
        
        drdtheta = drdtheta - ...
                   (dadtheta(:,end).' .* int_cta/(a(end)+1)^2 - ...
                    a(end)/(a(end)+1) .* (stheta./a + ctheta .* diag(dadtheta).'./a./a)).' * a;
        
        drdtheta = drdtheta * TR/T1(n);
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT1dtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * dadT1 + ...
                       sin(2*theta.') * TR * (-1/T1(n)^2)         *  a;
        dadT1dtheta = triu(dadT1dtheta);
        
        drdT1dtheta = - (-stheta./a.^2 .* dadT1 + diag(dadT1dtheta).'./a./a .* ctheta - 2 .* dadT1 .* diag(dadtheta).'./a./a./a .* ctheta).' * a - ...
                        (stheta./a + diag(dadtheta).'./a./a .* ctheta).' * dadT1;
        drdT1dtheta = triu(drdT1dtheta);
        
        drdT1dtheta = drdT1dtheta - ...
                   (dadT1dtheta(:,end).' .* int_cta/(a(end)+1)^2 + dadtheta(:,end).' .* (dint_dta_dT1/(a(end)+1)^2 - 2 .* dadT1(end) .* int_cta/(a(end)+1)^3) - ...
                    dadT1(end)/(a(end)+1)^2 .* ( stheta./a           + ctheta .* diag(dadtheta).'./a./a) - ...
                    a(end)    /(a(end)+1)   .* (-stheta.*dadT1./a./a + ctheta .* diag(dadT1dtheta).'./a./a - 2 .* ctheta .* diag(dadtheta).' .* dadT1./a./a./a)).' * a - ...
                   (dadtheta(:,end).' .* int_cta/(a(end)+1)^2 - ...
                    a(end)/(a(end)+1) .* (stheta./a + ctheta .* diag(dadtheta).'./a./a)).' * dadT1;
        
        drdT1dtheta = drdT1dtheta + bsxfun(@times, dadT1dtheta,      int_cta.' - a(end)/(a(end)+1) * int_cta(end));
        drdT1dtheta = drdT1dtheta + bsxfun(@times,    dadtheta, dint_dta_dT1.' - dadT1(end)/(a(end)+1)^2 * int_cta(end) - a(end)/(a(end)+1) * dint_dta_dT1(end));
        
        drdT1dtheta = drdT1dtheta * TR/T1(n) - drdtheta/T1(n);
        
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT2dtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * dadT2 + ...
                       sin(2*theta.') * TR * (1/T2(n)^2)          *  a;
        dadT2dtheta = triu(dadT2dtheta);
        
        drdT2dtheta = - (-stheta./a.^2 .* dadT2 + diag(dadT2dtheta).'./a./a .* ctheta - 2 .* dadT2 .* diag(dadtheta).'./a./a./a .* ctheta).' * a - ...
                        (stheta./a + diag(dadtheta).'./a./a .* ctheta).' * dadT2;
        drdT2dtheta = triu(drdT2dtheta);
        
        drdT2dtheta = drdT2dtheta - ...
                   (dadT2dtheta(:,end).' .* int_cta/(a(end)+1)^2 + dadtheta(:,end).' .* (dint_dta_dT2/(a(end)+1)^2 - 2 .* dadT2(end) .* int_cta/(a(end)+1)^3) - ...
                    dadT2(end)/(a(end)+1)^2 .* ( stheta./a           + ctheta .* diag(dadtheta).'./a./a) - ...
                    a(end)    /(a(end)+1)   .* (-stheta.*dadT2./a./a + ctheta .* diag(dadT2dtheta).'./a./a - 2 .* ctheta .* diag(dadtheta).' .* dadT2./a./a./a)).' * a - ...
                   (dadtheta(:,end).' .* int_cta/(a(end)+1)^2 - ...
                    a(end)/(a(end)+1) .* (stheta./a + ctheta .* diag(dadtheta).'./a./a)).' * dadT2;
        
        drdT2dtheta = drdT2dtheta + bsxfun(@times, dadT2dtheta,      int_cta.' - a(end)/(a(end)+1) * int_cta(end));
        drdT2dtheta = drdT2dtheta + bsxfun(@times,    dadtheta, dint_dta_dT2.' - dadT2(end)/(a(end)+1)^2 * int_cta(end) - a(end)/(a(end)+1) * dint_dta_dT2(end));
        
        drdT2dtheta = drdT2dtheta * TR/T1(n);

        
        %% covert from spherical to Cartesian coordinates
        dy(:,:,3,n) = bsxfun(@times, drdT2dtheta, stheta).' + diag(drdT2 .* ctheta);
        dy(:,:,2,n) = bsxfun(@times, drdT1dtheta, stheta).' + diag(drdT1 .* ctheta);
        dy(:,:,1,n) = bsxfun(@times, drdtheta,    stheta).' + diag( r    .* ctheta);
        dz(:,:,  n) = bsxfun(@times, drdtheta,    ctheta).' - diag( r    .* stheta);
    end
end