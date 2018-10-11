function [y,z,dy,dz] = radial_relaxation_simulator_closed_form(theta, TR, T1, T2, r0, B0, B1)


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

if nargout > 2
    dy = zeros(Npulse,Npulse,3,Nspin);
end

%% Loop over different T1, T2 values
for n = Nspin:-1:1
    %% Calculate the spin evolution
    a = exp(- TR * cumsum(stheta.^2/T2(n) + ctheta.^2/T1(n)));   
    int_cta  = cumsum(ctheta./a);
    r = a .* (r0 + TR/T1(n) * int_cta);
    
    %% Calculate the derivative of the spin evolution wrt. T1
    dadT1 =  a    .*       TR        .* cumsum(ctheta.^2/T1(n)^2);
    drdT1 = dadT1 .* (r0 + TR/T1(n)   * int_cta) + ...
             a    .* (   - TR/T1(n)^2 * int_cta + TR/T1(n) * cumsum(-dadT1./a./a .* ctheta));
         
    %% Calculate the derivative of the spin evolution wrt. T2
    dadT2 =  a    .*       TR      .* cumsum(stheta.^2/T2(n)^2);
    drdT2 = dadT2 .* (r0 + TR/T1(n) * int_cta) + ...
             a    .* (     TR/T1(n) * cumsum(-dadT2./a./a .* ctheta));
    
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
        
        drdtheta = - (stheta./a + diag(dadtheta).'./a./a .* ctheta).' * (a .* TR/T1(n));
        drdtheta = triu(drdtheta) + bsxfun(@times, dadtheta, r0 + TR/T1(n) * int_cta.');
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT1dtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * dadT1 + ...
                       sin(2*theta.') * TR * (-1/T1(n)^2)         *  a;
        dadT1dtheta = triu(dadT1dtheta);
        
        drdT1dtheta = - (stheta./a + diag(dadtheta).'./a./a .* ctheta).' * (TR/T1(n) * (dadT1 - a/T1(n))) + ...
                        (stheta.*dadT1./a./a - diag(dadT1dtheta).'./a./a .* ctheta + 2*diag(dadtheta).'.*ctheta .* dadT1./a./a./a).' * (TR/T1(n) * a); 
        drdT1dtheta = triu(drdT1dtheta) + bsxfun(@times, dadT1dtheta, r0 + TR/T1(n) * int_cta.') + ...
                                          bsxfun(@times, dadtheta,    (- TR/T1(n)^2 * int_cta.'  - TR/T1(n) * cumsum(ctheta.*dadT1./a./a).'));
        
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT2dtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * dadT2 + ...
                       sin(2*theta.') * TR * (1/T2(n)^2)          *  a;
        dadT2dtheta = triu(dadT2dtheta);
        
        drdT2dtheta = - (stheta./a + diag(dadtheta).' .* ctheta./a.^2).' * (TR/T1(n) * dadT2) + ...
                        (stheta.*dadT2./a./a - diag(dadT2dtheta).'./a./a .* ctheta + 2*diag(dadtheta).'.*ctheta .* dadT2./a./a./a).' * (TR/T1(n) * a); 
        drdT2dtheta = triu(drdT2dtheta) + bsxfun(@times, dadT2dtheta, r0 + TR/T1(n) * int_cta.') + ...
                                          bsxfun(@times, dadtheta,    (- TR/T1(n) * cumsum(ctheta.*dadT2./a./a).'));

        
        %% covert from spherical to Cartesian coordinates
        dy(:,:,3,n) = bsxfun(@times, drdT2dtheta, stheta).' + diag(drdT2 .* ctheta);
        dy(:,:,2,n) = bsxfun(@times, drdT1dtheta, stheta).' + diag(drdT1 .* ctheta);
        dy(:,:,1,n) = bsxfun(@times, drdtheta,    stheta).' + diag( r    .* ctheta);
        dz(:,:,  n) = bsxfun(@times, drdtheta,    ctheta).' - diag( r    .* stheta);
    end
end