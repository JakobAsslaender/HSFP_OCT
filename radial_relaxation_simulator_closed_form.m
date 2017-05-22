function [y,z,dy,dz] = radial_relaxation_simulator_analytical(theta, TR, T1, T2, r0)


%% Dimensions
Nspin  = length(T1);
Npulse = length(theta);

%% Some precalculations
stheta = sin(theta);
ctheta = cos(theta);

%% Loop over different T1, T2 values
dy = zeros(Npulse,Npulse,3,Nspin);
for n = Nspin:-1:1
    %% Calculate the spin evolution
    a = exp(- TR * cumsum(stheta.^2/T2(n) + ctheta.^2/T1(n)));   
    int_cta  = cumsum(ctheta./a);
    r = a .* (r0 + TR/T1(n) * int_cta);
    
    %% Calculate the derivative of the spin evolution wrt. T1
    dadT1 =  a    .*       TR        .* cumsum(ctheta.^2/T1(n)^2);
    drdT1 = dadT1 .* (r0 + TR/T1(n)   * int_cta) + ...
             a    .* (   - TR/T1(n)^2 * int_cta + TR/T1(n) * cumsum(-dadT1 .* ctheta./a.^2));
         
    %% Calculate the derivative of the spin evolution wrt. T2
    dadT2 =  a    .*       TR      .* cumsum(stheta.^2/T2(n)^2);
    drdT2 = dadT2 .* (r0 + TR/T1(n) * int_cta) + ...
             a    .* (     TR/T1(n) * cumsum(-dadT2 .* ctheta./a.^2));
    
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
        
        drdtheta = - (stheta./a + diag(dadtheta).' .* ctheta./a.^2).' * (a .* TR/T1(n));
        drdtheta = triu(drdtheta) + bsxfun(@times, dadtheta, r0 + TR/T1(n) * int_cta.');
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT1dtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * dadT1 + ...
                       sin(2*theta.') * TR * (-1/T1(n)^2)         *  a;
        dadT1dtheta = triu(dadT1dtheta);
        
        drdT1dtheta = - (stheta./a + diag(dadtheta).' .* ctheta./a.^2).' * (TR/T1(n) * (dadT1 - a/T1(n))) + ...
                        (stheta./a.^2.*dadT1 - diag(dadT1dtheta).' .* ctheta./a.^2 + 2*diag(dadtheta).'.*ctheta./a.^3 .* dadT1).' * (TR/T1(n) * a); 
        drdT1dtheta = triu(drdT1dtheta) + bsxfun(@times, dadT1dtheta, r0 + TR/T1(n) * int_cta.') + ...
                                          bsxfun(@times, dadtheta,    (- TR/T1(n)^2 * int_cta.'  - TR/T1(n) * cumsum(ctheta./a.^2.*dadT1).'));
        
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT2dtheta = (sin(2*theta.') * TR * (1/T1(n) - 1/T2(n))) * dadT2 + ...
                       sin(2*theta.') * TR * (1/T2(n)^2)          *  a;
        dadT2dtheta = triu(dadT2dtheta);
        
        drdT2dtheta = - (stheta./a + diag(dadtheta).' .* ctheta./a.^2).' * (TR/T1(n) * dadT2) + ...
                        (stheta./a.^2.*dadT2 - diag(dadT2dtheta).' .* ctheta./a.^2 + 2*diag(dadtheta).'.*ctheta./a.^3 .* dadT2).' * (TR/T1(n) * a); 
        drdT2dtheta = triu(drdT2dtheta) + bsxfun(@times, dadT2dtheta, r0 + TR/T1(n) * int_cta.') + ...
                                          bsxfun(@times, dadtheta,    (- TR/T1(n) * cumsum(ctheta./a.^2.*dadT2).'));

        
        %% covert from spherical to Cartesian coordinates
        dy(:,:,3,n) = bsxfun(@times, drdT2dtheta, stheta).' + diag(drdT2 .* ctheta);
        dy(:,:,2,n) = bsxfun(@times, drdT1dtheta, stheta).' + diag(drdT1 .* ctheta);
        dy(:,:,1,n) = bsxfun(@times, drdtheta,    stheta).' + diag( r    .* ctheta);
        dz(:,:,  n) = bsxfun(@times, drdtheta,    ctheta).' - diag( r    .* stheta);
    end
end