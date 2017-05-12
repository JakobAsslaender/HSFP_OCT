function [y,z,dy,dz] = radial_relaxation_simulator_analytical(theta, TR, T1, T2, r0)


%% Dimensions
Nspin  = length(T1);
Npulse = length(theta);

%% Loop over different T1, T2 values
for n = Nspin:-1:1
    a = exp(- TR * cumsum(sin(theta).^2/T2(n) + cos(theta).^2/T1(n)));
    r = a .* (r0 + TR/T1(n) * cumsum(cos(theta)./a));
    
    dadT1 =  a    .*       TR        .* cumsum(cos(theta).^2/T1(n)^2);
    drdT1 = dadT1 .* (r0 + TR/T1(n)   * cumsum(cos(theta)./a)) + ...
             a    .* (   - TR/T1(n)^2 * cumsum(cos(theta)./a) + TR/T1(n) * cumsum(-dadT1 .* cos(theta)./a.^2));
         
    dadT2 =  a    .*       TR      .* cumsum(sin(theta).^2/T2(n)^2);
    drdT2 = dadT2 .* (r0 + TR/T1(n) * cumsum(cos(theta)./a)) + ...
             a    .* (     TR/T1(n) * cumsum(-dadT2 .* cos(theta)./a.^2));
    
    y(:,1,n) =  r    .* sin(theta);
    y(:,2,n) = drdT1 .* sin(theta);
    y(:,3,n) = drdT2 .* sin(theta);
    z = (r .* cos(theta)).';
    
    if nargout > 2
        dadtheta    = zeros(1,Npulse);
        dadT1dtheta = zeros(1,Npulse);
        dadT2dtheta = zeros(1,Npulse);
        
        % changed pulse,affected signal
        drdtheta    = zeros(Npulse);
        drdT1dtheta = zeros(Npulse);
        drdT2dtheta = zeros(Npulse);
        
        int_cta  = cumsum(cos(theta)./a);
        int_ctadT1 = cumsum(cos(theta)./a.^2.*dadT1);
        int_ctadT2 = cumsum(cos(theta)./a.^2.*dadT2);
        
        for k = Npulse:-1:1
            dadtheta(  k:end) = a(k:end) * TR * sin(2*theta(k)) * (1/T1(n) - 1/T2(n));
            drdtheta(k, :   ) = dadtheta .* (r0 + TR/T1(n) * int_cta(k));
            drdtheta(k,k:end) = drdtheta(k,k:end) - a(k:end) .* TR/T1(n) .* (sin(theta(k))/a(k) + dadtheta(k)*cos(theta(k))/a(k)^2);
            
            dadT1dtheta(  k:end) = dadT1(k:end) * TR * sin(2*theta(k)) * (1/T1(n) - 1/T2(n)) + ...
                                    a   (k:end) * TR * sin(2*theta(k)) * (-1/T1(n)^2);
            drdT1dtheta(k, :   ) = dadT1dtheta .* (r0 + TR/T1(n)   * int_cta(k)) + ...
                                   dadtheta    .* (   - TR/T1(n)^2 * int_cta(k)  - ...
                                                        TR/T1(n)   * int_ctadT1(k));
            drdT1dtheta(k,k:end) = drdT1dtheta(k,k:end) - dadT1(k:end) .* TR/T1(n)    .* (sin(theta(k))/a(k) + dadtheta(k)*cos(theta(k))/a(k)^2) + ...
                                                           a   (k:end) .* (TR/T1(n).^2 .* (sin(theta(k))/a(k) + dadtheta(k)*cos(theta(k))/a(k)^2) + ...
                                                                           TR/T1(n)    .* (sin(theta(k))/a(k).^2.*dadT1(k) - dadT1dtheta(k)*cos(theta(k))/a(k)^2 + 2 * dadtheta(k)*cos(theta(k))/a(k)^3 * dadT1(k)));
            
            dadT2dtheta(  k:end) = dadT2(k:end) * TR * sin(2*theta(k)) * (1/T1(n) - 1/T2(n)) + ...
                                    a   (k:end) * TR * sin(2*theta(k)) * (1/T2(n)^2);
            drdT2dtheta(k, :   ) = dadT2dtheta .* (r0 + TR/T1(n)   * int_cta(k)) - ...
                                   dadtheta    .*       TR/T1(n)   * int_ctadT2(k);
            drdT2dtheta(k,k:end) = drdT2dtheta(k,k:end) - dadT2(k:end) .* TR/T1(n)    .* (sin(theta(k))/a(k) + dadtheta(k)*cos(theta(k))/a(k)^2) + ...
                                                           a   (k:end) .* TR/T1(n)    .* (sin(theta(k))/a(k).^2.*dadT2(k) - dadT2dtheta(k)*cos(theta(k))/a(k)^2 + 2 * dadtheta(k)*cos(theta(k))/a(k)^3 * dadT2(k));
        end
        
        dy(:,:,3,n) = (drdT2dtheta .* repmat(sin(theta), [Npulse 1])).' + diag(drdT2 .* cos(theta));
        dy(:,:,2,n) = (drdT1dtheta .* repmat(sin(theta), [Npulse 1])).' + diag(drdT1 .* cos(theta));
        dy(:,:,1,n) = (drdtheta    .* repmat(sin(theta), [Npulse 1])).' + diag(r .* cos(theta));
        dz(:,:,  n) = (drdtheta    .* repmat(cos(theta), [Npulse 1])).' - diag(r .* sin(theta));
    end
end