function [y,z,dy,dz] = look_locker_simulator(theta, TR, T1, r0)


%% Dimensions
Nspin  = length(T1);
Npulse = length(theta);

if nargout > 2
    dy = zeros(Npulse,Npulse,3,Nspin);
    dz = dy;
end

%% Loop over different T1, T2 values
for n = Nspin:-1:1
    %% Calculate the spin evolution
    T1s = (1/T1(n) - 1/TR * log(cos(theta))).^(-1);
    a   = exp(- TR * cumsum(1./T1s));
    r   = a .* (r0 + TR/T1(n) * cumsum(1./a));
    
    %% Calculate the derivative of the spin evolution wrt. T1
    dT1sdT1 = (1/T1(n) - 1/TR * log(cos(theta))).^(-2) .* T1(n)^(-2);
    dadT1 =  a    .*       TR        .* cumsum(dT1sdT1./T1s.^2);
    drdT1 = dadT1 .* (r0 + TR/T1(n) * cumsum(1./a))  + ...
             a    .* (   - TR/T1(n)^2 * cumsum(1./a) - TR/T1(n) * cumsum(dadT1./a.^2));
         
    
    %% covert from spherical to Cartesian coordinates
    y(:,2,n) = drdT1 .* sin(theta);
    y(:,1,n) =  r    .* sin(theta);
    z = (r .* cos(theta)).';
    
    % Calculate the derivatives wrt. theta only if we need it
    if nargout > 2        
        %% Calculate the derivative of the spin evolution wrt. theta
        dT1sdtheta = - (1/T1(n) - 1/TR * log(cos(theta))).^(-2) .* 1/TR .* tan(theta);
        dadtheta   = - (dT1sdtheta./T1s.^2).' * a .* TR;
        dadtheta = triu(dadtheta);
        
        drdtheta   = (diag(dadtheta)./a.'./a.') * (a * (- TR/T1(n)));
        drdtheta = triu(drdtheta) + bsxfun(@times, dadtheta, (r0 + TR/T1(n) * cumsum(1./a)).');
        
        %% Calculate the derivative of the spin evolution wrt. T1 and theta
        dadT1dtheta   = -(dT1sdtheta./T1s.^2).' * dadT1 .* TR;                          
        dadT1dtheta = triu(dadT1dtheta);
        
        drdT1dtheta   = -   (diag(dadT1dtheta)./a.'./a.')              * ( a    * (TR/T1(n))) ...
                        + 2*(diag(dadtheta)./a.' .* dadT1.'./a.'./a.') * ( a    * (TR/T1(n))) ...
                        -   (diag(dadtheta)./a.'./a.')                 * (dadT1 * (TR/T1(n))) ...
                        +   (diag(dadtheta)./a.'./a.')                 * ( a    * (TR/T1(n)^2));
        drdT1dtheta   = triu(drdT1dtheta) + bsxfun(@times, dadT1dtheta, (r0 + TR/T1(n) * cumsum(1./a   )).') ...
                                          + bsxfun(@times,    dadtheta, ( - TR/T1(n)^2 * cumsum(1./a   )).') ...
                                          + bsxfun(@times,    dadtheta, ( - TR/T1(n)   * cumsum(dadT1./a.^2)).');
        
%         %% covert from spherical to Cartesian coordinates
        dy(:,:,2,n) = -bsxfun(@times, drdT1dtheta, sin(theta)).' + diag(drdT1 .* cos(theta));
        dy(:,:,1,n) = -bsxfun(@times, drdtheta,    sin(theta)).' + diag( r    .* cos(theta));
        dz(:,:,  n) = -bsxfun(@times, drdtheta,    cos(theta)).' - diag( r    .* sin(theta));
    end
end