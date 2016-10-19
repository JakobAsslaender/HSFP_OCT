function [y,z,dy,dz] = radial_relaxation_simulator(phi, TR, T1, T2, r0)

%% Bloch Simulation
Nspin  = length(T1);
Npulse = length(phi);


D = @(T,T1,T2,phi)[
    1,... 
    0;...
    cos(phi)/(sin(phi)^2*T1/T2+cos(phi)^2) * (1-exp(-T*(sin(phi)^2/T2 + cos(phi)^2/T1))),... 
    exp(-T*(sin(phi)^2/T2 + cos(phi)^2/T1))];

dDdphi = @(T,T1,T2,phi)[
    0,... 
    0;...
    (sin(phi)*(exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2)) - 1))/(cos(phi)^2 + (T1*sin(phi)^2)/T2) - (cos(phi)*(exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2)) - 1)*(2*cos(phi)*sin(phi) - (2*T1*cos(phi)*sin(phi))/T2))/(cos(phi)^2 + (T1*sin(phi)^2)/T2)^2 - (T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)*((2*cos(phi)*sin(phi))/T1 - (2*cos(phi)*sin(phi))/T2))/(cos(phi)^2 + (T1*sin(phi)^2)/T2),... 
    T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*((2*cos(phi)*sin(phi))/T1 - (2*cos(phi)*sin(phi))/T2)];

dDdT1 = @(T,T1,T2,phi)[
 0,...
 0;...
 (cos(phi)*sin(phi)^2*(exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2)) - 1))/(T2*(cos(phi)^2 + (T1*sin(phi)^2)/T2)^2) - (T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^3)/(T1^2*(cos(phi)^2 + (T1*sin(phi)^2)/T2)),...
 (T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^2)/T1^2];

dDdT1dphi = @(T,T1,T2,phi)[
 0,...
 0;...
 (exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*sin(phi)*(2*T1^3*T2^2*cos(phi)^4 - T1^4*T2*sin(phi)^4 - 2*T^2*T2^3*cos(phi)^8 - 2*T1^3*T2^2*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^4 + 3*T1^3*T2^2*cos(phi)^2*sin(phi)^2 + 2*T^2*T1^3*cos(phi)^4*sin(phi)^4 + 2*T*T1^2*T2^2*cos(phi)^6 + 2*T^2*T1*T2^2*cos(phi)^8 + T1^4*T2*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*sin(phi)^4 - 2*T*T1^4*cos(phi)^2*sin(phi)^4 - 2*T1^4*T2*cos(phi)^2*sin(phi)^2 + T*T1*T2^3*cos(phi)^6 + 5*T*T1^3*T2*cos(phi)^2*sin(phi)^4 - 3*T1^3*T2^2*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^2*sin(phi)^2 + 6*T*T1^2*T2^2*cos(phi)^4*sin(phi)^2 - 4*T^2*T1*T2^2*cos(phi)^6*sin(phi)^2 - 2*T^2*T1^2*T2*cos(phi)^4*sin(phi)^4 + 4*T^2*T1^2*T2*cos(phi)^6*sin(phi)^2 + 2*T1^4*T2*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^2*sin(phi)^2))/(T1^3*(T2*cos(phi)^2 + T1*sin(phi)^2)^3),...
 -(2*T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)*sin(phi)*(T1*T2 + T*T1*cos(phi)^2 - T*T2*cos(phi)^2))/(T1^3*T2)];

dDdT2 = @(T,T1,T2,phi)[
 0,...
 0;...
 - (T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)*sin(phi)^2)/(T2^2*(cos(phi)^2 + (T1*sin(phi)^2)/T2)) - (T1*cos(phi)*sin(phi)^2*(exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2)) - 1))/(T2^2*(cos(phi)^2 + (T1*sin(phi)^2)/T2)^2),...
 (T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*sin(phi)^2)/T2^2];

dDdT2dphi = @(T,T1,T2,phi)[
 0,...
 0;...
 -(exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*sin(phi)*(2*T1^2*T2^3*cos(phi)^4 - T1^3*T2^2*sin(phi)^4 - 2*T1^2*T2^3*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^4 + T1^3*T2^2*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*sin(phi)^4 + 3*T1^2*T2^3*cos(phi)^2*sin(phi)^2 - 2*T1^3*T2^2*cos(phi)^2*sin(phi)^2 - 2*T^2*T1^3*cos(phi)^2*sin(phi)^6 + 2*T^2*T2^3*cos(phi)^6*sin(phi)^2 + 2*T*T1*T2^3*cos(phi)^6 - T*T1^3*T2*sin(phi)^6 + 3*T*T1*T2^3*cos(phi)^4*sin(phi)^2 - 2*T*T1^3*T2*cos(phi)^2*sin(phi)^4 - 3*T1^2*T2^3*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^2*sin(phi)^2 + 2*T1^3*T2^2*exp(T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)^2*sin(phi)^2 + 2*T*T1^2*T2^2*cos(phi)^2*sin(phi)^4 + 4*T^2*T1*T2^2*cos(phi)^4*sin(phi)^4 - 2*T^2*T1*T2^2*cos(phi)^6*sin(phi)^2 + 2*T^2*T1^2*T2*cos(phi)^2*sin(phi)^6 - 4*T^2*T1^2*T2*cos(phi)^4*sin(phi)^4))/(T1*T2^2*(T2*cos(phi)^2 + T1*sin(phi)^2)^3),...
 (2*T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)*sin(phi)*(T1*T2 - T*T1*sin(phi)^2 + T*T2*sin(phi)^2))/(T1*T2^3)];



y  = zeros(Npulse,3,Nspin);        % 'y' component of the magnetization
z  = zeros(Npulse  ,Nspin);        % 'z' component of the magnetization
if nargout > 2
    dy = zeros(Npulse,Npulse,3,Nspin); % derivative of 'y' relatively to phi(i)
    dz = zeros(Npulse,Npulse  ,Nspin); % derivative of 'z' relatively to phi(i)
end

for n = Nspin:-1:1
    r     = [1;r0];
    drdT1 = [0; 0];
    drdT2 = [0; 0];
    
    %            changed pulse,affected signal
    drdT1dphi = zeros(2,Npulse,Npulse);
    drdT2dphi = zeros(2,Npulse,Npulse);
    
    if nargout > 2
        dr      = zeros(Npulse,Npulse);
        drdphi  = ones(2,Npulse);
        drdphik = ones(2,Npulse);
        U       = zeros(Npulse+1,2,2,Npulse);
        for k = 1:Npulse+1
            U(k,:,:,k) = eye(2);
        end
    end
    
    for k = 1:Npulse
        if nargout > 2            
            if k>1
                drdT1dphi(:,1:k,  k) = D(TR,T1(n),T2(n),phi(k)) * drdT1dphi(:,1:k,k-1);
                drdT2dphi(:,1:k,  k) = D(TR,T1(n),T2(n),phi(k)) * drdT2dphi(:,1:k,k-1);
                drdT1dphi(:,1:k-1,k) = drdT1dphi(:,1:k-1,k) + dDdT1(TR,T1(n),T2(n),phi(k)) * drdphi(:,1:k-1);
                drdT2dphi(:,1:k-1,k) = drdT2dphi(:,1:k-1,k) + dDdT2(TR,T1(n),T2(n),phi(k)) * drdphi(:,1:k-1);
            end            
            drdT1dphi(:,k,k) = drdT1dphi(:,k,k) + dDdT1dphi(TR,T1(n),T2(n),phi(k)) * r + dDdphi(TR,T1(n),T2(n),phi(k)) * drdT1;
            drdT2dphi(:,k,k) = drdT2dphi(:,k,k) + dDdT2dphi(TR,T1(n),T2(n),phi(k)) * r + dDdphi(TR,T1(n),T2(n),phi(k)) * drdT2;
            
            drdphi (:,1:k-1) =      D(TR,T1(n),T2(n),phi(k)) * drdphi(:,1:k-1);
            drdphi (:,  k)   = dDdphi(TR,T1(n),T2(n),phi(k)) * r;
            drdphik(:,  k)   = drdphi(:,k);
            
            U(k+1,:,:,1:k) = reshape(D(TR,T1(n),T2(n),phi(k)) * reshape(U(k,:,:,1:k), [2 2*k]), [1 2 2 k]);
        end
        
        drdT1 = dDdT1(TR,T1(n),T2(n),phi(k)) * r + D(TR,T1(n),T2(n),phi(k)) * drdT1;
        drdT2 = dDdT2(TR,T1(n),T2(n),phi(k)) * r + D(TR,T1(n),T2(n),phi(k)) * drdT2;
        r     =     D(TR,T1(n),T2(n),phi(k)) * r;
        
        y(k,1,n) = r    (2) .* sin(phi(k));
        y(k,2,n) = drdT1(2) .* sin(phi(k));
        y(k,3,n) = drdT2(2) .* sin(phi(k));
        z(k,  n) = r    (2) .* cos(phi(k));
    end

    %% If nargout > 2, the derivative of y relatively to u is also computed
    if nargout > 2
        for k = 1:Npulse % index of the changed pulse
            dM = reshape(U(k+1:Npulse+1,:,:,k+1), [2*(Npulse-k+1) 2]) * drdphik(:,k);
            dr(k:Npulse,k) = reshape(dM, [(Npulse-k+1) 2])*[0; 1];
        end
        dz(:,:,  n) = dr(:,:)                     .* repmat(cos(phi).' , [1 Npulse]) - diag(z(:,  n) .* tan(phi).');
        dy(:,:,1,n) = dr(:,:)                     .* repmat(sin(phi).' , [1 Npulse]) + diag(y(:,1,n) .* cot(phi).');
        dy(:,:,2,n) = squeeze(drdT1dphi(2,:,:)).' .* repmat(sin(phi).' , [1 Npulse]) + diag(y(:,2,n) .* cot(phi).');
        dy(:,:,3,n) = squeeze(drdT2dphi(2,:,:)).' .* repmat(sin(phi).' , [1 Npulse]) + diag(y(:,3,n) .* cot(phi).');
    end
end


end