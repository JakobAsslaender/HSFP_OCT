function [y,z,dy,dz] = radial_relaxation_simulator(phi, TR, T1, T2)

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

dDdT2 = @(T,T1,T2,phi)[
 0,...
 0;...
 - (T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*cos(phi)*sin(phi)^2)/(T2^2*(cos(phi)^2 + (T1*sin(phi)^2)/T2)) - (T1*cos(phi)*sin(phi)^2*(exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2)) - 1))/(T2^2*(cos(phi)^2 + (T1*sin(phi)^2)/T2)^2),...
(T*exp(-T*(cos(phi)^2/T1 + sin(phi)^2/T2))*sin(phi)^2)/T2^2];


y  = zeros(Npulse,3,Nspin);        % 'y' component of the magnetization
z  = zeros(Npulse  ,Nspin);        % 'z' component of the magnetization
if nargout > 2
    dr = zeros(Npulse,Npulse,3,Nspin); % derivative of 'r' relatively to phi(i)
    dy = zeros(Npulse,Npulse,3,Nspin); % derivative of 'y' relatively to phi(i)
    dz = zeros(Npulse,Npulse  ,Nspin); % derivative of 'z' relatively to phi(i)
end

for n = Nspin:-1:1
    r     = [1; 1];
    drdT1 = [0; 0];
    drdT2 = [0; 0];
    
    dMT1da = zeros(2,Npulse,Npulse);
    dMT2da = zeros(2,Npulse,Npulse);
    
    if nargout > 2
        drdphi  = ones(2,Npulse+1);
        drdphik = ones(2,Npulse);
        U       = zeros(Npulse+1,2,2,Npulse+1);
        for k = 1:Npulse+1
            U(k,:,:,k) = eye(2);
        end
    end
    
    for k = 1:Npulse
        if nargout > 2
            drdphi (:,2:k) =      D(TR,T1(n),T2(n),phi(k)) * drdphi(:,2:k);
            drdphi (:,k+1) = dDdphi(TR,T1(n),T2(n),phi(k)) * r;
            drdphik(:,  k) = drdphi(:,k+1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dMT1da(:,1:k,k) = D(TR,T1(n),T2(n),phi(k)) * dMT1da(:,1:k,k-1) + dDdT1dphi(TR,T1(n),T2(n),phi(k)) * drdphi(:,1:k);
            dMT2da(:,1:k,k) = D(TR,T1(n),T2(n),phi(k)) * dMT2da(:,1:k,k-1) + dDdT2dphi(TR,T1(n),T2(n),phi(k)) * drdphi(:,1:k);
            
            dMT1da(:,  k,k) =                  dMT1da(:,  k,k)    + Dn * (dR(fa(k)) * dMdT1);
            dMT2da(:,  k,k) =                  dMT2da(:,  k,k)    + Dn * (dR(fa(k)) * dMdT2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            U(k+1,:,:,1:k) = reshape(D(TR,T1(n),T2(n),phi(k)) * reshape(U(k,:,:,1:k), [2 2*k]), [1 2 2 k]);
%             dMdT1 = dDdT1n * RM + Dn * (R(fa(k)) * dMdT1);
%             dMdT2 = dDdT2n * RM + Dn * (R(fa(k)) * dMdT2);
            r     = squeeze(U(k+1,:,:,1))*[1; 1];
        else
            drdT1 = dDdT1(TR,T1(n),T2(n),phi(k)) * r + D(TR,T1(n),T2(n),phi(k)) * drdT1;
            drdT2 = dDdT2(TR,T1(n),T2(n),phi(k)) * r + D(TR,T1(n),T2(n),phi(k)) * drdT2;
            r     =     D(TR,T1(n),T2(n),phi(k)) * r;
        end
        
        y(k,1,n) = r    (2) .* sin(phi(k));
        y(k,2,n) = drdT1(2) .* sin(phi(k));
        y(k,3,n) = drdT2(2) .* sin(phi(k));
        z(k,  n) = r    (2) .* cos(phi(k));
    end

    %% If nargout > 2, the derivative of y relatively to u is also computed
    if nargout > 2
        for k = 1:Npulse % index of the changed pulse
            dM = reshape(U(k+1:Npulse+1,:,:,k+1), [2*(Npulse-k+1) 2]) * drdphik(:,k);
            dr(k:Npulse,k,1,n) = reshape(dM, [(Npulse-k+1) 2])*[0; 1];
        end
        dy(:,:,1,n) = dr(:,:,1,n)   .* repmat(sin(phi).' , [1 Npulse 1 Nspin])...
                    + repmat(diag(y(:,1,n) .* cot(phi).'), [1   1    1 Nspin]);
%         dy(:,:,2,n) = squeeze(dMT1da(2,:,:)).';
%         dy(:,:,3,n) = squeeze(dMT2da(2,:,:)).';
        dy = dy(:,:,1); %dy(:,:,2) = 0;
    end
        y = y(:,1); %y(:,2) = 1;
end


end