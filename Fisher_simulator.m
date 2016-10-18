function [y,z,dy,dz] = Fisher_simulator(fa, TR, T1, T2, M0)

%% Bloch Simulation
%  using the evolution operator for a magnetisation vector of 2+1D dimension

% output : y(k,n) and dy(k,i,n)
%    k : time parameter
%    n : number of the spin
%    i : time parameter for the derivative at this time.

Nspin  = length(T1);
Npulse = length(fa);

R = @(u) [
    1    0        0    ;...
    0  cos(u)  -sin(u) ;...
    0  sin(u)   cos(u)];

dR = @(u) [
    0     0        0    ;...
    0  -sin(u)  -cos(u) ;...
    0   cos(u)  -sin(u)];

D = @(T,T1,T2)[
    1                  0          0      ;...
    0             exp(-T/T2)      0      ;...
    (1-exp(-T/T1))     0      exp(-T/T1)];

dDdT1 = @(T,T1,T2)[
    0                      0             0;...
    0                      0             0;...
    -exp(-T/T1)*(T/T1^2)   0   exp(-T/T1)*(T/T1^2)];

dDdT2 = @(T,T1,T2)[
    0                0               0;...
    0      exp(-T/T2)*(T/T2^2)       0;...
    0                0               0];


y  = zeros(Npulse,3,Nspin);        % 'y' component of the magnetization
z  = zeros(Npulse  ,Nspin);        % 'z' component of the magnetization
if nargout > 2
    dy = zeros(Npulse,Npulse,3,Nspin); % derivative of 'y' relatively to u(i)
    dz = zeros(Npulse,Npulse  ,Nspin); % derivative of 'z' relatively to u(i)
end

for n = Nspin:-1:1
    Dn     =     D(TR,T1(n),T2(n));
    dDdT1n = dDdT1(TR,T1(n),T2(n));
    dDdT2n = dDdT2(TR,T1(n),T2(n));
    
    dMT1da = zeros(3,Npulse,Npulse);
    dMT2da = zeros(3,Npulse,Npulse);
    dRM    = zeros(3,Npulse);
    
    if nargout > 2
        dMk   = zeros(3, Npulse);
        U     = zeros(Npulse+1,3,3,Npulse+1);
        for k = 1:Npulse+1
            U(k,:,:,k) = eye(3);
        end
    end
    
    dMdT1 = [0; 0; 0];
    dMdT2 = [0; 0; 0];
    M     =  M0;
    for k = 1:Npulse
        RM = R(fa(k)) * M;
        if nargout > 2
            dRM(:,1:k-1) = (R(fa(k)) * Dn) * dRM(:,1:k-1);
            dRM(:,  k)   = dR(fa(k))       * M;
            dMk(:,  k)   =             Dn  * dRM(:,k);
            if k==1
                dMT1da(:,1,1) = dDdT1n * dRM(:,1);
                dMT2da(:,1,1) = dDdT2n * dRM(:,1);
            else
                dMT1da(:,1:k,k) = Dn * (R(fa(k)) * dMT1da(:,1:k,k-1)) + dDdT1n * dRM(:,1:k);
                dMT2da(:,1:k,k) = Dn * (R(fa(k)) * dMT2da(:,1:k,k-1)) + dDdT2n * dRM(:,1:k);
                dMT1da(:,  k,k) =                  dMT1da(:,  k,k)    + Dn * (dR(fa(k)) * dMdT1);
                dMT2da(:,  k,k) =                  dMT2da(:,  k,k)    + Dn * (dR(fa(k)) * dMdT2);
            end
            
            U(k+1,:,:,1:k) = reshape(Dn * R(fa(k)) * reshape(U(k,:,:,1:k), [3 3*k]), [1 3 3 k]);
            dMdT1 = dDdT1n * RM + Dn * (R(fa(k)) * dMdT1);
            dMdT2 = dDdT2n * RM + Dn * (R(fa(k)) * dMdT2);
            M     = squeeze(U(k+1,:,:,1))*M0;
        else
            dMdT1 = dDdT1n * RM + Dn * (R(fa(k)) * dMdT1);
            dMdT2 = dDdT2n * RM + Dn * (R(fa(k)) * dMdT2);
            M     =     Dn * RM;
        end
        y(k,1,n) = [0 1 0]*M;
        y(k,2,n) = [0 1 0]*dMdT1;
        y(k,3,n) = [0 1 0]*dMdT2;
        z(k,1,n) = [0 0 1]*M;
    end
    
    %% If nargout > 2, the derivative of y relatively to u is also computed
    if nargout > 2
        for k = 1:Npulse % index of the changed pulse
            dM = reshape(U(k+1:Npulse+1,:,:,k+1), [3*(Npulse-k+1) 3]) * dMk(:,k);
            dy(k:Npulse,k,1,n) = reshape(dM, [(Npulse-k+1) 3])*[0; 1; 0];
            dz(k:Npulse,k,1,n) = reshape(dM, [(Npulse-k+1) 3])*[0; 0; 1];
        end
        dy(:,:,2,n) = squeeze(dMT1da(2,:,:)).';
        dy(:,:,3,n) = squeeze(dMT2da(2,:,:)).';
    end
    
end
end