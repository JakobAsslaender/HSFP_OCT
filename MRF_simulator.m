function [y,z,dy,dz] = MRF_simulator(fa, TR, T1, T2, M0)

%% Bloch Simulation
%  using the evolution operator for a magnetisation vector of 2+1D dimension

% output : y(k,n) and dy(k,i,n)
%    k : time parameter
%    n : number of the spin
%    i : time parameter for the derivative at this time.

Nspin  = length(T1);
Npulse = length(fa);

% if (isfield(GLOBAL, 'y') && isfield(GLOBAL, 'u') && all(u == GLOBAL.u) && (nargout < 3 || (isfield(GLOBAL, 'dy') && ~isempty(GLOBAL.dy))))
% if 0
%     y = GLOBAL.y;
%     z = GLOBAL.z;
%     if nargout > 2
%         dy = GLOBAL.dy;
%         dz = GLOBAL.dz;
%     end
% else
    R = @(u) [
        1    0        0    ;...
        0  cos(u)  -sin(u) ;...
        0  sin(u)   cos(u)];
    
    dR = @(u) [
        0     0        0    ;...
        0  -sin(u)  -cos(u) ;...
        0   cos(u)  -sin(u)];
    
    D = @(T,T1,T2)[      
        1            0          0      ;...
        0        exp(-T/T2)     0      ;...
        (1-exp(-T/T1))     0      exp(-T/T1)];
    
    
    y  = zeros(Npulse,Nspin);        % 'y' component of the magnetization
    z  = zeros(Npulse,Nspin);        % 'z' component of the magnetization
    if nargout > 2
        dy = zeros(Npulse,Npulse,Nspin); % derivative of 'y' relatively to u(i)
        dz = zeros(Npulse,Npulse,Nspin); % derivative of 'z' relatively to u(i)
    end
    
    for n = Nspin:-1:1
        Dn = D(TR/2,T1(n),T2(n));
        
        if nargout > 2
            dMk   = zeros(3, Npulse);            
            U = zeros(Npulse+1,3,3,Npulse+1);
            for k = 1:Npulse+1
                U(k,:,:,k) = eye(3);
            end
        end
        
        M = M0;
        for k = 1:Npulse
            % Compute the state at the time 'k'
            if nargout > 2
                dMk(:,k) = Dn * (dR(fa(k)) * (Dn * M));
                
                U(k+1,:,:,1:k) = reshape(Dn *  R(fa(k)) * Dn * reshape(U(k,:,:,1:k), [3 3*k]), [1 3 3 k]);
                M = squeeze(U(k+1,:,:,1))*M0;
            else
                M = Dn *  (R(fa(k)) * (Dn * M));
            end
            y(k,n) = [0 1 0]*M;
            z(k,n) = [0 0 1]*M;
        end
        
        %% If nargout > 1, the derivative of y relatively to u is also computed
        if nargout > 2
            for i = 1:Npulse % index of the changed pulse
                dM = reshape(U(i+1:Npulse+1,:,:,i+1), [3*(Npulse-i+1) 3])*dMk(:,i);
                dy(i:Npulse,i,n) = reshape(dM, [(Npulse-i+1) 3])*[0; 1; 0];
                dz(i:Npulse,i,n) = reshape(dM, [(Npulse-i+1) 3])*[0; 0; 1];
            end
        end
        
    end
    
%     GLOBAL.u = u;
%     GLOBAL.y = y;
%     GLOBAL.z = z;
%     if nargout > 2
%         GLOBAL.dy = dy;
%         GLOBAL.dz = dz;
%     else
%         GLOBAL.dy = [];
%         GLOBAL.dz = [];
%     end
end