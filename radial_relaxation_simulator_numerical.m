function [y,z,dy,dz] = radial_relaxation_simulator_numerical(theta, TR, T1, T2, r0)


%% Dimensions
Nspin  = length(T1);
Npulse = length(theta);
B1     = ones(1, length(T1));

%% Define Spin Evolution Matrices
D = @(T,T1,T2,B1,theta)[
    1,... 
    0;...
    cos(B1*theta)/(sin(B1*theta)^2*T1/T2+cos(B1*theta)^2) * (1-exp(-T*(sin(B1*theta)^2/T2 + cos(B1*theta)^2/T1))),... 
    exp(-T*(sin(B1*theta)^2/T2 + cos(B1*theta)^2/T1))];

dDdtheta = @(T,T1,T2,B1,theta)[
    0,... 
    0;...
    (B1*sin(B1*theta)*(exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2)) - 1))/(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2) - (cos(B1*theta)*(2*B1*cos(B1*theta)*sin(B1*theta) - (2*B1*T1*cos(B1*theta)*sin(B1*theta))/T2)*(exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2)) - 1))/(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2)^2 - (T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)*((2*B1*cos(B1*theta)*sin(B1*theta))/T1 - (2*B1*cos(B1*theta)*sin(B1*theta))/T2))/(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2),...
    T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*((2*B1*cos(B1*theta)*sin(B1*theta))/T1 - (2*B1*cos(B1*theta)*sin(B1*theta))/T2)];

dDdT1 = @(T,T1,T2,B1,theta)[
    0,...
    0;...
    (cos(B1*theta)*sin(B1*theta)^2*(exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2)) - 1))/(T2*(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2)^2) - (T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^3)/(T1^2*(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2)),...
    (T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^2)/T1^2];

dDdT1dtheta = @(T,T1,T2,B1,theta)[
    0,...
    0;...
    (B1*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*sin(B1*theta)*(2*T1^3*T2^2*cos(B1*theta)^4 - T1^4*T2*sin(B1*theta)^4 - 2*T^2*T2^3*cos(B1*theta)^8 - 2*T*T1^4*cos(B1*theta)^2*sin(B1*theta)^4 - 2*T1^4*T2*cos(B1*theta)^2*sin(B1*theta)^2 + 3*T1^3*T2^2*cos(B1*theta)^2*sin(B1*theta)^2 + 2*T^2*T1^3*cos(B1*theta)^4*sin(B1*theta)^4 + T*T1*T2^3*cos(B1*theta)^6 + T1^4*T2*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*sin(B1*theta)^4 + 2*T*T1^2*T2^2*cos(B1*theta)^6 + 2*T^2*T1*T2^2*cos(B1*theta)^8 - 2*T1^3*T2^2*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^4 + 5*T*T1^3*T2*cos(B1*theta)^2*sin(B1*theta)^4 + 2*T1^4*T2*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^2*sin(B1*theta)^2 + 6*T*T1^2*T2^2*cos(B1*theta)^4*sin(B1*theta)^2 - 4*T^2*T1*T2^2*cos(B1*theta)^6*sin(B1*theta)^2 - 2*T^2*T1^2*T2*cos(B1*theta)^4*sin(B1*theta)^4 + 4*T^2*T1^2*T2*cos(B1*theta)^6*sin(B1*theta)^2 - 3*T1^3*T2^2*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^2*sin(B1*theta)^2))/(T1^3*(T2*cos(B1*theta)^2 + T1*sin(B1*theta)^2)^3),...
    -(2*B1*T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)*sin(B1*theta)*(T1*T2 + T*T1*cos(B1*theta)^2 - T*T2*cos(B1*theta)^2))/(T1^3*T2)];
 
dDdT2 = @(T,T1,T2,B1,theta)[
    0,...
    0;...
    -(T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)*sin(B1*theta)^2)/(T2^2*(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2)) - (T1*cos(B1*theta)*sin(B1*theta)^2*(exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2)) - 1))/(T2^2*(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2)^2),...
    (T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*sin(B1*theta)^2)/T2^2];

dDdT2dtheta = @(T,T1,T2,B1,theta)[
    0,...
    0;...
    -(B1*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*sin(B1*theta)*(2*T1^2*T2^3*cos(B1*theta)^4 - T1^3*T2^2*sin(B1*theta)^4 + T1^3*T2^2*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*sin(B1*theta)^4 + 3*T1^2*T2^3*cos(B1*theta)^2*sin(B1*theta)^2 - 2*T1^3*T2^2*cos(B1*theta)^2*sin(B1*theta)^2 - 2*T^2*T1^3*cos(B1*theta)^2*sin(B1*theta)^6 + 2*T^2*T2^3*cos(B1*theta)^6*sin(B1*theta)^2 + 2*T*T1*T2^3*cos(B1*theta)^6 - T*T1^3*T2*sin(B1*theta)^6 - 2*T1^2*T2^3*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^4 + 3*T*T1*T2^3*cos(B1*theta)^4*sin(B1*theta)^2 - 2*T*T1^3*T2*cos(B1*theta)^2*sin(B1*theta)^4 + 2*T*T1^2*T2^2*cos(B1*theta)^2*sin(B1*theta)^4 + 4*T^2*T1*T2^2*cos(B1*theta)^4*sin(B1*theta)^4 - 2*T^2*T1*T2^2*cos(B1*theta)^6*sin(B1*theta)^2 + 2*T^2*T1^2*T2*cos(B1*theta)^2*sin(B1*theta)^6 - 4*T^2*T1^2*T2*cos(B1*theta)^4*sin(B1*theta)^4 - 3*T1^2*T2^3*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^2*sin(B1*theta)^2 + 2*T1^3*T2^2*exp(T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)^2*sin(B1*theta)^2))/(T1*T2^2*(T2*cos(B1*theta)^2 + T1*sin(B1*theta)^2)^3),...
    (2*B1*T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)*sin(B1*theta)*(T1*T2 - T*T1*sin(B1*theta)^2 + T*T2*sin(B1*theta)^2))/(T1*T2^3)];

dDdB1 = @(T,T1,T2,B1,theta)[
    0,... 
    0;...
    (theta*sin(B1*theta)*(exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2)) - 1))/(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2) - (cos(B1*theta)*(2*theta*cos(B1*theta)*sin(B1*theta) - (2*T1*theta*cos(B1*theta)*sin(B1*theta))/T2)*(exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2)) - 1))/(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2)^2 - (T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*cos(B1*theta)*((2*theta*cos(B1*theta)*sin(B1*theta))/T1 - (2*theta*cos(B1*theta)*sin(B1*theta))/T2))/(cos(B1*theta)^2 + (T1*sin(B1*theta)^2)/T2),...
    T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*((2*theta*cos(B1*theta)*sin(B1*theta))/T1 - (2*theta*cos(B1*theta)*sin(B1*theta))/T2)];
 
dDdB1dtheta = @(T,T1,T2,B1,theta)[
    0,... 
    0;...
    (T2*sin(B1*theta)*(exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2)) - 1))/(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2) + (T2*cos(B1*theta)*(T1 - T2)*(sin(2*B1*theta) + 2*B1*theta*cos(2*B1*theta))*(exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2)) - 1))/(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)^2 + (B1*T2*theta*cos(B1*theta)*(exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2)) - 1))/(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2) + (T*exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2))*cos(B1*theta)*(T1 - T2)*(sin(2*B1*theta) + 2*B1*theta*cos(2*B1*theta)))/(T1*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)) - (2*B1*T2*theta*cos(B1*theta)*sin(2*B1*theta)^2*(T1 - T2)^2*(exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2)) - 1))/(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)^3 - (2*B1*T2*theta*sin(B1*theta)*sin(2*B1*theta)*(T1 - T2)*(exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2)) - 1))/(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)^2 - (2*B1*T*theta*exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2))*sin(B1*theta)*sin(2*B1*theta)*(T1 - T2))/(T1*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)) - (2*B1*T*theta*exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2))*cos(B1*theta)*sin(2*B1*theta)^2*(T1 - T2)^2)/(T1*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)^2) - (B1*T^2*theta*exp(-(T*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2))/(T1*T2))*cos(B1*theta)*sin(2*B1*theta)^2*(T1 - T2)^2)/(T1^2*T2*(T2 + T1*sin(B1*theta)^2 - T2*sin(B1*theta)^2)),...
    -(2*T*exp(-T*(cos(B1*theta)^2/T1 + sin(B1*theta)^2/T2))*(T1 - T2)*(T1*T2*cos(B1*theta)*sin(B1*theta) + B1*T1*T2*theta*cos(B1*theta)^2 - B1*T1*T2*theta*sin(B1*theta)^2 - 2*B1*T*T1*theta*cos(B1*theta)^2*sin(B1*theta)^2 + 2*B1*T*T2*theta*cos(B1*theta)^2*sin(B1*theta)^2))/(T1^2*T2^2)];


%% Allocate Memory
r  = zeros(Npulse,3,Nspin);        % radial component of the magnetization and its T1, T2, B1 derivatives
if nargout > 2
    dy = zeros(Npulse,Npulse,3,Nspin); % derivative of 'y' relatively to theta(i)
    dz = zeros(Npulse,Npulse  ,Nspin); % derivative of 'z' relatively to theta(i)
end

%% Loop over different T1, T2, B1 values, if desired
for n = Nspin:-1:1
    ri    = [1;r0];
    drdT1 = [0; 0];
    drdT2 = [0; 0];
    drdB1 = [0; 0];
    
    %            changed pulse,affected signal
    drdT1dtheta = zeros(2,Npulse,Npulse);
    drdT2dtheta = zeros(2,Npulse,Npulse);
    drdB1dtheta = zeros(2,Npulse,Npulse);
    
    if nargout > 2
        drdtheta = zeros(Npulse,Npulse);
        dri      = ones(2,Npulse);
        drk      = ones(2,Npulse);
        U        = zeros(Npulse+1,2,2,Npulse);
        for k = 1:Npulse+1
            U(k,:,:,k) = eye(2);
        end
    end
    
    %% Main loop: Calculate spin evolution
    for k = 1:Npulse
        if nargout > 2            
            if k>1
                drdT1dtheta(:,1:k,  k) = D(TR,T1(n),T2(n),B1(n),theta(k)) * drdT1dtheta(:,1:k,k-1);
                drdT2dtheta(:,1:k,  k) = D(TR,T1(n),T2(n),B1(n),theta(k)) * drdT2dtheta(:,1:k,k-1);
                drdB1dtheta(:,1:k,  k) = D(TR,T1(n),T2(n),B1(n),theta(k)) * drdB1dtheta(:,1:k,k-1);
                
                drdT1dtheta(:,1:k-1,k) = drdT1dtheta(:,1:k-1,k) + dDdT1(TR,T1(n),T2(n),B1(n),theta(k)) * dri(:,1:k-1);
                drdT2dtheta(:,1:k-1,k) = drdT2dtheta(:,1:k-1,k) + dDdT2(TR,T1(n),T2(n),B1(n),theta(k)) * dri(:,1:k-1);
                drdB1dtheta(:,1:k-1,k) = drdB1dtheta(:,1:k-1,k) + dDdB1(TR,T1(n),T2(n),B1(n),theta(k)) * dri(:,1:k-1);
            end            
            drdT1dtheta(:,k,k) = drdT1dtheta(:,k,k) + dDdT1dtheta(TR,T1(n),T2(n),B1(n),theta(k)) * ri + dDdtheta(TR,T1(n),T2(n),B1(n),theta(k)) * drdT1;
            drdT2dtheta(:,k,k) = drdT2dtheta(:,k,k) + dDdT2dtheta(TR,T1(n),T2(n),B1(n),theta(k)) * ri + dDdtheta(TR,T1(n),T2(n),B1(n),theta(k)) * drdT2;
            drdB1dtheta(:,k,k) = drdB1dtheta(:,k,k) + dDdB1dtheta(TR,T1(n),T2(n),B1(n),theta(k)) * ri + dDdtheta(TR,T1(n),T2(n),B1(n),theta(k)) * drdB1;
            
            dri(:,1:k-1) =        D(TR,T1(n),T2(n),B1(n),theta(k)) * dri(:,1:k-1);
            dri(:,  k)   = dDdtheta(TR,T1(n),T2(n),B1(n),theta(k)) * ri;
            drk(:,  k)   = dri(:,k);
            
            U(k+1,:,:,1:k) = reshape(D(TR,T1(n),T2(n),B1(n),theta(k)) * reshape(U(k,:,:,1:k), [2 2*k]), [1 2 2 k]);
        end
        
        drdT1 = dDdT1(TR,T1(n),T2(n),B1(n),theta(k)) * ri + D(TR,T1(n),T2(n),B1(n),theta(k)) * drdT1;
        drdT2 = dDdT2(TR,T1(n),T2(n),B1(n),theta(k)) * ri + D(TR,T1(n),T2(n),B1(n),theta(k)) * drdT2;
        drdB1 = dDdB1(TR,T1(n),T2(n),B1(n),theta(k)) * ri + D(TR,T1(n),T2(n),B1(n),theta(k)) * drdB1;
        ri    =     D(TR,T1(n),T2(n),B1(n),theta(k)) * ri;
        
     
        r(k,1,n) = ri   (2);
        r(k,2,n) = drdT1(2);
        r(k,3,n) = drdT2(2);
        r(k,4,n) = drdB1(2);
    end
    

    %% If nargout > 2, loop once more to calculate all dy/dtheta_k 
    if nargout > 2
        for k = 1:Npulse % index of the changed pulse
            dM = reshape(U(k+1:Npulse+1,:,:,k+1), [2*(Npulse-k+1) 2]) * drk(:,k);
            drdtheta(k:Npulse,k) = reshape(dM, [(Npulse-k+1) 2])*[0; 1];
        end
        dz(:,:,  n) =         drdtheta              .* repmat(cos(B1(n)*theta.'), [1 Npulse]) - diag(r(:,1,n) .* B1(n) .* sin(B1(n)*theta.'));
        dy(:,:,1,n) =         drdtheta              .* repmat(sin(B1(n)*theta.'), [1 Npulse]) + diag(r(:,1,n) .* B1(n) .* cos(B1(n)*theta.'));
        dy(:,:,2,n) = squeeze(drdT1dtheta(2,:,:)).' .* repmat(sin(B1(n)*theta.'), [1 Npulse]) + diag(r(:,2,n) .* B1(n) .* cos(B1(n)*theta.'));
        dy(:,:,3,n) = squeeze(drdT2dtheta(2,:,:)).' .* repmat(sin(B1(n)*theta.'), [1 Npulse]) + diag(r(:,3,n) .* B1(n) .* cos(B1(n)*theta.'));
        dy(:,:,4,n) = squeeze(drdB1dtheta(2,:,:)).' .* repmat(sin(B1(n)*theta.'), [1 Npulse]) + diag(r(:,4,n) .* B1(n) .* cos(B1(n)*theta.')) ...
                    +         drdtheta   .* repmat(theta.' .* cos(B1(n)*theta.'), [1 Npulse]) + diag(r(:,1,n) .*         (cos(B1(n)*theta.') - B1(n) .* theta.' .* sin(B1(n)*theta.')));
    end
end

z = squeeze(r(:,1,:)) .* repmat(cos(B1(n)*theta.'), [1   Nspin]);
y =         r         .* repmat(sin(B1(n)*theta.'), [1 4 Nspin]);

% dy/dB1 has a second term:
y(:,4,:) = y(:,4,:) + r(:,1,:) .* repmat(theta.' .* cos(B1(n)*theta.'), [1 1 Nspin]);

end