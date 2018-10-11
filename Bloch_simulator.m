function [xy,z] = Bloch_simulator(alpha, TR, T1, T2, r0, B0, B1)

TE = TR/2;
Rot = @(alpha) [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

xy = zeros(length(alpha),length(T1));
z  = zeros(length(alpha),length(T1));

M = complex(zeros(2, length(T1)));
M(2,:) = r0;  % initialization

% M = Rot(B1*abs(alpha(1))) * M;
% M(1,:) = M(1,:) .* exp(-1i * (angle(alpha(1)) + B0*TE) - TE./T2);
% M(2,:) = 1 + (M(2,:)-1) .* exp(-TE./T1);

for ip=1:length(alpha)
    % Apply RF-pulse (that now acts only on the real part)
    M(1,:) = M(1,:) .* exp(1i*angle(alpha(ip)));  
    Mtmp = Rot(B1*abs(alpha(ip))) * real(M);
    M = Mtmp + 1i*imag(M);
        
    % Free precession and T2 relaxation
    M(1,:) = M(1,:) .* exp(-1i * (angle(alpha(ip)) + B0*TE) - TE./T2);    
    % T1 relaxation
    M(2,:) = 1 + (M(2,:)-1) .* exp(-TE./T1);
    
    % store the signal at the echo time
    xy(ip,:) = M(1,:);
     z(ip,:) = M(2,:);
    
    % Free precession and T2 relaxation 
    M(1,:) = M(1,:) .* exp(-1i*B0*TE - TE./T2);    
    % T1 relaxation
    M(2,:) = 1 + (M(2,:)-1) .* exp(-TE./T1);
end
end