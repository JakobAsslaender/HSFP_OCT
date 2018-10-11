function [y,z] = Bloch_simulator_MRF(alpha, TR, T1, T2, r0, B0, B1)

Rot = @(alpha) [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

y = zeros(length(alpha),length(T1));
z = zeros(length(alpha),length(T1));

M = complex(zeros(2, length(T1)));
M(2,:) = r0;  % initialization

for ip=1:length(alpha)
    % Apply RF-pulse (that now acts only on the real part)
    Mtmp = Rot(B1*alpha(ip)*round(cos(ip*pi))) * real(M);
    M = Mtmp + 1i*imag(M);
    
    % store the signal at the echo time
    y(ip,:) = M(1,:)*cos(pi*ip);
    z(ip,:) = M(2,:);
    
    % Free precession and T2 relaxation
    M(1,:) = M(1,:) .* exp(-1i*B0*TR(ip) - TR(ip)./T2);    
    % T1 relaxation
    M(2,:) = 1 + (M(2,:)-1) .* exp(-TR(ip)./T1);
end
end