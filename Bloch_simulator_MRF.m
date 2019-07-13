function [y,z] = Bloch_simulator_MRF(alpha, TR, T1, T2, r0, B0, B1)

TR = TR/2; % actually TE now

y = zeros(length(alpha),length(T1));
z = zeros(length(alpha),length(T1));

M = complex(zeros(2, length(T1)));
M(2,:) = r0;  % initialization

ca = reshape(cos(B1*alpha), 1,1,[]);
sa = reshape(sin(B1*alpha), 1,1,[]);
R = [ca -sa; sa ca];

Inv = [-sin(B1 * pi/2).^2 0; 0 cos(B1 * pi)];

for ip=1:length(alpha)
    % Apply RF-pulse (that now acts only on the real part)
    if abs(alpha(ip)) > 0.9 * pi
        Mtmp = Inv * real(M);
    else
%     Mtmp = Rot(B1*alpha(ip)*round(cos(ip*pi))) * real(M);
        Mtmp = R(:,:,ip) * real(M);
    end
    M = Mtmp + 1i*imag(M);
    
    % Free precession and T2 relaxation
    M(1,:) = M(1,:) .* exp(-TR(ip) * (1i*B0 + 1./T2));   
    % T1 relaxation
    M(2,:) = 1 + (M(2,:)-1) .* exp(-TR(ip)./T1);
    
    % store the signal at the echo time
%     y(ip,:) = M(1,:)*cos(pi*ip);
    y(ip,:) = M(1,:);
    z(ip,:) = M(2,:);
    
    % Free precession and T2 relaxation
    M(1,:) = M(1,:) .* exp(-TR(ip) * (1i*B0 + 1./T2));
    % T1 relaxation
    M(2,:) = 1 + (M(2,:)-1) .* exp(-TR(ip)./T1);
end
end