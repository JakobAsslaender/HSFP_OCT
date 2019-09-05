function [s, b, ds, db] = simulate_MT_ODE(x, TR, t, m0s, T1, T2f, R, T2s, bDerivatives, bfinite_pulse_correction)

x(2,:) = x(2,:) * 1e-3; % convert ms to s
Tmax = t(end);

if bfinite_pulse_correction
    xfun = @(t) finite_pulse_correction(x, t, TR, Tmax);
    TR = 'finite_pulse_correction';
else
    xfun = @(t) hann_interpolation(t, Tmax, x);
end
xt = xfun(t);

if nargout > 2
    r0 = [m0s-1; m0s; 1; 1; 1; 0; zeros(12+18*length(x(:)),1)];
elseif bDerivatives
    r0 = [m0s-1; m0s; 1; 1; 1; 0; zeros(12,1)];
else
    r0 = [m0s-1; m0s; 1];
end

f = @(t,r) radial_MT_ODE(t,r, xfun, TR, Tmax, m0s, T1, T2f, R, T2s);
% f = @(t,r) radial_MT_ODE_XTrans(t,r, xfun, TR, Tmax, m0s, T1, T2f, R, T2s);

% options = odeset('RelTol', 1e-6);
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
% options = odeset('RelTol', 1e5, 'AbsTol', 1e5, 'MaxStep', TR/10, 'InitialStep', TR/10);
% options = odeset('RelTol', 1e5, 'AbsTol', 1e5, 'MaxStep', 4.5e-4, 'InitialStep', 4.5e-4);

for ir = 1:2
    [~, b] = ode45(f,t,r0,options);
    r0 = b(end,:);
    r0(1:3:end) = -r0(1:3:end); % anti periodic boundary conditions for the free pool
    r0(2:3:end) =  r0(2:3:end) * (1 - pi^2 * T2s/1e3/max(x(2,:))); % periodic boundary conditions for the semi-solid pool, attenuated by the inversion pulse
end

s = b(:,1:3:end); % extract the free pool only

if nargout > 2
    ds = s(:,7:end);
    s = s(:,1:6); % 1st is the signal, 2:6 are the derivatives wrt. T1 etc.
    
    ds = reshape(ds, size(ds,1), [], 2, length(x));
    ds = permute(ds, [1 4 3 2]);
    ds = bsxfun(@times, ds, sin(xt(1,:).'));
    
    Nb = size(ds, 2)-1;
    for ix = 1:6
        for ti = 0:Nb
            w  = (1-cos(pi * (1 - t/Tmax*Nb + ti)))/2;
            w(abs(t/Tmax*Nb - ti) >= 1) = 0;
            ds(:,ti+1,1,ix) = ds(:,ti+1,1,ix) + (s(:,ix) .* cos(xt(1,:).') .* w.');
            ds(:,ti+1,2,ix) = ds(:,ti+1,2,ix) + (s(:,ix) .* cos(xt(1,:).') .* w.');
        end
    end
    ds(:,:,2,:) = ds(:,:,2,:) * 1e-3;
    ds = permute(ds, [1 3 2 4]);
    ds = reshape(ds, size(ds,1), [], 6);
    db = [];
end

s = bsxfun(@times, s, sin(xt(1,:).'));

end
