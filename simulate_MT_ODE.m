function [s, b, ds, db] = simulate_MT_ODE(x, TR, t, m0s, T1, T2f, R, T2s, bDerivatives, bfinite_pulse_correction)

% x(2,:) = x(2,:) * 1e-3; % convert ms to s
Tmax = t(end);

if bfinite_pulse_correction
    if size(x, 2) == length(TR:TR:Tmax)
        xtmp = x.';
    else
        xtmp = hann_interpolation(TR:TR:Tmax, Tmax, x).';
    end
    xfun = @(t) finite_pulse_correction(xtmp, t, TR, Tmax);
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



N=floor(Tmax/TR)-1;
fpgrid=TR/2+TR*(0:N).';
rfbeg=fpgrid(2:end)-TR/2-TRF/2;
rfend=fpgrid(2:end)-TR/2+TRF/2;
span=[0;sort([rfbeg;rfend]);Tmax];
span=cat(1,span(1:end-1).',span(2:end).');

f = @(t,r) radial_MT_ODE(t,r, xfun, TR, Tmax, m0s, T1, T2f, R, T2s);
% f = @(t,r) radial_MT_ODE_XTrans(t,r, xfun, TR, Tmax, m0s, T1, T2f, R, T2s);

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
% options = odeset('RelTol', 1e5, 'AbsTol', 1e5, 'MaxStep', TR/10, 'InitialStep', TR/10);
% options = odeset('RelTol', 1e5, 'AbsTol', 1e5, 'MaxStep', 4.5e-4, 'InitialStep', 4.5e-4);

for ir = 1:3
    for j=1:size(span,2)
        switch mod(j,2)
                case 1
                    f = @(t,r) radial_MT_ODE(t,r, xfun, 'free_precession', Tmax, m0s, T1, T2f, R, T2s);
                case 0
                    f = @(t,r) radial_MT_ODE(t,r, xfun, 'rf_pulse', Tmax, m0s, T1, T2f, R, T2s);
        end
        sol = ode45(f,[span(1,j),span(2,j)],r0,options);
        r0 = deval(sol,span(2,j));
        locs=find( t<=span(2,j) & t>=span(1,j));
        b(locs,:)=deval(sol,t(locs));
    end
    
    r0(1:3:end) = -r0(1:3:end); % anti periodic boundary conditions for the free pool
%     r0(2:3:end) =  r0(2:3:end) * (1 - pi^2 * T2s/1e3/max(x(2,:))); % periodic boundary conditions for the semi-solid pool, attenuated by the inversion pulse
    r0(2:3:end) =  r0(2:3:end) * exp(- pi^2 * T2s/max(x(2,:))); % periodic boundary conditions for the semi-solid pool, attenuated by the inversion pulse
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
%     ds(:,:,2,:) = ds(:,:,2,:) * 1e-3;
    ds = permute(ds, [1 3 2 4]);
    ds = reshape(ds, size(ds,1), [], 6);
    db = [];
end

s = bsxfun(@times, s, sin(xt(1,:).'));

end
