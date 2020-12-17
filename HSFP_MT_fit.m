function [F, J] = HSFP_MT_fit(param,t,x,TR,u, bfinite_pulse_correction)

if length(param) > 6
    % complex valued M0...
    M0 = param(1) + 1i * param(2);
    ds = simulate_MT_ODE(x, TR, t, param(3), param(4), param(5), param(6), param(7), true, bfinite_pulse_correction);
    
    s = M0 * ds(:,1);
    F = u*s;
    F = [real(F); imag(F)];
    
    ds(:,2:end) = M0 * ds(:,2:end);
    Jc = u*ds;
    J(:,1) = [real(Jc(:,1)); imag(Jc(:,1))];
    J(:,2) = [-imag(Jc(:,1)); real(Jc(:,1))];
    J(:,3:7) = [real(Jc(:,2:end)); imag(Jc(:,2:end))];    
else
    % real valued M0...
    s = simulate_MT_ODE(x, TR, t, param(2), param(3), param(4), param(5), param(6), true, bfinite_pulse_correction);
    
    J = [s(:,1), param(1) * s(:,2:end)];
    s = param(1) * s(:,1);
    
    F = u*s;
    F = real(F);
    
    J = u*J;
    J = real(J);
end

end
