function [F, J] = HSFP_MT_fit(param,t,x,TR,u)


% s = simulate_MT_ODE_sep_T1(x, TR, t, param(3), param(4), param(5), param(6), param(7), param(8));
s = simulate_MT_ODE(x, TR, t, param(3), param(4), param(5), param(6), param(7));

J = [s(:,1), s(:,1) * 1i, (param(1) + 1i * param(2)) * s(:,2:end)];
s = (param(1) + 1i * param(2)) * s(:,1);

F = u*s;
F = [real(F); imag(F)];

J = u*J;
J = [real(J); imag(J)];

end
