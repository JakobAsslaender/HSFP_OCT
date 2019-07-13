function [y] = radial_relaxation_simulator_analytical(c, T1, T2, r0)

syms t x x2
theta(t) = poly2sym(c, t);

a(t) = exp(-vpaintegral(sin(theta(x))^2 /T2 + cos(theta(x))^2 /T1, x, 0, t));

r(t) = (a(t)*(T1*r0 + vpaintegral(cos(theta(x2))*exp(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t)))/T1;
drdT1(t) = (exp(-vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(vpaintegral(cos(theta(x2))*exp(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2))*(T1*vpaintegral(1 - cos(theta(x))^2, x, 0, x2) - vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)), x2, 0, t) + T1^2*T2*r0))/(T1^3*T2) - (exp(-vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + vpaintegral(cos(theta(x2))*exp(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t)))/T1^2 - (exp(-vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + vpaintegral(cos(theta(x2))*exp(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t))*(vpaintegral(1 - cos(theta(x))^2, x, 0, t)/(T1*T2) - vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1^2*T2)))/T1;
drdT2(t) = (exp(-vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + vpaintegral(cos(theta(x2))*exp(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t))*(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2^2) - vpaintegral(cos(theta(x))^2, x, 0, t)/(T1*T2)))/T1 + (exp(-vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*vpaintegral(cos(theta(x2))*exp(vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2))*(T2*vpaintegral(cos(theta(x))^2, x, 0, x2) - vpaintegral(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)), x2, 0, t))/(T1^2*T2^2);

% r(t) = (exp(-int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + int(cos(theta(x2))*exp(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t)))/T1;
% drdT1(t) = (exp(-int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(int(cos(theta(x2))*exp(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2))*(T1*int(1 - cos(theta(x))^2, x, 0, x2) - int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)), x2, 0, t) + T1^2*T2*r0))/(T1^3*T2) - (exp(-int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + int(cos(theta(x2))*exp(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t)))/T1^2 - (exp(-int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + int(cos(theta(x2))*exp(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t))*(int(1 - cos(theta(x))^2, x, 0, t)/(T1*T2) - int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1^2*T2)))/T1;
% drdT2(t) = (exp(-int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*(T1*r0 + int(cos(theta(x2))*exp(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2)), x2, 0, t))*(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2^2) - int(cos(theta(x))^2, x, 0, t)/(T1*T2)))/T1 + (exp(-int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, t)/(T1*T2))*int(cos(theta(x2))*exp(int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)/(T1*T2))*(T2*int(cos(theta(x))^2, x, 0, x2) - int(T1 - T1*cos(theta(x))^2 + T2*cos(theta(x))^2, x, 0, x2)), x2, 0, t))/(T1^2*T2^2);

y(1) = r(t) * sin(theta(t));
y(2) = drdT1(t) * sin(theta(t));
y(3) = drdT2(t) * sin(theta(t));

end