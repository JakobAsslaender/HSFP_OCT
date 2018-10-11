function drdt = radial_Bloch_ODE(t,r, theta_t, T1, T2)

    A     = [-cos(theta_t(t))^2/T1 - sin(theta_t(t))^2/T2, cos(theta_t(t))/T1;
                                                        0,                  0];
    
    dAdT1 = [                      cos(theta_t(t))^2/T1^2,-cos(theta_t(t))/T1^2;
                                                        0,                  0];
                                                    
    dAdT2 = [                      sin(theta_t(t))^2/T2^2,                  0;
                                                        0,                  0];
                                                    
    E = [ A   , zeros(2), zeros(2);...
         dAdT1,   A     , zeros(2);...
         dAdT2, zeros(2),   A     ];
     
    drdt = E * r;
%     drdt = cos(theta_t(t))/T1 - r * (cos(theta_t(t))^2/T1 + sin(theta_t(t))^2/T2);
end