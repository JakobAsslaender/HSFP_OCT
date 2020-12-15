function dydt = radial_MT_ODE(t,y, xfun, TR, Tmax, m0s, T1, T2f, R, T2s)

x = xfun(t);
theta = x(1,:);

if ischar(TR)
%     nbar = x(2,:)*1e-3; % The 1e-3 is for converting T2s from ms to s
    if strcmp(TR,'free_precession')
        nbar = x(2,:).*0;
    else
        nbar = x(2,:);
    end
    if length(y) > 18
        error('Calculating the derivative wrt. theta is only implemented for x(2,:) being TRF rather than the attenuation rate');
    end
else
    TRF  = x(2,:);
    nbar = 4*theta^2 / TRF / TR;
%     nbar = 4*theta^2 * 1e-3 / TRF / TR;
end

ct = cos(theta);
st = sin(theta);


A    = [-ct^2/T1 - st^2/T2f - ct^2*R*m0s,               R*ct*(1-m0s)     , (1-m0s)*ct/T1;...
                R*ct*m0s                , - (1/T1 + R*(1-m0s) + nbar*T2s),     m0s/T1   ;...
                   0                    ,                   0            ,        0     ];

dydt = A * reshape(y, 3, []);
dydt = dydt(:);

if length(y) > 3
    dAdx = [-ct^2*R   ,     -R*ct   ,     -ct/T1   ;... % dAdm0s
              R*ct    ,        R    ,       1/T1   ;...
                0     ,        0    ,        0     ;...
             ct^2/T1^2,        0    , -(1-m0s)*ct/T1^2;... % dAdT1
                0     ,     1/T1^2  ,    -m0s/T1^2 ;...
                0     ,        0    ,        0     ;...
            st^2/T2f^2,        0    ,        0     ;... % dAdT2f
                0     ,        0    ,        0     ;...
                0     ,        0    ,        0     ;...
            -ct^2*m0s ,   ct*(1-m0s),        0     ;... % dAdR
             ct*m0s   ,    -(1-m0s) ,        0     ;...
                0     ,        0    ,        0     ;...
                0     ,        0    ,        0     ;... % dAdT2s
                0     ,      -nbar  ,        0     ;...
                0     ,        0    ,        0     ];
    
    
    yidx = false(size(y));
    yidx([1:18:length(y), 2:18:length(y), 3:18:length(y)]) = true;
    
    dydt(~yidx) = dydt(~yidx) + reshape(dAdx * reshape(y(yidx), 3, []), [], 1);
end

if length(y) > 18
    % Derivative wrt. theta
    dAdth = [2*ct*st * (1/T1 - 1/T2f +R*m0s),   -R*st*(1-m0s)  , -(1-m0s)*st/T1;...
                    -R*st*m0s               , -2*T2s*nbar/theta,       0       ;...
                       0                    ,        0         ,       0       ];
    
    dAdxth = [ 2*ct*st * R  ,     R*st        ,     st/T1    ;... % dAdm0s
                   -R*st    ,        0        ,       0      ;...
                   0        ,        0        ,       0      ;...
              -2*ct*st/T1^2 ,        0        , (1-m0s)*st/T1^2;... % dAdT1
                   0        ,        0        ,       0       ;...
                   0        ,        0        ,       0       ;...
              2*st*ct/T2f^2 ,        0        ,       0       ;... % dAdT2f
                   0        ,        0        ,       0       ;...
                   0        ,        0        ,       0       ;...
                 2*ct*st*m0s,    -st*(1-m0s)  ,       0       ;... % dAdR
                 -st*m0s    ,        0        ,       0       ;...
                   0        ,        0        ,       0       ;...
                   0        ,        0        ,       0       ;... % dAdT2s
                   0        ,  -2*nbar/theta  ,       0       ;...
                   0        ,        0        ,       0       ];
    
    
    % Derivative wrt. TRF
    dAdTRF = [0,      0       ,0;...
              0, T2s*nbar/TRF ,0;...
              0,      0       ,0];
    
    dAdxTRF = [0,     0    , 0;... % dAdm0s
               0,     0    , 0;...
               0,     0    , 0;...
               0,     0    , 0;... % dAdT1
               0,     0    , 0;...
               0,     0    , 0;...
               0,     0    , 0;... % dAdT2f
               0,     0    , 0;...
               0,     0    , 0;...
               0,     0    , 0;... % dAdR
               0,     0    , 0;...
               0,     0    , 0;...
               0,     0    , 0;... % dAdT2s
               0, nbar/TRF , 0;...
               0,     0    , 0];
    
    Nb = length(y)/36-1.5; % -1 since y contains 21 non-derivative entries and -1 for support points vs. gaps
    for iy = 19:36:length(y)
        ti = ((iy-1)/36 - 1/2);
        if abs(t/Tmax*Nb - ti) < 1
            w  = (1-cos(pi * (1 - t/Tmax*Nb + ti)))/2;
            dydt(iy   :iy+17) = dydt(iy   :iy+17) + w * reshape(dAdth * reshape(y(1:18), 3, []), [], 1);
            dydt(iy+3 :iy+17) = dydt(iy+3 :iy+17) + w * dAdxth * y(1:3);
            dydt(iy+18:iy+35) = dydt(iy+18:iy+35) + w * reshape(dAdTRF * reshape(y(1:18), 3, []), [], 1);
            dydt(iy+21:iy+35) = dydt(iy+21:iy+35) + w * dAdxTRF * y(1:3);
        end
    end
end