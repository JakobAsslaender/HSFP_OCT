function dydt = radial_MT_ODE_no_m0f(t,y, theta_t, TRF_t, TR, Tmax, m0s, T1f, T2f, R, T1s, T2s)

theta = theta_t(t);
TRF   = TRF_t(t);
nbar = T2s / 1e3 / TRF / TR;

ct = cos(theta);
st = sin(theta);

% A    = [-ct^2/T1f - st^2/T2f - R*m0s,                R*ct*(1-m0s)           , (1-m0s)*ct/T1f;...
%                 R*ct*m0s            , - (1/T1s + R*(1-m0s) + nbar*4*theta^2),     m0s/T1s   ;...
%                    0                ,                   0                   ,        0      ];
% 
% dAdx = [   -R     ,     -R*ct       ,    -ct/T1f   ;... % dAdm0s
%           R*ct    ,        R        ,      1/T1s   ;...
%             0     ,        0        ,        0     ;...
%         ct^2/T1f^2,        0        , -(1-m0s)*ct/T1f^2;... % dAdT1f
%             0     ,        0        ,        0     ;...
%             0     ,        0        ,        0     ;...
%         st^2/T2f^2,        0        ,        0     ;... % dAdT2f
%             0     ,        0        ,        0     ;...
%             0     ,        0        ,        0     ;...
%           -m0s    ,   ct*(1-m0s)    ,        0     ;... % dAdR
%          ct*m0s   ,    -(1-m0s)     ,        0     ;...
%             0     ,        0        ,        0     ;...
%             0     ,        0        ,        0     ;... % dAdT1s
%             0     ,     1/T1s^2     ,   -m0s/T1s^2 ;...
%             0     ,        0        ,        0     ;...
%             0     ,        0        ,        0     ;... % dAdT2s
%             0     , -4e-3*theta^2/TRF/TR,    0     ;...
%             0     ,        0        ,        0     ];
        

A    = [-ct^2/T1f - st^2/T2f - ct^2*R*m0s,                R*ct*(1-m0s)           , (1-m0s)*ct/T1f;...
                R*ct*m0s                 , - (1/T1s + R*(1-m0s) + nbar*4*theta^2),     m0s/T1s   ;...
                   0                     ,                   0                   ,        0      ];

dAdx = [-ct^2*R   ,     -R*ct       ,    -ct/T1f   ;... % dAdm0s
          R*ct    ,        R        ,      1/T1s   ;...
            0     ,        0        ,        0     ;...
        ct^2/T1f^2,        0        , -(1-m0s)*ct/T1f^2;... % dAdT1f
            0     ,        0        ,        0     ;...
            0     ,        0        ,        0     ;...
        st^2/T2f^2,        0        ,        0     ;... % dAdT2f
            0     ,        0        ,        0     ;...
            0     ,        0        ,        0     ;...
        -ct^2*m0s ,   ct*(1-m0s)    ,        0     ;... % dAdR
         ct*m0s   ,    -(1-m0s)     ,        0     ;...
            0     ,        0        ,        0     ;...
            0     ,        0        ,        0     ;... % dAdT1s
            0     ,     1/T1s^2     ,   -m0s/T1s^2 ;...
            0     ,        0        ,        0     ;...
            0     ,        0        ,        0     ;... % dAdT2s
            0     , -4e-3*theta^2/TRF/TR,    0     ;...
            0     ,        0        ,        0     ];

dydt = A * reshape(y, 3, []);
dydt = dydt(:);

yidx = false(size(y));
yidx([1:21:length(y), 2:21:length(y), 3:21:length(y)]) = true;

dydt(~yidx) = dydt(~yidx) + reshape(dAdx * reshape(y(yidx), 3, []), [], 1);

% Nb = length(y)/42-1.5; % -1 since y contains 21 non-derivative entries and -1 for support points vs. gaps
% for iy = 22:42:length(y)
%     ti = ((iy-1)/42 - 1/2);
%     if abs(t/Tmax*Nb - ti) < 1
%         w  = (1-cos(pi * (1 - t/Tmax*Nb + ti)))/2;
%         
% % Nb = length(y)/21-2; % -1 since y contains 21 non-derivative entries and -1 for support points vs. gaps
% % for iy = 22:21:length(y)
% %     ti = ((iy-1)/21 - 1);
% %     if abs(t/Tmax*Nb - ti) < 1
% %         w  = (1-cos(pi * (1 - t/Tmax*Nb + ti)))/2;
%         
% %         delta = zeros(1,Nb+1); delta(min(max(floor(t/Tmax*Nb),0),Nb)+1) = 1;
% %         w2 = hann_interpolation(t, Tmax, delta);
% %         if (t/Tmax*Nb - ti) < 0
% %             w2 = 1 - w2;
% %         end         
% %         if abs(w - w2) > 1e-10
% %             error;
% %         end
%         
%         
%         % Derivative wrt. theta
%         dAdth = [2*ct*st*w * (1/T1f - 1/T2f),   -R*st*w*(1-m0s)  , -(1-m0s)*st*w/T1f;...
%                         -R*st*w*m0s         , -nbar*8*theta*w,       0      ;...
%                              0              ,        0       ,       0      ];
%     
%         dAdxth = [       0        ,     R*st*w      ,   st*w/T1f    ;... % dAdm0s
%                        -R*st*w    ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;...
%                   -2*ct*st*w/T1f^2,        0        , (1-m0s)*st*w/T1f^2;... % dAdT1f
%                          0        ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;...
%                   2*st*ct*w/T2f^2 ,        0        ,       0       ;... % dAdT2f
%                          0        ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;...
%                          0        ,    -st*w*(1-m0s)    ,       0       ;... % dAdR
%                      -st*w*m0s    ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;... % dAdT1s
%                          0        ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;...
%                          0        ,        0        ,       0       ;... % dAdT2s
%                          0        , 8e-3*theta*w/TRF/TR,       0       ;...
%                          0        ,        0        ,       0       ];
%     
%         dydt(iy  :iy+20) = dydt(iy  :iy+20) + reshape(dAdth * reshape(y(1:21), 3, []), [], 1);
%         dydt(iy+3:iy+20) = dydt(iy+3:iy+20) + dAdxth * y(1:3);
%         
%         % Derivative wrt. TRF aka. TP
%         dAdTP = [0,             0           ,0;...
%                  0, T2s*4e-3*theta^2*w/TRF^2/TR,0 ;...
%                  0,             0           ,0];
% 
%         dAdxTP = [0,          0           , 0;... % dAdm0s
%                   0,          0           , 0;...
%                   0,          0           , 0;...
%                   0,          0           , 0;... % dAdT1f
%                   0,          0           , 0;...
%                   0,          0           , 0;...
%                   0,          0           , 0;... % dAdT2f
%                   0,          0           , 0;...
%                   0,          0           , 0;...
%                   0,          0           , 0;... % dAdR
%                   0,          0           , 0;...
%                   0,          0           , 0;...
%                   0,          0           , 0;... % dAdT1s
%                   0,          0           , 0;...
%                   0,          0           , 0;...
%                   0,          0           , 0;... % dAdT2s
%                   0, -4e-3*theta^2*w/TRF^2/TR, 0;...
%                   0,          0           , 0];
%               
%         dydt(iy+21:iy+41) = dydt(iy+21:iy+41) + reshape(dAdTP * reshape(y(1:21), 3, []), [], 1);
%         dydt(iy+24:iy+41) = dydt(iy+24:iy+41) + dAdxTP * y(1:3);
%     end    
% end
end