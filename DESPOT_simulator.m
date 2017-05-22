function [y,z,dy,dz] = DESPOT_simulator(alpha, TR, T1, T2, Npulse_spgr)


%% Dimensions
Npulse = length(alpha);

[y1,~,dy1,~] = spgr_simulator(alpha(1:Npulse_spgr), TR, T1);
[y2,~,dy2,~] = ssfp_simulator(alpha(Npulse_spgr+1:end), T1, T2);

y = zeros(Npulse,3);
y(1:Npulse_spgr,1:2) = y1;
y(Npulse_spgr+1:end,:) = y2;

dy = zeros(Npulse,Npulse,3);
dy(1:Npulse_spgr,1:Npulse_spgr,1:2) = dy1;
dy(Npulse_spgr+1:end,Npulse_spgr+1:end,:) = dy2;

z = 0;
dz = 0;

end