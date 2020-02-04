function [y,z,dy,dz] = DESPOT_simulator(alpha, TR, T1, T2, Npulse_spgr, B0, B1, FWHM)

if nargin < 8 || isempty(FWHM)
    FWHM = 0;
end

%% Dimensions
Nspin  = length(T1);
Npulse = length(alpha);

y = zeros(Npulse,3,Nspin);
z = zeros(Npulse,  Nspin);

if nargout > 2
    dy = zeros(Npulse,Npulse,3);
    dz = zeros(Npulse,Npulse,3);
end

%% Loop over different T1, T2 values
for n = Nspin:-1:1
    if nargout > 2
%         [y1,z1,dy1,~] = spgr_simulator(alpha(1:Npulse_spgr), TR, T1(n));
        [y1,z1,dy1,~] = spgr_simulator(B1 * alpha(1:Npulse_spgr), TR, T1(n));
        [y2,z2,dy2,~] = ssfp_simulator( alpha(Npulse_spgr+1:end), TR, T1(n), T2(n), B0, B1);
        
        dy(1:Npulse_spgr,1:Npulse_spgr,        1:2,n) = dy1;
        dy(Npulse_spgr+1:end,Npulse_spgr+1:end,  :,n) = dy2;
    else
%         [y1, z1] = spgr_simulator(alpha(1:Npulse_spgr), TR, T1(n));
        [y1, z1] = spgr_simulator(B1 * alpha(1:Npulse_spgr), TR, T1(n));
        [y2, z2] = ssfp_simulator( alpha(Npulse_spgr+1:end), TR, T1(n), T2(n), B0, B1);
    end
    
    y(1:Npulse_spgr,    1:2,n) = y1 * exp(-TR/2*pi*FWHM);
    y(Npulse_spgr+1:end,  :,n) = y2;
    
    z(1:Npulse_spgr,    n) = z1;
    z(Npulse_spgr+1:end,n) = z2;
end
end