function [c,ceq,Dc,Dceq]  = constr_pSSFP(fa, TR, T1, T2, M0)

ceq     = [];                 % Nonlinear equalities at x
Dceq    = [];


if nargout > 2
    [y,~,Dy,~] = MRF_simulator(fa, TR, T1, T2, M0);
else
    y = MRF_simulator(fa, TR, T1, T2, M0);
end

%% Calculate Constraints
if nargout > 2
    y(1:2:end,:) = -y(1:2:end,:);       
    Dy(1:2:end,:,:) = -Dy(1:2:end,:,:);    
    Dc = Dy;
    Dc = permute(Dc, [2 1 3]);
    c = y;
    
%     c = [c, (abs(u').^2 - pi.^2)];
%     Dc = cat(3, Dc, diag(2*u));
    
    c = col(c);
    Dc = Dc(:,:);
else
    y(1:2:end,:) = -y(1:2:end,:);
    c = y;
%     c = [c, (abs(u').^2 - pi.^2)];
    c = col(c);
end