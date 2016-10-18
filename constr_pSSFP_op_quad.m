function [c,ceq,Dc,Dceq]  = constr_pSSFP_op_quad(fa, TR, T1, T2, M0)

ceq     = [];                 % Nonlinear equalities at x
Dceq    = [];


if nargout > 2
    [y,z,Dy,Dz] = MRF_simulator(fa, TR, T1, T2, M0);
else
    [y,z] = MRF_simulator(fa, TR, T1, T2, M0);
end

%% Calculate Constraints
if nargout > 2
    y(1:2:end,:) = -y(1:2:end,:);       
    Dy(1:2:end,:,:) = -Dy(1:2:end,:,:);
    
    Dy = Dy .* repmat(reshape(z, [size(z,1) 1 size(z,2)]), [1 size(Dy,2) 1]);
    Dc = Dy + Dz .* repmat(reshape(y, [size(z,1) 1 size(z,2)]), [1 size(Dy,2) 1]);

    Dc = permute(Dc, [2 1 3]);
    
    c = y .* z;
    
    c = col(c);
    Dc = Dc(:,:);
else
    y(1:2:end,:) = -y(1:2:end,:);
    c = y .* z;
    c = col(c);
end