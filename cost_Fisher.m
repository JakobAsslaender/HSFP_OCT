function [C, gr] = cost_Fisher(fa, simulator)


C = 0;
if nargout > 1
    [y,~,dy,~] = simulator(fa);
else
    y = simulator(fa);
end

Nspin  = size(y,3);
Npulse = length(fa);

for n=1:Nspin
    I = y(:,:,n)'*y(:,:,n);
    Im1 = inv(I);
    
    % Optimize for the average of all parameters and all spins
%     C = C + sum(diag(Im1));
    
    % Optimize only for the worst parameter of the worst spin
    %     [Cn, param_idx] = max(diag(Im1));
    %     disp(param_idx);
    param_idx = 1;
    Cn = Im1(param_idx, param_idx);
    if Cn > C
        C = Cn;
        n_idx = n;
    end
end

% Normalize the cost
C = C/(Nspin^2);


%% Calculate the gradient
if nargout > 1
    gr = zeros(1,Npulse);
    
    % Optimize for the average of all parameters and all spins
%     for k = 1:Npulse
%         for n = 1:Nspin
%             for m = size(dy,3):-1:1
%                 for l = size(y,2):-1:1
%                     dIda(l,m) = 2 * real(dy(:,k,m,n)'*y(:,l,n));
%                 end
%             end
%             gr(k) = gr(k) - sum(diag(Im1 * dIda * Im1));
%         end
%     end
%     gr = gr /(Nspin^2) ;
    
    
    % Optimize only for the worst parameter of the worst spin
    for k = 1:Npulse
        for n = 1:Nspin
            for m = size(dy,3):-1:1
                for l = size(y,2):-1:1
                    dIda(l,m) = 2 * real(dy(:,k,m,n)'*y(:,l,n));
                end
            end
            crb = diag(Im1 * dIda * Im1);
            gr(k) = gr(k) - crb(param_idx);
        end
    end
    gr = gr /(Nspin^2);
end

end