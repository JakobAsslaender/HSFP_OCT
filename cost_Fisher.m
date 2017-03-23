function [C, gr] = cost_Fisher(fa, simulator, weights, correlate)

if nargin < 4 || isempty(correlate)
    correlate = ones(size(weights));
end

C = 0;
if nargout > 1
    [y,z,dy,~] = simulator(fa);
else
    y = simulator(fa);
end

Nspin  = size(y,3);
Npulse = length(fa);

for n=1:Nspin
    I = y(:,correlate,n)'*y(:,correlate,n);
    Im1 = inv(I);
    
    % Optimize for the average of all parameters and all spins
    C = C + sum(diag(Im1) .* weights(:,n));
    
    % Optimize only for the worst parameter of the worst spin
    %     [Cn, param_idx] = max(diag(Im1));
    %     disp(param_idx);
%     Cn = Im1(param_idx, param_idx);
%     if Cn > C
%         C = Cn;
%         n_idx = n;
%     end
end

% Normalize the cost
C = C/Nspin;


%% Calculate the gradient
if nargout > 1
    gr = zeros(1,Npulse);
    
    % Optimize for the average of all parameters and all spins
    for k = 1:Npulse
        for n = 1:Nspin
            for m = length(correlate):-1:1
                for l = length(correlate):-1:1
                    dIda(l,m) = 2 * real(dy(:,k,correlate(m),n)'*y(:,correlate(l),n));
                end
            end
            gr(k) = gr(k) - sum(diag(Im1 * dIda * Im1) .* weights(:,n));
        end
    end
    gr = gr /(Nspin^2) ;
    
    
    % Optimize only for the worst parameter of the worst spin
%     for k = 1:Npulse
%         for n = 1:Nspin
%             for m = size(dy,3):-1:1
%                 for l = size(y,2):-1:1
%                     dIda(l,m) = 2 * real(dy(:,k,m,n)'*y(:,l,n));
%                 end
%             end
%             crb = diag(Im1 * dIda * Im1);
%             gr(k) = gr(k) - crb(param_idx);
%         end
%     end
%     gr = gr /(Nspin^2);
end


%% Plot result
% y(:,2,:) = -y(:,2,:) .* repmat(reshape(T1,1,1,[]), [size(y,1) 1 1]);
% y(:,3,:) =  y(:,3,:) .* repmat(reshape(T2,1,1,[]), [size(y,1) 1 1]);

% figure(1); 
% subplot(2,2,2); hold off;
% plot(fa    /pi);
% xlabel('t (s)'); ylabel('\theta/\pi'); %legend('initial','optimized');
% 
% subplot(2,2, [1 3]);
% hold off; 
% plot(squeeze(y(:,1,:)),squeeze(z), 'o-');
% hold all; 
% plot(sin(0:.01:pi), cos(0:.01:pi), 'black');
% plot(-sin(0:.01:pi), cos(0:.01:pi), 'black');
% plot([0 0], [-1 1], 'black');
% plot([-1 1], [0 0], 'black');
% % plot( sqrt(T2/T1 * (1/4 - ((0:.01:1) - .5).^2)), 0:.01:1, 'red');
% % plot(-sqrt(T2/T1 * (1/4 - ((0:.01:1) - .5).^2)), 0:.01:1, 'red');
% xlabel('y'); ylabel('z'); axis equal;
% 
% subplot(2,2,4);
% hold off;
% plot(squeeze(y(:,1,:))); 
% xlabel('t (s)'); ylabel('Fisher Information'); 
% drawnow;

end