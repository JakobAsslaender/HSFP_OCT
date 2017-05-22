function [C, gr] = CRB(theta, simulator, weights, correlate)

if nargin < 4 || isempty(correlate)
    correlate = ones(size(weights));
end

C = 0;
if nargout > 1
    [y,~,dy,~] = simulator(theta);
else
    y = simulator(theta);
end

Nspin  = size(y,3);
Npulse = length(theta);

for n=1:Nspin
    I = y(:,correlate,n)'*y(:,correlate,n);
    Im1 = inv(I);
    
    % Optimize for the average of all parameters and all spins
    C = C + sum(diag(Im1) .* weights(correlate,n));
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
            gr(k) = gr(k) - sum(diag(Im1 * dIda * Im1) .* weights(correlate,n));
        end
    end
    gr = gr /(Nspin^2);
end


%% Plot result
% figure(1); 
% subplot(2,2,2); hold off;
% plot(theta/pi);
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

%% export dynamics
% global idx TR T1 T2 
% if idx < 1000
%     ID = fopen(['~/Documents/Output/Talks/2017_04_06_radial_relaxation/Figures/OCT_dynamics_T1_pi2_iter_', num2str(idx), '.txt'], 'w');
%     fprintf(ID, 't_s theta z y dydT1 dydT2 \n');
%     for itheta = 1:length(theta)
%         fprintf(ID, '%f %f %f %f %f %f \n', itheta*TR, theta(itheta)/pi, z(itheta), y(itheta,1), -y(itheta,2) * T1, y(itheta,3) * T2);
%     end
%     fclose(ID);
% end
% idx = idx + 1;

end