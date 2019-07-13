function [C, gr] = CRB(x, simulator, weights, correlate)

if nargin < 4 || isempty(correlate)
    correlate = 1:length(weights);
end

if nargout > 1
    [y,~,dy,~] = simulator(x);
else
    y = simulator(x);
end

I = y(:,correlate)'*y(:,correlate);

Im1 = inv(I);

% Optimize for the average of all parameters
C = sum(diag(Im1) .* weights(correlate));

% Normalize the cost; the abs is just in case the inversion failed
C = abs(C);



%% Calculate the gradient
if nargout > 1
    Nt = size(dy,2);
    gr = zeros(1,Nt);
    
    % Optimize for the average of all parameters
    for k = 1:Nt
        for m = length(correlate):-1:1
            for l = length(correlate):-1:1
                dIdx(l,m) = 2 * real(dy(:,k,correlate(m))'*y(:,correlate(l)));
            end
        end
        gr(k) = - sum(diag(I \ dIdx / I) .* weights(correlate));
    end
    
end

end