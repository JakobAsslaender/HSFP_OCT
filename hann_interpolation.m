function xi = hann_interpolation(t, Tmax, x)
    Nb = length(x)-1;

    x1 = x(:,min(max(floor(t/Tmax*Nb),0),Nb)+1);
    x2 = x(:,min(max(ceil (t/Tmax*Nb),0),Nb)+1);
    xi = (x1 .* (1-cos(pi * (1 - mod(t/Tmax*Nb, 1)))) + ...
          x2 .* (1-cos(pi *      mod(t/Tmax*Nb, 1))))/2;
      
% nearest neighbor interpolation
%     xi = x(min(max(round (t/Tmax*Nb),0),length(x)-1)+1);
end