function xc = finite_pulse_correction(x, t, TR, Tmax)

xc = hann_interpolation(t, Tmax, x);

theta1 = hann_interpolation(floor(t/TR)*TR, Tmax, x(1,:));
theta2 = hann_interpolation( ceil(t/TR)*TR, Tmax, x(1,:));
theta_corr = abs(theta1 - (mod(t,TR) - (TR - xc(2,:))/2)./xc(2,:) .* (theta2 + theta1));

xc(1,:) = min(xc(1,:), theta_corr);

idx = (mod(t, TR) >= TR/2 - xc(2,:)/2) .* (mod(t, TR) <= TR/2 + xc(2,:)/2);
xc(2,:) = (theta1 + theta2).^2 ./ xc(2,:).^2;
xc(2,~idx) = 0;

end