function fpt = finite_pulse_correction(x,t, TR,Tmax)

TH=x(:,1);
TRF=x(:,2);
t=mod(t,Tmax);
t=reshape(t,[],1);
T=t/TR;
N=length(TH);
I=mod(floor(T),N)+1;
J=round(T);
s=t-J*TR;
J=mod(J-1,N)+1;
TRFJ=TRF(J);
THJ=TH(J);
THJ1=TH(mod(J,N)+1);
fpt=zeros(length(t),2);
trans=(abs(s)<=TRFJ/2);
% fpt(:,2)= trans.*((THJ+THJ1)./TRFJ).^2;
fpt(:,2)= ((THJ+THJ1).^2./TRFJ/TR);
THJ=THJ.*(-1).^(J+1);
THJ1=THJ1.*(-1).^(J);
THI=TH(I).*(-1).^(I+1);
sat=s./TRFJ+1/2;
trans=trans&(mod(J,N)~=0);
fpt(:,1)= trans.*(THJ.*(1-sat)+THJ1.*sat) + (1-trans).*THI;
fpt=fpt.';
fpt(1,:) = abs(fpt(1,:));
end

% function xc = finite_pulse_correction(x, t, TR, Tmax)
% 
% xc = hann_interpolation(t, Tmax, x);
% 
% theta1 = hann_interpolation(floor(t/TR)*TR, Tmax, x(1,:));
% theta2 = hann_interpolation( ceil(t/TR)*TR, Tmax, x(1,:));
% 
% % theta_corr = abs(theta1 - (mod(t,TR) - (TR - xc(2,:))/2)./xc(2,:) .* (theta2 + theta1));
% % xc(1,:) = min(xc(1,:), theta_corr);
% 
% theta_corr = abs(theta1 - (mod(t-TR/2,TR) - (TR - xc(2,:))/2)./xc(2,:) .* (theta2 + theta1));
% xc(1,:) = min(theta1, theta_corr);
% 
% idx = (mod(t, TR) >= TR/2 - xc(2,:)/2) .* (mod(t, TR) <= TR/2 + xc(2,:)/2);
% xc(2,:) = (theta1 + theta2).^2 ./ xc(2,:).^2;
% xc(2,~idx) = 0;
% 
% end