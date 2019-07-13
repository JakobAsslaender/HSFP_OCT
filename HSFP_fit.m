function [F, J] = HSFP_fit(param,t,theta,TR,B0,B1,u,model,TRF, TBP)


if nargin<10 || isempty(TBP)
    slice_profile = 1;
elseif TBP == 2
    slice_profile = [1., 0.987207, 0.949667, 0.889817, 0.811461, 0.719419, 0.619089, 0.515981, 0.415261, 0.321371, 0.237735, 0.166607, 0.10903, 0.0649297, 0.0332976, 0.0124423];
    ta = [   0         0         0         0         0         0         0         0         0         0         0         0;...
        0.0014    0.0014    0.0014    0.0014    0.0014    0.0014    0.0013    0.0013    0.0013    0.0013    0.0013    0.0013;...
        0.0225    0.0225    0.0224    0.0222    0.0220    0.0218    0.0215    0.0211    0.0207    0.0203    0.0198    0.0193;...
        0.1034    0.1031    0.1022    0.1006    0.0985    0.0959    0.0927    0.0892    0.0853    0.0812    0.0769    0.0725;...
        0.2715    0.2699    0.2653    0.2577    0.2475    0.2350    0.2207    0.2051    0.1887    0.1722    0.1560    0.1406;...
        0.5048    0.4998    0.4852    0.4618    0.4309    0.3943    0.3540    0.3122    0.2709    0.2323    0.1981    0.1696;...
        0.7366    0.7252    0.6922    0.6402    0.5736    0.4974    0.4174    0.3392    0.2679    0.2083    0.1633    0.1339;...
        0.9013    0.8815    0.8247    0.7372    0.6283    0.5086    0.3889    0.2788    0.1864    0.1186    0.0806    0.0689;...
        0.9791    0.9523    0.8759    0.7607    0.6215    0.4743    0.3336    0.2110    0.1142    0.0467    0.0147    0.0222;...
        0.9988    0.9690    0.8847    0.7593    0.6105    0.4569    0.3146    0.1950    0.1046    0.0440    0.0096    0.0052];
else
    error('Only TBP = 2 implemented');
end

switch model
    case {'HSFP_IR', 'HSFP_anti'}
        if TRF == 350e-6 % 350us rect. pulse
            theta = repmat(theta, [90 1]);
            for ii = 1:3
                theta(ii,:) = theta(ii,:) * (ii-0.5)/3.5;
                theta(end-ii+1,:) = theta(end-ii+1,:) * (ii-0.5)/3.5;
            end
            theta = theta(:);
            TR = TR/90;
            B0 = B0 * 90;
        elseif TRF == 1e-3 % 1000us rect. pulse
            theta = repmat(theta, [45 1]);
            for ii = 1:5
                theta(ii,:) = theta(ii,:) * ii/5;
                theta(end-ii+1,:) = theta(end-ii+1,:) * (ii-1)/5;
            end
            theta = theta(:);
            TR = TR/45;
            B0 = B0 * 45;
        elseif TRF ~= 0 % assume instant pulse
            error('Pulse duration not implemented');
        end
        
        s = 0;
        for is = 1:length(slice_profile)
            switch model
                case 'HSFP_IR'
                    s = s + radial_relaxation_simulator_closed_form(theta, TR, param(3), param(4), -1, B0, B1*slice_profile(is));
                case 'HSFP_anti'
                    s = s + radial_relaxation_simulator_anti_periodic(theta, TR, param(3), param(4), B0, B1*slice_profile(is));
            end                
%             s = s + param(1) * exp(1i*param(2)) * hybrid_state_simulator(theta, TR, param(3), param(4), 'anti', B0, B1*slice_profile(is));
%             s = s + param(1) * exp(1i*param(2)) * hybrid_state_simulator(theta, TR, param(3), param(4), -1, B0, B1*slice_profile(is));
            
        end
        s = s / length(slice_profile);
        
        if TRF == 350e-6 % 350us rect. pulse
            s = s(45:90:end,:);
        elseif TRF == 1e-3 % 1000us rect. pulse
            s = s(23:45:end,:);
        end
        
        J = [s(:,1), 1i*s(:,1), (param(1) + 1i * param(2)) * s(:,2:end)];
        s = (param(1) + 1i * param(2)) * s(:,1);        
    case 'HSFP_SS'
        s = 0;
        TR = TR/180;
        B0 = B0 * 180;
        for is = 1:size(ta,2)
            ts = repmat(theta, [180 1]);
            for ii = 1:9 % 1ms sinc pulse, TBP 2
                ts(ii,:) = ts(ii,:) * ta(ii,is);
                ts(end-ii+1,:) = ts(end-ii+1,:) * ta(ii,is);
            end
            ts(10:end-9,:) = ts(10:end-9,:) * ta(end,is);
            ts = ts(:);
            
            s = s + radial_relaxation_simulator_closed_form(ts, TR, param(3), param(4), -1, B0, B1);
        end
        s = s / size(ta,2);
        
        s = s(90:180:end,:);
        
        J = [s(:,1), 1i*s(:,1), (param(1) + 1i * param(2)) * s(:,2:end)];
        s = (param(1) + 1i * param(2)) * s(:,1);
        
    case 'Bloch'
        alpha = theta(2:end)' + theta(1:end-1)';
        alpha = [theta(end); alpha; pi; -alpha;];
        alpha(2:2:end) = -alpha(2:2:end);
        idx_alpha = [false(length(theta), 1); true(length(theta), 1)];
        
        if TRF == 1e-3 % 1000us rect. pulse
            alpha = repmat(alpha', [10 1]);
            alpha = alpha / 10;
            alpha(:,852) = 0;
            alpha(6,852) = pi;
            alpha = alpha(:);
            
            TRs = TR/45*ones(size(alpha));
            TRs(10:10:end) = TR - 9*TR/45;
        elseif TRF == 0
            TRs = TR*ones(size(alpha));
        else
            error('Pulse duration not implemented');
        end
        
        s = 0;
        for is = 1:length(slice_profile)
            s = (param(1) + 1i * param(2)) * Bloch_simulator_MRF(alpha, TRs, param(3), param(4), -1, B0, B1*slice_profile(is));
        end
        s = s / length(slice_profile);
        
        if TRF == 1e-3 % 1000us rect. pulse
            s = s(10:10:end);
        end
        
        s(2:2:end) = -s(2:2:end);
        s = s(idx_alpha);
        
    case 'LL'
        s = 0;
        for is = 1:length(slice_profile)
            st = (param(1) + 1i * param(2)) * look_locker_simulator(theta * B1 * slice_profile(is), TR, param(3), -1);
            s = s + st(:,1);
        end
        s = s / length(slice_profile);
end

F = u*s;
F = [real(F); imag(F)];

if nargout > 1
    J = u*J;
    J = [real(J); imag(J)];
end
end
