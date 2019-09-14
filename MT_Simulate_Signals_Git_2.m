%%
ijob = round(str2double(getenv('SLURM_ARRAY_TASK_ID')));
disp(ijob)
rng(ijob)
if ijob <= 400
    uniform_sampling = false;
else
    uniform_sampling = true;
end
%%
bDerivatives = true;
bfinite_pulse_correction = true;
load('/scratch/kl3141/MRF_matlab/MRF_CRB/data/OCT_resutl_MT_100_3s_low_SAR.mat', 'x');
TR = 4.5e-3;
t = TR/2:TR:3;
% m0s = .2;
% T1 = 1;
% T2f = .07;
% R = 50;
% T2s = 40e-3;
N = 5;

if uniform_sampling
    m0s = rand(N,1) * .7;
    T1 = rand(N,1) * 2.8 + .2;
    T2f = T1 .* (rand(N,1) * .5 + .05);
    R = 490 * rand(N,1) + 10;
    T2s = .2e-3 + rand(N,1) * 150e-3;
    
else
    m0s = zeros(N,1);
    T1  = zeros(N,1);
    T2f = zeros(N,1);
    R   = zeros(N,1);
    T2s = zeros(N,1);
    
    for n = N:-1:1
        while m0s(n) <= 0
            m0s(n) = 0.1806 + randn() * 0.1189;
        end
        while T1(n) < .1
            T1(n) = 1.7251 + randn() * 1.1767;
        end
        while (T2f(n) < .01 || T2f(n) > T1(n))
            T2f(n) = T1(n) .* (0.0689 + randn() * 0.0897);
        end
        while R(n) < 10
            R(n) = 113.9322 + randn() * 101.7831;
        end
        while T2s(n) < .2e-3
            T2s(n) = 0.0815 + randn() * 0.0428;
        end
    end
end

%%
for n = N:-1:1
    y = simulate_MT_ODE(x, TR, t, m0s(n), T1(n), T2f(n), R(n), T2s(n), bDerivatives, bfinite_pulse_correction);
    s(:,n) = y(:,1);
    
    I = y(:,1:4)'*y(:,1:4);
    % Contains the CRB in this order: PD, m0s, T1, T2f, R, T2s
    CRB_all(:,n) = diag(inv(I));
    
%     progresscounter(N);
end

% Just in case something went wrong...
CRB_all = abs(CRB_all);

%%
if uniform_sampling
    save(['/scratch/kl3141/MRF_matlab/MRF_CRB/Git_results/Uniform_2/ijob_', num2str(ijob)], 'm0s', 'T1', 'T2f', 'R', 'T2s', 's', 'CRB_all');
else
    save(['/scratch/kl3141/MRF_matlab/MRF_CRB/Git_results/Gaussian_2/ijob_', num2str(ijob)], 'm0s', 'T1', 'T2f', 'R', 'T2s', 's', 'CRB_all');
end