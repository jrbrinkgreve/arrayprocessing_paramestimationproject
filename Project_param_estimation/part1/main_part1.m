%main file:
%array and measurement settings
close all

M = 3; %antennas
N = 20; %time samples
Delta = 0.5; %array spacing in wavelengths
m = 5; %smoothing degree
num_loops = 1000; %number of loops for statistical averaging
SNRS = [0 4 8 12 16 20];


%signal specifications
d = 2;
theta = [-20 30] / 180 * pi; %need to pass it to the array response in radians
f = [0.10 0.12];



 %preallocate
    theta_f_storage = zeros(num_loops, 2*d, length(SNRS));
    theta_f_storage_joint = zeros(num_loops, 2*d, length(SNRS));
    


%loop over snrs
parfor j = 1:length(SNRS) %run (in parallel) for every snr:
    
    SNR = SNRS(j);
    
    
    %generate new data and estimate
    for i = 1:num_loops
        
        %generate data
        [X, A, S] = gendata(M,N,Delta,theta,f,SNR);
        
        
        theta_estimate = sort(esprit(X,d, Delta)) ;   %sort to prevent arbitrary ordering, but this is a problem in practice!!!
        f_estimate = sort(espritfreq(X, d));         %of course not in general the lowest angles and
                                                      %frequencies
                                                      %correspond. the
                                                      %whole point for jd
        
        [theta_estimate_joint, f_estimate_joint] = joint_esprit_smoothed(X, d, Delta, m);
    


        %store
        theta_f_storage(i, :, j) = [theta_estimate.', f_estimate.'];
        theta_f_storage_joint(i,:,j) = [theta_estimate_joint.', f_estimate_joint.'];
     
    end
end




%noiseless test:
[X, A, S] = gendata(M,N,Delta,theta,f,1e9);  
theta_estimate = sort(esprit(X,d, Delta)) ;   %sort to prevent arbitrary ordering, but this is a problem in practice!!!
f_estimate = sort(espritfreq(X, d));         %of course not in general the lowest angles and
                                                      %frequencies
                                                      %correspond. the
                                                      %whole point for jd

%compute beamformers:
[spatial_beamformer, temporal_beamformer] =  calc_beamformers(M, Delta, N, theta_estimate, f_estimate);

S_est = spatial_beamformer' * X;
A_est = X * temporal_beamformer';

%print errors
A_zf_error  = norm(A - A_est, 'fro') / numel(A)
S_zf_error = norm(S - S_est, 'fro') / numel(S)





npoints = 1000;
angles = linspace(-pi/2, pi/2, npoints);   % radians
freqs = linspace(0, 1, npoints);          % normalized frequency

%testdata
[X, A, S] = gendata(M, N, Delta, angles, freqs, 10);  %to scan over all angles and freqs with 1000 npoints points



%need to figure out this part below still AAAAAAAAA
%{
theta_estimate = sort(esprit(X,d, Delta)) ;   %sort to prevent arbitrary ordering, but this is a problem in practice!!!
f_estimate = sort(espritfreq(X, d));         %of course not in general the lowest angles and
                                                      %frequencies
                                                      %correspond. the
                                                      %whole point for jd

%compute beamformers at 10 snr now:
[spatial_beamformer, temporal_beamformer] =  calc_beamformers(M, Delta, N, theta_estimate, f_estimate);

S_est = spatial_beamformer' * X;
A_est = X * temporal_beamformer';
%}




%get responses
angular_response  = abs(spatial_beamformer' * A);         % d x npoints
temporal_response = abs(S * temporal_beamformer');        % d x npoints

%in dB
angular_response_db  = 20 * log10(angular_response ./ max(angular_response, [], 2));
temporal_response_db = 20 * log10(temporal_response ./ max(temporal_response, [], 2));

%angular response plot
figure;
plot(angles * 180 / pi, angular_response_db.', 'LineWidth', 1.2);
xlabel('Angle (degrees)', 'FontSize', 12);
ylabel('Gain (dB)', 'FontSize', 12);
title('Angular Response of Beamformers', 'FontSize', 14);
grid on;
%ylim([-60 5]);
xlim([-90 90]);
legend(arrayfun(@(i) sprintf('Source %d', i), 1:size(spatial_beamformer, 2), 'UniformOutput', false), ...
       'Location', 'best');
set(gca, 'FontSize', 20);

%freq response plot
figure;
plot(freqs, temporal_response_db.', 'LineWidth', 1.2);
xlabel('Normalized Frequency (cycles/sample)', 'FontSize', 12);
ylabel('Gain (dB)', 'FontSize', 12);
title('Frequency Response of Beamformers', 'FontSize', 14);
grid on;
ylim([-60 5]);
xlim([0 1]);
legend(arrayfun(@(i) sprintf('Source %d', i), 1:size(temporal_beamformer, 2), 'UniformOutput', false), ...
       'Location', 'best');
set(gca, 'FontSize', 20);






%postprocessing

%preallocate result arrays
mean_errors_esprit = zeros(length(SNRS), 2*d); 
var_errors_esprit = zeros(length(SNRS), 2*d);  

mean_errors_joint = zeros(length(SNRS), 2*d);
var_errors_joint = zeros(length(SNRS), 2*d);

%estimate mean / var
parfor j = 1:length(SNRS)
    data_esprit = theta_f_storage(:, :, j);
    data_joint = theta_f_storage_joint(:, :, j);
    
    %mean
    mean_errors_esprit(j, :) = mean(data_esprit, 1);
    mean_errors_joint(j, :) = mean(data_joint, 1);
    
    %var
    var_errors_esprit(j, :) = var(data_esprit, 0, 1);  %0 means normalize by N-1, unbiased estimator
    var_errors_joint(j, :) = var(data_joint, 0, 1);
end





% Plot 1: Variance of ANGLE estimates
figure;
plot(SNRS, var_errors_esprit(:,1:d), '-o', 'LineWidth', 1.5); hold on;
plot(SNRS, var_errors_joint(:,1:d), '--s', 'LineWidth', 1.5);
title('Variance of Angle Estimates (θ) vs SNR');
xlabel('SNR (dB)');
ylabel('Variance (degrees²)');
legend([compose("ESPRIT Var(θ_%d)", 1:d), compose("Joint ESPRIT Var(θ_%d)", 1:d)], 'Location', 'northeast');
set(gca, 'FontSize', 20);
grid on;

% Plot 2: Variance of FREQUENCY estimates
figure;
plot(SNRS, var_errors_esprit(:,d+1:end), '-o', 'LineWidth', 1.5); hold on;
plot(SNRS, var_errors_joint(:,d+1:end), '--s', 'LineWidth', 1.5);
title('Variance of Frequency Estimates (f) vs SNR');
xlabel('SNR (dB)');
ylabel('Variance');
legend([compose("ESPRIT Var(f_%d)", 1:d), compose("Joint ESPRIT Var(f_%d)", 1:d)], 'Location', 'northeast');
set(gca, 'FontSize', 20);
grid on;

% Plot 3: Mean of ANGLE estimates
figure;
plot(SNRS, mean_errors_esprit(:,1:d), '-o', 'LineWidth', 1.5); hold on;
plot(SNRS, mean_errors_joint(:,1:d), '--s', 'LineWidth', 1.5);
title('Mean of Angle Estimates (θ) vs SNR');
xlabel('SNR (dB)');
ylabel('Mean (degrees)');
legend([compose("ESPRIT Mean(θ_%d)", 1:d), compose("Joint ESPRIT Mean(θ_%d)", 1:d)], 'Location', 'northeast');
set(gca, 'FontSize', 20);
grid on;

% Plot 4: Mean of FREQUENCY estimates
figure;
plot(SNRS, mean_errors_esprit(:,d+1:end), '-o', 'LineWidth', 1.5); hold on;
plot(SNRS, mean_errors_joint(:,d+1:end), '--s', 'LineWidth', 1.5);
title('Mean of Frequency Estimates (f) vs SNR');
xlabel('SNR (dB)');
ylabel('Mean (normalized frequency)');
legend([compose("ESPRIT Mean(f_%d)", 1:d), compose("Joint ESPRIT Mean(f_%d)", 1:d)], 'Location', 'northeast');
set(gca, 'FontSize', 20);
grid on;


%note: plotting scripts with some help of an LLM




%computing zero-forcing beamformers:








%analysis:

%--------------
