%Demo step segmentation
%JJD 2022-02-19
tic;
run_data = load('sample_80min_uninterrupted_run_01_14_2019.mat');
run_data = run_data.myRun;
%Load actigraph data from an uninterrupted run
% formatted as matrix with x-y-z-r axes as columns in that order

n_demo = round(max(size(run_data))*0.2); %Only demo first ~16min of run, for speed
%Works fine on the whole run though

x = run_data(1:n_demo,1);
y = run_data(1:n_demo,2);
z = run_data(1:n_demo,3);
r = run_data(1:n_demo,4);


%Params for step segmentation
fs = 100; %Hz
interp_type = 'linear'; %or cubic
n_interp = 101; %0-100% of step
segment_location = 'valley'; %to grab in flight phase vs. max GRF
segment_axis = 'r'; 
cluster_axis = 'all'; %All works best if device orientation is unknown

%Step segmentation happens in this function
[Xi, Yi, Zi, Ri, LR_cluster, step_dur] = RSTAR_segment_steps(x, y, z, r, ...
                                            fs, interp_type, n_interp, ...
                                            segment_location, segment_axis, cluster_axis);

%% Plot all steps and average step

figure('units', 'normalized', 'position', [0.1 0.1 0.5 0.6]);

alf = 0.05;

%Mind the transpose
subplot(2,2,1);
hold on;
plot(Xi(LR_cluster == 1,:)', 'color', [0 0 0.5 alf]);
plot(Xi(LR_cluster == 2,:)', 'color', [0 0.3 0.9 alf]);

plot(mean(Xi(LR_cluster == 1,:)), 'color', [1 1 1], 'LineWidth',2);
plot(mean(Xi(LR_cluster == 2,:)), 'color', [1 1 1], 'LineWidth',2);
xlim([0 100]);
title('X axis')


subplot(2,2,2);
hold on;
plot(Yi', 'color', [0.5 0 0 0.01]);
plot(mean(Yi),'color', [1 1 1], 'LineWidth',2);
xlim([0 100]);
title('Y axis')


subplot(2,2,3);
hold on;
plot(Zi', 'color', [0.15 0.8 0.15 alf]);
plot(mean(Zi),'color', [1 1 1], 'LineWidth',2);
xlim([0 100]);
title('Z axis')


subplot(2,2,4);
hold on;
plot(Ri', 'color', [0.9 0 0.5 alf]);
plot(mean(Ri),'color', [1 1 1], 'LineWidth',2);
xlim([0 100]);
title('Resultant');

toc;