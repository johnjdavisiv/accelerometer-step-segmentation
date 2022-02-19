function [Xi, Yi, Zi, Ri, LR_cluster, step_dur, varargout] = RSTAR_segment_steps(x, y, z, r, ...
                                            fs, interp_type, n_interp, ...
                                            segment_location, segment_axis, cluster_axis)
%RSTAR_segment_steps 
%John J. Davis IV


%RSTAR: Rapid Step Assignment for Running
%Despite the name this also works for walking and other cyclical activities.

% Step segmentation for accelerometer data in running

%Applies a modified/improved version of the step segmentation algorithm of
%Urbanek et al. 2017, along with left/right step discrimination via PCA and 
%k-means clustering presented by JJD at ASB/ISB 2019

%Cite: 
% Urbanek, J.K., Harezlak, J., Glynn, N.W., Harris, T., Crainiceanu, C. and Zipunnikov, V., 
% 2017. Stride variability measures derived from wrist-and hip-worn accelerometers. 
% Gait & posture, 52, pp.217-223.

% Davis et al. 2019 Assessing Left-right Asymmetry in Running Using Wearable Accelerometery 
% and Automated Step Segmentation. XXVII Congress of the International Society of Biomechanics, 
% Calgary, Canada, August 2019.


%Input: 

%  x, y, z, and r   -   The x, y, z, and resultant (vector magnitude) axis
%components of an accelerometer signal. The signal is assumed to be
%*continuous* - please apply a running or walking detection algorithm first
%before dumping data into this algorithm. See CARLclassifier for usage
%in the case of detecting running. Ensure that this signal is long enough
%to get meaningful results from an FFT. I recommend at least 8-10 seconds 
%of continuous accelerometer data. Longer is better for slow walking

%  fs  -  sample frequency, in hz

%  interp_type  -  'linear' or 'cubic'.  Dictates whether the interpolated
%  steps are interpolated using linear or cubic spline interpolations.
%  Splines look prettier but might not be a truthful representation of the
%  signal for very fast running with low sample rates.

%  n_interp  -  number of interpolation points for each step (goes from
%  zero to n_interp). Conventionally 101 but be careful with interplating a
%  short step (e.g. fast running - can only be ~30 samples per step)

%  segment_location - 'peak' or 'valley' - indicates where in the gait
%  cycle you want to "chop" steps: the global maxima within the step, or
%  the global minima. Most signal processing papers use peaks, but for
%  biomechanists, valleys are more relevant.

%  segment_axis - 'x', 'y', 'z', or 'r' - indicates which axis to attempt
%  to segment steps on. R should be your default choice for running

%  cluster_axis  -  'x', 'y', 'z', 'all', or 'none'. Dictates which axis is 
%  used to cluster left and right steps. At the waist, 'x' is the best 
%  choice for walking and running but at the wrist 'all' is better because
%  the orientation of the device is unknown. If 'none' is specified, the LR
%  cluster is just a column of 1s

%Output:

%  Xi, Yi, Zi, Ri  -  A matrix of interpolated steps on each axis. The
%  matrix will have n_steps rows and n_interp columns.

%  LR_cluster  -  a vector indicating which cluster each step belongs to.
%  Eventually will be 1 for left and 2 for right. For now, only does
%  arbitrary clustering (need objective data on which step is which!)

%  step_duration  -  a vector of length n_steps that specifies the
%  duration, in seconds, of each step.  

%Optional output

%  x_raw, y_raw, z_raw, r_raw  -  Cell arrays of the raw data from each
%  step. There will be n_steps different cell arrays, and their length will
%  be variable (because the step duration for each is variable)

%  n_steps - integer number of steps detected

%Known limitations:

% Very slow walking (< 0.9 m/s or < 2 mph ) may not be appropriately
% segmented because noise from postural sway pollutes the frequency
% domain so much that it is hard to resolve the dominant frequency. In
% faster walking and especially running, the dominant FFT peak always
% corresponds to cadence so the algorithm works quite well. More
% sophisticated methods may be needed in this case to define the
% appropriate filter parameters. Sometimes, using a longer sample (30sec or
% more) of continuous walking can help resolve the dominant FFT peak better


%% segmentation
switch lower(segment_axis)
    case 'x'
        detrend_v = detrend(x, 'constant');
    case 'y'
        detrend_v = detrend(y, 'constant');
    case 'z'
        detrend_v = detrend(z, 'constant');
    case 'r'
        detrend_v = detrend(r, 'constant');
    otherwise
        error('''%s'' is not a valid segmentation axis. Specify ''x'', ''y'', ''z'', or ''r''.');
end
        
%Uncomment to view frequency spectrum
%plot(freqHz, myFFTmag);

lp_cutoff = 5; %Lowpass filter cutoff - 5 Hz works just fine...
norm_cutoff_freq = lp_cutoff/(fs/2); %set normalized cutoff frequency
filter_order = 4; %fourth-order

%Perform lowpass filter
[bee,ayy] = butter(filter_order,norm_cutoff_freq);  
filt_v = filtfilt(bee,ayy, detrend_v);
% 
% plot(r); hold on; plot(filt_v);

%% ----------------   Find the peaks and chop into steps    ----------

%use findpeaks() to find the valleys (multiply signal by -1)
%then "grab" these and normalize them, storing them in a matrix

if strcmpi(segment_location, 'valley')
    [~, locs] = findpeaks(-1*filt_v);
elseif strcmpi(segment_location, 'peak')
    [~, locs] = findpeaks(filt_v);
else 
    error('''%s'' is not a valid segmentation location. Specify ''peak'' or ''valley''.', segment_location);
end

%We'll always get one fewer full step than number of valleys
n_steps = length(locs) -1;

Xi = zeros(n_steps, n_interp);
Yi = zeros(n_steps, n_interp);
Zi = zeros(n_steps, n_interp);
Ri = zeros(n_steps, n_interp);

x_raw = cell(n_steps,1);
y_raw = cell(n_steps,1);
z_raw = cell(n_steps,1);
r_raw = cell(n_steps,1);

step_dur = zeros(n_steps,1);

%For each step...
for a=1:n_steps

    %Select the data, then normalize it 0-100, then store it
    working_r = r(locs(a):(locs(a+1)-1));
    
    %Need to be a "half open interval"  [hence the minus one]
    %otherwise step time will always be 100th of a second too long.
    working_x = x(locs(a):(locs(a+1)-1));
    working_y = y(locs(a):(locs(a+1)-1));
    working_z = z(locs(a):(locs(a+1)-1));
    
    %Save step duration
    step_dur(a) = length(working_r)/fs; %step duration, in seconds
    
    %Save raw step data
    x_raw{a} = working_x;
    y_raw{a} = working_y;
    z_raw{a} = working_z;
    r_raw{a} = working_r;
    
    t_raw = 1:length(working_r);
    t_interp = linspace(1,length(working_r), n_interp);
    
    switch lower(interp_type)
        case 'linear'
            %Interpolating using linear interpolation
            interp_r = interp1(t_raw, working_r, t_interp);
            interp_x = interp1(t_raw, working_x, t_interp);
            interp_y = interp1(t_raw, working_y, t_interp);
            interp_z = interp1(t_raw, working_z, t_interp);           

        case 'cubic'
            %Interpolating using cubic splines
            interp_r = spline(t_raw, working_r, t_interp);
            interp_x = spline(t_raw, working_x, t_interp);
            interp_y = spline(t_raw, working_y, t_interp);
            interp_z = spline(t_raw, working_z, t_interp);
        otherwise
            error('%s is not a valid interpolation method. Use ''linear'' or ''cubic''.', ...
                interp_type);
    end
    
    %Save interpolated steps
    Ri(a,:) = interp_r;
    Xi(a,:) = interp_x;
    Yi(a,:) = interp_y;
    Zi(a,:) = interp_z;
        
end

if nargout > 6
    varargout{1} = x_raw;
    varargout{2} = y_raw;
    varargout{3} = z_raw;
    varargout{4} = r_raw;
    varargout{5} = n_steps;
end


%% ----------     PCA and k-means clustering     ------
%Separates left and right steps

switch lower(cluster_axis)
    case 'x'
        M_pca =  Xi;
    case 'y'
        M_pca = Yi;
    case 'z'
        M_pca = Zi;
    case 'all'
        M_pca = [Xi, Yi, Zi];
    case 'none'
        LR_cluster = ones(size(Xi,1),1);
        return  %Just skip clustering if 'none' is specified
        %Slightly bad programming practice but saves computation time
    otherwise
        error('''%s'' is not a valid axis for clustering. Use ''x'', ''y'', ''z'', ''all'', or ''none''.',...
            cluster_axis);
end

%PCA, keep scores - uses SVD by default.
[~,score,~] = pca(M_pca);

%uncomment to view PCA clusters if desired
%plot(score(:,1), score(:,2), 'k.', 'MarkerSize', 12);

%Perform k-means clustering - defaults to 3 PCs, looks for two clusters
%Defaults to 3 PCs and 5 replicates
LR_cluster = kmeans([score(:,1), score(:,2), score(:,3)], 2, 'Replicates', 10);
%Franti and Sieranoja 2019 say that # of replicates in PCA is not such a
%big deal when you have really well separated clusters, which is what we
%will have here.

%Later I will find some rule to determine which is left and right...

end

