% Stefan Haufe, 2019
%
% required toolboxes: 
% MVGC: https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html
% export_fig: https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig

%% load data and choose settings

% initialize and save random number generator seed
rng('shuffle');
s = rng;
save('results/rng_seed', 's');

% % load random number generator seed
% load('results/rng_seed')
% rng(s)

% load New York Head model for a reduced set of 57 EEG channels
load sa_nyhead_57channel

% colormap
load cm17

% 3D leadfield for a subset of 1K cortical locations
L_3D = sa.cortex75K.V_fem(:, sa.cortex1K.in_from_cortex75K, :);

% 1D leadfield at the same locations, orientation is perpendicular to
% cortical surface
L_normal = sa.cortex75K.V_fem_normal(:, sa.cortex1K.in_from_cortex75K);

% number of EEG channels
M = length(sa.clab_electrodes);

% number of epochs
Nepo = 200;

% length of epochs in samples
Lepo = 100;

% overall number of samples
T = Nepo*Lepo;

% true lag in samples, if directed interaction is modeled as delay between
% band-limited time series
lag = 10;

% True AR model order, if directed interaction is modeled as bivariate AR
% model
morder_true = 5;

% use same order in estimation, choice of this parameter not being the
% point here
morder_est = morder_true;

% number of sources
% only the first two are interacting, the remaining ones are all mutually
% independent
Nsources = 20;

% choose whether interaction should be modeled as AR process or through a
% simple delay
interact_type = 'ar'; % 'lag' or 'ar'

% number of repetitions of the experiment
nrep = 10;

%% 
for irep = 1:nrep

  s = rng;
  save(['results/rep' num2str(irep) '_rng_seed'], 's');

  %% generate the source time series

  if isequal(interact_type, 'lag');
    % model interaction through lag

    % construct filtered noise
    sources = randn(T+lag, Nsources);
    [b a] = butter(2, [0.1 0.3]);
    sources = filtfilt(b, a, sources);

    % second time series is just lagged version of first one
    sources(:, 2) = circshift(sources(:, 1), lag);
    sources = sources(lag+1:end, :);
    sources = sources';
  else % ('ar')
    % model interaction betwenn first and second source through bivariate AR model
    sources = zeros(Nsources, T);
    sources(1:2, :) = gen_ar_biv(T, morder_true);

    % remaining sources are independent, modeled as univariate AR
    for isource = 3:Nsources
      sources(isource, :) = gen_ar_uni(T, morder_true);
    end

    % apply highpass to suppress fluctuations slower than epoch length
    [b a] = butter(2, 0.02, 'high');
    sources = filtfilt(b, a, sources')';
  end

  % unify scale of all sources
  sources = zscore(sources')';

  % calculate connectivity on ground-truth data between first source (sender) and all others using Granger
  % causality (GC) and time-reversed Granger causality (TRGC)
  % Here, the MVG toolbox (Barnett/Seth) is used. We follow the procedure 
  % time series -> VAR -> AUTOCOV -> GC recommended by Barnett/Seth 
  TRGC_true_sender = [];
  GC_true_sender = [];
  for isource = 2:Nsources
    % convert time series into VAR representation
    [A, SIG, E] = tsdata_to_var(reshape(sources([1 isource], :), 2, Lepo, Nepo), morder_est);

    % convert VAR into autocov sequence
    Gorig = var_to_autocov(A, SIG, 100);

    % convert autocov sequence into GC scores
    GC_true_sender(isource, 1) = autocov_to_mvgc(Gorig, 2, 1) - autocov_to_mvgc(Gorig, 1, 2);

    % transpose autocov sequence to obtain autocov seq. of time-reversed data
    Grev = permute(Gorig, [2 1 3]);

    % compute GC on time-reversed data
    TRGC_true_sender(isource, 1) = autocov_to_mvgc(Grev, 1, 2) - autocov_to_mvgc(Grev, 2, 1);
  end
  % TRGC is simply the difference between GC on original and time-reversed
  % data
  TRGC_true_sender = GC_true_sender - TRGC_true_sender;

  % calculate connectivity on ground-truth data  between second source (receiver) and all others
  TRGC_true_receiver = [];
  GC_true_receiver = [];
  for isource = [1 3:Nsources]
    [A, SIG, E] = tsdata_to_var(reshape(sources([2 isource], :), 2, Lepo, Nepo), morder_est);
    Gorig = var_to_autocov(A, SIG, 100);
    GC_true_receiver(isource, 1) = autocov_to_mvgc(Gorig, 2, 1) - autocov_to_mvgc(Gorig, 1, 2);

    Grev = permute(Gorig, [2 1 3]);
    TRGC_true_receiver(isource, 1) = autocov_to_mvgc(Grev, 1, 2) - autocov_to_mvgc(Grev, 2, 1);
  end
  TRGC_true_receiver = GC_true_receiver - TRGC_true_receiver;

  % plot net outflow for both measures for the original source time series
  % this looks pretty good for both: source one sends to source two, and
  % source two receives from source one, nothing else is seen
  figure; plot([GC_true_sender TRGC_true_sender]); legend('GC', 'TRGC'); grid on
  figure; plot([GC_true_receiver TRGC_true_receiver]); legend('GC', 'TRGC'); grid on

  %% generate pseudo sensor data 

  % number of voxels in source space (1K)
  Nvox = length(sa.cortex1K.in_from_cortex75K);

  % choose locations for voxels. Here, the first source is constrained to be
  % in the left hemisphere, and the second source is constrained to be in the
  % right hemisphere. The other (independent brain noise) sources can be
  % anywhere
  ind_sources = randi(Nvox, Nsources, 1);
  ind_sources(1) = randi(150, 1, 1);
  ind_sources(2) = 503+randi(150, 1, 1);

  % plot the ground truth: the sender (cyan ball) sends to the receiver
  % (magenta ball). The remaining sources (white balls) do not participate in any
  % interaction.
  allplots_cortex_nyhead_leftright(sa.cortex75K, 0*sa.cortex1K.in_to_cortex75K_geod, [0 0], cm17, '', 1, ['figures/rep' num2str(irep) '/ground_truth'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % plot TRGC outflow from sending source location (blue ball) to all other 
  % voxels. The correct result is a red patch around the receiver
  % (magenta ball) and zero connectivity everywhere else including the
  % vicinity of the blue ball. This is achieved by both GC and TRGC
  pl = zeros(Nvox, 1);
  pl(ind_sources) = TRGC_true_sender;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/ground_truth_TRGC_sender'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % same for GC
  pl = zeros(Nvox, 1);
  pl(ind_sources) = GC_true_sender;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/ground_truth_GC_sender'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});


  % plot TRGC outflow (negative inflow) from receiving source location (magenta ball) to all other 
  % voxels. The correct result is a blue patch (negative outflow =
  % inflow) around the sender (blue ball) and zero connectivity everywhere else including the
  % vicinity of the magenta ball. 
  pl = zeros(Nvox, 1);
  pl(ind_sources) = TRGC_true_receiver;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/ground_truth_TRGC_receiver'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % same for GC
  pl = zeros(Nvox, 1);
  pl(ind_sources) = GC_true_receiver;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/ground_truth_GC_receiver'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});



  % project sources to EEG sensors using 1D leadfield (assuming perpendicular
  % orientation.
  sensors = L_normal(:, ind_sources)*sources;
  sensors = sensors ./ norm(sensors, 'fro');

  % generate some white sensor noise
  noise = randn(size(sensors));
  noise = noise ./ norm(noise, 'fro');

  % add a small amount (5%) of sensor noise. This makes the data
  % full rank, which is useful for beamformer (lcmv) analysis
  eeg = sensors + 0.05*noise;

  %% eloreta source analysis

  % calculate 3D eLORETA inverse filter with a small amount of regularization
  P_eloreta = mkfilt_eloreta_v2(L_3D, 0.05);

  % apply eLORETA filter to obtain source estimates 
  est_sources_eloreta = reshape(P_eloreta(:, :)'*eeg, Nvox, 3, Nepo*Lepo);
  % est_sources_eloreta = tprod(eeg, [-1 3], P_eloreta, [-1 1 2]);

  % plot voxel-wise source power
  po_eloreta = sum(est_sources_eloreta(:, :).^2, 2);
  pl = po_eloreta;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [0 max(abs(pl))], cm17a, '', 1, ['figures/rep' num2str(irep) '/eloreta_power'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});


  % calculate directed connectivity on source-reconstruted data between sending voxel and all other
  % voxels. To this end, GC and TRGC are calculated between 3D spaces
  TRGC_est_sender_eloreta = [];
  GC_est_sender_eloreta = [];
  disp('Loop over voxels. This will take ~3 minutes')
  for ivox = setdiff(1:Nvox, ind_sources(1))
  %   ivox
    est_sources_ = reshape(permute(reshape(est_sources_eloreta([ind_sources(1) ivox], :, :), 2, 3, Lepo, Nepo), [2 1 3 4]), 3*2, Lepo, Nepo);

    [A, SIG, E] = tsdata_to_var(est_sources_, morder_est);
    Gorig = var_to_autocov(A, SIG, 100);
    GC_est_sender_eloreta(ivox, 1) = autocov_to_mvgc(Gorig, [4 5 6], [1 2 3]) - autocov_to_mvgc(Gorig, [1 2 3], [4 5 6]);

    Grev = permute(Gorig, [2 1 3]);
    TRGC_est_sender_eloreta(ivox, 1) = autocov_to_mvgc(Grev, [4 5 6], [1 2 3]) - autocov_to_mvgc(Grev, [1 2 3], [4 5 6]);
  end
  TRGC_est_sender_eloreta = GC_est_sender_eloreta - TRGC_est_sender_eloreta;

  % plot TRGC outflow from sending source location (blue ball) to all other 
  % voxels. The correct result would be a red patch around the receiver
  % (magenta ball) and zero connectivity everywhere else including the
  % vicinity of the blue ball. This is usually achieved for TRGC.
  pl = TRGC_est_sender_eloreta;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/eloreta_TRGC_sender'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % plot the same for GC. Here, additional spurious connectivity patterns are
  % frequentlyl observed, for example in the vicinity of the sender (cyan ball). 
  pl = GC_est_sender_eloreta;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/eloreta_GC_sender'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % calculate directed connectivity on source-reconstruted data between receiving voxel and all other
  % voxels. 
  TRGC_est_receiver_eloreta = [];
  GC_est_receiver_eloreta = [];
  disp('Loop over voxels. This will take ~3 minutes')
  for ivox = setdiff(1:Nvox, ind_sources(2))
  %   ivox
    est_sources_ = reshape(permute(reshape(est_sources_eloreta([ind_sources(2) ivox], :, :), 2, 3, Lepo, Nepo), [2 1 3 4]), 3*2, Lepo, Nepo);

    [A, SIG, E] = tsdata_to_var(est_sources_, morder_est);
    Gorig = var_to_autocov(A, SIG, 100);
    GC_est_receiver_eloreta(ivox, 1) = autocov_to_mvgc(Gorig, [4 5 6], [1 2 3]) - autocov_to_mvgc(Gorig, [1 2 3], [4 5 6]);

    Grev = permute(Gorig, [2 1 3]);
    TRGC_est_receiver_eloreta(ivox, 1) = autocov_to_mvgc(Grev, [4 5 6], [1 2 3]) - autocov_to_mvgc(Grev, [1 2 3], [4 5 6]);
  end
  TRGC_est_receiver_eloreta = GC_est_receiver_eloreta - TRGC_est_receiver_eloreta;


  % plot TRGC outflow (negative inflow) from receiving source location (magenta ball) to all other 
  % voxels. The correct result would be a blue patch (negative outflow =
  % inflow) around the sender (blue ball) and zero connectivity everywhere else including the
  % vicinity of the magenta ball. 
  pl = TRGC_est_receiver_eloreta;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/eloreta_TRGC_receiver'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % plot the same for GC.
  pl = GC_est_receiver_eloreta;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/eloreta_GC_receiver'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  %% LCMV source analysis 

  [P_lcmv, ~, po_lcmv] = mkfilt_lcmv(L_3D, cov(eeg'));
  [Pn, ~, pon] = mkfilt_lcmv(L_3D, eye(M));
  est_sources_lcmv = tprod(eeg, [-1 3], P_lcmv, [-1 1 2]);

  % plot source power 
  po_lcmv = po_lcmv ./ pon;
  pl = po_lcmv;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [0 max(abs(pl))], cm17a, '', 1, ['figures/rep' num2str(irep) '/lcmv_power'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});


  % calculate directed connectivity on source-reconstruted data between sending voxel and all other
  % voxels. To this end, GC and TRGC are calculated between 3D spaces
  TRGC_est_sender_lcmv = [];
  GC_est_sender_lcmv = [];
  disp('Loop over voxels. This will take ~3 minutes')
  for ivox = setdiff(1:Nvox, ind_sources(1))
  %   ivox
    est_sources_ = reshape(permute(reshape(est_sources_lcmv([ind_sources(1) ivox], :, :), 2, 3, Lepo, Nepo), [2 1 3 4]), 3*2, Lepo, Nepo);

    [A, SIG, E] = tsdata_to_var(est_sources_, morder_est);
    Gorig = var_to_autocov(A, SIG, 100);
    GC_est_sender_lcmv(ivox, 1) = autocov_to_mvgc(Gorig, [4 5 6], [1 2 3]) - autocov_to_mvgc(Gorig, [1 2 3], [4 5 6]);

    Grev = permute(Gorig, [2 1 3]);
    TRGC_est_sender_lcmv(ivox, 1) = autocov_to_mvgc(Grev, [4 5 6], [1 2 3]) - autocov_to_mvgc(Grev, [1 2 3], [4 5 6]);
  end
  TRGC_est_sender_lcmv = GC_est_sender_lcmv - TRGC_est_sender_lcmv;

  % plot TRGC outflow from sending source location (blue ball) to all other 
  % voxels. The correct result would be a red patch around the receiver
  % (magenta ball) and zero connectivity everywhere else including the
  % vicinity of the blue ball. This is usually achieved for TRGC.
  pl = TRGC_est_sender_lcmv;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/lcmv_TRGC_sender'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % plot the same for GC. Here, additional spurious connectivity patterns are
  % frequentlyl observed, for example in the vicinity of the sender (cyan ball). 
  pl = GC_est_sender_lcmv;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/lcmv_GC_sender'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % calculate directed connectivity on source-reconstruted data between receiving voxel and all other
  % voxels. 
  TRGC_est_receiver_lcmv = [];
  GC_est_receiver_lcmv = [];
  disp('Loop over voxels. This will take ~3 minutes')
  for ivox = setdiff(1:Nvox, ind_sources(2))
  %   ivox
    est_sources_ = reshape(permute(reshape(est_sources_lcmv([ind_sources(2) ivox], :, :), 2, 3, Lepo, Nepo), [2 1 3 4]), 3*2, Lepo, Nepo);

    [A, SIG, E] = tsdata_to_var(est_sources_, morder_est);
    Gorig = var_to_autocov(A, SIG, 100);
    GC_est_receiver_lcmv(ivox, 1) = autocov_to_mvgc(Gorig, [4 5 6], [1 2 3]) - autocov_to_mvgc(Gorig, [1 2 3], [4 5 6]);

    Grev = permute(Gorig, [2 1 3]);
    TRGC_est_receiver_lcmv(ivox, 1) = autocov_to_mvgc(Grev, [4 5 6], [1 2 3]) - autocov_to_mvgc(Grev, [1 2 3], [4 5 6]);
  end
  TRGC_est_receiver_lcmv = GC_est_receiver_lcmv - TRGC_est_receiver_lcmv;


  % plot TRGC outflow (negative inflow) from receiving source location (magenta ball) to all other 
  % voxels. The correct result would be a blue patch (negative outflow =
  % inflow) around the sender (blue ball) and zero connectivity everywhere else including the
  % vicinity of the magenta ball. 
  pl = TRGC_est_receiver_lcmv;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/lcmv_TRGC_receiver'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  % plot the same for GC.
  pl = GC_est_receiver_lcmv;
  allplots_cortex_nyhead_leftright(sa.cortex75K, pl(sa.cortex1K.in_to_cortex75K_geod), [-max(abs(pl)) max(abs(pl))], cm17, '', 1, ['figures/rep' num2str(irep) '/lcmv_GCreceiver'], ...
    {sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(1)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(2)), :), ...
    sa.cortex75K.vc_smooth(sa.cortex1K.in_from_cortex75K(ind_sources(3:end)), :)});

  close all
end
