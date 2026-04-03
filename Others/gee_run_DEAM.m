function [] = gee_deam_analysis( ...
  deam_data_dir, ...
  results_save_dir, ...
  ix_file)

  tools_lib_dir = [ ...
    '/Users', ...
    '/sresarkar', ...
    '/matlabLibrarys'];

  addpath(genpath([tools_lib_dir, '/SparseReg']) );
  addpath(genpath([tools_lib_dir, '/TensorReg']) );
  addpath(genpath([tools_lib_dir, '/tensor_toolbox']) );
  addpath(genpath([tools_lib_dir, '/GEE']) );

  deam_data_dir = char(deam_data_dir);
  results_save_dir = char(results_save_dir);

  assert( ...
    isfolder(deam_data_dir), ...
    ['The provided directory ', deam_data_dir, ' was not found.']);

  if ~isfolder(results_save_dir)
    mkdir(results_save_dir);

    assert(isfolder(results_save_dir) );
  end

  if isstring(ix_file) || ischar(ix_file)
    ix_file = str2double(ix_file);
  end


  tStart = tic;

  fname = sprintf("DEAM_split%d.mat", ix_file);
  deam_data_full_path_fx = fullfile(deam_data_dir, fname);
  assert(isfile(deam_data_full_path_fx), "File not found: %s", deam_data_full_path_fx);

  disp("Loading file:");
  disp(fname);

  deam_data = load(deam_data_full_path_fx);

  lambda_range = [1e-4, 10];
  n_lambdas = 10;

  if (n_lambdas > 1)
    log_lambda_range = log(lambda_range);
    log_lam_step = abs(diff(log_lambda_range) ) / (n_lambdas - 1);
    log_lambdas = log_lambda_range(1):log_lam_step:log_lambda_range(2);
    lambdas = exp(log_lambdas);
  else
    lambdas = lambda_range(1);
  end

  dist = 'normal';
  covar_type = 'equicorr';
  pentype = 'enet';
  penparam = 1;

  t_dim = 2;
  m_dim_train = size(deam_data.X_train, t_dim + 1);
  m_dim_test = size(deam_data.X_test, t_dim + 1);
  P_dim = prod(size(deam_data.X_train, 1:t_dim) );

  X_train = reshape(deam_data.X_train, P_dim, m_dim_train)';
  X_test = reshape(deam_data.X_test, P_dim, m_dim_test)';

  y_train = deam_data.y_train; %.(attribute_of_interest);
  y_test = deam_data.y_test; %.(attribute_of_interest);

  group_labels_train = deam_data.group_labels_train;
  group_labels_test = deam_data.group_labels_test;


  Gtrain = findgroups(group_labels_train);
  time_train = zeros(m_dim_train,1);
  for g = 1:max(Gtrain)
      idx = find(Gtrain == g);
      time_train(idx) = (1:numel(idx))';
  end


  if isrow(time_train)
    time_train = time_train';
  end


  Gtest = findgroups(group_labels_test); % Gtest = 1..n*, one per observation
  time_test = zeros(m_dim_test,1);
  for g = 1:max(Gtest)
      idx = find(Gtest == g);
      time_test(idx) = (1:numel(idx))';
  end

  if isrow(time_test)
    time_test = time_test';
  end

  assert(length(y_train) == m_dim_train);
  assert(length(y_test) == m_dim_test);
  assert(length(group_labels_train) == m_dim_train);
  assert(length(group_labels_test) == m_dim_test);
  assert(length(time_train) == m_dim_train);
  assert(length(time_test) == m_dim_test);

  err_mse = gee_cv( ...
    X_train, ...
    y_train, ...
    group_labels_train, ...
    time_train, ...
    dist, ...
    covar_type, ...
    lambdas, ...
    pentype, ...
    penparam);

  ix_min_mse = max(find(err_mse == min(err_mse) ));
  cv_lambda = lambdas(ix_min_mse);

  [betahat, alphahat, stats] = gee_sparsereg( ...
    findgroups(group_labels_train), ...
    time_train, ...
    X_train, ...
    y_train, ...
    dist, ...
    covar_type, ...
    cv_lambda, ...
    'penalty', pentype, ...
    'penparam', penparam);

  resids = y_train - (X_train * betahat);

  y_preds = gee_equicorr_predict( ...
    resids,  ...
    group_labels_train, ...
    betahat, ...
    alphahat, ...
    X_test, ...
    group_labels_test);

  
  e2 = (y_test - y_preds).^2;
  mse_per_subject = splitapply(@mean, e2, Gtest);      % length n*, each is (1/m_i) sum_j e2_ij
  RMSPE = sqrt(mean(mse_per_subject));                 % (1/n*) sum_i mse_i, then sqrt
  

  runtime_seconds = toc(tStart);
  disp(['Total runtime (seconds): ', num2str(runtime_seconds)]);


  % ---- Construct filename ----
  save_filename = fullfile( ...
      results_save_dir, ...
      sprintf('gee_DEAM_Rep_%d.mat', ix_file) );

  save(save_filename, ...
      'y_preds', ...
      'y_test', ...
      'RMSPE', ...
      'cv_lambda', ...
      'betahat', ...
      'alphahat', ...
      'runtime_seconds');

  disp(['Saved results to: ', save_filename]);

  end