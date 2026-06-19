function y_preds = gee_equicorr_predict( ...
  resids_train, ...
  group_ids_train, ...
  b_hat, ...
  alpha_hat, ...
  X_test, ...
  group_ids_test)

  m_dim_train = length(group_ids_train);
  m_dim_test = length(group_ids_test);

  assert(length(resids_train) == m_dim_train);
  assert(size(X_test, 1) == m_dim_test);

  resids_train = double(resids_train);
  X_test = double(X_test);

  if isrow(resids_train)
    resids_train = resids_train';
  end

  if isrow(group_ids_train)
    group_ids_train = group_ids_train';
  end

  if isrow(group_ids_test)
    group_ids_test = group_ids_test';
  end

  [~, train_order_ixs] = sort(group_ids_train);

  resids_train = resids_train(train_order_ixs,1);
  group_ids_train = group_ids_train(train_order_ixs,1);

  [~, test_order_ixs] = sort(group_ids_test);

  X_test = X_test(test_order_ixs,:);
  group_ids_test = group_ids_test(test_order_ixs,1);

  % group_ids_all = [group_ids_train; group_ids_test];
  % 
  % group_ixs = findgroups(group_ids_all);
  % train_group_ixs = group_ixs(1:m_dim_train,1);
  % test_group_ixs = group_ixs((m_dim_train + 1):end,1);
  % 
  % train_group_ixs_uniq = unique(train_group_ixs);
  % test_group_ixs_uniq = unique(test_group_ixs);
  % 
  % train_group_ixs_no_missing = findgroups(train_group_ixs);
  % train_groups_start_stop_ixs = [ ...
  %   splitapply(@min, (1:m_dim_train)', train_group_ixs_no_missing), ...
  %   splitapply(@max, (1:m_dim_train)', train_group_ixs_no_missing)]';
  % 
  % test_group_ixs_no_missing = findgroups(test_group_ixs);
  % test_groups_start_stop_ixs = [ ...
  %   splitapply(@min, (1:m_dim_test)', test_group_ixs_no_missing), ...
  %   splitapply(@max, (1:m_dim_test)', test_group_ixs_no_missing)]';

  % assume group ids are 1..n and consistent between train/test
  [train_ids_uniq, ~, train_g] = unique(group_ids_train, 'stable');
  [test_ids_uniq,  ~, test_g ] = unique(group_ids_test,  'stable');

  train_groups_start_stop_ixs = [ ...
   accumarray(train_g, (1:m_dim_train)', [], @min), ...
   accumarray(train_g, (1:m_dim_train)', [], @max) ]';

  test_groups_start_stop_ixs = [ ...
    accumarray(test_g, (1:m_dim_test)', [], @min), ...
    accumarray(test_g, (1:m_dim_test)', [], @max) ]';

  y_preds = X_test * b_hat;

  n_test_groups = length(test_ids_uniq);

  for gx = 1:n_test_groups
    %%ix_train_start_stop = find(train_group_ixs_uniq == test_group_ixs_uniq(gx) );
    ix_train_start_stop = find(train_ids_uniq == test_ids_uniq(gx));

    if length(ix_train_start_stop) == 1
      ax_train = train_groups_start_stop_ixs(1,ix_train_start_stop);
      bx_train = train_groups_start_stop_ixs(2,ix_train_start_stop);
      ax_test = test_groups_start_stop_ixs(1,gx);
      bx_test = test_groups_start_stop_ixs(2,gx);

      m_dim_train_gx = bx_train - ax_train + 1;

      pred_adjust_gx = (sum(resids_train(ax_train:bx_train) ) * alpha_hat) / (1 + ((m_dim_train_gx - 1) * alpha_hat) );
      y_preds(ax_test:bx_test) = y_preds(ax_test:bx_test) + pred_adjust_gx;
    end
  end

  [~, orig_test_order_ixs] = sort(test_order_ixs);
  y_preds = y_preds(orig_test_order_ixs,1);
end