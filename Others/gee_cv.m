function [err_mse, err_se] = gee_cv( ...
  X, y, ...
  group_ids, ...
  time, ...
  dist, ...
  covar_type, ...
  lambdas, ...
  pentype, ...
  penparam)
%{
DESCRIPTION:
  Get the 10-fold cross validation MSE's at the input lambda(s).

PARAMETERS:
            X:
            y:
    group_ids:
         time:
         dist:
   covar_type:
      lambdas:
      pentype:
     penparam:

RETURN VALUE:
  10-fold cross validation MSE's corresponding to the provided lambda value(s).
%}

  assert(length(size(X) ) == 2);

  m_dim = size(X, 1);

  assert(length(y) == m_dim);
  assert(length(group_ids) == m_dim);
  assert(length(time) == m_dim);

  X = double(X);
  y = double(y);

  if isrow(y)
    y = y';
  end

  if isrow(group_ids)
    group_ids = group_ids';
  end

  if isrow(time)
    time = time';
  end

  n_lams = length(lambdas);
  n_cv_grps = 10;

  %[~, time_order_ixs] = sort(time);
  %[~, group_order_ixs] = sort(group_ids(time_order_ixs) );
  %order_ixs = time_order_ixs(group_order_ixs);

  %X = X(order_ixs,:);
  %y = y(order_ixs,1);
  %group_ids = group_ids(order_ixs,1);
  %time = time(order_ixs,1);

  % --- Sort observations by (group, time) ---
  [~, order_ixs] = sortrows([group_ids(:), time(:)], [1 2]);

  X = X(order_ixs,:);
  y = y(order_ixs,1);
  group_ids = group_ids(order_ixs,1);
  time = time(order_ixs,1);


  group_ixs = findgroups(group_ids);
  group_start_stop_ixs = [ ...
    splitapply(@min, (1:m_dim)', group_ixs), ...
    splitapply(@max, (1:m_dim)', group_ixs)]';

  n_groups = size(group_start_stop_ixs, 2);
  n_cv_grps = min(n_cv_grps, n_groups);
  rndm_ixs = randsample(n_groups, n_groups);
  n_leftover = mod(n_groups, n_cv_grps);

  if (n_leftover > 0)
    rndm_ixs = [rndm_ixs; nan(n_cv_grps - n_leftover, 1)];
  end

  test_groups_ixs = reshape(rndm_ixs, n_cv_grps, [])';
  test_mses = repelem(nan, n_lams, n_cv_grps);

  for gx = 1:n_cv_grps
    test_group_gx_ixs = test_groups_ixs(:,gx);
    test_group_gx_ixs = sort(test_group_gx_ixs(~isnan(test_group_gx_ixs) ));

    test_ixs = arrayfun( ...
      @(grp_ix) group_start_stop_ixs(1,grp_ix):group_start_stop_ixs(2,grp_ix), ...
      test_group_gx_ixs, ...
      'UniformOutput', false);

    test_ixs = [test_ixs{:}];

    if isrow(test_ixs)
      test_ixs = test_ixs';
    end

    train_ixs = setdiff(1:m_dim, test_ixs);

    if isrow(train_ixs)
      train_ixs = train_ixs';
    end

    for lx = 1:n_lams
      train_group_labels = findgroups(group_ids(train_ixs,1) );
      [betahat, alphahat, stats] = gee_sparsereg( ...
        train_group_labels, ...
        time(train_ixs,1), ...
        X(train_ixs,:), ...
        y(train_ixs,1), ...
        dist, ...
        covar_type, ...
        lambdas(lx), ...
        'penalty', pentype, ...
        'penparam', penparam);

      y_hat = X(test_ixs,:) * betahat;

      %test_mses(lx,gx) = mse(y_hat, y(test_ixs,1) );
      test_mses(lx,gx) = mean((y(test_ixs,1) - y_hat).^2);
    end
  end

  err_mse = mean(test_mses, 2);
  err_se = std(test_mses')' / sqrt(n_cv_grps);
end