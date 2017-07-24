% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   See below.
% Notes:         Requires scripts "CGMoms_Kan.m"
%                                   "vectorial_nchoosek.m"

function [final_moms] = CGMoms(alpha_vals, MU, SIGMA, verbose)
  % Compute all moments up to a given degree of the multivariate Gaussian
  % distribution given by its mean vector and its covariance matrix.
  % The considered distribution is the classical n-dimensional normalised
  % Gaussian distribution:
  % (2 * pi)^(- 0.5 * n) * (det(SIGMA))^(- 0.5) * exp(- 0.5 * (x - MU)' * inv(SIGMA) * (x - MU))
  % References: Kan, R. (2008). From moments of sum to moments of product.
  %             Journal of Multivariate Analysis, 99(3):542 - 554.
  %             Willink, R. (2005). Normal moments and hermite
  %             polynomials. Statistics & Probability Letters,
  %             73(3):271 - 275.
  % @param alpha_vals full matrix of orders of the wanted moments (assumed
  %                   to be a matrix, size s * n where s is the number of
  %                   moments up to a maximum degree d)
  % @param MU mean vector of the multivariate Gaussian distribution
  %           (assumed to be a column, size n * 1)
  % @param SIGMA covariance matrix of the multivariate Gaussian
  %              distribution (assumed to be a positive definite matrix,
  %              size n * n)
  % @param (optional) verbose asks the script to be silent (if set to 0) or
  %                           talkative (default or if set anything not 0)
  % @return the list of all moments of the multivariate Gaussian
  %         distribution paramatrised by MU and SIGMA up to degree d,
  %         sorted according to the given matrix of orders

  % Check parameters. %%%%%%%%%
  if (~ exist('verbose', 'var') || isempty(verbose))
    verb = 1;
  else
    if verbose ~= 0
      verb = 1;
    else
      verb = 0;
    end
  end
  n = size(alpha_vals, 2); % Get the dimension of the problem.
  d = max(sum(alpha_vals, 2)); % Get the maximum order of moments needed.

  MU = reshape(MU, length(MU), 1); % Make MU a column.

  if (size(SIGMA, 1) ~= n) || (size(MU, 1) ~= n)
    error(['[', mfilename, '] dimension of orders, dimension of mean vector and dimension of covariance matrix are not coherent.']);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if verb == 1
    disp(['[', mfilename, '] Starting computation of ', num2str(size(alpha_vals, 1)), ' moments.'])
  end
  if ~ any(MU) && isdiag(SIGMA)
    % If the distribution is standard, use the closed form for %%
    % quickness and precision.                                 %%
    if verb == 1
      disp(['[', mfilename, '] Distribution is standard, we use the closed form.'])
    end
    d = max(sum(alpha_vals, 2));
    sigma = (2 * SIGMA(1, 1)) ^ 0.5; % convert to classic notations
    marg_moms = zeros(1 + d, 1); % Set all marginal moments to zeros (even and odd orders).
    K = 1 + (0:2:d); % Even orders indices (shifted of 1, because 0th order is at index 1 in Matlab).
    marg_moms(K, 1) = gamma(K / 2) .* realpow(sigma * ones(size(K)), K); % Formula for even order (K is already under the form (r + 1) where r is the order).
    final_moms = zeros(size(alpha_vals, 1), 1); % Prepare storage.
    is_odd = mod(sum(alpha_vals, 2), 2); % Vector which is 1 where total degree is odd (and 0 otherwise).
    final_moms(is_odd == 0, 1) = prod(marg_moms(alpha_vals(is_odd == 0, :) + 1), 2); % Set nonzero (even total degree) moments to their values (product of marginal moments).
    final_moms = final_moms * (2 * pi) ^ (- 0.5 * n) * det(SIGMA) ^ (- 0.5); % normalise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
    % If not, use the full formula. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verb == 1
      disp(['[', mfilename, '] Distribution is not standard, we use the full formulae.'])
    end
    
    % First moments. %%%%%%%%%%%%%%
    % Prepare pre-list for values of \alpha. This should contain all
    % combinations up to degree d, but containing at least one 0 and the
    % combination containing only ones. For example, for n=3 and d=4, this
    % should contain:
    % [000], [111], [200], [210], [400], [300], [310], ...
    % but not:
    % [211], [410], ...
 
    % Find row indexes in alpha_vals for which at least one coefficient is 0.
    at_least_one_zero = find(any(alpha_vals == 0, 2));
    % Select those rows.
    pre_alpha = alpha_vals(at_least_one_zero, :);
    % Add the row containing only ones if d>=n;
    if d >= n
      pre_alpha = [pre_alpha; ones(1, n)];
    end
 
    % Prepare corresponding list of moments.
    pre_moms = zeros(size(pre_alpha, 1), 1);
    % First n+1 moments are known.
    pre_moms(1:size(MU, 1) + 1) = [1; MU];
 
    % Use explicit formulas to compute first moments.
    for i = n + 2:size(pre_alpha, 1)
      alpha = pre_alpha(i, :);
      pre_moms(i) = CGMoms_Kan(alpha, MU, SIGMA);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % Save first moments computed %
    % and deduce which moments    %
    % still should be computed.   %
    final_moms = - 1 * ones(size(alpha_vals, 1), 1);
    already_computed = [];
    for k = 1:size(pre_moms, 1)
      index = find(ismember(alpha_vals, pre_alpha(k, :), 'rows'));
      already_computed = [already_computed, index];
      final_moms(index) = pre_moms(k);
    end
    to_be_computed = setdiff(1:size(final_moms, 1), already_computed);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % Remaining moments: use %%%%%%
    % Willink's recurrence   %%%%%%
    % formula.               %%%%%%
    for ind = to_be_computed
      alpha = alpha_vals(ind, :);
      k = find(alpha == max(alpha));
      k = k(1);
   
      % Construct the considered alpha_prev and check if it has indeed
      % already been computed.
      alpha_prev = alpha;
      alpha_prev(k) = alpha_prev(k) - 1;
      ind_alpha_prev = find(ismember(alpha_vals, alpha_prev, 'rows'));
      ind_prev = find(ismember(already_computed, ind_alpha_prev));
      if isempty(ind_prev)
        disp('[', mfilename, '] alpha_prev was not previously computed, we use Kan''s formula as emergency. Consider adding alpha_prev before current alpha in the list of moments to be computed for acceleration.');
        mu_prev = CGMoms_Kan(alpha_prev, MU, SIGMA);
      else
        mu_prev = final_moms(ind_alpha_prev);
      end
   
      tot = 0;
      for j = 1:n
        % Construct the considered alpha_prev_prev and check if it has
        % indeed already been computed.
        alpha_prev_prev = alpha_prev;
        alpha_prev_prev(j) = alpha_prev_prev(j) - 1;
        ind_alpha_prev_prev = find(ismember(alpha_vals, alpha_prev_prev, 'rows'));
        ind_prev_prev = find(ismember(already_computed, ind_alpha_prev_prev));
        if isempty(ind_prev_prev)
          disp('[', mfilename, '] alpha_prev_prev was not previously computed, we use Kan''s formula as emergency. Consider adding alpha_prev_prev before current alpha in the list of moments to be computed for acceleration.');
          mu_prev_prev = CGMoms_Kan(alpha_prev_prev, MU, SIGMA);
        else
          mu_prev_prev = final_moms(ind_alpha_prev_prev);
          %[tot, SIGMA(k, j), alpha_prev(j), mu_prev_prev]
        end
        tot = tot + SIGMA(k, j) * alpha_prev(j) * mu_prev_prev;
      end
      %[alpha, MU(k), mu_prev, tot, final_moms(ind), final_moms(ind)]
      final_moms(ind) = MU(k) * mu_prev + tot;
      already_computed = [already_computed, ind];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % Check if all moments have %%%
    % been computed.            %%%
    to_be_computed_final = setdiff(1:size(final_moms, 1), already_computed);
    if isempty(to_be_computed_final)
      disp(['[', mfilename, '] Full formulae: all moments computed.']);
    else
      to_be_computed_final
      error('[', mfilename, '] Full formulae: some moments still have to be computed.');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
