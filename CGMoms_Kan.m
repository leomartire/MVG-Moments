% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   See below.
% Notes:         Requires script "vectorial_nchoosek.m"

function mom = CGMoms_Kan(alpha, MU, SIGMA)
  % Compute the moment of given order of the multivariate Gaussian
  % distribution given by its mean vector and its covariance matrix.
  % Reference: Proposition 2 in [Kan, R. (2008). From moments of sum to
  %            moments of product. Journal of Multivariate Analysis,
  %            99(3):542 - 554].
  % @param alpha order of the wanted moment (assumed to be a line,
  %              size 1 * n)
  % @param MU mean vector of the multivariate Gaussian distribution
  %           (assumed to be a column, size n * 1)
  % @param SIGMA covariance matrix of the multivariate Gaussian
  %              distribution (assumed to be a positive definite matrix,
  %              size n * n)
  % @return the moment of order alpha of the multivariate Gaussian
  %         distribution given by MU and SIGMA

  % The first n sums of the formula are viewed as one sum on all possible
  % combinations. Example (n=2):
  % instead of
  %   \sum_{i=0}^{2} \sum_{j=0}^{3} X(i,j),
  % do
  %   \sum_{K\in\{(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3),(2,0),(2,1),(2,2),(2,3)\}} X(K).
  t = [];
  for i = 1:size(alpha, 2)
    t = [t, '0:', num2str(alpha(i)), ', '];
  end
  t = t(1:end - 2);
  beta_vals = eval(['combvec(', t, ')'])';

  % As introduced before, consider the sum on all possible combinations,
  % indexed by i. The sum indexed by r in the formula is considered as
  % another sum in itself.
  mom = 0;
  card_a = sum(alpha);
  semi_card_a = floor(card_a / 2);
  for i = 1:size(beta_vals, 1)
    beta = beta_vals(i, :);
    h = (alpha / 2 - beta)';
    card_b = sum(beta);
    for r = 0:semi_card_a
      % (-1)^{\sum_{i=1}^n(v_i)} from the formula.
      m_one_power = (- 1) ^ card_b;
      % (s_1 choose v_1) * ... * (s_n choose v_n) from the formula.
      nchoosek_terms = prod(vectorial_nchoosek(alpha, beta)); % vectorial_nchoosek computes the vector [nchoosek(alpha(1), beta(1)), nchoosek(alpha(2), beta(2)), ... nchoosek(alpha(n), beta(n))].
      % Numerator of the fraction in the formula.
      h_terms = ((0.5 * h'*SIGMA*h)^r) * ((h' * MU) ^ (card_a - 2 * r));
      % Denominator of the fraction in the formula.
      factorials = factorial(r) * factorial(card_a - 2 * r);
      % Add the product of the terms to the total.
      mom = mom + m_one_power * nchoosek_terms * h_terms / factorials;
    end
  end
end
