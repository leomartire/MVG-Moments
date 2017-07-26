% Author:        Léo Martire.
% Mail:          leo.martire@outlook.com
% Description:   See below.
% Notes:         N/A.

function [NCK] = vectorial_nchoosek(N, K)
  % Computes nchoosek for each pair of components in the given vectors.
  % Used in compute_Gaussian_moments_Kan.
  % @param N "top" vector
  % @param K "bottom" vector
  % @return the vector [nchoosek(N(1), K(1)), nchoosek(N(2), K(2)), ... nchoosek(N(n), K(n))]

  if size(N) ~= size(K)
    NCK = - 1;
    error('[vectorial_nchoosek] size mismatch.');
  else
    S = size(N);
    s = length(N);
    N = reshape(N, s, 1);
    K = reshape(K, s, 1);
    NCK = [];
    for i = 1:s
      NCK(i) = nchoosek(N(i), K(i));
    end
    NCK = reshape(NCK, S(1), S(2)); % Set result back to N original shape.
  end
end