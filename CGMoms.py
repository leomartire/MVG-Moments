import numpy as np
import math

###############################
# If SciPy cannot be found,   #
# find here re-codings of its #
# "comb" and "gamma"          #
# functions.                  #
###############################
scipy_found = False
try :
  import scipy
  scipy_found = True
except ImportError, e:
  scipy_found = False
  print("INFO: SciPy not found, functions will use self-made functions instead of SciPy's ones.")

if (scipy_found == False) :
  import math
  def scalar_comb(n, k):
    return (math.factorial(n) / (math.factorial(k) * math.factorial(n - k)))

  def comb_without_scipy(N, K):
    if (N.shape != K.shape) :
      NCK = -1
      print('ERROR: size mismatch.')
    else :
      S = N.shape
      s = N.size
      N = np.reshape(N, (1, s))# make as line vector
      K = np.reshape(K, (1, s))# make as line vector
      NCK = list()
      for i in range(0, s) :
        NCK.append(scalar_comb(N[0, i], K[0, i]))
      NCK = np.reshape(np.array(NCK), S)# Set result back to N original shape.
    return (NCK)

  def gamma_without_scipy(X):
    X = np.array(X)
    X = np.reshape(X, X.size)
    G = np.zeros(X.size)
    for i in range(0, X.size) :
      G[i] = math.gamma(X[i])
    return (G)
else :
  from scipy.special import comb # Imported for CGMoms_Kan.
  from scipy.special import gamma # Imported for CGMoms.
###############################

###############################
# Main functions.             #
###############################
def CGMoms_Kan(alpha = None, MU = None, SIGMA = None) :
  # Compute the moment of given order of the multivariate Gaussian distribution given by its mean vector and its covariance matrix.
  # Reference: Proposition 2 in [Kan, R. (2008). From moments of sum to moments of product. Journal of Multivariate Analysis, 99(3):542 - 554].
  # @param alpha order of the wanted moment
  # @param MU mean vector of the multivariate Gaussian distribution
  # @param SIGMA covariance matrix of the multivariate Gaussian distribution
  # @return the moment of order alpha of the multivariate Gaussian distribution given by MU and SIGMA

  alpha = np.array(alpha) # Ensure alpha is a np.array.
  alpha = np.reshape(alpha, alpha.size) # Ensure alpha is a vector.
  MU = np.array(MU) # Ensure MU is a np.array.
  MU = np.reshape(MU, MU.size) # Ensure MU is a vector.
  SIGMA = np.array(SIGMA) # Ensure SIGMA is a np.array.
  
  # The first n sums of the formula are viewed as one sum on all possible combinations. Example (n=2):
  # instead of
  #   \sum_{i=0**{2 \sum_{j=0**{3 X(i,j),
  # do
  #   \sum_{K\in\{(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3),(2,0),(2,1),(2,2),(2,3)\ X(K).
  list_of_indices=["np.arange(0," + str(alpha[i] + 1) + ")" for i in range(0, alpha.size)]
  exec("beta_vals=np.array(np.meshgrid(" + (", ".join(list_of_indices)) + ")).T.reshape(-1," + str(alpha.size) + ")")
  
  # As introduced before, consider the sum on all possible combinations, indexed by i. The sum indexed by r in the formula is considered as another sum in itself.
  mom = 0
  card_a = np.sum(alpha)
  semi_card_a = (int)(np.floor(card_a / 2))
  for i in range(0,beta_vals.shape[0]) :
    beta = beta_vals[i, :]
    h = alpha / 2.0 - beta
    card_b = np.sum(beta)
    for r in range(0,semi_card_a + 1) :
      # (-1)^{\sum_{i=1}^n(v_i)} from the formula.
      m_one_power = (-1) ** card_b
      # (s_1 choose v_1) * ...*(s_n choose v_n) from the formula.
      if (scipy_found == False) :
        nchoosek_terms = np.prod(comb_without_scipy(alpha, beta))
      else :
        # TODO: check if this works.
        print("WARNING: Test first.")
        # nchoosek_terms = np.prod(comb(alpha, beta, exact=True, repetition=False))
      # Numerator of the fraction in the formula.
      h_terms = ((0.5 * np.dot(np.dot(h, SIGMA), h)) ** r) * (np.dot(h, MU)) ** (card_a - 2 * r)
      # Denominator of the fraction in the formula.
      factorials = math.factorial(r) * math.factorial(card_a - 2 * r)
      # Add the product of the terms to the total.
      mom = mom + m_one_power * nchoosek_terms * h_terms / factorials
  return(mom)
  

def CGMoms(alpha_vals, MU, SIGMA, verbose = False) :
  # Compute all moments up to a given degree of the multivariate Gaussian distribution given by its mean vector and its covariance matrix. The considered distribution is the classical n-dimensional normalised Gaussian distribution:
  # (2 * pi)^(- 0.5 * n) * (det(SIGMA))^(- 0.5) * exp(- 0.5 * (x - MU)' * inv(SIGMA) * (x - MU))
  # References: Kan, R. (2008). From moments of sum to moments of product.
  #             Journal of Multivariate Analysis, 99(3):542 - 554.
  #             Willink, R. (2005). Normal moments and hermite
  #             polynomials. Statistics & Probability Letters,
  #             73(3):271 - 275.
  # @param alpha_vals full matrix of orders of the wanted moments (assumed to be a matrix, size s * n where s is the number of moments up to a maximum degree d)
  # @param MU mean vector of the multivariate Gaussian distribution (assumed to be a column, size n * 1)
  # @param SIGMA covariance matrix of the multivariate Gaussian distribution (assumed to be a positive definite matrix, size n * n)
  # @param (optional) verbose asks the script to be silent (default, or if set to False) or talkative (if set to True)
  # @return the list of all moments of the multivariate Gaussian distribution paramatrised by MU and SIGMA up to degree d, sorted according to the given matrix of orders
  
  alpha_vals = np.array(alpha_vals) # Ensure alpha_vals is a np.array.
  MU = np.array(MU) # Ensure MU is a np.array.
  MU = np.reshape(MU, MU.size) # Ensure MU is a vector.
  SIGMA = np.array(SIGMA) # Ensure SIGMA is a np.array.

  n = alpha_vals.shape[1] # Get the dimension of the problem.
  d = np.max(np.sum(alpha_vals, 1)) # Get the maximum order of moments needed.

  if(np.all(MU == np.zeros(n)) and np.all(SIGMA==SIGMA[1,1]*np.eye(n))):
    if (verbose) :
      print("INFO: Distribution is standard, we use the closed form.")
    sigma = (2 * SIGMA[1, 1]) ** 0.5; # convert to classic notations
    marg_moms = np.zeros(1 + d); # Set all marginal moments to zeros (even and odd orders).
    K = np.arange(0,d+1,2); # Even orders indices.
    if (scipy_found == False) :
      marg_moms[K] = sigma ** (K+1) * gamma_without_scipy((K+1)/2.0)
    else :
      # TODO: check if this works.
      print("WARNING: Test first.")
      # marg_moms[K] = sigma ** (K+1) * gamma((K+1)/2.0)
    final_moms = np.zeros(alpha_vals.shape[0]) # Prepare storage.
    is_odd = np.mod(np.sum(alpha_vals, 1), 2) # Vector which is 1 where total degree is odd (and 0 otherwise).
    final_moms[is_odd == 0] = np.prod(marg_moms[alpha_vals[is_odd == 0]], 1) # Set nonzero (even total degree) moments to their values (product of marginal moments).
    final_moms = final_moms * (2 * math.pi) ** (- 0.5 * n) * np.prod(np.diag(SIGMA)) ** (- 0.5) # normalise (note det(SIGMA) = np.prod(np.diag(SIGMA)) here)
  else:
    if (verbose) :
      print("INFO: Distribution is not standard, we use the full formulae.")
    
    # First moments. ##############
    # Prepare pre-list for values of \alpha. This should contain all
    # combinations up to degree d, but containing at least one 0 and the
    # combination containing only ones. For example, for n=3 and d=4, this
    # should contain:
    # [000], [111], [200], [210], [400], [300], [310], ...
    # but not:
    # [211], [410], ...

    # Find row indexes in alpha_vals for which at least one coefficient is 0.
    at_least_one_zero = np.array([0 in alpha_vals[i, :] for i in range(0, alpha_vals.shape[0])]);
    # Select those rows.
    pre_alpha = alpha_vals[at_least_one_zero, :];
    # Add the row containing only ones if d>=n;
    if(d >= n):
      pre_alpha = np.vstack((pre_alpha, np.ones(n))).astype(int)

    # Prepare corresponding list of moments.
    pre_moms = np.zeros(pre_alpha.shape[0]);
    # First n+1 moments are known.
    pre_moms[0:MU.size + 1] = np.hstack((1, MU));

    # Use explicit formulas to compute first moments.
    for i in range(n+1,pre_alpha.shape[0]):
      alpha = pre_alpha[i, :];
      pre_moms[i] = CGMoms_Kan(alpha, MU, SIGMA);
    ###############################

    # Save first moments computed #
    # and deduce which moments    #
    # still should be computed.   #
    final_moms = - 1 * np.ones(alpha_vals.shape[0])
    already_computed = list()
    for k in range(0, pre_moms.size):
      index = np.where(np.all(alpha_vals == pre_alpha[k, :], 1)==True)[0][0]
      already_computed.append(index)
      final_moms[index] = pre_moms[k]
    to_be_computed = np.setdiff1d(range(0, final_moms.size), already_computed)
    ###############################

    # Remaining moments: use ######
    # Willink's recurrence   ######
    # formula.               ######
    for ind in to_be_computed:
      alpha = np.copy(alpha_vals[ind, :])
      k = np.where(alpha == np.max(alpha))[0][0]

      # Construct the considered alpha_prev and check if it has indeed
      # already been computed.
      alpha_prev = np.copy(alpha)
      alpha_prev[k] -= 1
      ind_alpha_prev = np.where(np.all(alpha_vals == alpha_prev, 1)==True)[0][0]
      if (ind_alpha_prev in already_computed):
        mu_prev = final_moms[ind_alpha_prev]
      else:
        print('WARNING: alpha_prev was not previously computed, we use Kan''s formula as emergency. Consider adding alpha_prev before current alpha in the list of moments to be computed for acceleration.')
        mu_prev = CGMoms_Kan(alpha_prev, MU, SIGMA)

      tot = 0
      for j in range(0,n):
        # Construct the considered alpha_prev_prev and check if it has
        # indeed already been computed.
        alpha_prev_prev = np.copy(alpha_prev)
        alpha_prev_prev[j] -= 1
        ind_alpha_prev_prev = np.where(np.all(alpha_vals == alpha_prev_prev, 1) == True)[0][0]
        if (ind_alpha_prev_prev in already_computed):
          mu_prev_prev = final_moms[ind_alpha_prev_prev]
        else:
          disp('[', mfilename, '] alpha_prev_prev was not previously computed, we use Kan''s formula as emergency. Consider adding alpha_prev_prev before current alpha in the list of moments to be computed for acceleration.')
          mu_prev_prev = CGMoms_Kan(alpha_prev_prev, MU, SIGMA)
        tot = tot + SIGMA[k, j] * alpha_prev[j] * mu_prev_prev
      final_moms[ind] = MU[k] * mu_prev + tot
      already_computed.append(ind)
    ###############################

    # Check if all moments have ###
    # been computed.            ###
    to_be_computed_final = np.setdiff1d(range(0,final_moms.size), already_computed)
    if (to_be_computed_final.size == 0):
      print('INFO: Full formulae: all moments computed.')
      print(verbose)
    else:
      to_be_computed_final
      print('INFO: Full formulae: some moments still have to be computed.')
    ###############################

  return(final_moms)
