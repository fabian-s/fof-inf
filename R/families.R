# Family table shared by the generator and the fitter. For each family:
# - fit():   the mgcv family object pffr is fitted with;
# - linkinv: inverse link (response scale = linkinv(link scale));
# - draw(z, eta, sigma): response values from latent Gaussian values z. Gaussian: eta +
#   sigma z. Other families: quantiles at u = pnorm(z) (Gaussian copula; needs unit-variance
#   z), so every point has exactly the family's marginal law with mean linkinv(eta).
#   `sigma` is the error SD, used by gaussian and scat only.

scat_scale <- function(sigma) sigma / sqrt(SCAT_DF / (SCAT_DF - 2)) # Var = sigma^2

FAMILIES <- list(
  gaussian = list(
    fit = \() gaussian(),
    linkinv = identity,
    draw = \(z, eta, sigma) eta + sigma * z
  ),
  scat = list(
    fit = \() mgcv::scat(),
    linkinv = identity,
    draw = \(z, eta, sigma) eta + scat_scale(sigma) * qt(pnorm(z), df = SCAT_DF)
  ),
  poisson = list(
    fit = \() poisson(),
    linkinv = exp,
    draw = \(z, eta, sigma) qpois(pnorm(z), exp(eta))
  ),
  binomial = list(
    fit = \() binomial(),
    linkinv = plogis,
    draw = \(z, eta, sigma) qbinom(pnorm(z), 1, plogis(eta))
  ),
  betar = list(
    fit = \() mgcv::betar(link = "logit"),
    linkinv = plogis,
    draw = \(z, eta, sigma) {
      mu <- plogis(eta)
      qbeta(pnorm(z), mu * BETAR_PHI, (1 - mu) * BETAR_PHI)
    }
  ),
  nb = list(
    fit = \() mgcv::nb(),
    linkinv = exp,
    draw = \(z, eta, sigma) qnbinom(pnorm(z), size = NB_SIZE, mu = exp(eta))
  )
)
