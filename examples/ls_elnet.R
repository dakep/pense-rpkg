# Compute the LS-EN regularization path for Freeny's revenue data
# (see ?freeny)
data(freeny)
x <- as.matrix(freeny[ , 2:5])

regpath <- elnet(x, freeny$y, alpha = 0.75)
plot(regpath)

# Extract the coefficients at a certain penalization level
coef(regpath, lambda = regpath$lambda[5])

# What penalization level leads to good prediction performance?
cv_results <- elnet_cv(x, freeny$y, alpha = 0.75, cv_repl = 25,
                       cv_k = 4, cv_measure = 'tau')
plot(cv_results, se_mult = 1)
plot(cv_results, se_mult = 1, what = 'coef.path')

# Extract the coefficients at the penalization level with
# smallest prediction error ...
coef(cv_results)
# ... or at the penalization level with prediction error
# statistically indistinguishable from the minimum.
coef(cv_results, lambda = 'se')
