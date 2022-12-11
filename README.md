## Missing data methods
EM algorithm with simple & multivariate linear regression

### Notes on EM-Simple Linear Regression:
1. High missingness rates require more iterations of M-steps to achieve a convergence, it (also high noise) cause lower efficiency
2. applying EM over full likelihood performs better than apply it over observed likelihood: alway better coverage rate on X, (slightly) more efficient, less bias
3. EM is more efficient than CC when missingness rate increases, but slightly more biased under MCAR
4. EM-full tend to outperform CC when both noise and missingness rate are high; still CC is good due to MCAR
5. Misspecifying won't necessarily affect the estimation on beta
