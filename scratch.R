n <- 5000
beta0 <- 0
beta <- c(1, 1)
gamma0 <- -1
gamma <- c(1, 1)
phi <- 5
theta <- 0.4

# used to generate sigma
rSq_theoretical <- 0.8
# used to generate scale for logistic regression (control for missingness rate)
rSq_latent <- 0.80

sigma <- theory_rSq_to_sigma(beta, phi, theta, rSq_theoretical)
df <- generate_full_data(beta0, beta, phi, theta, sigma, n, misspecified = FALSE)
df$z <- as.character(df$z)
s <- latent_rSq_to_logistic_scale(gamma, rSq_latent, df)
df.completeness <- generate_completeness(gamma0, gamma, df, s)
df_incomplete <- generate_incomplete_data(df, df.completeness)

write.csv(df_incomplete, "./simlated_data.csv", row.names = FALSE)
# full
lm_model <- lm(y ~ ., data = df) 
lm_model
mean((df$y - lm_model$fitted.values)^2)
# missing
mean(is.na(df_incomplete$z))
lm_model <- lm(y ~ ., data = df_incomplete) 
lm_model
mean((na.omit(df_incomplete)$y - lm_model$fitted.values)^2)



