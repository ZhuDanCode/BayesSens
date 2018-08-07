# VAR_CEML <- function(X, Y, beta, S, M, n, K, b_0, B_0, S_0, nu_0, T0){
#     XTY <- crossprod(X, Y)
#     YTY <- crossprod(Y, Y)
#     XTX <- crossprod(X, X)
#     mbeta <- rowMeans(beta)
#     cbeta <- cov(t(beta))
#     Beta <- t(mvnrnd(mbeta,cbeta,M))
#     l <- numeric(M)
#     p <- numeric(M)
#     CE <- numeric(M)
#     S <- mean(S,3) / T0
#     invS <- solve(S)
#     for (i in 1:M) {
#       Sigma <- rWishart(1, T0, S)
#       B <- t(reshape(Beta[,i],K,n))
#       W <- YTY - B*XTY - t(B*XTY) + B * XTX * t(B)
#       l(i) <-  -T0 * n/2 * log(2*pi) + T0/2 * log(det(Sigma)) - 1/2 * trace(Sigma*W)
#       Sigma <- solve(Sigma)
#       p(i) <- linvwishpdf(Sigma, nu_0, S_0) + lmvnpdf(Beta[,i], b_0, B_0)
#       CE(i) <- linvwishpdf(Sigma, T0, invS) + lmvnpdf(Beta[,i], mbeta, cbeta)
#     }
#     mlog <- max(l+p-CE)
#     ML <- log(mean(exp(l+p-CE-mlog)))+mlog
# }
#
#
# VAR_Chibs_ML <- function(X, Y, beta, sigma, M,
#                          nu_1, invB, invVb,
#                          n, K, b_0, B_0, S_0, nu_0, T0){
#   d <- n * K
#   beta <- reshape(rowMeans(beta),K,n)
#   sigma <- mean(sigma,3)
#   Y_minus_X_beta <- Y - X %*% beta
#   W <- t(Y_minus_X_beta) %*% (Y_minus_X_beta)
#   S <- W+S_0
#   Sinv <- solve(S)
#   logf <- 0
#   b <- reshape(beta, d, 1)
#   for (g in 1:M) {
#     Sigma <- rWishart(1, nu_1, Sinv)  # matrix
#     Vg <- solve(invB + Sigma %x% XTX)
#     XTY_Sigma <- reshape(XTY %*% Sigma, d, 1)
#     bg <- Vg %*% (invVb + XTY_Sigma)
#     logf <- logf + lmvnpdf(b, bg, Vg)
#   }
#   l <- -T * n/2 * log(2*pi) - T/2 * log(det(sigma)) - 1/2 * trace(solve(sigma, W))
#   prior <- lmvnpdf(b,b_0,B_0) + linvwishpdf(sigma, nu_0, S_0)
#   post <- linvwishpdf(sigma,nu_1,S) + logf / M
#   ML <- l + prior - post
# }
#
#
# reshape <- function(data0, dim0) {
#   dim(data0) <- dim0
#   data0
# }
