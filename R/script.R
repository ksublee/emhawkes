
# mu <- c(0.1, 0.1)
# alpha <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=T)
# beta <- c(0.9, 0.9, 0.9, 0.9)
# mark_hawkes <- function(){rgeom(n=1, prob=0.5) + 1}
# h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta, mark_hawkes=mark_hawkes)
# res <- hsim(h1)


#
# mu <- matrix(c(0.08, 0.08, 0.05, 0.05), nrow = 4)
#
#
# alpha <- function(param = c(alpha11 = 0, alpha12 = 0.4, alpha33 = 0.5, alpha34 = 0.3)){
#   matrix(c(param["alpha11"], param["alpha12"], 0, 0,
#            param["alpha12"], param["alpha11"], 0, 0,
#            0, 0, param["alpha33"], param["alpha34"],
#            0, 0, param["alpha34"], param["alpha33"]), nrow = 4, byrow = T)
# }
#
#
# beta <- matrix(c(rep(0.6, 8), rep(1.2, 8)), nrow = 4, byrow = T)
#
# impact <- function(param = c(alpha1n=0, alpha1w=0.2, alpha2n=0.001, alpha2w=0.1),
#                    n=n, N=N, ...){
#
#   Psi <- matrix(c(0, 0, param['alpha1w'], param['alpha1n'],
#                   0, 0, param['alpha1n'], param['alpha1w'],
#                   param['alpha2w'], param['alpha2n'], 0, 0,
#                   param['alpha2n'], param['alpha2w'], 0, 0), nrow=4, byrow=T)
#
#   ind <- N[,"N1"][n] - N[,"N2"][n] > N[,"N3"][n] - N[,"N4"][n] + 0.5
#
#   km <- matrix(c(!ind, !ind, !ind, !ind,
#                  ind, ind, ind, ind,
#                  ind, ind, ind, ind,
#                  !ind, !ind, !ind, !ind), nrow = 4, byrow = T)
#
#   km * Psi
# }
#
# h <- new("hspec",
#          mu = mu, alpha = alpha, beta = beta, impact = impact)
#
# hr <- hsim(h, size=1000)
#
# #plot(hr$arrival, hr$N[,'N1'] - hr$N[,'N2'], type='s')
# #lines(hr$N[,'N3'] - hr$N[,'N4'], type='s', col='red')
#
# fit <- hfit(h, hr$inter_arrival, hr$type)
# summary(fit)


# mu <- c(0.15, 0.15)
# alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=T)
# beta <- c(2.6, 2.6, 2.6, 2.6)
# rmark <- function(param = c(p=0.65), ...){
#   rgeom(1, p=param[1]) + 1
# }
#
# impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
#   ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
#   alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
#   #alpha * ma * matrix( c(rep(param["eta1"], 2), rep(param["eta2"], 2)), nrow=2)
# }
# h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
#           rmark = rmark,
#           impact=impact)
# res <- hsim(h1, size=1000, lambda0 = matrix(rep(0.1,4), nrow=2))
#
# fit <- hfit(h1,
#             inter_arrival = res$inter_arrival,
#             type = res$type,
#             mark = res$mark,
#             lambda0 = matrix(rep(0.1,4), nrow=2))
#
# summary(fit)
#
#
# data <- read.csv("G:/Data/Tick/Hawkes_mark/simulation/full_3/1.csv", header = F)
# mu <- c(0.13, 0.11)
# alpha <- matrix(c(0.7, 0.6, 0.5, 0.4), nrow=2, byrow=T)
# beta <- c(2.4, 2.2, 2.42, 2.26)
# impact <- function(param = c(eta1=0.2), alpha, n, mark, ...){
#   ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
#   alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
#   #alpha * ma * matrix( c(rep(param["eta1"], 2), rep(param["eta2"], 2)), nrow=2)
# }
# h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
#           impact=impact)
#
# fit <- hfit(h1,
#             inter_arrival = c(0,data[[3]]),
#             type = c(0,data[[2]]),
#             mark = c(0,data[[4]]/0.005),
#             lambda0 = matrix(rep(0.1,4), nrow=2))
#
# summary(fit)


# cm <- matrix(rep(0, 5*1), ncol=5)
#
# for (j in 1:1){
#   mu <- c(0.15, 0.15)
#   alpha <- matrix(c(0.75, 0.6, 0.6, 0.75), nrow=2, byrow=T)
#   beta <- c(2.6, 2.6, 2.6, 2.6)
#   rmark <- function(param = c(p=0.65), ...){
#     rgeom(1, p=param[1]) + 1
#   }
#
#   impact <- function(param = c(eta1=0.2), alpha=alpha, n = n, mark = mark, ...){
#     ma <- matrix(rep(mark[n]-1, 4), nrow = 2)
#     alpha * ma * matrix( rep(param["eta1"], 4), nrow=2)
#
#   }
#   h1 <- new("hspec", mu=mu, alpha=alpha, beta=beta,
#             rmark = rmark,
#             impact=impact)
#   res <- hsim(h1, size=5000, lambda0 = matrix(rep(0.1,4), nrow=2))
#
#   fit <- hfit(h1,
#               inter_arrival = res$inter_arrival,
#               type = res$type,
#               mark = res$mark,
#               lambda0 = matrix(rep(0.1,4), nrow=2))
#
#   cm[j, ] <- coef(fit)
# }
#
# colMeans(cm[1:10,])


# mu <- function(param = c(mu1 = 0.08, eta1 = 1, eta2 =1), n=n, N=N, ...){
#
#
#   level <- N[,"N1"][n-1] - N[,"N2"][n-1] - (N[,"N3"][n-1] - N[,"N4"][n-1])
#   matrix( c(param["mu1"], param["eta1"]*level, param["eta2"]*level, param["mu1"]), nrow = 4)
#
#
# }
# alpha <- function(param = c(alpha11 = 0.75, alpha14 = 0.6, alpha41 = 0.6, alpha44 = 0.75)){
#   matrix(c(param["alpha11"], 0, 0, param["alpha14"],
#            0, 0, 0, 0,
#            0, 0, 0, 0,
#            param["alpha41"], 0, 0, param["alpha44"]), nrow = 4, byrow = T)
# }
#
# beta <- matrix(rep(2.6, 16), nrow=4, byrow=T)
# h <- new("hspec", mu, alpha, beta)
# hr <- hsim(h, size=1000)
#
# Nau <- hr$N[,"N1"]
# Nad <- hr$N[,"N2"]
# Nbu <- hr$N[,"N3"]
# Nbd <- hr$N[,"N4"]
#
# plot(Nau - Nad, type = "s")
# lines(Nbu - Nbd - 1, col="red", type="s")



