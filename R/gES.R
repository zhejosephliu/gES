# graphical exponential screening estimator
ges.mcmc <- function(x1, x2, hard.thres = 0.05, num.burn.in = 4000, num.estimation = 4000, plot.mcmc = TRUE) {
    n1 <- nrow(x1)
    n2 <- nrow(x2)

    if (ncol(x1) != ncol(x2)) {
        stop("ERROR: inconsistent number of columns of x1 and x2!")
    }
    p <- ncol(x1)

    # one-dimensional mapping for edges
    edge.mapping <- NULL
    k <- 0
    for(i in 1:(p - 1)) {
        for(j in (i + 1):p) {
            k <- k + 1
            edge.mapping <- rbind(edge.mapping, c(k, i, j))
        }
    }

    # missing edges
    edge.missing <- edge.mapping[, 2:3]

    # gES estimator of the precision matrix
    estimator <- matrix(0, p, p)

    # covariance matrix hard thresholding estimator
    S <- var(x2) * (abs(var(x2)) > hard.thres * sqrt(log(p) / n2))

    m <- 0
    count <- 0
    track <- NULL

    # MCMC
    while (1) {
        count <- count + 1
        edge.missing.backup <- edge.missing

        # number of existing edges before update
        num.existing.edge.backup <- p * (p - 1) / 2 - nrow(edge.missing)

        # edge to be updated
        r <- sample(1:(p * (p - 1) / 2), size = 1)

        if (num.existing.edge.backup < (p * (p - 1) / 2)) {
            index <- (2 * p - edge.missing[, 1]) * (edge.missing[, 1] - 1) / 2 +
                     edge.missing[, 2] - edge.missing[, 1]
            if (r %in% index) {
                id <- edge.missing[which(index == r), ]
                edge.missing <- matrix(as.matrix(edge.missing[-which(index == r), ]), ncol = 2)
            } else {
                id <- edge.mapping[r, 2:3]
                edge.missing <- rbind(edge.missing, edge.mapping[r, 2:3])
            }
        } else {
            id <- edge.mapping[r, 2:3]
            edge.missing <- rbind(edge.missing, edge.mapping[r, 2:3])
        }

        # number of existing edges after update
        num.existing.edge <- p * (p - 1) / 2 - nrow(edge.missing)

        w <- num.existing.edge - num.existing.edge.backup

        # glasso estimator before update
        res.glasso.backup <- glasso(var(x1), rho = 1e-10, zero = edge.missing.backup,
                                    penalize.diagonal = FALSE, start = "cold")$wi
        res.glasso.backup <- (res.glasso.backup + t(res.glasso.backup))/2

        # glasso estimator after update
        res.glasso <- glasso(var(x1), rho = 1e-10, zero = edge.missing,
                             penalize.diagonal = FALSE, start = "cold")$wi
        res.glasso <- (res.glasso + t(res.glasso))/2

        # ratio of prior distributions
        if (num.existing.edge == 0) {
            pi.ratio <- 1 / ((num.existing.edge.backup / (exp(1) * p * (p - 1))) ^ num.existing.edge.backup)
        } else {
            if (num.existing.edge.backup == 0) {
                pi.ratio <- (num.existing.edge / (exp(1) * p * (p - 1))) ^ num.existing.edge
            } else {
                pi.ratio <- ((1 + w / num.existing.edge.backup) ^ num.existing.edge) *
                            ((num.existing.edge.backup / (exp(1) * p * (p - 1))) ^ w)
            }
        }

        # ratio of weights
        const <- n2 / 2
        weight.backup <- sum(diag(res.glasso.backup %*% S))-log(det(res.glasso.backup))
        weight <- sum(diag(res.glasso %*% S)) - log(det(res.glasso))
        total.ratio <- min(1, pi.ratio * exp(const * (weight.backup - weight)))

        # acccept or reject the proposal
        rr <- runif(1)
        if (rr > total.ratio) {
            edge.missing <- edge.missing.backup
            current.estimator <- res.glasso.backup
        } else {
            current.estimator <- res.glasso
        }

        if (count > num.burn.in) {
            if (count > 0) {
                estimator <- estimator + current.estimator
                m <- m + 1
            }
        }

        track <- c(track, nrow(edge.missing))
        if (count >= num.burn.in + num.estimation) {
            break
        }
    }

    # convergence plot
    if (plot.mcmc == TRUE) {
        plot(p * (p - 1) / 2 - track, type = "l", xlab = "Number of MCMC iterations",
             ylab = "Number of selected edges")
    }

    # gES estimator of the precision matrix
    est.ges <- estimator / m

    return(est.ges)
}
