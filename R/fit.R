# Native R implementation. Oracle provenance and limitations are in VALIDATION.md.
fit <- function(data, y, time, unit=NULL, gvar=NULL, treatment=NULL,
                cluster=NULL, covariates=NULL, event_window=NULL,
                data_mode="auto", weighting="reweighted", covariate_lag=1L,
                inference="native", assignments=NULL, reps=99L, seed=0L, null=0, ...) {
  if (length(list(...))) stop("unsupported options")
  if (null!=0) stop("UNSUPPORTED: only null=0; nonzero sharp-null adjustment is not implemented")
  if (!data_mode %in% c("auto", "panel", "rcs")) stop("invalid data_mode")
  if (data_mode == "rcs" || is.null(unit)) stop("UNSUPPORTED: genuine panel required")
  if (is.null(gvar) && is.null(treatment)) stop("supply gvar or treatment")
  if (length(event_window)!=2 || any(!is.finite(event_window)) ||
      any(event_window!=floor(event_window)) || event_window[1]>event_window[2] ||
      event_window[2]<0) stop("explicit ordered integer event_window required")
  if (!weighting %in% c("variance", "reweighted")) stop("invalid weighting")
  if (!covariate_lag %in% c(0,1)) stop("covariate_lag must be 0 or 1")
  if (is.null(cluster)) cluster <- unit
  f <- data[order(data[[unit]], data[[time]]), , drop=FALSE]
  f$.u <- f[[unit]]; f$.t <- f[[time]]; f$.c <- f[[cluster]]
  if (anyNA(f[c(unit,time,cluster)]) || any(!is.finite(f$.t)) ||
      any(f$.t != floor(f$.t))) stop("nonmissing IDs and integer period time required")
  key <- function(u,t) paste(nchar(as.character(u)),u,t,sep=":")
  if (anyDuplicated(key(f$.u,f$.t))) stop("duplicate unit-time")
  if (any(vapply(split(f$.c,f$.u), function(v) length(unique(v))!=1, logical(1))))
    stop("cluster must be constant within unit")
  if (!is.null(gvar)) {
    f$.g <- as.numeric(f[[gvar]])
    f$.g[is.na(f$.g) | f$.g==0] <- Inf
    if (any(f$.g == -Inf) || any(is.finite(f$.g) & f$.g!=floor(f$.g))) stop("invalid cohort")
    if (any(vapply(split(f$.g,f$.u),function(v) length(unique(v))!=1,logical(1))))
      stop("gvar must be constant within unit")
  }
  at <- function(ff,col,h) ff[[col]][match(key(ff$.u,ff$.t+h),key(ff$.u,ff$.t))]
  if (!is.null(treatment)) {
    f$.d <- f[[treatment]]
    if (anyNA(f$.d) || any(!f$.d %in% c(0,1))) stop("binary nonmissing treatment required")
    for (ii in split(seq_len(nrow(f)),f$.u)) {
      if (any(diff(f$.d[ii])<0)) stop("nonabsorbing treatment unsupported")
      if (is.null(gvar) && any(diff(f$.d[ii])==1 & diff(f$.t[ii])!=1))
        stop("adoption in gap; supply known gvar")
    }
    if (!is.null(gvar) && any(f$.d != as.numeric(f$.t>=f$.g))) stop("inconsistent cohort/treatment")
    if (is.null(gvar)) {
      f$.g <- Inf
      for (ii in split(seq_len(nrow(f)),f$.u)) {
        z <- f$.t[ii][f$.d[ii]==1]
        if(length(z)) f$.g[ii] <- min(z)
      }
    }
  } else f$.d <- as.numeric(f$.t>=f$.g)
  f$.y <- as.numeric(f[[y]])
  regress <- function(ff, hh) {
    z <- rowMeans(vapply(hh,function(h) at(ff,".y",h),numeric(nrow(ff)))) - at(ff,".y",-1)
    d <- ff$.d-at(ff,".d",-1)
    clean <- (d==1 | at(ff,".d",max(0,hh))==0)
    x <- if(length(covariates)) sapply(covariates,function(v) at(ff,v,-covariate_lag)) else matrix(numeric(0),nrow(ff),0)
    x <- matrix(x,nrow=nrow(ff))
    keep <- which(clean & is.finite(d) & is.finite(z) & rowSums(!is.finite(x))==0)
    if (!length(keep)) stop("IDENTIFICATION_FAILURE: empty clean sample")
    n <- ave(d[keep], ff$.t[keep], FUN=length)
    nt <- ave(d[keep],ff$.t[keep],FUN=sum)
    keep <- keep[if(weighting=="reweighted" && length(covariates)) nt<n else if(weighting=="reweighted") nt>0 & nt<n else n>1]
    if(!length(keep) || !any(d[keep]==1) || !any(d[keep]==0)) stop("IDENTIFICATION_FAILURE: no overlap")
    z <- z[keep]; d <- d[keep]; x <- x[keep,,drop=FALSE]
    tt <- ff$.t[keep]
    fe <- outer(tt,sort(unique(tt)),"==")*1
    cl <- ff$.c[keep]; G <- length(unique(cl)); N <- length(d)
    if (weighting=="reweighted" && length(covariates)) {
      X <- cbind(fe,x); X0 <- X[d==0,,drop=FALSE]
      if(qr(X0)$rank<ncol(X0)) stop("IDENTIFICATION_FAILURE: singular control regression")
      B <- solve(crossprod(X0)); b <- B %*% crossprod(X0,z[d==0])
      r <- as.vector(z-X%*%b); est <- mean(r[d==1]); m <- colMeans(X[d==1,,drop=FALSE])
      score <- (d*(r-est)-(1-d)*r*as.vector(X%*%B%*%m)*sum(d))/sum(d)
      correction <- 1; df <- NA_real_; convention <- "RA unconditional CR0"
    } else {
      X <- cbind(d,fe,x); w <- rep(1,N)
      if(weighting=="reweighted") w <- ave(d,tt,FUN=length)/ave(1-d,tt,FUN=sum)
      K <- qr(X)$rank
      if(K<ncol(X)) stop("IDENTIFICATION_FAILURE: singular LP regression")
      B <- solve(crossprod(X,w*X)); b <- B%*%crossprod(X,w*z)
      est <- b[1]; score <- as.vector(X%*%B[,1])*w*as.vector(z-X%*%b)
      correction <- if(G>1 && N>K) G/(G-1)*(N-1)/(N-K) else NA_real_
      df <- G-1; convention <- "cluster CR1: G/(G-1)*(N-1)/(N-K), time FE counted"
    }
    cs <- rowsum(score,cl,reorder=FALSE)
    se <- if(G>1) sqrt(correction*sum(cs^2)) else NA_real_
    p <- if(!is.finite(se) || se==0) NA_real_ else if(is.na(df)) 2*pnorm(-abs(est/se)) else 2*pt(-abs(est/se),df)
    list(estimate=as.numeric(est),native_se=se,native_df=df,native_p=p,
         native_vce=convention,native_validation="official lpdid 1.1.0 deterministic fixture; not full certification",nobs=N,
         n_treated=sum(d),n_clusters=G,cohorts_used=sort(unique(tt[d==1])),sample_rows=keep-1L,
         support_times=sort(unique(tt)),estimation_clusters=sort(unique(cl)))
  }
  post <- seq(max(0,event_window[1]),event_window[2])
  point <- function(ff) regress(ff,post)
  out <- point(f)
  out$status <- "PASS"; out$estimator <- "lpdid"; out$data_mode <- "panel"
  out$weighting <- weighting; out$event_window <- event_window; out$pooled_horizons <- post
  out$target_label <- sprintf("%s pooled LP-DiD ATT [%s,%s]",weighting,min(post),max(post))
  out$covariates_used <- covariates; out$covariate_timing <- paste0("t-",covariate_lag)
  out$horizons <- lapply(seq(event_window[1],event_window[2]),function(h) {
    r <- if(h == -1) list(estimate=0,native_se=NA,status="REFERENCE") else
      tryCatch(c(regress(f,h),list(status="PASS")),error=function(e)
        list(estimate=NA,status="IDENTIFICATION_FAILURE",error=conditionMessage(e)))
    c(list(horizon=h),r)
  })
  if(any(!inference %in% c("point","native","jackknife","ri_beta","ri_t_native"))) stop("unknown inference")
  clusters <- sort(unique(f$.c))
  jk <- function(ff) {
    full <- point(ff)
    vals <- vapply(clusters,function(g) {
      reduced <- point(ff[ff$.c!=g,,drop=FALSE])
      reduced$estimate
    },numeric(1))
    G <- length(vals); se <- sqrt((G-1)/G*sum((vals-mean(vals))^2))
    tolerance <- 64*.Machine$double.eps*max(1,abs(full$estimate),abs(vals))
    if(!is.finite(se) || se<=tolerance) stop("IDENTIFICATION_FAILURE: zero jackknife variance")
    list(se=se,df=G-1,deletions=vals,p=2*pt(-abs((full$estimate-null)/se),G-1),
         cluster_count=G,cluster_scope="all input inference clusters",
         inactive_clusters=setdiff(clusters,full$estimation_clusters),
         support="surviving pooled cohort/time support")
  }
  out$inference <- list()
  if(any(c("jackknife") %in% inference)) out$inference$jackknife <- jk(f)
  if(any(c("ri_beta","ri_t_native") %in% inference)) {
    if(any(vapply(split(f$.g,f$.c),function(v) length(unique(v))!=1,logical(1)))) stop("one cohort per assignment cluster required")
    observed <- vapply(clusters,function(g) f$.g[match(g,f$.c)],numeric(1))
    if(is.null(assignments)) {
      had_seed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
      if(had_seed) old_seed <- get(".Random.seed",.GlobalEnv)
      on.exit(if(had_seed) assign(".Random.seed",old_seed,.GlobalEnv) else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
      set.seed(seed)
      total <- exp(lgamma(length(observed)+1)-sum(lgamma(as.numeric(table(observed))+1)))
      R <- min(reps,round(total)-1)
      if(R<1 || reps!=floor(reps)) stop("positive reps and nontrivial assignment space required")
      rows <- list(); seen <- paste(observed,collapse=",")
      for(i in seq_len(max(1000,100*R))) {
        a <- observed[sample.int(length(observed))]; k <- paste(a,collapse=",")
        if(!(k %in% seen)) {rows[[length(rows)+1]]<-a; seen<-c(seen,k)}
        if(length(rows)==R) break
      }
      if(length(rows)!=R) stop("unique sampling exhausted; supply assignments")
      assignments <- do.call(rbind,rows)
    }
    assignments <- as.matrix(assignments)
    assignments[is.na(assignments) | assignments==0] <- Inf
    if(ncol(assignments)!=length(clusters) || nrow(assignments)<1) stop("invalid assignment dimensions")
    keys <- apply(assignments,1,paste,collapse=",")
    if(anyDuplicated(keys) || paste(observed,collapse=",") %in% keys ||
       any(apply(assignments,1,function(a) !identical(sort(as.numeric(a)),sort(observed))))) stop("invalid cohort-preserving assignments")
    draws <- numeric(nrow(assignments)); tn <- draws
    for(i in seq_len(nrow(assignments))) {
      ff <- f; ff$.g <- assignments[i,match(ff$.c,clusters)]; ff$.d <- as.numeric(ff$.t>=ff$.g)
      pr <- point(ff); draws[i] <- pr$estimate-null
      if("ri_t_native" %in% inference) {
        if(!is.finite(pr$native_se) || pr$native_se<=.Machine$double.eps) stop("invalid randomized native variance")
        tn[i]<-draws[i]/pr$native_se
      }
    }
    tailp <- function(a,b) {
      tolerance <- 64*.Machine$double.eps*pmax(1,abs(a),abs(b))
      (1+sum(abs(a)>=abs(b)-tolerance))/(length(a)+1)
    }
    if("ri_beta" %in% inference) out$inference$ri_beta <- list(p=tailp(draws,out$estimate-null),statistic=out$estimate-null,draws=draws)
    if("ri_t_native" %in% inference) {
      if(!is.finite(out$native_se) || out$native_se<=.Machine$double.eps) stop("invalid native variance")
      stat<-(out$estimate-null)/out$native_se
      out$inference$ri_t_native<-list(p=tailp(tn,stat),statistic=stat,draws=tn)
    }
    assignments[is.infinite(assignments)] <- 0
    out$assignments <- assignments; out$assignment_clusters <- clusters
  }
  out
}
