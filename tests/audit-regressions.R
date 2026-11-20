library(fastlpdidR)
i <- rep(1:12,each=10);t <- rep(1:10,12)
g <- ifelse(i<=4,3,ifelse(i<=8,5,0));D <- g>0 & t>=g
d <- data.frame(i=i,t=t,g=g,y=i*.2+t*.3+sin(i*t*.7)+D*(2+.2*g+.1*(t-g)))
call <- function(data,...) fit(data,'y','t',unit='i',gvar='g',event_window=c(0,2),...)
singleton <- d[!d$i %in% 6:8,]
stopifnot(is.finite(call(singleton,inference='point')$estimate))
jk <- call(singleton,inference='jackknife')$inference$jackknife
physical <- vapply(sort(unique(singleton$i)),function(i) call(singleton[singleton$i!=i,],inference='point')$estimate,numeric(1))
stopifnot(isTRUE(all.equal(unname(jk$deletions),unname(physical),tolerance=1e-10)))
ok <- call(d,inference='jackknife')
stopifnot(length(ok$inference$jackknife$deletions)==12,is.finite(ok$inference$jackknife$se))
extra <- rbind(d,data.frame(i=99,t=1,g=0,y=1))
result <- call(extra,inference='jackknife')
stopifnot(result$n_clusters==12,result$inference$jackknife$cluster_count==13,
          identical(result$inference$jackknife$inactive_clusters,99))
cat('LPDID_R_AUDIT_PASS\n')
