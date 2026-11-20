# Independent lm oracle follows the frozen author regression equations.
# This is not a claim that official lpdid/fixest has been executed.
args <- commandArgs(trailingOnly=TRUE)
if(file.exists("R/fit.R")) source("R/fit.R") else library(fastlpdidR)
if(length(args)<4) {
  input <- if(file.exists("tests/panel.csv")) "tests/panel.csv" else "panel.csv"
  args <- c(input, tempfile(fileext=".csv"),tempfile(fileext=".csv"),tempfile(fileext=".csv"))
}
f <- read.csv(args[1]); out <- list(); cross <- list(); k<-0
at <- function(v,h) f[[v]][match(paste(f$i,f$t+h),paste(f$i,f$t))]
for(wgt in c("variance","reweighted")) for(cv in c(FALSE,TRUE)) {
  result <- fit(f,"y","t",unit="i",gvar="g",event_window=c(-2,2),
                weighting=wgt,covariates=if(cv) "x" else NULL)
  for(h in c(-2,0,1,2,999)) {
    hh <- if(h==999) 0:2 else h
    z <- rowMeans(sapply(hh,function(j) at("y",j)))-at("y",-1)
    d <- f$d-at("d",-1); x<-at("x",-1)
    ok <- (d==1 | at("d",max(0,hh))==0) & complete.cases(z,d)
    if(cv) ok<-ok & complete.cases(x)
    dat <- data.frame(z=z,d=d,x=x,t=f$t,i=f$i)[which(ok),]
    dat$n <- ave(dat$d,dat$t,FUN=length); dat$nt<-ave(dat$d,dat$t,FUN=sum)
    if(wgt=="reweighted" && cv) dat<-subset(dat,nt<n) else if(wgt=="reweighted") dat<-subset(dat,nt>0 & nt<n) else dat<-subset(dat,n>1)
    dat$w<-if(wgt=="reweighted") dat$n/(dat$n-dat$nt) else 1
    se<-NA
    if(wgt=="reweighted" && cv) {
      m<-lm(z~factor(t)+x,data=subset(dat,d==0))
      est<-mean(dat$z[dat$d==1]-predict(m,newdata=subset(dat,d==1)))
    } else {
      m<-lm(if(cv) z~d+factor(t)+x else z~d+factor(t),data=dat,weights=w)
      est<-coef(m)["d"]; X<-model.matrix(m); B<-solve(crossprod(X,dat$w*X))
      sc<-rowsum(X*as.numeric(resid(m)*dat$w),dat$i)
      G<-length(unique(dat$i)); N<-nrow(dat); K<-ncol(X)
      V<-B%*%crossprod(sc)%*%B*G/(G-1)*(N-1)/(N-K)
      se<-sqrt(V["d","d"])
    }
    got<-if(h==999) result else Filter(function(r) r$horizon==h,result$horizons)[[1]]
    stopifnot(abs(got$estimate-est)<1e-10)
    if(is.finite(se)) stopifnot(abs(got$native_se-se)<1e-10)
    k<-k+1; out[[k]]<-data.frame(weighting=wgt,covariates=as.integer(cv),horizon=h,estimate=est,se=se)
    cross[[k]]<-data.frame(weighting=wgt,covariates=as.integer(cv),horizon=h,estimate=got$estimate,se=got$native_se)
  }
  a<-fit(transform(f,y=1.31*y-.73*y2),"y","t",unit="i",gvar="g",event_window=c(-2,2),weighting=wgt,covariates=if(cv) "x" else NULL)
  b<-fit(transform(f,y=y2),"y","t",unit="i",gvar="g",event_window=c(-2,2),weighting=wgt,covariates=if(cv) "x" else NULL)
  stopifnot(abs(a$estimate-1.31*result$estimate+.73*b$estimate)<1e-10)
}
write.csv(do.call(rbind,out),args[2],row.names=FALSE)
write.csv(do.call(rbind,cross),args[3],row.names=FALSE)
obs<-f$g[!duplicated(f$i)]
ass<-rbind(c(tail(obs,1),head(obs,-1)),c(tail(obs,2),head(obs,-2)))
z<-fit(f,"y","t",unit="i",gvar="g",event_window=c(-2,2),inference=c("jackknife","ri_beta","ri_t_native"),assignments=ass)
write.csv(data.frame(estimate=z$estimate,native_se=z$native_se,jk_se=z$inference$jackknife$se,
                    ri_beta=z$inference$ri_beta$p,ri_t_native=z$inference$ri_t_native$p),args[4],row.names=FALSE)
stopifnot(z$inference$jackknife$df==11,length(z$inference$jackknife$deletions)==12)
bad<-try(fit(f,"y","t",unit="i",gvar="g",event_window=c(0,2),data_mode="rcs"),silent=TRUE)
stopifnot(inherits(bad,"try-error"))
cat("PASS: 20 independent lm point comparisons; 15 CR1 comparisons; 4 linearity checks; JK/RI smoke and RCS rejection\n")
