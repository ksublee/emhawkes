# Fixed fixtures for unequal-decay thinning; no fitting is timed.
source("tools/simulation-profile-cases.R")
equal_profile_case <- simulation_profile_case
simulation_profile_case <- function(id) {
  if(id %in% c("univariate","multi_equal","state_current")) return(equal_profile_case(id))
  d <- switch(id, small_unequal=1L, multi_unequal=4L, large_unequal=8L,
    rejection=2L, signed_unequal=2L, stop("Unknown fixture"))
  k <- 2L*d
  beta <- outer(seq(.8,1.6,length.out=d),seq(.6,2.4,length.out=k))
  alpha <- matrix(.35/k,d,k)
  mu <- rep(.5,d); initial <- matrix(.03/k,d,k)
  mapping <- lapply(seq_len(d),function(i) seq(i,k,by=d))
  if(id=="rejection") {
    mu <- rep(.01,d); beta <- beta*20
    alpha <- matrix(1,d,k); initial <- matrix(2,d,k)
  }
  if(id=="signed_unequal") {
    beta <- matrix(rep(c(1,2,1,2),d),d,byrow=TRUE)
    alpha <- matrix(rep(c(.08,-.02,.08,-.02),d),d,byrow=TRUE)
    initial <- alpha/2; mapping <- list(c(2L,1L),c(4L,3L))
  }
  list(object=new("hspec",mu=mu,alpha=alpha,beta=beta,type_col_map=mapping),
    initial=initial,seed=7361L,meta=data.frame(id=id,d=d,k=k,arrival="thinning"))
}
