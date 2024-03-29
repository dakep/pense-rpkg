% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_options.R
\name{en_dal_options}
\alias{en_dal_options}
\title{Use the DAL Elastic Net Algorithm}
\usage{
en_dal_options(
  max_it = 100,
  max_inner_it = 100,
  eta_multiplier = 2,
  eta_start_conservative = 0.01,
  eta_start_aggressive = 1,
  lambda_relchange_aggressive = 0.25
)
}
\arguments{
\item{max_it}{maximum number of (outer) iterations.}

\item{max_inner_it}{maximum number of (inner) iterations in each outer iteration.}

\item{eta_multiplier}{multiplier for the barrier parameter. In each iteration, the barrier must be more restrictive
(i.e., the multiplier must be > 1).}

\item{eta_start_conservative}{conservative initial barrier parameter. This is used if the previous penalty is
undefined or too far away.}

\item{eta_start_aggressive}{aggressive initial barrier parameter. This is used if the previous penalty is close.}

\item{lambda_relchange_aggressive}{how close must the lambda parameter from the previous penalty term be to use
an aggressive initial barrier parameter (i.e., what constitutes "too far").}
}
\value{
options for the DAL EN algorithm.
}
\description{
Use the DAL Elastic Net Algorithm
}
\seealso{
Other EN algorithms: 
\code{\link{en_admm_options}()},
\code{\link{en_cd_options}()},
\code{\link{en_lars_options}()}
}
\concept{EN algorithms}
