% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaD.R
\name{metaG}
\alias{metaG}
\title{Creat demographic tables}
\usage{
metaG(df, trt_var)
}
\arguments{
\item{df, }{input clinical data frame,}

\item{trt_var, }{treatment variable}
}
\value{
A data farme showing summary for clinical variables for subjects invovled in the study
}
\description{
This function create a demographic table. It assumes that the first column
contains the subjectId, the subsequent columns are clinical variables:diabetes,currentsmoker,africanamerican,
statintreatment,women,hbpmed, totcvd,stroke, statin, chd, mi, totalchol,trig,hdlc,ldlc, age, bmi, hscrp, tnfr2,
il6, glyca17, sicam
}
