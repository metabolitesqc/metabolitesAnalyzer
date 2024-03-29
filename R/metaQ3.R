#' Checking metabolites quality3
#'
#' This function check the metabolites data quality. It assumes that the first column of df_sample
#' contains the subjectId, with or without year, plate_well in the second and third column.  in the subsequent columns
#' are metabolites measurements. It assumes the first column of df_pp contains plate_well position for the pooled plasma
#' samples.
#
#' @param df_sample, experimental metabolites data
#' @param df_pp, pooled plasma data
#' @param plot_missing, defaults to true
#' @param output_missing, defaults to true
#' @param plot_CV, defaults to true
#' @param output_CV, defaults to true
#' @param output_meta_summary, defaults to true
#' @param plot_well_median, defaults to true
#'
#' @return list showing missing, cv, meta_summary
#' @export
metaQ3 = function(df_sample,
                 df_pp,
                 plot_missing= TRUE,
                 output_missing = TRUE,
                 plot_CV = TRUE,
                 output_CV = TRUE,
                 output_meta_summary = TRUE,
                 plot_well_median = TRUE

                 ){

  pp_median = get_moving_median(df_pp)
  plot_moving_pp = plot_moving_median(pp_median)
  plot_moving_pp =plot_moving_pp +
    theme(legend.position = "none")

  sample_median = get_moving_median(df_sample)
  readr::write_csv(sample_median, "sample_median.csv" )
  plot_moving_sample = plot_moving_median(sample_median)

  pp_sample_median = get_merged_moving_median(sample_median, pp_median)
  plot_moving_combined = get_median_plots_for_both(pp_sample_median)


  na_count = get_na_count(df_sample)
  plot_na = get_na_count_plot(na_count, df_sample)

  #get coefficient of variation
  cv_res = get_cv(df_pp)
  plot_cv = get_pp_cv_plot(df_pp)

  psych_summary = get_psych_summary(df_sample)


  pdf('well_median.pdf', width = 6, height = 9)
  egg::ggarrange(plot_moving_pp, plot_moving_sample, plot_moving_combined, heights = c(1.5,1.5, 1.5))
  dev.off()

  pdf('cv_and_missing.pdf', width = 8, height = 12)
  egg::ggarrange(plot_na, plot_cv, heights = c(1.5,1.5))
  dev.off()


  return(list(na_count = na_count,  cv_res = cv_res, summary_report = psych_summary, plot_list=list(
    plot_moving_pp, plot_moving_sample, plot_moving_combined, plot_na, plot_cv)))

}






get_na_count = function(df_sample){

  meta = names(df_sample)[grepl("[0-9]", names(df_sample))]

  sample_missing_df =df_sample[ names(df_sample) %in% c( meta)]


  na_count_samples  = as.data.frame(sapply(sample_missing_df,
                                           function(y) sum(is.na(y))))



  names(na_count_samples) ="NA_count"

  na_count = na_count_samples %>%
    dplyr::mutate(meta_name = rownames(.)) %>%
    dplyr::arrange(NA_count) %>%
    dplyr::mutate(rank= 1: dim(.)[1],
                  percentage = 100*round(NA_count/dim(sample_missing_df)[1], 3))

  return(na_count)

}




get_na_count_plot = function(na_count, df_sample){

  index5 = which.min(abs(na_count$percentage -5))
  index10 = which.min(abs(na_count$percentage -10))
  index20 = which.min(abs(na_count$percentage -20))
  index30 = which.min(abs(na_count$percentage -30))

  na_count_imp= na_count[c(index5, index10, index20, index30),]%>%
    dplyr::select(rank, percentage,NA_count)%>%
    dplyr::mutate(percentage2 = round(percentage, 0))


  p=na_count %>%
    ggplot2::ggplot(aes(rank, NA_count)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("missing count and missing percentage by metabolite")+
    ggplot2::xlab("metabolite's missing count rank")+

    theme_bw()+
    geom_vline(xintercept=na_count_imp$rank[[1]] , linetype = "dashed", color = "red") +
    geom_vline(xintercept=na_count_imp$rank[[2]] , linetype = "dashed", color = "red") +
    geom_vline(xintercept=na_count_imp$rank[[3]] , linetype = "dashed", color = "red") +
    geom_vline(xintercept=na_count_imp$rank[[4]] , linetype = "dashed", color = "red") +
    geom_hline(yintercept=na_count_imp$NA_count[[1]], linetype = "dashed", color = "red") +
    geom_hline(yintercept=na_count_imp$NA_count[[2]] , linetype = "dashed", color = "red")+
    geom_hline(yintercept=na_count_imp$NA_count[[3]] , linetype = "dashed", color = "red")+
    geom_hline(yintercept=na_count_imp$NA_count[[4]] , linetype = "dashed", color = "red")+

    ggplot2::scale_y_continuous(limits=c(-10, nrow(df_sample)+10), expand = c(0, 0), breaks =
                                  sort(c(seq(0, round (nrow(df_sample), digits=-2 ),
                                             length.out= 5),
                                         round(na_count_imp$NA_count[[1]],0),
                                         round(na_count_imp$NA_count[[2]],0),
                                         round(na_count_imp$NA_count[[3]],0),
                                         round(na_count_imp$NA_count[[4]],0))),
                                name = expression("missing count out of total"),
                                sec.axis = sec_axis(~ . * 100 / nrow(df_sample) ,
                                                    name = "missing percentage out of total",
                                                    breaks = sort(c(seq(0, 100, length.out=5),
                                                                    na_count_imp$percentage2[[1]],
                                                                    na_count_imp$percentage2[[2]],
                                                                    na_count_imp$percentage2[[3]],
                                                                    na_count_imp$percentage2[[4]])))

    )+
    ggplot2::scale_x_continuous(limits=c(-15, nrow(na_count)+15), expand = c(0, 0),
                                breaks = sort(c(seq(0, round (nrow(na_count), digits=-3 ),
                                                    length.out= 5),
                                                round(na_count_imp$rank[[1]],0),
                                                round(na_count_imp$rank[[2]],0 ),
                                                round(na_count_imp$rank[[3]],0),
                                                round(na_count_imp$rank[[4]],0))))+
    ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
                   axis.text.x = element_text(angle = 90, hjust = 1))


  return(p)


}




get_cv = function(df_pp){
  pp_meta_df = df_pp %>%
    dplyr::select( -plate_well )
  cv_res = as.data.frame(sapply(pp_meta_df, raster::cv, na.rm = TRUE))
  names(cv_res) = c("CV")
  cv_res = cv_res %>%
    tibble::rownames_to_column(var = "meta_name")
  return(cv_res)
}



get_pp_cv_plot = function(df_pp){

  pp_meta_df = df_pp %>%
    dplyr::select( -plate_well )


  cv_res = as.data.frame(sapply(pp_meta_df, raster::cv))
  names(cv_res) = c("CV")

  cv_res = cv_res[order(cv_res$CV),,drop=FALSE]

  cv_res= cv_res %>%
    dplyr::filter(!is.na(CV)) %>%
    dplyr::mutate( percentage = round((1:dim(.)[1])*100/dim(.)[1], 4)) %>%
    dplyr::filter( percentage < 99.5) %>%
    dplyr::mutate( percentage = round((1:dim(.)[1])*100/dim(.)[1], 4))

  cv_res = cv_res %>%
    dplyr::mutate(CV = round(CV, 1),
                  percentage = round(percentage, 1))

  p5<-round(100*sum(as.numeric(cv_res$CV<=5.0), na.rm = TRUE)/sum(!is.na(cv_res$CV), na.rm = TRUE), 0)
  p10<-round(100*sum(as.numeric(cv_res$CV<=10.0), na.rm = TRUE)/sum(!is.na(cv_res$CV), na.rm = TRUE), 0)
  p20<-round(100*sum(as.numeric(cv_res$CV<=20.0), na.rm = TRUE)/sum(!is.na(cv_res$CV), na.rm = TRUE), 0)
  p30<-round(100*sum(as.numeric(cv_res$CV<=30.0), na.rm = TRUE)/sum(!is.na(cv_res$CV), na.rm = TRUE), 0)

  xmax<-round(max(cv_res$CV, na.rm=T), 0)

  p<-cv_res %>%
    ggplot2::ggplot(aes(CV, percentage))+labs(x="CV %", y="cumulative % of metabolites") + xlim(-15, xmax+15)+
    ggplot2::geom_point() +
    theme_bw() +
    ggplot2::ggtitle("% of PP metabolites with CVs less than given threshold") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))+
    ggplot2::geom_segment(aes(x = -15, y = p5, xend = 5.0, yend = p5), linetype = "dashed", color = "blue") +
    ggplot2::geom_segment(aes(x = -15, y = p10, xend = 10.0, yend = p10), linetype = "dashed", color = "blue") +
    ggplot2::geom_segment(aes(x = -15, y = p20, xend = 20.0, yend = p20), linetype = "dashed", color = "blue") +
    ggplot2::geom_segment(aes(x = -15, y = p30, xend = 30.0, yend = p30), linetype = "dashed", color = "blue") +
    ggplot2::geom_segment(aes(x = 5.0, y = -7, xend = 5.0, yend = p5), linetype = "dashed", color = "red") +
    ggplot2::geom_segment(aes(x = 10.0, y = -7, xend = 10.0, yend = p10), linetype = "dashed", color = "red") +
    ggplot2::geom_segment(aes(x = 20.0, y = -7, xend = 20.0, yend = p20), linetype = "dashed", color = "red") +
    ggplot2::geom_segment(aes(x = 30.0, y = -7, xend = 30.0, yend = p30), linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    ggplot2::scale_y_continuous(limits=c(-7, 107), expand = c(0, 0), breaks = sort(c(seq(0, 100, length.out=5), p5, p10, p20, p30)))+
    ggplot2::scale_x_continuous(limits=c(-15, xmax+15), expand = c(0, 0), breaks = sort(c(seq(0, trunc(xmax/100)*100, length.out=5), 5, 10, 20, 30, trunc(xmax/100)*100, 0)))

  return(p)


}






get_imputated_sample = function(df_sample){

  #remove more than 20% missing
  df_sample <- df_sample[,colSums(is.na(df_sample))< nrow(df_sample)*.2]

  if("year" %in% names(df_sample)){

    df_sample = df_sample %>%
      tidyr::separate(plate_well, c("plate","well"), remove =FALSE) %>%
      dplyr::select(plate_well,  year, plate, well, everything()) %>%
      dplyr::mutate(plate = as.numeric(plate))

    # replace NA with 0.25 of min
    NA2mean <- function(x) replace(x, is.na(x), min(x, na.rm = TRUE)*0.25)
    df_sample[-c(1:5)][]<- lapply(df_sample[-c(1:5)], NA2mean)

  }else{

    df_sample = df_sample %>%
      tidyr::separate(plate_well, c("plate","well"), remove =FALSE) %>%
      dplyr::select(plate_well,  plate, well, everything()) %>%
      dplyr::mutate(plate = as.numeric(plate))

    # replace NA with 0.25 of min
    NA2mean <- function(x) replace(x, is.na(x), min(x, na.rm = TRUE)*0.25)
    df_sample[-c(1:4)][]<- lapply(df_sample[-c(1:4)], NA2mean)


  }

  return(df_sample)

}




get_moving_median = function(df){
  meta = names(df)[grepl("[0-9]", names(df))]

  df =df[ names(df) %in% c("plate_well", meta)]
  df= as.data.frame(df)

  row.names(df) =df$plate_well
  df$plate_well= NULL

  df = as.data.frame(t(df))

  row_median = as.data.frame( apply(df, 2,median, na.rm = TRUE))

  names(row_median) ="row_median"

  median_df = row_median %>%
    tibble::rownames_to_column(var ="plate_well") %>%
    tidyr::separate(plate_well, c("plate", "well")) %>%
    dplyr::mutate(plate = as.numeric(plate),
                  well = as.numeric(well)) %>%
    dplyr::mutate(order =(plate-1)*96+well )%>%
    dplyr::mutate(plate = as.factor(plate))

  return(median_df)


}

plot_moving_median = function(median_df){

  p =median_df %>%
    ggplot2::ggplot(aes(order, row_median, color = plate)) +
    ggplot2::geom_point()+
    ggplot2::ylab("median") +
    ggplot2::ggtitle("median vs analysis order")+
    ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = 0.5)) +
    theme_bw()

  return(p)
}



get_merged_moving_median =function(sample_median, pp_median){

  sample_median = sample_median %>%
    dplyr::mutate(group ="sample")

  pp_median = pp_median %>%
    dplyr::mutate(group ="pp")

  median_df = dplyr::bind_rows(sample_median,  pp_median)

  return(median_df)
}


get_median_plots_for_both = function(median_df){


  p =median_df %>%
    ggplot2::ggplot(aes(order, row_median, color = group)) +
    ggplot2::geom_point()+
    ggplot2::ylab("median") +
    ggplot2::ggtitle("median vs analysis order")+
    ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))+
    theme_bw()

  return(p)
}

get_psych_summary = function(df){
  df_sum = psych::describe(df, fast= FALSE)

  df_sum =df_sum %>%
    dplyr::add_rownames(var = "meta")

  meta_name = names(df)[grepl("[0-9]", names(df))]

  df_sum =df_sum %>%
    dplyr::filter( meta %in% meta_name ) %>%
    dplyr::mutate( NA_count = nrow(df) - n) %>%
    dplyr::select( -vars, - mad, -trimmed, -n) %>%
    dplyr::select(meta, NA_count, everything())

  return(df_sum)
}




