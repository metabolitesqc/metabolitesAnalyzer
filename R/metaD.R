#' Creat demographic tables
#'
#' This function create a demographic table. It assumes that the first column
#' contains the subjectId, the subsequent columns are clinical variables:diabetes,currentsmoker,africanamerican,
#' statintreatment,women,hbpmed, totcvd,stroke, statin, chd, mi, totalchol,trig,hdlc,ldlc, age, bmi, hscrp, tnfr2,
#' il6, glyca17, sicam
#'
#' @param df,  input clinical data frame,
#' @param trt_var,  treatment variable
#' @return A data farme showing summary for clinical variables for subjects invovled in the study
#' @export
metaG = function(df, trt_var){

  names(df) = stringr::str_replace_all(names(df), "_", "")
  names(df) = stringr::str_to_lower(names(df))

  if(trt_var %in% names(df)){


    df_placebo =df %>%
      dplyr::filter( (!!sym(trt_var))==0)

    df_trt =df %>%
      dplyr::filter( (!!sym(trt_var))==1)


    df_total = df

    cat_trt = get_cat_var_summary(df_trt, "Treatment/case")
    cat_placebo = get_cat_var_summary(df_placebo, "Placebo/control")
    cat_total = get_cat_var_summary(df_total, "Total")


    con_trt = get_con_var_summary(df_trt, "Treatment/case")
    con_placebo = get_con_var_summary(df_placebo, "Placebo/control")
    con_total = get_con_var_summary(df_total, "Total")


    trt_df = dplyr::bind_rows( cat_trt, con_trt)
    placebo_df = dplyr::bind_rows(cat_placebo, con_placebo)
    total_df = dplyr::bind_rows( cat_total, con_total )

    final_res = trt_df %>%
      dplyr::left_join(placebo_df) %>%
      dplyr::left_join(total_df)

    names(final_res) = c("clin", paste("Treatment/case",nrow(df_trt)), paste("Placebo/control", nrow(df_placebo)) , paste("Total", nrow(df_total)) )
    return(final_res)
  }
  else{
    stop("Error, make sure the treatment variable exists in your clinical dataset")
  }
}


cat_var = c("diabetes","currentsmoker","africanamerican","statintreatment", "women", "hbpmed", "totcvd", "stroke", "statin", "chd", "mi")
con_var = c("totalchol","trig","hdlc","ldlc","age","bmi", "hscrp","tnfr2", "il6" ,"glyca17", "sicam" )


get_cat_var_summary = function(df, name){

  cat_var_df = as.data.frame(cat_var,
                             stringsAsFactors = FALSE)

  names(cat_var_df) = "cat_var"

  cat_var = cat_var_df %>%
    dplyr::mutate(cat_flag = ifelse(cat_var  %in% names(df), 1, 0)) %>%
    dplyr::filter(cat_flag ==1) %>%
    dplyr::pull(cat_var)


  cat_df = df %>%
    dplyr::select(subjectid, cat_var)

  res = as.data.frame( sapply(cat_df[-1], mean, na.rm =T))
  names(res) = "mean"
  mean_res = res %>%
    tibble::rownames_to_column(var ="clin") %>%
    dplyr::mutate( mean = round(mean*100,1) )

  res = as.data.frame( sapply(cat_df[-1], sum, na.rm =T))
  names(res) = "sum"
  sum_res = res %>%
    dplyr::add_rownames(var ="clin")

  results = mean_res %>%
    dplyr::left_join(sum_res, by ="clin") %>%
    dplyr::mutate( col_name = paste(sum, " (", mean, "%)", sep="")) %>%
    dplyr::select(-mean, -sum)

  names(results)[2] = name

  return(results)
}



get_con_var_summary = function(df, name){

  # for continious var
  con_var_df = as.data.frame(con_var,
                             stringsAsFactors = FALSE)

  names(con_var_df) = "con_var"

  con_var = con_var_df %>%
    dplyr::mutate(con_flag = ifelse(con_var  %in% names(df), 1, 0)) %>%
    dplyr::filter(con_flag ==1) %>%
    dplyr::pull(con_var)

  con_df = df %>%
    dplyr::select(subjectid, con_var)

  res = as.data.frame( sapply(con_df[-1], median, na.rm =T))
  names(res) = "median"
  median_res = res %>%
    dplyr::add_rownames(var ="clin") %>%
    dplyr::mutate(median = round(median,1))

  res = as.data.frame( sapply(con_df[-1], quantile, 0.25, na.rm =T))
  names(res) = "low"
  low_res = res %>%
    dplyr::mutate(low = round(low,1))
  rownames(low_res) = rownames(median_res)

  res = as.data.frame( sapply(con_df[-1], quantile, 0.75, na.rm =T))
  names(res) = "high"
  high_res = res %>%
    dplyr::mutate(high = round(high,1))
  rownames(high_res) = rownames(median_res)

  results = bind_cols(median_res, low_res, high_res) %>%
    dplyr::mutate( col_name = paste(median, " (", low,"-", high,  ")", sep="")) %>%
    dplyr::select(-median, -low, -high )

  names(results)[2] = name

  return(results)
}



