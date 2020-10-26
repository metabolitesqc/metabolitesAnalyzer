#' Metabolites analysis
#'
#' This function doing metabolites vs treatment association, metabolites vs clinical variable association.
#' It assumes that the first column of df_meta
#' contains the subjectId, with or without year, plate_well in the second and third column.  in the subsequent columns
#' are metabolites measurements. It assumes the first column of df_clinical contains subjectId for patients invovled
#' in the study, it then contains year if it has measurements in two years. trt_var is a character, showign the treatment variable, and clinical_var is a character vector showing
#' clinical variable used to plot the rainplot to show the clinical varialbe and metabolites association.
#' adjust_var shows the clinical variables used in the model for adjustment, it is a character vector.
#
#' @param df_meta, experimental metabolites data
#' @param df_clinical, clinical data
#' @param trt_var, treatment variable
#' @param clinical_var, clinical variables in the rainplot
#' @param adjustment_var, adjustment variable in the linear model
#' @param no_signal_correction, defaults to true for data with two years
#' @param plot_volcano, defaults to true
#' @param output_trt_association, defaults to true
#' @param output_clin_rainplot, defaults to true
#'
#' @return list showing linear model output
#' @export
metaA = function(df_meta,
                 df_clinical,
                 trt_var,
                 clinical_var,  #in rainplot
                 adjustment_var =c(),
                 no_signal_correction = TRUE, #means it will use delta
                 plot_volcano = TRUE,
                 output_trt_association = TRUE,
                 output_clin_rainplot = TRUE){

  linear_model_res = get_linear_model_res( df_meta,
                                           df_clinical,
                                           trt_var,
                                           adjustment_var = c(),
                                           no_signal_correction)

  waterfall_plots = get_waterfall_plots(linear_model_res)
  ggsave("walterfall_plots.pdf", width =18, height =18)


  if(output_clin_rainplot == TRUE){

    rplot = get_rainplot( linear_model_res,
                          df_meta,
                          df_clinical, #either bl or with both years
                          clinical_var,
                          no_signal_correction )
    ggsave("rplot.pdf", width =18, height =18)
  }


  if(no_signal_correction == TRUE & plot_volcano == TRUE){
    df_fold = get_fold_change(df_meta, df_clinical, trt_var)
    volcano_plot = get_volcano_plots(linear_model_res,df_fold)
    ggsave("volcano_plot.pdf", width=13, height =10,limitsize = FALSE)
  }

  return( list(linear_model_res =linear_model_res))
}



get_linear_model_res = function(df_meta,
                                 df_clinical,
                                 trt_var,
                                 adjustment_var =c(),
                                 no_signal_correction){

  if(no_signal_correction == TRUE){
    df_meta = get_sample_delta(df_meta)
  }else{
    df_meta = get_combat_corrected_data(df_meta)
  }


  df = df_meta %>%
    dplyr::inner_join(df_clinical, by = "subjectId")

  meta = names(df)[grepl("[0-9]", names(df))]

  if(  length(adjustment_var) == 0){

    df = df[ names(df) %in% c(meta, trt_var) ]
    variables = c(trt_var)

    f = as.formula(
      paste("meta_reading",
            paste(variables),
            sep = " ~ "))

    df_long = df %>%
      dplyr::select(trt_var, everything()) %>%
      tidyr::gather(key="meta", value="meta_reading", -variables) %>%
      tidyr::nest(data = c(variables, meta_reading))

  }else{
    adjustment_var = unlist(stringr::str_split(adjustment_var, " "))

    df =df[ names(df) %in% c(meta, trt_var, adjustment_var) ]
    variables = c(trt_var, adjustment_var)

    f = as.formula(
      paste("meta_reading",
            paste(variables, collapse = " + "),
            sep = " ~ "))

    df_long = df %>%
      dplyr::select(trt_var, adjustment_var, everything()) %>%
      tidyr::gather(key="meta", value="meta_reading", -variables) %>%
      tidyr::nest(data = c(variables, meta_reading))
  }

  f2 = as.formula(
    paste("meta_reading",
          paste(trt_var, collapse = " + "),
          sep = " ~ ")
  )

  lm_report = df_long%>%
    dplyr::mutate( model = purrr::map(.$data, ~lm(f, data =.x))) %>%
    dplyr::mutate(glance = purrr::map(model, broom::tidy, conf.int = TRUE, conf.level = 0.95)) %>%
    dplyr::select(meta, glance) %>%
    tidyr::unnest(cols= c(glance)) %>%
    dplyr::filter(term == trt_var) %>%
    dplyr::select(meta, term, estimate,conf.low, conf.high, p.value) %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate( estimate = round(estimate, 3),
                   conf.low = round(conf.low, 3),
                   conf.high = round(conf.high, 3))

  if(length( unique(df_clinical[[trt_var]])) ==2 ){
    wilcox_report = df_long %>%
      dplyr::mutate( test = purrr::map(.$data, ~wilcox.test(f2, data =.x ))) %>%
      dplyr::mutate( wilcox_p_value = purrr::map_dbl(.$test, ~.x$p.value)) %>%
      dplyr::mutate( p.adjust =p.adjust(wilcox_p_value, method = "fdr") ) %>%
      dplyr::select(meta, wilcox_p_value, p.adjust)

    report = wilcox_report %>%
      dplyr::left_join(lm_report, by ="meta") %>%
      dplyr::arrange( p.adjust)

    readr::write_csv( report, "lm_report.csv")

    return(report)
  }else{

    readr::write_csv( lm_report, "lm_report.csv")
    return(lm_report)

  }
}





get_waterfall_plots = function(linear_model_res){


  index_max = which.max(linear_model_res$estimate)
  linear_model_res = linear_model_res[-index_max,]

  df_w = linear_model_res %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate( legend = ifelse(estimate >=0, "promoted","suppressed")) %>%
    dplyr::mutate( legend = as.factor(legend)) %>%
    dplyr::mutate(log10_p = -log10(wilcox_p_value)) %>%
    dplyr::mutate( log10_p = round(log10_p,2))

  scale.min=0
  scale.max= max(df_w$log10_p, na.rm = T)
  middle.pt = 1.3

  by_scale = round(scale.max/8,2)

  p= df_w %>%
    ggplot(aes( reorder(meta, -estimate ), estimate, fill =log10_p))+
    geom_bar(stat ="identity", color=NA )+
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_fill_gradient2( low = "black", mid = "white",
                          high = "red", midpoint = middle.pt,
                          breaks=seq(from=scale.min, to=scale.max, by=by_scale),
                          name = "p.value")+
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    ylab("coefficients") +
    xlab("Individual metabolite") +
    ggtitle("waterfall plot for treatment") +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

  return(p)
}



get_volcano_plots = function(linear_model_res,
                             df_fold
){

  df_p = linear_model_res %>%
    dplyr::select(meta,  wilcox_p_value, p.adjust)

  df_vol = df_fold %>%
    dplyr::left_join(df_p, by ="meta") %>%
    dplyr::select(wilcox_p_value, p.adjust, fold_change) %>%
    dplyr::mutate(log2_fold = log2(fold_change),
                  `-log10_p` = -log10(wilcox_p_value)) %>%
    dplyr::mutate(size = round(`-log10_p`,0)) %>%
    dplyr::mutate(fdr = ifelse(p.adjust<0.05, 1, 0) ) %>%
    dplyr::mutate( fdr = as.factor(fdr))

  vsize = expression( paste("-", log[10], " p-value"))
  vx = expression(paste(log[2]," fold change Size of rectangles is proportional to -", log[10]," p-value"))
  vy = expression( paste("-", log[10], " p-value"))

  p= df_vol%>%
    ggplot(aes(log2_fold, `-log10_p`, shape = fdr) )+
    geom_point(aes(size = size))+
    xlab(vx)+
    xlim(-1,1.8) +
    ylab(vy)+
    geom_hline(yintercept=1.3 ,linetype="dotted") +
    geom_vline(xintercept=log2(1.1487),linetype="dotted") +
    geom_vline(xintercept=-log2(1.1487) , linetype="dotted")+
    theme_bw()+
    theme(
      legend.title = element_text(color = "black", size = 17),
      legend.text = element_text(color = "black", size = 17),
      axis.text.x = element_text(colour="black"),
      axis.text.y = element_text(colour="black"),
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17)


    )+
    labs(shape="FDR significance",  size= vsize)
  return(p)
}



get_fold_change = function(df_meta, df_clinical, trt_var){

  df_imputated_sample = get_imputated_sample(df_meta)

  meta = names(df_imputated_sample)[grepl("[0-9]", names(df_imputated_sample))]

  df_imputated_sample = df_imputated_sample[ names(df_imputated_sample) %in%
                                               c("subjectId","year", meta)]

  df_imputated_sample = df_imputated_sample %>%
    dplyr::arrange(subjectId) %>%
    dplyr::add_count(subjectId) %>%
    dplyr::filter(n==2) %>%
    dplyr::select(-n)

  df_clinical = df_clinical[ names(df_clinical)  %in% c("subjectId", trt_var)] %>%
    dplyr::distinct()

  df_fish_trt_sub = df_clinical %>%
    dplyr::filter( !!as.symbol(trt_var) ==1) %>%
    dplyr::pull(subjectId)

  df_fish_ctl_sub = df_clinical %>%
    dplyr::filter( !!as.symbol(trt_var)  ==0) %>%
    dplyr::pull(subjectId)

  df_bl = df_imputated_sample %>%
    dplyr::filter( year ==0) %>%
    dplyr::select(-year)

  df_yr1 = df_imputated_sample %>%
    dplyr::filter( year ==1)%>%
    dplyr::select(-year)

  df_year1_trt = df_yr1 %>%
    dplyr::filter(subjectId %in% df_fish_trt_sub)

  df_year1_no_trt = df_yr1 %>%
    dplyr::filter(subjectId %in% df_fish_ctl_sub)

  df_bl_trt = df_bl %>%
    dplyr::filter(subjectId %in% df_fish_trt_sub)

  df_bl_no_trt =df_bl %>%
    dplyr::filter(subjectId %in% df_fish_ctl_sub)

  #calculate (df_year1_trt/df_bl_trt)/(df_year1_no_trt/df_bl_no_trt)
  a=(df_year1_trt[-c(1,2)]/df_bl_trt[-c(1,2)])

  b=(df_year1_no_trt[-c(1,2)]/df_bl_no_trt[-c(1,2)])

  c= as.data.frame (sapply(a, mean))

  d = as.data.frame (sapply(b, mean))

  e = c/d

  names(e) = "fold_change"

  fold_df = e %>%
    tibble::rownames_to_column( var = "meta") %>%
    dplyr::arrange(-fold_change)

  return(fold_df)

}






get_imputated_sample= function(df_meta){

  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]

  #remove more than 20% missing
  df_data <- df_data[,colSums(is.na(df_data))< nrow(df_data)*.2]
  df_data[] <- lapply(df_data, NA2mean)

  res = dplyr::bind_cols( df_info, df_data)
  return(res)

}


#log transform, centered to median, sd =1
get_normalized_data = function(df_meta){

  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]

  #remove more than 20% missing
  df_data <- df_data[,colSums(is.na(df_data))< nrow(df_data)*.2]

  df_data = as.data.frame(apply( df_data,
                                 2,
                                 function(y) ( log(y) - median(log(y), na.rm = T))/ sd( log(y), na.rm = T)))

  res = dplyr::bind_cols( df_info, df_data)
  return(res)

}


get_normalized_data_no_log = function(df_meta){

  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]

  #remove more than 20% missing
  df_data <- df_data[,colSums(is.na(df_data))< nrow(df_data)*.2]

  df_data = as.data.frame(apply( df_data,
                                 2,
                                 function(y) ( y - median(y, na.rm = T))/ sd( y, na.rm = T)))

  res = dplyr::bind_cols( df_info, df_data)
  return(res)

}



get_sample_delta = function(df_meta){

  df_meta = get_imputated_sample(df_meta)

  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]

  df_info = df_info %>%
    dplyr::select(subjectId, year)

  df_meta = dplyr::bind_cols( df_info, df_data)

  if("year" %in% names(df_meta)){

    df_for_test = df_meta %>%
      dplyr::add_count(subjectId) %>%
      dplyr::filter(n==2) %>%
      dplyr::select(-n) %>%
      dplyr::arrange(subjectId)

    df_bl = df_for_test %>%
      dplyr::filter(year ==0 ) %>%
      dplyr::select(-year)

    df_y1 = df_for_test %>%
      dplyr::filter(year ==1 ) %>%
      dplyr::select(-year)

    df_bl[-1] = log(df_bl[-1] )
    df_y1[-1] = log( df_y1[-1])

    df_delta = df_bl
    df_delta[-1] = df_y1[-1] - df_bl[-1]

    df_delta = get_normalized_data_no_log(df_delta)

  }else(
    stop()
  )

  return(df_delta)
}




get_combat_corrected_data= function(df_meta){

  df_meta = get_imputated_sample_1(df_meta)
  df_meta = get_normalized_data(df_meta)

  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]
  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  df_info = df_info %>%
    tidyr::separate( plate_well, c("plate", "well"), sep="_", remove = FALSE) %>%
    dplyr::mutate(plate = as.numeric(plate))

  batch = df_info$plate

  edata = df_meta[names(df_meta) %in% c('subjectId', names(df_data))]

  rownames(edata) = edata$subjectId
  edata$subjectId = NULL

  edata = t(edata)

  combat_edata = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  combat_edata = as.data.frame(t(combat_edata))

  res = dplyr::bind_cols(df_info, combat_edata)

  res =res %>%
    dplyr::select(-plate, -well, -plate_well) %>%
    dplyr::select(subjectId, everything())

  return(res)

}


NA2mean <- function(x) replace(x, is.na(x), min(x, na.rm = TRUE)*0.25)

get_rainplot= function(linear_model_res,
                       df_meta,
                       df_clinical, #either bl or both years
                       clinical_var,
                       no_signal_correction = TRUE){

  top_30_meta = linear_model_res %>%
    dplyr::arrange(wilcox_p_value) %>%
    dplyr::slice_min( wilcox_p_value, n=30) %>%
    dplyr::pull(meta)

  if(no_signal_correction == TRUE){
    df_meta = get_sample_delta(df_meta)
    readr::write_csv(df_meta, "df_meta_from_get_rainplot_function.csv")
  }else{
    df_meta = get_combat_corrected_data(df_meta)
  }

  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]
  df_data = df_meta[ top_30_meta ]

  meta_name = names(df_data)

  df_meta = dplyr::bind_cols(df_info, df_data)
  df_clinical = get_clinical_df_for_rplt( df_clinical, clinical_var)

  plot_data = get_clin_association( df_clinical, df_meta, top_30_meta, clinical_var)

  p = produce_rplot(plot_data)

  return(p)
}


get_clinical_df_for_rplt = function(df_clinical,
                                    clinical_var,
                                    no_signal_correction = TRUE
                                    ){


  if(no_signal_correction == TRUE){
    df_clinical = df_clinical[names(df_clinical) %in% c("subjectId","year", clinical_var)]

    df_clinical =df_clinical %>%
      dplyr::arrange(subjectId)

    clin_y1 = df_clinical %>%
      dplyr::filter(year ==1) %>%
      dplyr::select( - year)

    clin_bl = df_clinical %>%
      dplyr::filter(year ==0) %>%
      dplyr::select( - year)

    delta_clinical = clin_y1 #temp assign, will adjust below
    delta_clinical[-1] = clin_y1[-1] - clin_bl[-1]

    df_clinical = delta_clinical

  }else{
    df_clinical = df_clinical[names(df_clinical) %in% c("subjectId", clinical_var)]
  }


  df_clinical[c(2: ncol(df_clinical))] = as.data.frame(apply( df_clinical[c(2: ncol(df_clinical))],
                                                                    2,
                                                                    function(y) (y - mean(y, na.rm=TRUE)) / sd(y, na.rm = TRUE) ))
  return(df_clinical)

}





get_clin_association = function(df_clinical, df_meta, top_30_meta, clinical_var){

  df_data = df_meta[, top_30_meta ]


  meta_name = names(df_data)
  name1 = meta_name[1]
  name2 = meta_name[ length(meta_name)]

  clin_name1 = clinical_var[1]
  clin_name2 = clinical_var[ length(clinical_var)]

  df_clin = df_clinical[ c("subjectId",clinical_var )]
  df_data = df_meta[ c("subjectId", top_30_meta)]


  df = df_data %>%
    dplyr::inner_join (df_clin, by ="subjectId")


  nested_long = df %>%
    tidyr::gather( name1:name2, key ="meta", value="meta_reading" ) %>%
    tidyr::gather(clin_name1:clin_name2, key = "clin", value="clinical_value") %>%
    dplyr::group_by( meta, clin) %>%
    tidyr::nest()

  lm_res=nested_long %>%
    dplyr::mutate( model = purrr::map(data, clin_model )) %>%
    dplyr::mutate( glance = purrr::map(model, broom::tidy)) %>%
    dplyr::select(clin, meta, glance) %>%
    tidyr::unnest(cols= c(glance)) %>%
    dplyr::filter(term =="meta_reading") %>%
    dplyr::select(meta, clin, estimate, p.value)

  return(lm_res)
}


clin_model = function(df){

  df = df %>%
    dplyr::select(-subjectId)
  lm(clinical_value~meta_reading , data= df)
}


produce_rplot = function(plot_data){

  plot_data = plot_data %>%
    dplyr::rename(response = clin,
                  term = meta)%>%
    dplyr::mutate(estimate = round(estimate,2))

  # term as mz and rt --------------------------------------------------
  plot_data = plot_data %>%
    tidyr::separate(term, c("a", "b","c"), sep="_") %>%
    dplyr::mutate(a = as.numeric(a),
                  b = as.numeric(b)) %>%
    dplyr::mutate(a = round(a, 4),
                  b = round(b, 2)) %>%
    dplyr::mutate( term = paste(a, b, sep="  ")) %>%
    dplyr::select(-a, -b, -c)

  # P-Value transformation --------------------------------------------------
  plot_data <-
    plot_data %>%
    dplyr::mutate(p.value = -1 * log10(p.value))

  # Calculate symmetric limits based on most extreme value
  max_abs_estimate <- max(abs(plot_data$estimate))
  max_lim <- max_abs_estimate
  min_lim = -1 * max_lim


  # Ordering by Cluster -----------------------------------------------------

  # Convert to matrix and reshape for clustering.
  cluster_data <-
    plot_data %>%
    dplyr::select(response, term, estimate) %>%
    tidyr::spread(response, estimate)

  rnms <-
    cluster_data$term

  cluster_data <-
    cluster_data %>%
    dplyr::select(-term) %>%
    as.matrix()

  rownames(cluster_data) <- rnms

  # cluster dependent variable terms
  clust <- hclust(dist(cluster_data), method = 'ward.D2')

  # `clust$order` orders `term` into clusters
  term_order <-
    clust$labels[clust$order]

  # Convert term to a factor, ordered by `term_order`
  plot_data_clo <-
    plot_data %>%
    dplyr::mutate(term = factor(term, levels = term_order))

  # Theme + Palette ---------------------------------------------------------
  ## Palette

  palette <-
    # Blue
    c("#053061",
      "#313695",
      "#4575b4",
      "#74add1",
      "#abd9e9",
      "#e0f3f8",
      "#fee090",
      "#fdae61",
      "#f46d43",
      "#d73027",
      "#a50026",
      '#67001f')
  # Red

  ## theme
  thm <-
    # Good starting theme + set text size
    theme_light(base_size = 18) +
    theme(
      # Remove axis ticks and titles
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),

      # Remove gridlines and boxes
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      legend.key = element_blank(),

      # White backgrounds
      panel.background = element_rect(fill = 'white'),
      plot.background = element_rect(fill = 'white'),
      legend.background = element_rect(fill = 'white'),

      # Angle text
      axis.text.x.top  = element_text(angle = 45, hjust = 0)
    )

  rainplot <-
    # Use cluter ordered data
    ggplot(plot_data_clo) +
    geom_point(aes(x = response, y = term, colour = estimate, size = p.value)) +
    scale_x_discrete(position = 'top') +
    scale_size_area(expression(paste(-log[10]('P-value'))), max_size = 12) +
    scale_color_gradientn(
      'Effect Size Estimate',
      colors = palette,
      limits = c(min_lim, max_lim),
      breaks = c(min_lim, min_lim / 2, 0 , max_lim / 2, max_lim)
    ) +
    thm

}

