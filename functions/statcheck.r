# statcheck function adapted from the statcheck package by Michele Nuijten
# see github.com/michelenuijten/statcheck
statcheck <- function(
  x
){
  # Create empty objects
  raw_result = NULL
  genotype = NULL
  snp = NULL
  dv = NULL
  or_comparison = NULL
  or_result = NULL
  ci_confidence = NULL
  ci_lb = NULL
  ci_ub = NULL
  p_comparison = NULL
  p_result  =  NULL
  p_recalc  =  NULL
  test_stat_recalc <- NULL
  se_recalc = NULL
  
  if (length(x)==0) return(NULL)
  
  if (is.null(names(x))) names(x) <-  1:length(x)
  
  # create counter for loops
  index <- 1
  
  message("Extracting statistics...")
  pb <- txtProgressBar(max=length(x),style=3)
  for (i in 1:length(x)){
    
    txt <- x[i]
    
    # identify sequence of results and extract their text
    resLoc <- gregexpr("(genotype).*?rs[0-9]{1,10}.*?(associate|relate|correlate)[a-z]{1,}?\\s(with|to).*?(odds\\sratio|\\(?OR\\)?).*?[0-9]{2}\\%\\s(confidence\\sinterval|\\(?ci\\)?).*?(p.*?\\s?[0-9]?.[0-9]{1,5})",
                       txt,
                       ignore.case = TRUE)[[1]]
    resContext <- substring(txt,
                            resLoc,
                            resLoc + attr(resLoc, "match.length") + 3)
    
    # locate data in context strings
    # genotype
    locator_genotype <- gregexpr("genotype.?[A-Za-z]{2}", resContext, ignore.case = TRUE)
    
    # SNP
    locator_snp <- gregexpr("\\srs[0-9]{1,10}\\s", resContext, ignore.case = TRUE)
    
    # DV
    locator_ass <- gregexpr("(associat[a-z]{1,}|
                           correlat[a-z]{1,}|
                           relat[a-z]{1,})\\s(with|to)\\s", resContext, ignore.case = TRUE)
    
    locator_dv <- gregexpr("(associat[a-z]{1,}|
                           correlat[a-z]{1,}|
                           relat[a-z]{1,})\\s(with|to)\\s[a-z]{1,}", resContext, ignore.case = TRUE)
    
    # Odds ratio
    locator_or <- gregexpr("(odds ratio.*|or.*?).*?[<>=]", resContext, ignore.case = TRUE)
    
    # Confidence interval
    locator_ci <- gregexpr("[0-9]{2}\\%.*?(confidence interval.*|ci.*?)[=:;].*?([0-9]{1,}?[.]?[0-9]{1,}.*?[0-9]{1,}?[.]?[0-9]{1,})", resContext, ignore.case = TRUE)
    
    # P-value (if present)
    locator_p <- gregexpr("p.?[<>=]", resContext, ignore.case = TRUE)
    
    for (j in 1:length(resContext)){
      # extract data from each context with use of locators
      # genotype
      genotype_ind <- substr(resContext[[j]], locator_genotype[[j]] + 10 - 1, locator_genotype[[j]] + 10)
      # SNP
      snp_temp <- stringr::str_match_all(resContext[[j]],
                                "\\srs[0-9]{1,10}\\s")[[1]][1]
      snp_ind <- gsub(pattern = "\\s", x = snp_temp, replacement = "")
      # DV
      dv_ind <- gsub(pattern = "\\s",
                     x = substr(resContext[[j]], locator_ass[[j]] + attr(locator_ass[[j]], "match.length"), locator_dv[[j]] + attr(locator_dv[[j]], "match.length")),
                     "")
      # OR
      or_temp <- substr(resContext[[j]], locator_or[[j]], locator_or[[j]] + attr(locator_or[[j]], "match.length") + 5)
      or_comp <- str_match_all(or_temp,
                               "[<>=]")[[1]][1]
      or_ind <- as.numeric(str_match_all(or_temp,
                              "[0-9]{1,}?[.][0-9]{1,3}")[[1]][1])
      # CI
      ci_temp <- substr(resContext[[j]], locator_ci[[j]], locator_ci[[j]] + attr(locator_ci[[j]], "match.length") + 5)
      ci_conf <- as.numeric(substr(str_match_all(ci_temp, "[0-9]{1,2}\\%\\s?[cC]")[[1]][,1], 0, 2))
      ci_conc <- str_match_all(ci_temp, "[0-9]{1,}?[.][0-9]{1,3}")[[1]][,1]
      ci_lb_ind <- as.numeric(min(ci_conc))
      ci_ub_ind <- as.numeric(max(ci_conc))
      # P-value
      p_temp <- substr(resContext[[j]], locator_p[[j]], locator_p[[j]] + attr(locator_p[[j]], "match.length") + 5)
      p_comp <- str_match_all(p_temp,
                              "[<>=]")[[1]][1]
      p_ind <- as.numeric(str_match_all(p_temp,
                             "[0-9]{1,}?[.][0-9]{1,3}")[[1]][1])
      
      # recalculate p-value
      # from Altman, D. G., & Bland, J. M. (2011). How to obtain the P value from a confidence interval. BMJ , 343, d2304.
      se_ind <- (ci_ub_ind - ci_lb_ind) / (2 * qnorm((1 - (ci_conf / 100)) / 2, lower.tail = FALSE))
      z_ind <- or_ind / se_ind
      p_recalc_ind <- pnorm(z_ind, lower.tail = FALSE) * 2
      
      
      # write back results into main objects
      raw_result[index] = resContext[[j]]
      genotype[index] = genotype_ind
      snp[index] = snp_ind
      dv[index] = dv_ind
      or_comparison[index] = or_comp
      or_result[index] = or_ind
      ci_confidence[index] = ci_conf
      ci_lb[index] = ci_lb_ind
      ci_ub[index] = ci_ub_ind
      p_comparison[index] = p_comp
      p_result [index] =  p_ind
      p_recalc [index] =  p_recalc_ind
      test_stat_recalc[index] = z_ind
      se_recalc[index] = se_ind
      
      index = index + 1
    }
  }
  
  
  # data frame
  Res <- data.frame(raw_result = raw_result,
                    genotype = genotype,
                    snp = snp,
                    dv = dv,
                    or_comparison = or_comparison,
                    or_result = or_result,
                    ci_confidence = ci_confidence,
                    ci_lb = ci_lb,
                    ci_ub = ci_ub,
                    p_comparison = p_comparison,
                    p_result  = p_result,
                    p_recalc  = p_recalc,
                    test_stat_recalc = test_stat_recalc,
                    se_recalc = se_recalc)

class(Res) <- c("statcheck","data.frame")

###--------------------------------------------------------------------- 
# Return message when there are no results
if(nrow(Res)>0){
  
  
  return(Res) 
} else {
  Res <- cat("statcheck did not find any results\n")
}
}