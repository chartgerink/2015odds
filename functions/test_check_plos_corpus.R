source('functions/statcheck.r')

filename <- list.files(recursive = TRUE)[grep('test_corpus_res_and_disc_', list.files(recursive = TRUE))]

full <- data.frame(NULL)

for(j in 1:length(filename)){
  
  txt <- readChar(filename[j], file.info(filename[j])$size)
  
  loc <- gregexpr('rs[0-9]{1,10}', text = txt)
  
  pre <- NULL
  value <- NULL
  post <- NULL
  
  for(i in 1:length(loc[[1]])){
    value[i] <- substr(txt, start = loc[[1]][i], stop = loc[[1]][i] + attr(loc[[1]], "match.length")[1] - 1)
    # value[i] <- substr(temp, start = 1, stop = tail(gregexpr(pattern = "[0-9]", temp)[[1]], 1))
    
    pre[i] <- substr(txt, start = loc[[1]][i] - 100, stop = loc[[1]][i])
    
    post[i] <- substr(txt, start = loc[[1]][i] + attr(loc[[1]], "match.length")[1],
                      stop = loc[[1]][i] + attr(loc[[1]], "match.length")[1] + 100)
  }
  
  temp <- data.frame(file = filename[j],
                     value = as.character(value),
                     pre = pre,
                     post = post)
  
  full <- rbind(full, temp)
}



write.csv(full, "test.csv", row.names = FALSE)