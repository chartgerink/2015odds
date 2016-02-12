if(!require(rplos)){install.packages('rplos')}

temp <- searchplos('rs12979860', fl = c('id'), limit = 1)

test_corpus <- searchplos('rs12979860', fl = c('id', 'results_and_discussion'), limit = temp$meta$numFound)

for(i in 1:dim(test_corpus$data)[1]){
  write.table(test_corpus$data$id[i],
              row.names = FALSE, col.names = FALSE,
              sprintf('data/test_corpus/test_corpus_id_%s',
                      ifelse(i < 10, paste0('00', i),
                             ifelse(i < 100, paste0('0', i),
                                    i))))
  write.table(test_corpus$data$results_and_discussion[i],
              row.names = FALSE, col.names = FALSE,
              sprintf('data/test_corpus/test_corpus_res_and_disc_%s',
                      ifelse(i < 10, paste0('00', i), 
                             ifelse(i < 100, paste0('0', i),
                                    i))))
}