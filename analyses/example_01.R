source('functions/statcheck.r')

txt <- "esults: In this Chinese elderly population, prevalence of overweight, central obesity, diabetes, dyslipidemia, hypertension,
and MetS were 48.3%, 71.0%, 32.4%, 75.7%, 68.3% and 54.5%, respectively. In the cross-sectional analyses, no SNP was
found to be associated with MetS. Genotype TT of SNP rs4402960 within the gene IGF2BP2 was associated with overweight
(odds ratio (OR) = 0.479, 95% confidence interval (CI): 0.316-0.724, p = 0.001) and genotype CA of SNP rs1801131 within the
gene MTHFR was associated with hypertension (OR = 1.560, 95% CI: 1.194–2.240, p = 0.001). However, these associations
were not observed in the longitudinal analyses"
x <- statcheck(txt)

write.csv(x, 'data/example_01.csv', row.names = FALSE)
