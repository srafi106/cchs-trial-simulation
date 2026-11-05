# Reading and recoding of the CCHS 2022 PUMF.
# Variable codes follow the 2022 PUMF data dictionary.

read_pumf = function(path) {
  as.data.frame(data.table::fread(path))
}

recode_pumf = function(d) {
  d$cvd = ifelse(d$CCCDGCAR == 1, 1L, ifelse(d$CCCDGCAR == 2, 0L, NA_integer_))
  d$age_grp = factor(d$DHHGAGE, levels = 3:5,
                     labels = c("A35to49", "A50to64", "A65plus"))
  d$sex = factor(d$DHH_SEX, levels = c(1, 2), labels = c("Male", "Female"))
  d$inc_q = factor(d$INCDGRCA, levels = 1:5, labels = paste0("Q", 1:5))
  d$hr = factor(d$GEODGHR4)
  d$prov = factor(d$GEOGPRV)
  d$analytic = !is.na(d$DHHGAGE) & d$DHHGAGE >= 3 & !is.na(d$cvd) &
    !is.na(d$inc_q) & !is.na(d$sex)
  d
}

# Respondents aged 35 or older with a valid outcome, sex and income quintile
analytic_sample = function(d) {
  dat = d[d$analytic, ]
  dat$hr = droplevels(dat$hr)
  dat$prov = droplevels(dat$prov)
  dat
}
