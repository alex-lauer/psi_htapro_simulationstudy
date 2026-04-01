number_of_events <- list()
for (i in 1:200) {
number_of_events[[i]] <-
sim_out[["AFTER"]][[i]][["stop"]] %>%
  group_by(arm) %>%
  summarise(sum = sum(PRO_status))
}

number_of_events_old <- list()
for (i in 1:200) {
  number_of_events_old[[i]] <-
    old_sim[["AFTER"]][[i]][["stop"]] %>%
    group_by(arm) %>%
    summarise(sum = sum(status, na.rm = TRUE))
}

a<- bind_rows(number_of_events)

b<- bind_rows(number_of_events_old)

cbind(a,b) %>% View()


number_of_events_full <- list()
for (i in 1:200) {
  number_of_events_full[[i]] <- sim_out[["AFTER"]][[i]][["full"]] %>%
    group_by(arm) %>%
    summarise(sum = sum(PRO_status))
}


c<- bind_rows(number_of_events_full)
cbind(c,a) %>% View()
