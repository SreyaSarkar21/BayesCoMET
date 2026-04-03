#### Consider dynamic annotations of 200 songs (song id 1001-1200) ####

rm(list=ls())
library(tidyverse)

dynamicValence <- read.csv("~/DEAM_data_paper/annotations/annotations averaged per song/dynamic (per second annotations)/valence.csv",
                           header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
## Considering the first 200 songs for application
dynamicValence <- dynamicValence[dynamicValence$song_id > 1000 & dynamicValence$song_id <= 1200, -1]
# Extract numeric timestamps from colnames
timestamps <- as.numeric(gsub("sample_(\\d+)ms", "\\1", colnames(dynamicValence)))
keep_cols <- which(timestamps >= 15000 & timestamps <= 43500)
dynamicValence <- dynamicValence[, keep_cols]
time_sec <- (as.numeric(str_extract(colnames(dynamicValence), "\\d+")) / 1000)
colnames(dynamicValence) <- as.character(time_sec)

dynamicArousal <- read.csv("~/DEAM_data_paper/annotations/annotations averaged per song/dynamic (per second annotations)/arousal.csv",
                           header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
dynamicArousal <- dynamicArousal[dynamicArousal$song_id > 1000 & dynamicArousal$song_id <= 1200, -1]

# Extract numeric timestamps from colnames
timestamps <- as.numeric(gsub("sample_(\\d+)ms", "\\1", colnames(dynamicArousal)))
keep_cols <- which(timestamps >= 15000 & timestamps <= 43500)
dynamicArousal <- dynamicArousal[, keep_cols]
time_sec <- (as.numeric(str_extract(colnames(dynamicArousal), "\\d+")) / 1000)
colnames(dynamicArousal) <- as.character(time_sec)


avgOverTime <- function(dat, k, byrow = FALSE) {
    if(byrow) {
        nvar <- nrow(dat)
        # Spliting into consecutive blocks of size k
        idx <- split(seq_len(nvar),
                     ceiling(seq_len(nvar) / k))
        
        # Keeping only full blocks of size k
        idx <- idx[sapply(idx, length) == k]
        
        # Averaging within each block (song-wise)
        avg_dat <- sapply(idx, function(rows) {
            colMeans(dat[rows, , drop = FALSE])
        })
    } else {
        nvar <- ncol(dat)
        # Spliting into consecutive blocks of size k
        idx <- split(seq_len(nvar), ceiling(seq_len(nvar) / k))
        
        # Keeping only full blocks of size k
        idx <- idx[sapply(idx, length) == k]
        
        # Averaging within each block (song-wise)
        avg_dat <- sapply(idx, function(cols) {
            rowMeans(dat[, cols, drop = FALSE])
        })
    }
    
    list(avg_dat = as.data.frame(avg_dat), idx = idx)
    
}


### averaging valence over k consecutive timestamps
### k = 2 means averaging observations at t-th and (t+0.5)-th second timestamps
### t = 15, ..., 43
res <- avgOverTime(dynamicValence, k = 2, byrow = FALSE)
avg_valence_k2 <- res$avg_dat
# Labeling columns by mean of each block
new_cols <- lapply(res$idx, function(cols) mean(time_sec[cols]))
colnames(avg_valence_k2) <- as.character(new_cols)

### averaging arousal over k consecutive timestamps
res <- avgOverTime(dynamicArousal, k = 2, byrow = FALSE)
avg_arousal_k2 <- res$avg_dat
# Labeling columns by mean of each block
new_cols <- lapply(res$idx, function(cols) mean(time_sec[cols]))
colnames(avg_arousal_k2) <- as.character(new_cols)


n <- nrow(avg_arousal_k2); mis <- rep(ncol(avg_arousal_k2), n)
N <- sum(mis)
song_ids <- 1001:1200
y_all <- matrix(NA, N, 2)
y_all[, 1] <- as.vector(t(avg_valence_k2)) 
y_all[, 2] <- as.vector(t(avg_arousal_k2))
colnames(y_all) <- c("valence", "arousal")
yij_pc1 <- prcomp(y_all, scale. = TRUE)$x[, 1] ## first pc as response

############### 12 spectral features x their 4 summary averaged over k consecutive timestamps #################
Xijlist <- vector("list", N)
mis_cumsum <- cumsum(mis)
mis_starts <- c(1, mis_cumsum[-length(mis)] + 1)
for(i in seq_along(song_ids)) {
    feature_file <- readr::read_delim(paste0("~/DEAM_data_paper/features/", song_ids[i], ".csv"), delim = ";")
    feature_file <- feature_file[feature_file$frameTime >=15, ]
    feature_file <- column_to_rownames(feature_file, 'frameTime')
    select_feature_file <- feature_file %>% select(matches("spectral"))
    ordered_cols <- colnames(select_feature_file) %>% tibble(col = .) %>%
        mutate(feature = str_extract(col, "spectral[^_]+"),
               stat = case_when(str_detect(col, "_sma_amean$") ~ "amean",
                                str_detect(col, "_sma_de_amean$") ~ "de_amean",
                                str_detect(col, "_sma_stddev$") ~ "stddev",
                                str_detect(col, "_sma_de_stddev$") ~ "de_stddev"),
               stat = factor(stat, levels = c("amean", "de_amean", "stddev", "de_stddev"))) %>% arrange(feature, stat)
    select_feature_file <- select_feature_file[, ordered_cols$col]
    res <- avgOverTime(dat = select_feature_file, k = 2, byrow = TRUE) ## True because timestamps are rows
    avg_select_feature_file <- t(res$avg_dat)
    rownames(avg_select_feature_file) <- colnames(avg_arousal_k2)
    Xijlist[mis_starts[i]:mis_cumsum[i]] <- lapply(seq_len(nrow(avg_select_feature_file)), function(foo) {
        array(data = as.numeric(avg_select_feature_file[foo, ]),
              dim = c(12, 4),
              dimnames = list(feature = unique(ordered_cols$feature),
                              stat = c("amean", "de_amean", "stddev", "de_stddev")))
    })
}

dat_spectral_songitimej <- list(y = yij_pc1,
                                xlist = Xijlist,
                                zlist = Xijlist)

