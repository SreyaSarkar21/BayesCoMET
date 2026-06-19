rm(list=ls())
library(readr)
library(tidyverse)
lfw_attributes <- readr::read_csv(file = "~/lfw_attributes.csv")
lfw_attributes$folder_name <- gsub(" ", "_", trimws(lfw_attributes$person))
#View(lfw_attributes)
unique_folders <- unique(lfw_attributes$folder_name)
base_dir <- "~/lfw-deepfunneled"
valid_folder_indx <- unname(which(sapply(file.path(base_dir, unique_folders),
                                         function(f) {length(list.files(f)) >= 4 & length(list.files(f)) <= 10})))

rowindx <- which(lfw_attributes$folder_name %in% unique_folders[valid_folder_indx])

lfw_attributes <- lfw_attributes[rowindx, ] ## attributes with at least 4 and at most 10 images per individual

attr_counts <- table(lfw_attributes$folder_name)
attr_df <- as.data.frame(attr_counts)
names(attr_df) <- c("folder_name", "n_csv")
img_df <- data.frame(
    folder_name = unique(lfw_attributes$folder_name),
    n_img = sapply(
        unique(lfw_attributes$folder_name),
        function(fn) length(list.files(file.path(base_dir, fn)))
    )
)
merged <- left_join(attr_df, img_df, by = "folder_name")
merged$match <- merged$n_csv == merged$n_img
table(merged$match)
#merged[merged$match == FALSE, ]

valid_folders <- merged$folder_name[merged$match]
lfw_attributes_clean <- lfw_attributes %>%
    filter(folder_name %in% valid_folders)

yall <- as.matrix(lfw_attributes_clean[, 3:75])
yijs <- prcomp(yall, scale. = TRUE)$x[, 1] ## first pc as response
n <- length(unique(lfw_attributes_clean$folder_name))
mi <- as.numeric(table(lfw_attributes_clean$folder_name))

library(imager)
library(dplyr)
library(stringr)
lfw_data_dir = "~/lfw-deepfunneled"
resized_images_save_dir = "~/resized_images_lfw_3D_mi4and10"

ppl_names = unique(lfw_attributes_clean$folder_name)

for (person_px in ppl_names) {
    from_dir_px = paste0(lfw_data_dir, "/", person_px)
    to_dir_px = paste0(resized_images_save_dir, "/", person_px)

    if (!dir.exists(to_dir_px)) {
        dir.create(to_dir_px)
        stopifnot(dir.exists(to_dir_px))
    }

    img_filenames_px = list.files(from_dir_px)
    n_images_px = length(img_filenames_px)

    for (ix in 1:n_images_px) {
        img <- imager::load.image(paste0(from_dir_px, "/", img_filenames_px[ix]))
        saveRDS(as.array(imager::resize(img, size_x=32, size_y=32))[,,1,], file = sprintf(
            "%s/%s.rds",
            to_dir_px,
            tools::file_path_sans_ext(img_filenames_px[ix])))
    }
}

# Load all .rds file paths
resized_images_save_dir = "~/resized_images_lfw_3D_mi4and10/"
rds_files <- list.files(resized_images_save_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

# Get the folder (person) name for each file
folders_X <- basename(dirname(rds_files))

# Split file paths by person
rds_split <- split(rds_files, folders_X)

# Match the number of images expected per person (from metadata)
names(mi) <- names(table(lfw_attributes_clean$folder_name))

# For each person, keep the first mi[person] files
rds_valid <- unlist(mapply(function(files, m) head(files, m),
                           rds_split[names(mi)], mi, SIMPLIFY = FALSE))

Xijlist <- lapply(rds_valid, readRDS)

dat <- list(yijs = yijs, Xijlist = Xijlist, Zijlist = Xijlist,
            mi = mi, n = n)

