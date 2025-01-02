library(lubridate)


#"data_before.csv" is the raw dataset
data_financial = read.csv2("data_before.csv")  # read function from package readxl
for (i in seq_along(data_financial$bankrupt_day)) {
  if ( is.na(data_financial$bankrupt_day[i]) ) {
    data_financial$bankrupt_day[i] = as.Date('2020-12-31')
  }
}
ind = which(data_financial$status == 1)
data_financial$start_day = sapply(data_financial$start_day, convert_to_ymd)
data_financial$bankrupt_day = sapply(data_financial$bankrupt_day, convert_to_ymd)
# Step 1: Count the occurrences of each year in both vectors
# Example vectors of years
years_vector1 <- year(data_financial$start_day)
years_vector2 <- year(data_financial$bankrupt_day[ind])
years_vector3 <- year(data_financial$bankrupt_day[-ind])

# Step 1: Count the occurrences of each year in all three vectors
year_counts1 <- table(years_vector1)
year_counts2 <- table(years_vector2)
year_counts3 <- table(years_vector3)

# Step 2: Get the union of all years (ensure all vectors have the same years)
all_years <- sort(union(union(names(year_counts1), names(year_counts2)), names(year_counts3)))

# Step 3: Convert year counts into numeric vectors for all three sets
year_counts1_aligned <- as.numeric(year_counts1[match(all_years, names(year_counts1))])
year_counts2_aligned <- as.numeric(year_counts2[match(all_years, names(year_counts2))])
year_counts3_aligned <- as.numeric(year_counts3[match(all_years, names(year_counts3))])

# Replace NAs with 0 for any missing years
year_counts1_aligned[is.na(year_counts1_aligned)] <- 0
year_counts2_aligned[is.na(year_counts2_aligned)] <- 0
year_counts3_aligned[is.na(year_counts3_aligned)] <- 0

# Step 4: Combine the counts into a matrix
counts_matrix <- rbind(year_counts1_aligned, year_counts2_aligned, year_counts3_aligned)

# Step 5: Create a grouped bar plot
barplot(counts_matrix,
        beside = TRUE,                     # Group the bars beside each other
        xlab = "Year",
        ylab = "Number of firms",
        names.arg = all_years,              # Use the year names for the X-axis
        col = c("green", "red", "orange"),  # Different colors for each vector
        legend.text = c("Listed firms", "Financial distress firms", "Censored firms"),  # Add a legend
        args.legend = list(x = "topleft"))  # Legend position

