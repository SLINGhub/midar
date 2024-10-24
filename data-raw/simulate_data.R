NULL
# #################################################################################
# ## Random dataset with 7 batches with batch effects
# # used ChatGTP package to generate the data
# ################################################################################
#
#
# library(dplyr)  # Load dplyr for tibble conversion
#
# set.seed(123)  # For reproducibility
#
# # Parameters
# n <- 700  # Total number of data points
# segments <- 7  # Number of batches (including the new batch)
# mean_main <- 1000  # Average for batch 1
#
# # Generate main batches with specified variances and means
# batch_size <- n / segments
#
# # Main segments with specified variances and means
# y_main <- c(
#   rnorm(batch_size, mean = mean_main, sd = 100),          # Segment 1: mean = 1000, sd = 100
#   rnorm(batch_size, mean = mean_main * 2, sd = 50),      # Segment 2: mean = 2000, sd = 50
#   rnorm(batch_size, mean = mean_main, sd = 100),          # Segment 3: mean = 1000, sd = 100
#   rnorm(batch_size, mean = mean_main * 2, sd = 200),     # Segment 4: mean = 2000, sd = 200
#   rnorm(batch_size, mean = mean_main, sd = 100),          # Segment 5: mean = 1000, sd = 100
#   rnorm(batch_size, mean = mean_main / 2, sd = 100 * 1.5),  # Segment 6: mean = 500, sd = 150
#   rnorm(batch_size, mean = mean_main, sd = 100)           # Segment 7: mean = 1000, sd = 100 (same as Segment 1)
# )
#
# # Create x values (linear increase)
# x_values <- seq(1, n)
#
# # Combine x and y into a data frame
# data <- data.frame(x = x_values, y = y_main)
#
# # Modify every 10th sample to have half the variance
# for (i in seq(10, n, by = 10)) {
#   # Determine the current segment
#   segment_index <- (i - 1) %/% batch_size + 1
#
#   # Get the mean and standard deviation for that segment
#   if (segment_index == 1) {
#     mean_value <- mean_main
#     sd_value <- 100
#   } else if (segment_index == 2) {
#     mean_value <- mean_main * 2
#     sd_value <- 50
#   } else if (segment_index == 3) {
#     mean_value <- mean_main
#     sd_value <- 100
#   } else if (segment_index == 4) {
#     mean_value <- mean_main * 2
#     sd_value <- 200
#   } else if (segment_index == 5) {
#     mean_value <- mean_main
#     sd_value <- 100
#   } else if (segment_index == 6) {
#     mean_value <- mean_main / 2
#     sd_value <- 150  # For Segment 6
#   } else {
#     mean_value <- mean_main
#     sd_value <- 100  # For Segment 7 (same as Segment 1)
#   }
#
#   # Generate a new value with half the variance
#   new_value <- rnorm(1, mean = mean_value, sd = sd_value / sqrt(2))
#   data$y[i] <- new_value
# }
#
# # Create a color vector and batch numbers
# colors <- rep("blue", n)
# colors[seq(10, n, by = 10)] <- "red"  # Set every 10th point to red
# data$batch <- rep(1:segments, each = batch_size)
#
# # Add color and type information to the data frame
# data$color <- ifelse(colors == "red", "red", "blue")
# data$type <- ifelse(data$color == "red", "BQC", "SPL")
#
# # Convert to a tibble
# data_tibble <- as_tibble(data)
# data_tibble_new <- data_tibble |> select(analysis_id = x, qc_type = type, batch_id = batch, conc = y)
# readr::write_excel_csv(data_tibble_new, file = "batch_effect-simdata-u1000-sd100_7batches.csv")
# #
# # # Plot the data with custom colors
# # plot(data_tibble$x, data_tibble$y, main = "Random Linear Data with Batch and Type Information", xlab = "X-axis", ylab = "Y-axis", pch = 19, col = colors)
# # abline(lm(y ~ x, data = data_tibble), col = "black")  # Add linear regression line
# #
# # # View the tibble
# # print(data_tibble)
#
#
#
# #################################################################################
# ## Random dataset with 7 batches with batch effects and within-batch nonlinear drift
# # used ChatGTP package to generate the data
# ################################################################################
#
# library(dplyr)  # Load dplyr for tibble conversion
#
# set.seed(123)  # For reproducibility
#
# # Parameters
# n <- 1400  # Total number of data points
# segments <- 7  # Number of batches
# mean_main <- 1000  # Average for batch 1
#
# # Generate main batches with specified variances and means
# batch_size <- n / segments
#
# # Main segments with specified variances and means
# y_main <- c(
#   rnorm(batch_size, mean = mean_main, sd = 100),          # Segment 1: mean = 1000, sd = 100
#   rnorm(batch_size, mean = mean_main * 2, sd = 50),      # Segment 2: mean = 2000, sd = 50
#   rnorm(batch_size, mean = mean_main, sd = 100),          # Segment 3: mean = 1000, sd = 100
#   rnorm(batch_size, mean = mean_main * 2, sd = 200),     # Segment 4: mean = 2000, sd = 200
#   rnorm(batch_size, mean = mean_main, sd = 100),          # Segment 5: mean = 1000, sd = 100
#   rnorm(batch_size, mean = mean_main / 2, sd = 100 * 1.5),  # Segment 6: mean = 500, sd = 150
#   rnorm(batch_size, mean = mean_main, sd = 100)           # Segment 7: mean = 1000, sd = 100
# )
#
# # Create x values (linear increase)
# x_values <- seq(1, n)
#
# # Define frequencies for each segment and decrease them by half
# frequencies <- c(0.5, 1, 1.5, 0.25, 1.25, 0.75, 1)  # Halved frequencies
# shifts <- runif(segments, min = 50, max = 150)  # Random shifts between 50 and 150
#
# # Add strong nonlinear drift for each batch with different frequencies and random shifts
# for (segment in 1:segments) {
#   start_index <- (segment - 1) * batch_size + 1
#   end_index <- segment * batch_size
#   # Stronger nonlinear drift using sine function with decreased frequencies and random shifts
#   drift <- 200 * sin(seq(0, frequencies[segment] * 2 * pi, length.out = batch_size)) + shifts[segment]  # Added random shift
#   y_main[start_index:end_index] <- y_main[start_index:end_index] + drift
# }
#
# # Combine x and y into a data frame
# data <- data.frame(x = x_values, y = y_main)
#
# # Modify every 10th sample to have half the variance
# for (i in seq(10, n, by = 10)) {
#   # Determine the current segment
#   segment_index <- (i - 1) %/% batch_size + 1
#
#   # Get the mean and standard deviation for that segment
#   if (segment_index == 1) {
#     mean_value <- mean_main
#     sd_value <- 100
#   } else if (segment_index == 2) {
#     mean_value <- mean_main * 2
#     sd_value <- 50;
#   } else if (segment_index == 3) {
#     mean_value <- mean_main;
#     sd_value <- 100;
#   } else if (segment_index == 4) {
#     mean_value <- mean_main * 2;
#     sd_value <- 200;
#   } else if (segment_index == 5) {
#     mean_value <- mean_main;
#     sd_value <- 100;
#   } else if (segment_index == 6) {
#     mean_value <- mean_main / 2;
#     sd_value <- 150;  # For Segment 6
#   } else {
#     mean_value <- mean_main;
#     sd_value <- 100;  # For Segment 7 (same as Segment 1)
#   }
#
#   # Generate a new value with half the variance
#   new_value <- rnorm(1, mean = mean_value, sd = sd_value / sqrt(2))
#   data$y[i] <- new_value
# }
#
# # Create a color vector and batch numbers
# colors <- rep("blue", n)
# colors[seq(10, n, by = 10)] <- "red"  # Set every 10th point to red
# data$batch <- rep(1:segments, each = batch_size)
#
# # Add color and type information to the data frame
# data$color <- ifelse(colors == "red", "red", "blue")
# data$type <- ifelse(data$color == "red", "BQC", "SPL")
#
# # Convert to a tibble
# data_tibble <- as_tibble(data)
#
# # Plot the data with custom colors
# plot(data_tibble$x, data_tibble$y, main = "Random Linear Data with Strong Nonlinear Drift (Random Shifts)", xlab = "X-axis", ylab = "Y-axis", pch = 19, col = colors)
# abline(lm(y ~ x, data = data_tibble), col = "black")  # Add linear regression line
#
# # Convert to a tibble
# data_tibble <- as_tibble(data)
# data_tibble_new <- data_tibble |> select(analysis_id = x, qc_type = type, batch_id = batch, conc = y)
# readr::write_excel_csv(data_tibble_new, file = "drift_batch_effect-simdata-u1000-sd100_7batches.csv")
