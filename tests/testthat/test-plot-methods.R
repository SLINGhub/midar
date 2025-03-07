#
# library(testthat)
# library(ggplot2)
# library(ggrepel)
# library(cli)
#
# test_that("custom ggplot print method handles overlap warnings correctly", {
#
#   # Sample data with potential overlapping points
#   overlap_data <- data.frame(
#     x = rep(1, 300),  # All points stacked vertically
#     y = 1:300,
#     label = "long_labellllllllllllllll"
#   )
#
#   no_overlap_data <- data.frame(
#     x = rep(1, 30),  # All points stacked vertically
#     y = 1:30,
#     label = "a"
#   )
#
#   # Create a plot that might trigger the "too many overlaps" warning
#   plot_with_overlaps <- ggplot(overlap_data, aes(x = x, y = y, label = label)) +
#     geom_point(size = 5) +
#     ggrepel::geom_text_repel(aes(label = .data$label), size = 15, na.rm = TRUE, max.overlaps = 1)
#
#   # Create plot with geom_text_repel
#   plot_without_overlaps <- ggplot(overlap_data, aes(x = x, y = y, label = label)) +
#     geom_point() +
#     geom_text_repel(size = 0.1)
#
#
#   # Create plot that might trigger overlap warnings using geom_text
#   plot_with_overlaps <- ggplot(plot_with_overlaps, aes(x = x, y = y, label = label)) +
#     geom_point(size = 15) +
#     ggrepel::geom_text_repel(aes(label = .data$label), size = 115, na.rm = TRUE, max.overlaps = 1)
#
#   # Use a PDF device to control plot size
#   temp_file <- tempfile(fileext = ".pdf")
#   pdf(temp_file, width = 1, height = 1)  # Setting small size to induce overlaps
#   warn_msg <- capture_message({print(plot_with_overlaps)})
#   dev.off()  # Close the device
#
#   # Verify if the warning message contains the expected custom warning
#   expect_true(grepl("data points are unlabeled \\(too many overlaps\\)\\. Please adjust corresponding label", warn_msg$args$text$str))
#
#   # Clean up: remove the temporary file
#   on.exit(unlink(temp_file))
# })
#
# test_that("no warning is produced when using geom_text_repel", {
#
#
#   # Use a PDF device to control plot size
#   temp_file <- tempfile(fileext = ".pdf")
#   pdf(temp_file, width = 111, height = 111)  # Small size, but using geom_text_repel
#   warn_msg <- capture_message(print(plot_without_overlaps))
#   dev.off()  # Close the device
#
#   # Ensure no warnings are captured
#   expect_true(length(warn_msg) == 0)
#
#   # Clean up: remove the temporary file
#   on.exit(unlink(temp_file))
# })
