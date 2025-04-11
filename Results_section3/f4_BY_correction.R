### Load f4-statistics results from ADMIXTOOLS output and correct Z-scores using the Benjamini-Yekutieli procedure

#######################

# Load the openxlsx package
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)

# Read 'f4..out' file
  f4_out_file <- 'f4..out'
  if (!file.exists(f4_out_file)) {
    warning(paste("File not found:", f4_out_file))
    next
  }

  # Read the data and assign column names
  data <- read.table(f4_out_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(data) < 10) {
    warning(paste("Data in", f4_out_file, "has fewer than 10 columns"))
    next
  }
  colnames(data) <- c('Outgroup', 'Pop_1', 'Pop_2', 'Pop_3', 'f4', 'SE', 'Z', 'n_ABBAs', 'n_BABAs', 'n_SNPs')

  # Convert relevant columns to appropriate data types
  data$f4 <- as.numeric(data$f4)
  data$SE <- as.numeric(data$SE)
  data$Z <- as.numeric(data$Z)
  data$n_ABBAs <- as.integer(data$n_ABBAs)
  data$n_BABAs <- as.integer(data$n_BABAs)
  data$n_SNPs <- as.integer(data$n_SNPs)

  # Convert Z-scores to two-tailed p-values using log.p = TRUE
  log_p_values <- pnorm(-abs(data$Z), log.p = TRUE) + log(2)

  # Apply Benjamini-Yekutieli correction manually on the log scale
  n_tests <- length(log_p_values)
  harmonic_number <- sum(1 / seq_len(n_tests))

  # Convert log p-values to p-values for adjustment
  p_values <- exp(log_p_values)

  # Apply BY correction using p.adjust
  adjusted_p_values <- p.adjust(p_values, method = "BY")

  # Compute log of adjusted p-values
  log_adjusted_p_values <- log(adjusted_p_values)

  # Convert adjusted log p-values back to adjusted Z-scores
  z_by <- qnorm(log_adjusted_p_values - log(2), lower.tail = FALSE, log.p = TRUE)
  z_by <- z_by * sign(data$Z)  # Preserve the sign

  # Add the adjusted Z-scores to the data frame
  data$Z_BY <- z_by

  # Write the results to {region}_results_ancients.tsv
  output_file <- "f4.tsv"
  data <- data[order(data$Z_BY), ]  # Sort data by the Z_BY column in ascending order
  write.table(data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

  output_file <- "f4.xlsx"
  write.xlsx(data, file = output_file, rowNames = FALSE)

  cat("Data has been saved as both a tsv and an Excel table.\n")

  # Check number of lines against 'f4.list' minus 3
  f4_list_file <- 'f4.list'
  if (file.exists(f4_list_file)) {
    f4_list_lines <- length(readLines(f4_list_file))
    expected_lines <- f4_list_lines - 4
    result_lines <- nrow(data)
    if (result_lines != expected_lines) {
      warning(paste("Number of lines in", output_file, "(", result_lines,
                    ") does not match expected number of lines (", expected_lines,
                    ") in"))
    }
  } else {
    warning(paste("File 'f4.list' not found. Cannot check line counts."))
  }
