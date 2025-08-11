data<-read.table('raw-data/UV-Vis/Pyrrolidine-DMBO/PD5.csv', header=T, skip=15, sep=',')

plot(data$Time.s., data$Absorbance, type='l', ylab='Absorbance (307 nm)', xlab='Time (sec.)', las=1,lwd=2)

# Estimate initial parameters
A0_start <- max(data$Absorbance) - min(data$Absorbance)
k_start <- 0.001  # guess based on timescale
Ainf_start <- min(data$Absorbance)

# Fit pseudo-first-order kinetics
fit <- nls(Absorbance ~ A0 * exp(-k * Time.s.) + Ainf,
           data = data,
           start = list(A0 = A0_start, k = k_start, Ainf = Ainf_start))

# Summary of fit
summary(fit)

# Add fitted curve to the plot
lines(data$Time.s., predict(fit), col='red', lwd=2)
legend("topright", legend=c("Data", "Fit"), col=c("black", "red"), lwd=2)

summary(fit)$coefficients

k <- summary(fit)$coefficients["k", "Estimate"]
print(k)



# Assuming you have a folder with all the data.


#' Extract Pseudo-First-Order Rate Constants from UV-Vis Data
#'
#' This function loops through all `.csv` files in a given folder,
#' fits a pseudo-first-order exponential decay model to UV-Vis absorbance data,
#' and extracts the fitted rate constant \eqn{k} from each file.
#'
#' The model fitted is: \deqn{A(t) = A_0 \cdot e^{-k t} + A_{\infty}}
#'
#' The data files are assumed to have a common format:
#' - CSV file
#' - First 15 lines skipped
#' - Columns named `Time.s.` and `Absorbance`
#'
#' @param folder_path Character string. Path to the folder containing the `.csv` data files.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{File}{Name of the file}
#'   \item{k}{Estimated pseudo-first-order rate constant}
#' }
#'
#' @examples
#' \dontrun{
#' # Extract rate constants from all UV-Vis traces in folder
#' k_results <- extract_k_values("raw-data/UV-Vis/")
#' print(k_results)
#' }
#'
#' @importFrom stats nls coef
#' @export
extract_k_values <- function(folder_path, drift = FALSE, plot_results = TRUE) {
  # List all CSV files in the folder
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

  # Initialize results list
  results <- data.frame(File = character(), k = numeric(), stringsAsFactors = FALSE)

  for (file in files) {
    # Read the data (skip header lines and assume same format)
    data <- tryCatch({
      read.table(file, header = TRUE, skip = 15, sep = ",")
    }, error = function(e) return(NULL))

    if (is.null(data) || !"Time.s." %in% names(data) || !"Absorbance" %in% names(data)) next

    # Estimate initial parameters
    A0_start <- max(data$Absorbance) - min(data$Absorbance)
    k_start <- 0.005
    Ainf_start <- min(data$Absorbance)
    m_start <- 0  # baseline slope for drift model

    # Choose formula based on drift argument
    if (drift) {
      formula_fit <- Absorbance ~ A0 * exp(-k * Time.s.) + Ainf + m * Time.s.
      start_list <- list(A0 = A0_start, k = k_start, Ainf = Ainf_start, m = m_start)
    } else {
      formula_fit <- Absorbance ~ A0 * exp(-k * Time.s.) + Ainf
      start_list <- list(A0 = A0_start, k = k_start, Ainf = Ainf_start)
    }

    # Fit the model
    fit <- tryCatch({
      nls(formula_fit,
          data = data,
          start = start_list,
          control = nls.control(warnOnly = TRUE))
    }, error = function(e) return(NULL))

    # Extract k value if fitting succeeded
    if (!is.null(fit)) {
      k_value <- summary(fit)$coefficients["k", "Estimate"]
      results <- rbind(results, data.frame(File = basename(file), k = k_value))

      # Plot if requested
      if (plot_results) {
        plot(data$Time.s., data$Absorbance, type = 'l',
             ylab = 'Absorbance (307 nm)', xlab = 'Time (sec.)',
             las = 1, lwd = 2,
             main = paste("Fit for", basename(file)))
        lines(data$Time.s., predict(fit), col = 'red', lwd = 2)
        legend("topright", legend = c("Data", "Fit"), col = c("black", "red"), lwd = 2)


      }
    }
  }

  return(results)
}

par(mfrow=c(3,3))
par(mfrow=c(3,2))
k_results <- extract_k_values("raw-data/UV-Vis/DBCO_Biodrop_250630/", drift=FALSE, plot_results = TRUE)
molarity_dbco <- c(0.03, 0.04, 0.05, 0.07, 0.10)  # in M
k_values_dbco <- k_results$k


k_results <- extract_k_values("raw-data/UV-Vis/Pyrrolidine-DMBO/", drift=FALSE, plot_results = TRUE)
print(k_results)

# Molarity of azide across conditions
molarity <- c(0.03, 0.03, 0.03, 0.04, 0.04, 0.05, 0.07, 0.10)  # in M

# Ensure k_results is a numeric vector of the same length
k_values <- k_results$k

k_values<- tapply(k_values, molarity, mean)
molarity <- unique(molarity)

# Simulate pseudo-first-order rate constants with some noise
#set.seed(42)  # for reproducibility
#true_k2 <- 25  # M^-1 s^-1
#k_values <- true_k2 * molarity + rnorm(length(molarity), mean = 0, sd = 0.5)


# Plot pseudo-first-order rate constants vs azide concentration
plot(molarity, k_values,
     xlab = "[Azide] (M)", ylab = expression(k[obs]~"(" * s^{-1} * ")"),
     pch = 16, las = 1, col = "blue", main = "Second-order kinetics fit")

# Fit linear model: k_obs = k2 * [azide]
fit_lm <- lm(k_values ~ molarity + 0)

# Add regression line
abline(fit_lm, col = "red", lwd = 2)

# Extract second-order rate constant (slope)
k2 <- coef(fit_lm)[["molarity"]]
print(paste("Second-order rate constant k2 =", signif(k2, 4), "M^-1 s^-1"))

# Optional: Add R^2 to plot
r2 <- summary(fit_lm)$r.squared
# Legend with math-style units
legend("topleft", legend = c(
  bquote(k[2] == .(signif(k2, 4)) ~ M^{-1} %.% s^{-1}),
  bquote(R^2 == .(signif(r2, 3)))
), bty = "n", text.col = "darkred")


normalize_abs <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


pd_2<-read.table("raw-data/UV-Vis/Pyrrolidine-DMBO/PD2.csv", skip = 15, sep = ",", header=TRUE)
dbco_2<-read.table("raw-data/UV-Vis/DBCO_Biodrop_250630/03_DBCO273-01_0701123614.csv", skip = 15, sep = ",", header=TRUE)


# Main plot
par(yaxs='i', xaxs='i')
plot(pd_2$Time.s., normalize_abs(pd_2$Absorbance),
     type='l', lwd=2, las=1,
     ylab='Absorbance (A.U)', xlab='Time (sec.)')
lines(dbco_2$Time.s.[1:102], normalize_abs(dbco_2$Absorbance)[1:102], lty=2)

# Add legend above the plot
legend("top", inset = -0.15,              # inset negative to move outside
       legend = c("Pyrrolidine-DMBO", "ADIBO-amine"),
       pch = c(16, 21),
       lty = c(1, 2),
       xpd = TRUE,                       # allow drawing outside plot
       horiz = TRUE,                     # horizontal legend
       bty = "n")                        # no box
par(xpd=FALSE)

# Inset plot (top-right)
par(fig=c(0.35, 0.95, 0.35, 0.95), new=TRUE, yaxs='i', xaxs='i')
plot(molarity_dbco, k_values_dbco, ylim=c(0, 0.04), xlim=c(0, 0.13),
     las=1, type='n', xlab="[Azide] (M)",
     ylab=expression(k[obs]~"(" * s^{-1} * ")"))
abline(fit_lm)
fit_lm_dbco <- lm(k_values_dbco ~ molarity_dbco + 0)
abline(fit_lm_dbco, lty=2)
points(molarity, k_values, pch=16)
points(molarity_dbco, k_values_dbco, pch=21, bg='white')

# Reset plotting region
par(fig=c(0, 1, 0, 1))

df<-data.frame(y = c(k_values, k_values_dbco), x = c(molarity, molarity_dbco), group = c(rep(1, length(molarity)), rep(0, length(molarity_dbco))))
model <- lm(y ~ x * group + 0, data = df)
summary(model)

# Extract summary
sm <- summary(model)

# Detect interaction term automatically
interaction_name <- grep(":", rownames(sm$coefficients), value = TRUE)

# Extract interaction coefficient info
coef_row <- sm$coefficients[interaction_name, ]

# Format p-value according to APA
p_val <- coef_row["Pr(>|t|)"]
p_str <- if (p_val < .001) {
  "p < .001"
} else {
  sprintf("p = %.3f", p_val) |> sub("^0", "", .)
}

# Print APA-style
cat(sprintf(
  "The interaction between x and group was significant, b = %.2f, SE = %.2f, t(%d) = %.2f, %s.\n",
  coef_row["Estimate"],
  coef_row["Std. Error"],
  df.residual(model),
  coef_row["t value"],
  p_str
))
