# XAD Automatic Differentiation Demo
# Run from odelia package root: Rscript examples/parabola/demo.R
#
# Objective: f(x,y) = -(x-2)² - (y-3)²  (maximum at (2, 3))

library(tidyverse)

# devtools::load_all(".", quiet = TRUE)
library(odelia)

# Link against the installed odelia library
so_name <- paste0("odelia", .Platform$dynlib.ext)
lib_path <- system.file("libs", package = "odelia")
lib_file <- list.files(lib_path, pattern = so_name, recursive = TRUE, full.names = TRUE)[1]
# lib_dir <- dirname(lib_file)

# Use full path to .so because it might not be named libodelia.so
Sys.setenv(PKG_LIBS = paste0(lib_file, " -Wl,-rpath,", dirname(lib_file)))

Rcpp::sourceCpp("inst/examples/parabola/src/parabola_interface.cpp", verbose = TRUE)

cat("XAD Automatic Differentiation Demo\n\n")

# Capture optimization path by wrapping the objective
path <- tibble(x = numeric(), y = numeric(), value = numeric())
objective <- function(p) {
    result <- parabola_eval(p)
    path <<- path |> add_row(x = p[1], y = p[2], value = result$value)
    -result$value  # negate for minimization
}
gradient <- function(p) -parabola_eval(p)$gradient

# Run L-BFGS-B
initial_guess <- c(-1, -1)
opt_result <- optim(
    par = initial_guess,
    fn = objective,
    gr = gradient,
    method = "L-BFGS-B"
)

cat("Optimization path:\n")
path <- path |> mutate(step = row_number() - 1)
print(path, n = Inf)

# Create grid for contour plot
grid <- expand_grid(
    x = seq(-2, 5, length.out = 100),
    y = seq(-2, 6, length.out = 100)
) |>
    rowwise() |>
    mutate(z = parabola_eval(c(x, y))$value) |>
    ungroup()

# Plot
p <- ggplot() +
    geom_contour_filled(data = grid, aes(x = x, y = y, z = z), bins = 15, alpha = 0.7) +
    geom_path(data = path, aes(x = x, y = y), color = "white", linewidth = 1,
              arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    geom_point(data = path, aes(x = x, y = y), color = "white", size = 2) +
    geom_point(data = path |> slice(1), aes(x = x, y = y), color = "red", size = 4) +
    annotate("text", x = path$x[1] - 0.3, y = path$y[1] - 0.5, 
             label = "Start", color = "red", fontface = "bold") +
    geom_point(data = path |> slice(n()), aes(x = x, y = y), color = "chartreuse", size = 4) +
    annotate("text", x = 2, y = 3.6, label = "Maximum\n(2, 3)", 
             color = "chartreuse", fontface = "bold") +
    geom_point(aes(x = 2, y = 3), shape = 4, color = "white", size = 5, stroke = 2) +
    labs(title = "Gradient Ascent with Automatic Differentiation (XAD)",
         subtitle = expression(f(x,y) == -(x-2)^2 - (y-3)^2),
         x = "x", y = "y", fill = "f(x,y)") +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

output_file <- "inst/examples/parabola/gradient_ascent.png"
ggsave(output_file, p, width = 8, height = 6, dpi = 150)

cat("\n✓ Plot saved to:", output_file, "\n")
cat("  Initial: (", initial_guess[1], ", ", initial_guess[2], ")\n", sep = "")
cat("  Final: (", round(opt_result$par[1], 4), ", ", round(opt_result$par[2], 4), ")\n", sep = "")
cat("  Steps:", nrow(path), "\n")
