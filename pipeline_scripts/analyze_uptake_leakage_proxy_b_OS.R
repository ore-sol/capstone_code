#to run: first: source("analyze_uptake_leakage_proxy_b_OS.R") then: (if files in working dir): run_batch()

run_interaction_indices <- function(uptake_file, leakage_file, output_prefix, exclude_same_metabolite = TRUE) {
# Calculate species-level effective facilitation, uptake competition, and a
# model-free competition-cooperation proxy b from experimental U and L data.

exclude_same_metabolite <- TRUE

clip01 <- function(x) {
  pmin(1, pmax(0, x))
}

minmax01 <- function(x) {
  x_range <- range(x, finite = TRUE)
  if (!all(is.finite(x_range)) || diff(x_range) <= .Machine$double.eps) {
    # No among-species variation means no relative evidence in either direction.
    return(rep(0.5, length(x)))
  }
  clip01((x - x_range[[1]]) / diff(x_range))
}

standardize_z <- function(x) {
  x_sd <- sd(x, na.rm = TRUE)
  if (!is.finite(x_sd) || x_sd <= .Machine$double.eps) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / x_sd
}

cat("Reading RDS files...\n")
uptake <- readRDS(uptake_file)
leakage <- readRDS(leakage_file)

if (length(dim(uptake)) != 2) {
  stop("The uptake object must be a species x metabolite matrix.")
}
if (length(dim(leakage)) != 3) {
  stop("The leakage object must be a species x reactant x product array.")
}
if (is.null(rownames(uptake)) || is.null(colnames(uptake))) {
  stop("The uptake matrix must have species row names and metabolite column names.")
}
if (any(vapply(dimnames(leakage), is.null, logical(1)))) {
  stop("All three leakage dimensions must be named: species, reactant, product.")
}
if (anyNA(uptake) || anyNA(leakage)) {
  stop("Uptake and leakage data cannot contain missing values.")
}
if (min(uptake) < 0 || min(leakage) < 0) {
  stop("Uptake and leakage values must be non-negative.")
}

# Align species by name rather than assuming that the two files use the same order.
leakage_species <- dimnames(leakage)[[1]]
common_species <- rownames(uptake)[rownames(uptake) %in% leakage_species]

if (length(common_species) < 2) {
  stop("At least two shared species are required to calculate pairwise indices.")
}

if (!identical(rownames(uptake), common_species)) {
  uptake <- uptake[common_species, , drop = FALSE]
}
if (!identical(leakage_species, common_species)) {
  leakage <- leakage[common_species, , , drop = FALSE]
}

n_species <- nrow(uptake)
uptake_metabolites <- colnames(uptake)
reactants <- dimnames(leakage)[[2]]
products <- dimnames(leakage)[[3]]

# Match L reactants/products to U columns. A product absent from U is treated as
# unavailable for uptake and therefore contributes zero effective facilitation.
reactant_in_u <- match(reactants, uptake_metabolites)
product_in_u <- match(products, uptake_metabolites)

cat(sprintf("Shared species: %d\n", n_species))
cat(sprintf(
  "Leakage reactants represented in U: %d/%d\n",
  sum(!is.na(reactant_in_u)), length(reactants)
))
cat(sprintf(
  "Leakage products represented in U: %d/%d\n",
  sum(!is.na(product_in_u)), length(products)
))

# The simulated uptake matrix has row sum 1. Row normalization makes the
# experimental calculation comparable and prevents broad-uptake species from
# receiving a larger score merely because they contain more nonzero entries.
uptake_row_sum <- rowSums(uptake)
if (any(uptake_row_sum <= 0)) {
  stop("Every species must have at least one positive uptake value.")
}
uptake_probability <- uptake / uptake_row_sum

# -----------------------------------------------------------------------------
# Competition index: mean pairwise cosine similarity of each species' U row.
# Larger values indicate greater uptake overlap and stronger potential competition.
# -----------------------------------------------------------------------------
uptake_norm <- sqrt(rowSums(uptake_probability^2))
cosine_matrix <- tcrossprod(uptake_probability) / (uptake_norm %o% uptake_norm)
diag(cosine_matrix) <- NA_real_
competition_index <- rowMeans(cosine_matrix, na.rm = TRUE)

# -----------------------------------------------------------------------------
# Effective facilitation efficiency for donor i:
#
# Q_i(a,b) = L_i(a,b) / sum_b L_i(a,b)
# A_j(b) = 1 when species j can take up product b, otherwise 0
# F_i = sum_{a,b} p_i(a) Q_i(a,b) mean_{j != i}[A_j(b)]
#       / sum_{a,b} p_i(a) Q_i(a,b)
#
# Thus F_i is the fraction of donor-weighted potential leakage that can be used
# by other community members. It lies in [0, 1] and is not inflated merely by
# having more annotated leakage edges. Terms with a == b are excluded.
# -----------------------------------------------------------------------------
uptake_presence_at_products <- matrix(
  0,
  nrow = n_species,
  ncol = length(products),
  dimnames = list(common_species, products)
)
matched_products <- !is.na(product_in_u)
uptake_presence_at_products[, matched_products] <- 1 * (
  uptake[, product_in_u[matched_products], drop = FALSE] > 0
)

total_product_accessibility <- colSums(uptake_presence_at_products)
same_metabolite_row <- match(products, reactants)
same_metabolite_product <- which(!is.na(same_metabolite_row))

effective_facilitation <- numeric(n_species)
potential_leakage_weight <- numeric(n_species)

cat("Calculating effective facilitation...\n")
for (i in seq_len(n_species)) {
  donor_reactant_uptake <- numeric(length(reactants))
  matched_reactants <- !is.na(reactant_in_u)
  donor_reactant_uptake[matched_reactants] <- uptake_probability[
    i, reactant_in_u[matched_reactants]
  ]
  
  leakage_distribution <- 1 * leakage[i, , , drop = TRUE]
  
  if (exclude_same_metabolite && length(same_metabolite_product) > 0) {
    leakage_distribution[cbind(
      same_metabolite_row[same_metabolite_product],
      same_metabolite_product
    )] <- 0
  }
  
  # Normalize the products of each reactant to a leakage distribution Q_i(a,b).
  leakage_row_sum <- rowSums(leakage_distribution)
  active_leakage_row <- leakage_row_sum > 0
  leakage_distribution[active_leakage_row, ] <- (
    leakage_distribution[active_leakage_row, , drop = FALSE] /
      leakage_row_sum[active_leakage_row]
  )
  
  # Weight every reactant -> product leakage probability by donor uptake.
  weighted_leakage <- sweep(
    leakage_distribution,
    MARGIN = 1,
    STATS = donor_reactant_uptake,
    FUN = "*"
  )
  
  leakage_to_product <- colSums(weighted_leakage)
  mean_other_accessibility <- (
    total_product_accessibility - uptake_presence_at_products[i, ]
  ) / (n_species - 1)
  
  potential_leakage_weight[[i]] <- sum(leakage_to_product)
  effective_gain <- sum(leakage_to_product * mean_other_accessibility)
  effective_facilitation[[i]] <- if (potential_leakage_weight[[i]] > 0) {
    effective_gain / potential_leakage_weight[[i]]
  } else {
    0
  }
}

# -----------------------------------------------------------------------------
# Model-free proxy b.
# Put C and F on equal standard-deviation scales, take their contrast, and map
# the contrast to [0, 1]:
#
# D_i = z(C_i) - z(F_i)
# b_i = [D_i - min(D)] / [max(D) - min(D)]
#
# This is the balanced linear contrast: neither C nor F is favoured merely
# because it has a different numerical range. It guarantees 0 <= b <= 1,
# increases with C holding F fixed, and decreases with F holding C fixed.
# -----------------------------------------------------------------------------
competition_scaled <- minmax01(competition_index)
facilitation_scaled <- minmax01(effective_facilitation)
competition_z <- standardize_z(competition_index)
facilitation_z <- standardize_z(effective_facilitation)
balance_contrast <- competition_z - facilitation_z
proxy_b <- minmax01(balance_contrast)

species_metrics <- data.frame(
  species = common_species,
  effective_facilitation = effective_facilitation,
  competition_cosine = competition_index,
  facilitation_scaled = facilitation_scaled,
  competition_scaled = competition_scaled,
  facilitation_z = facilitation_z,
  competition_z = competition_z,
  balance_contrast = balance_contrast,
  proxy_b = proxy_b,
  potential_leakage_weight = potential_leakage_weight,
  stringsAsFactors = FALSE
)
species_metrics <- species_metrics[order(species_metrics$proxy_b), ]

correlation_C_F <- cor(competition_index, effective_facilitation)
balanced_correlation_bound <- sqrt(max(0, (1 - correlation_C_F) / 2))

community_metrics <- data.frame(
  n_species = n_species,
  n_uptake_metabolites = ncol(uptake),
  n_leakage_reactants = length(reactants),
  n_leakage_products = length(products),
  mean_effective_facilitation = mean(effective_facilitation),
  mean_competition_cosine = mean(competition_index),
  mean_proxy_b = mean(proxy_b),
  correlation_competition_facilitation = correlation_C_F,
  correlation_b_facilitation = cor(proxy_b, effective_facilitation),
  correlation_b_competition = cor(proxy_b, competition_index),
  balanced_correlation_upper_bound = balanced_correlation_bound
)

species_csv <- paste0(output_prefix, "_species.csv")
community_csv <- paste0(output_prefix, "_community.csv")
plot_png <- paste0(output_prefix, "_dual_axis.png")
plot_pdf <- paste0(output_prefix, "_dual_axis.pdf")

write.csv(species_metrics, species_csv, row.names = FALSE)
write.csv(community_metrics, community_csv, row.names = FALSE)

draw_dual_axis_plot <- function(metrics) {
  cooperation_colour <- "#2166AC"
  competition_colour <- "#B2182B"
  point_order <- order(metrics$proxy_b)
  x <- metrics$proxy_b[point_order]
  cooperation <- metrics$effective_facilitation[point_order]
  competition <- metrics$competition_cosine[point_order]
  
  cooperation_ylim <- extendrange(cooperation, f = 0.08)
  competition_ylim <- extendrange(competition, f = 0.08)
  
  par(mar = c(5.2, 6.3, 2.0, 6.3), las = 1)
  plot(
    x, cooperation,
    xlim = c(0, 1), ylim = cooperation_ylim,
    pch = 16, cex = 0.9,
    col = adjustcolor(cooperation_colour, alpha.f = 0.72),
    xlab = "", ylab = "",
    axes = FALSE,
    bty = "n"
  )
  axis(1)
  axis(2, col.axis = cooperation_colour, col.ticks = cooperation_colour)
  mtext(
    expression("Competition-cooperation proxy " * b),
    side = 1, line = 3.2
  )
  mtext(
    "Cooperation index (effective facilitation)",
    side = 2, line = 4.5, col = cooperation_colour, las = 0
  )
  grid(col = "grey88", lty = 1)
  points(
    x, cooperation,
    pch = 16, cex = 0.9,
    col = adjustcolor(cooperation_colour, alpha.f = 0.72)
  )
  
  par(new = TRUE)
  plot(
    x, competition,
    xlim = c(0, 1), ylim = competition_ylim,
    type = "p", pch = 17, cex = 0.9,
    col = adjustcolor(competition_colour, alpha.f = 0.72),
    axes = FALSE, xlab = "", ylab = ""
  )
  axis(4, col.axis = competition_colour, col.ticks = competition_colour)
  mtext(
    "Competition index (mean cosine similarity)",
    side = 4, line = 4.5, col = competition_colour, las = 0
  )
  box()
  legend(
    "topright",
    legend = c("Effective facilitation", "Cosine similarity"),
    pch = c(16, 17),
    col = c(cooperation_colour, competition_colour),
    bty = "n"
  )
}

png(plot_png, width = 2400, height = 1800, res = 300)
draw_dual_axis_plot(species_metrics)
dev.off()

pdf(plot_pdf, width = 8, height = 6)
draw_dual_axis_plot(species_metrics)
dev.off()

cat("\nCommunity summary:\n")
print(community_metrics, row.names = FALSE, digits = 6)
cat("\nCreated files:\n")
cat(sprintf("  %s\n", species_csv))
cat(sprintf("  %s\n", community_csv))
cat(sprintf("  %s\n", plot_png))
cat(sprintf("  %s\n", plot_pdf))
}

run_batch <- function(data_dir = ".", out_dir = ".") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  uptake_files <- list.files(data_dir, pattern = "^uptake_matrix_.+\\.rds$", full.names = TRUE)
  if (length(uptake_files) == 0) stop("No uptake_matrix_*.rds files found in: ", data_dir)
  for (i in seq_along(uptake_files)) {
    uf  <- uptake_files[[i]]
    key <- sub("^uptake_matrix_(.+)\\.rds$", "\\1", basename(uf))
    lf  <- file.path(data_dir, paste0("leakage_tensor_", key, ".rds"))
    pfx <- file.path(out_dir,  paste0(key, "_interaction_indices"))
    cat(sprintf("\n[%d/%d] Sample: %s\n", i, length(uptake_files), key))
    if (!file.exists(lf)) { message("  SKIP — no matching leakage tensor: ", lf); next }
    tryCatch(
      run_interaction_indices(uf, lf, pfx),
      error = function(e) message("  FAIL — ", conditionMessage(e))
    )
  }
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  run_batch(
    data_dir = if (length(args) >= 1) args[[1]] else ".",
    out_dir  = if (length(args) >= 2) args[[2]] else "."
  )
}

