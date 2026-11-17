library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
setwd(args[[1]])
results <- read.csv('summary.csv', stringsAsFactors = FALSE)
runs <- read.csv('runs.csv', stringsAsFactors = FALSE)
results <- results[order(results$proteins, results$clusters), ]
datasets <- unique(results$dataset)
modes <- c("python", "rust", "rust-kmer", "wasm", "wasm-kmer")
mode_labels <- c("Python", "Rust", "Rust + k-mer", "WebAssembly", "WebAssembly + k-mer")
baseline <- results[results$mode == "rust", c("dataset", "median_retained_links")]
runs$links_vs_full_rust_percent <- 100 * runs$retained_links /
  baseline$median_retained_links[match(runs$dataset, baseline$dataset)]
for (name in c("results", "runs")) {
  data <- get(name)
  data$dataset <- factor(data$dataset, levels = datasets)
  data$mode <- factor(data$mode, levels = modes, labels = mode_labels)
  assign(name, data)
}

colours <- c("Python" = "#4477AA", "Rust" = "#228833",
             "Rust + k-mer" = "#CCBB44", "WebAssembly" = "#AA3377",
             "WebAssembly + k-mer" = "#66CCEE")
base_plot <- ggplot(results, aes(x = dataset, colour = mode, shape = mode)) +
  scale_colour_manual(values = colours, drop = FALSE) +
  labs(x = "Fixture", colour = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

runtime <- base_plot +
  geom_point(data = runs, aes(y = elapsed_seconds),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5, seed = 1),
             size = 1.5, alpha = 0.45) +
  geom_point(aes(y = median_elapsed_seconds), position = position_dodge(width = 0.5), size = 3) +
  scale_y_log10() +
  labs(y = "Time to figure JSON (seconds, log scale)")

links <- base_plot +
  geom_point(data = runs, aes(y = links_vs_full_rust_percent),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5, seed = 1),
             size = 1.5, alpha = 0.45, na.rm = TRUE) +
  geom_point(aes(y = links_vs_full_rust_percent),
             position = position_dodge(width = 0.5), size = 3, na.rm = TRUE) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey65") +
  labs(y = "Link count (% of unfiltered Rust)")

ggsave("benchmark-runtime.svg", runtime, device = grDevices::svg, width = 8, height = 4.5)
ggsave("benchmark-links.svg", links, device = grDevices::svg, width = 8, height = 4.5)

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
columns <- c("dataset", "clusters", "proteins", "mode", "median_elapsed_seconds",
             "speedup_vs_python", "median_retained_links", "median_links_shared_with_rust",
             "median_groups", "median_groups_shared_with_rust", "links_vs_full_rust_percent")
headers <- c("Fixture", "Clusters", "Proteins", "Mode", "Time (s)",
             "Speedup vs Python", "Retained links", "Shared links", "Groups", "Shared groups", "Links vs Rust (%)")
table <- results[columns]
for (column in columns) {
  if (is.numeric(table[[column]])) table[[column]] <- round(table[[column]], 2)
  table[[column]] <- escape_html(ifelse(is.na(table[[column]]), "—", as.character(table[[column]])))
}
rows <- apply(table, 1, function(row) paste0("<tr>", paste0("<td>", row, "</td>", collapse = ""), "</tr>"))
writeLines(c(
  "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>Clinker benchmark results</title>",
  "<style>body{font:16px sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;color:#222}table{border-collapse:collapse}th,td{padding:.5em;border-bottom:1px solid #ddd;text-align:left}img{max-width:100%}</style></head><body>",
  "<h1>Clinker benchmark results</h1><p>The table shows medians; plots show each repeat and its median. Link percentages compare counts with unfiltered Rust. Shared links and groups compare gene identities with unfiltered Rust.</p>",
  paste0("<table><thead><tr>", paste0("<th>", headers, "</th>", collapse = ""), "</tr></thead><tbody>"),
  rows, "</tbody></table>",
  "<h2>Elapsed time</h2><img src=\"benchmark-runtime.svg\" alt=\"Median elapsed time by fixture and mode\">",
  "<h2>Retained links</h2><img src=\"benchmark-links.svg\" alt=\"Retained links relative to unfiltered Rust\">",
  "<p>CLI time includes process startup and JSON transfer. Browser time starts with Chromium running and includes worker startup and transfers. Saving results and browser rendering are excluded. See runs.csv for comparison limitations.</p></body></html>"
), "benchmark-results.html")
