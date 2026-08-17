############################################################
# Demo analysis
# STEM121+ fibers at 12, 16, and 20 mm caudal
# TP vs Combination
# Mann-Whitney U test with Holm correction
############################################################

# Read demo data
df <- read.csv("demo/STEM121_demo_data.csv")

# Distances
distances <- c(
  "12 mm caudal",
  "16 mm caudal",
  "20 mm caudal"
)

# Mann-Whitney U test at each distance
results <- lapply(distances, function(d) {
  
  sub <- df[df$distance == d, ]
  
  tp <- sub$value[sub$group == "TP"]
  comb <- sub$value[sub$group == "Combination"]
  
  test <- wilcox.test(
    tp,
    comb,
    paired = FALSE,
    exact = TRUE,
    alternative = "two.sided"
  )
  
  data.frame(
    distance = d,
    n_TP = length(tp),
    n_Combination = length(comb),
    mean_TP = mean(tp),
    mean_Combination = mean(comb),
    p_value = test$p.value
  )
})

results <- do.call(rbind, results)

# Holm correction
results$p_Holm <- p.adjust(
  results$p_value,
  method = "holm"
)

# Display results
print(results)

# Save results
write.csv(
  results,
  "demo/STEM121_demo_results.csv",
  row.names = FALSE
)