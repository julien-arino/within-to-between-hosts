library(ggplot2)
df <- data.frame(x = 1:3, y = 1:3, cat = c("Cell", "Immune", "Viral"))
p <- ggplot(df, aes(x, y, fill = cat)) +
  geom_rect(aes(xmin=x-0.4, xmax=x+0.4, ymin=y-0.4, ymax=y+0.4), alpha = 0.15, color = NA) +
  scale_fill_manual(values = c(Cell = "#800080", Viral = "#A9A9A9", Immune = "#F5F5F5")) +
  guides(fill = guide_legend(override.aes = list(fill = c("#ECD9EC", "#FEFEFE", "#F2F2F2"), color = "black", linewidth = 1, size = 1)))
ggsave("test_legend.png", p, width=5, height=5)
