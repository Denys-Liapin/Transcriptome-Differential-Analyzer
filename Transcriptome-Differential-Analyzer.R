library(readxl)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggvenn)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)


# завантажили таблицю
rowData <- as.data.frame(read_excel("DenysData.xlsx"))
# видалили перші 2 рядки з технічною інформацією
data <- rowData[-(1:2), ]
# назва рядків = тексту в першому стовпчику
rownames(data) <- make.unique(as.character(data[, 1]))
# прибрали перший стовпчик
data <- data[ ,-1]

# виявлення дуплікатів
#duplicated_genes <- rowData[duplicated(rowData[, 1]), 1]
#unique_duplicates <- unique(duplicated_genes)
#writeLines(as.character(unique_duplicates), "duplicated_genes_list.txt")

# підготовка метаданних
coldata <- data.frame(
  row.names = c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"),
  Genotype = factor(c("Het","Het","Het","Het", "Hom","Hom","Hom","Hom", "Het","Het","Het","Het", "Hom","Hom","Hom","Hom")),
  Diet     = factor(c("Complete","Complete","Complete","Complete", "Complete","Complete","Complete","Complete", "HalfTyr","HalfTyr","HalfTyr","HalfTyr", "HalfTyr","HalfTyr","HalfTyr","HalfTyr"))
)

# перевірка чи співпадають колонки в data й coldata
all(colnames(data) == rownames(coldata))

# перевели таблицю з текстового формати у цифровий
data_numeric <- matrix(as.numeric(as.matrix(data)), 
                       nrow = nrow(data), 
                       dimnames = dimnames(data))

# валідацій ти синхронізація countData з colData, ініціалізація матриць результатів, Запис матриці дизайну (~ Genotype + Diet + Genotype:Diet)
dds <- DESeqDataSetFromMatrix(
  countData = data_numeric,
  colData   = coldata,
  design    = ~ Genotype + Diet + Genotype:Diet
)

# Оцінка факторів розміру, оцінка дисперсії, побудова GLM та тест Вальда
dds <- DESeq(dds)

# можна переглянути проміжні результати
# resultsNames(dds)

# зберігаємо проміжні результати Hom vs Het (100% Tyr) 
res_disease_100 <- results(dds, contrast = c("Genotype", "Hom", "Het"), alpha = 0.05)
# зберігаємо проміжні результати Hom vs Het (50% Tyr) 
res_disease_50 <- results(dds, contrast = list(c("Genotype_Hom_vs_Het", "GenotypeHom.DietHalfTyr")), alpha = 0.05)
# зберігаємо проміжні результати Hom 50% Tyr vs Hom 100% Tyr
res_diet_hom <- results(dds, contrast = list(c("Diet_HalfTyr_vs_Complete", "GenotypeHom.DietHalfTyr")), alpha = 0.05)
# зберігаємо проміжні результати Het 50% Tyr vs Het 100% Tyr
res_diet_het <- results(dds, contrast = c("Diet", "HalfTyr", "Complete"), alpha = 0.05)



# ____________________________Volcano plot______________________________________
# оголосили змінну
make_volcano <- NA
# Призначили функцію
make_volcano <- function(res_data, plot_title) {

# записуємо 10 найбільш статистично значуших генів
top_genes <- head(rownames(res_data[order(res_data$padj), ]), 10)

# підготовка Volcano plot
p <- EnhancedVolcano(
  res_data,
  lab = rownames(res_data),              # Назви генів беремо з назв рядків
  x = 'log2FoldChange',                  # Вісь X — величина ефекту
  y = 'padj',                            # Вісь Y — скоригований p-value (FDR)
  selectLab = top_genes,
  pCutoff = 0.05,                        # Горизонтальна лінія значущості
  FCcutoff = 1.0,                        # Вертикальні лінії сили ефекту (вдвічі)
  pointSize = 0.5,                       # Розмір точок-генів
  labSize = 2.0,                         # Розмір тексту з назвами топ-генів
  title = plot_title,
  subtitle = 'Volcano plot',
  legendPosition = 'right'
)

p_final <- p + 
  coord_cartesian(
    xlim = c(-4, 4),   # наближаємо вісь X: покаже тільки від -4 до +4
    ylim = c(0, 30)    # наближаємо вісь Y: покаже висоту від 0 до 30
  ) +
  # робимо дрібнішу сітку для зручності
  scale_x_continuous(breaks = seq(-4, 4, by = 1)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +

annotate(
  "text", 
  x = 0, y = 29,                        # Центр по X, але майже самий пік по Y (під залізом)
  label = "CONFIDENTIAL - DRAFT FOR REVIEW ONLY", 
  cex = 4.5,                            # Зменшив розмір шрифту, щоб ліг в один красивий рядок
  fontface = "bold.italic", 
  col = "darkgray", 
  alpha = 0.25                          # Зробив трохи чіткішим (25%), бо зверху він нічого не перекриває
)

return(p_final)
}

# зберігаємо варіант Volcano plot
# Disease effect 100% Tyr
p1 <- make_volcano(res_disease_100, 'Disease Effect: Hom vs Het (100% Tyr)')
ggsave("Volcano_Disease_100.png", plot = p1, width = 12, height = 9, dpi = 300)
# Disease effect 50% Tyr
p2 <- make_volcano(res_disease_50, 'Disease Effect: Hom vs Het (50% Tyr)')
ggsave("Volcano_Disease_50.png", plot = p2, width = 12, height = 9, dpi = 300)
# Diet effect on Het
p3 <- make_volcano(res_diet_het, 'Diet Effect: 50% vs 100% Tyr (Het)')
ggsave("Volcano_Diet_Het.png", plot = p3, width = 12, height = 9, dpi = 300)
# Diet effect on Hom
p4 <- make_volcano(res_diet_hom, 'Diet Effect: 50% vs 100% Tyr (Hom)')
ggsave("Volcano_Diet_Hom.png", plot = p4, width = 12, height = 9, dpi = 300)




# ____________________________MA plot___________________________________________
# Стабілізація дисперсії (Variance Stabilizing Transformation)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd) # Чиста матриця нормалізованих значень

# оголосили функцію
make_ma_plot <- function(res_data, plot_title) {
  
  # Перетворюємо результати в датафрейм
  df <- as.data.frame(res_data)
  
  # Створюємо колонку для кольору (Up, Down, NS)
  df$Significance <- "NS"
  df$Significance[df$padj < 0.05 & df$log2FoldChange > 1] <- "Up"
  df$Significance[df$padj < 0.05 & df$log2FoldChange < -1] <- "Down"
  
  # Захист від NA у padj
  df$Significance[is.na(df$Significance)] <- "NS"
  
  # Рахуємо точну кількість для легенди
  n_up <- sum(df$Significance == "Up")
  n_down <- sum(df$Significance == "Down")
  n_ns <- sum(df$Significance == "NS")
  
  # Перетворюємо у фактор із чіткими рівнями
  df$Significance <- factor(df$Significance, levels = c("Up", "Down", "NS"))
  df$GeneName <- rownames(df)
  
  # Знаходимо ТОП-10 генів для підписів
  top_genes_df <- df[!is.na(df$padj), ]
  top_genes <- head(rownames(top_genes_df[order(top_genes_df$padj), ]), 10)
  df$GeneName <- rownames(df)
  
  # Будуємо графік через базовий ggplot2
  # логарифмуємо вісь X 
  p <- ggplot(df, aes(x = baseMean, y = log2FoldChange, color = Significance)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_x_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      labels = c("1", "10", "100", "1K", "10K", "100K")
    ) +
    coord_cartesian(ylim = c(-4, 4)) + # Обмежуємо вісь Y, як і хотіли
    scale_color_manual(
      values = c("Up" = "#B31B21", "Down" = "#1465B2", "NS" = "darkgray"),
      labels = c("Up" = paste0("Up (", n_up, ")"), 
                 "Down" = paste0("Down (", n_down, ")"), 
                 "NS" = paste0("NS (", n_ns, ")"))
    ) +
    
    # Пунктирні лінії порогів
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, color = "black") +
    
    # Підписи осей та заголовки
    labs(
      title = plot_title,
      subtitle = paste0("Up: ", sum(df$Significance == "Up"), 
                        " | Down: ", sum(df$Significance == "Down")),
      x = "Mean of Normalized Counts (Log10 scale)",
      y = "Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )

  # Додаємо підписи ТОП-10 генів
  library(ggrepel)
  df_top <- df[df$GeneName %in% top_genes & df$Significance != "NS", ]
  if(nrow(df_top) > 0) {
    p <- p + geom_label_repel(
      data = df_top,
      aes(label = GeneName),
      size = 2.5,
      label.padding = unit(0.15, "lines"),
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.2, "lines"),
      color = "black", fill = "white", fontface = "bold",
      show.legend = FALSE
    )
  }
  
  # водяний знак
  p_final <- p + 
    annotate(
      "text", 
      x = max(df$baseMean, na.rm = TRUE) / 100, y = 3.7, 
      label = "CONFIDENTIAL - DRAFT FOR REVIEW ONLY", 
      cex = 4.5, fontface = "bold.italic", col = "darkgray", alpha = 0.25
    )
  
  return(p_final)
}


# MA plot для Hom vs Het (100% Tyr)
ma1 <- make_ma_plot(res_disease_100, 'MA Plot: Disease Effect - Hom vs Het (100% Tyr)')
ggsave("MA_Disease_100.png", plot = ma1, width = 11, height = 8, dpi = 300, bg = "white")

# MA plot для Hom vs Het (50% Tyr)
ma2 <- make_ma_plot(res_disease_50, 'MA Plot: Disease Effect - Hom vs Het (50% Tyr)')
ggsave("MA_Disease_50.png", plot = ma2, width = 11, height = 8, dpi = 300, bg = "white")

# MA plot для Het 50% Tyr vs Het 100% Tyr
ma3 <- make_ma_plot(res_diet_het, 'MA Plot: Diet Effect - Het 50% vs 100% Tyr')
ggsave("MA_Diet_Het.png", plot = ma3, width = 11, height = 8, dpi = 300, bg = "white")

# MA plot для Hom 50% Tyr vs Hom 100% Tyr
ma4 <- make_ma_plot(res_diet_hom, 'MA Plot: Diet Effect - Hom 50% vs 100% Tyr')
ggsave("MA_Diet_Hom.png", plot = ma4, width = 11, height = 8, dpi = 300, bg = "white")





# ____________________Heatmap of significant genes______________________________
vst_mat <- assay(vst(dds, blind = FALSE))















# _______________________________________________________________________
# СТВОРЮЄМО МЕТОД
make_venn_plot <- function(venn_data, colors, file_name) {
  
  # Будуємо базову діаграму Венна
  p <- ggvenn(
    venn_data, 
    fill_color = colors,       # Кольори передаємо як аргумент
    stroke_size = 0.5,         # Товщина ліній кругів
    set_name_size = 4,         # Розмір шрифту для назв груп
    text_size = 4.5            # Розмір цифр всередині кругів
  )
  
  # Додаємо водяний знак
  p_final <- p + 
    annotate(
      "text", 
      x = 0, y = 1.8,          # Координати верхньої межі для ggvenn
      label = "CONFIDENTIAL - DRAFT FOR REVIEW ONLY", 
      cex = 4, 
      fontface = "bold.italic", 
      col = "darkgray", 
      alpha = 0.25
    )
  
  # зберігаємо файл у високій якості
  ggsave(
    filename = file_name, 
    plot = p_final, 
    width = 8, 
    height = 8, 
    dpi = 300, 
    bg = "white"
  )
  
  return(p_final)
}


# витягуємо значущі гени
g_dis_100 <- rownames(subset(res_disease_100, padj < 0.05 & abs(log2FoldChange) > 1))
g_dis_50  <- rownames(subset(res_disease_50,  padj < 0.05 & abs(log2FoldChange) > 1))
g_diet_het <- rownames(subset(res_diet_het,   padj < 0.05 & abs(log2FoldChange) > 1))
g_diet_hom <- rownames(subset(res_diet_hom,   padj < 0.05 & abs(log2FoldChange) > 1))

# Створюємо список для порівняння генотипів
all_effects <- list(
  "Hom_vs_Het_100" = g_dis_100,
  "Hom_vs_Het_50"  = g_dis_50,
  "Diet_Effect_Het" = g_diet_het,
  "Diet_Effect_Hom" = g_diet_hom
)

# Генеруємо матрицю всіх можливих унікальних пар
combinations <- t(combn(names(all_effects), 2))
pairs_matrix <- t(combn(names(all_effects), 2))

# Палітра для графіків
venn_colors <- c("#1465B2", "#B31B21", "#2CA02C", "#FF7F0E", "#9467BD", "#8C564B")

# Запускаємо цикл, який автоматично все порівняє і відмалює
for (i in 1:nrow(pairs_matrix)) {
  name1 <- pairs_matrix[i, 1]
  name2 <- pairs_matrix[i, 2]
  
  # Формуємо пару для поточної діаграми
  current_pair <- list()
  current_pair[[name1]] <- all_effects[[name1]]
  current_pair[[name2]] <- all_effects[[name2]]
  
  # Складаємо назву файлу та заголовок графіка
  file_title <- paste0("Venn_", name1, "_vs_", name2, ".png")
  plot_title <- paste0("Venn: ", name1, " vs ", name2)
  
  # Вибираємо пару кольорів з палітри
  col_pair <- c(venn_colors[i], venn_colors[(i + 2) %% 6 + 1])
  
  # Викликаємо метод
  make_venn_plot(current_pair, col_pair, file_title)
}

