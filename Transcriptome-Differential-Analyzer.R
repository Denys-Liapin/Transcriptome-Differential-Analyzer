library(readxl)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)

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

make_volcano <- NA
make_volcano <- function(res_data, plot_title) {

# записуємо 10 найбільш статистично значуших генів (для подальшого підпису їх на графіках)
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
  scale_y_continuous(breaks = seq(0, 30, by = 5))

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