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

# зберігаємо проміжні результати 
res_disease_100 <- results(dds, contrast = c("Genotype", "Hom", "Het"), alpha = 0.05)

# записуємо 10 найбільш статистично значуших генів (для подальшого підпису їх на графіках)
top_genes <- head(rownames(res_disease_100[order(res_disease_100$padj), ]), 10)

# підготовка Volcano plot
p <- EnhancedVolcano(
  res_disease_100,
  lab = rownames(res_disease_100),       # Назви генів беремо з назв рядків
  x = 'log2FoldChange',                  # Вісь X — величина ефекту
  y = 'padj',                            # Вісь Y — скоригований p-value (FDR)
  selectLab = top_genes,
  pCutoff = 0.05,                        # Горизонтальна лінія значущості
  FCcutoff = 1.0,                        # Вертикальні лінії сили ефекту (вдвічі)
  pointSize = 0.5,                       # Розмір точок-генів
  labSize = 2.0,                         # Розмір тексту з назвами топ-генів
  title = 'Disease Effect: Hom vs Het (100% Tyr)',
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

# зберігаємо варіант Volcano plot
ggsave(
  filename = "Volcano_Disease.png", 
  plot = p_final, 
  width = 12, 
  height = 9, 
  dpi = 300
)