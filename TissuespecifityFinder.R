BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(GenomicRanges)
library(IRanges)
library(dplyr)

#### E
peaks <- read.table("/path/to/your/output_peaks.narrowPeak", header = FALSE) # REPLACE WITH YOUR ACTUAL PATH !!!!!

accessible_regions <- read.table("/path/to/your/*_multi_tissue_accessible_regions.bed", header = FALSE) # REPLACE WITH YOUR ACTUAL PATH !!!!!
print(accessible_regions)
# Load tissue-specific accessible chromatin regions as a BED file
accessible_chromatin <- "/path/to/your/*_multi_tissue_accessible_regions.bed" # REPLACE WITH YOUR ACTUAL PATH !!!!!
accessible_regions <- readPeakFile(accessible_chromatin)

# Convertir los picos a un objeto GRanges
peaks_gr <- GRanges(
  seqnames = peaks$V1,  # Cromosoma
  ranges = IRanges(start = peaks$V2, end = peaks$V3)  # Inicio y Fin
)

# Comprobacion
seqnames(accessible_regions)  # Acceder a las secuencias (cromosomas)
ranges(accessible_regions)    # Acceder a los rangos (start y end)
strand(accessible_regions)    # Acceder a la orientación (strand)

tissues <- accessible_regions$V11 #Add la columna de tejidos de mi GRange 

# Crear un nuevo GRanges con los valores extraídos más la columna de tejidos
accessible_gr <- GRanges(
  seqnames = seqnames(accessible_regions),  # Secuencias (cromosomas)
  ranges = ranges(accessible_regions),      # Rangos (start y end)
  strand = strand(accessible_regions),      # Cadena
  tissue = tissues                          # Columna de tejido
)

# Verificar el nuevo GRanges
print(accessible_gr)
print(length(accessible_gr))

# Crear una lista de tejidos para procesar
tissue_list <- unique(accessible_gr$tissue)

# Crear un data.frame vacío para almacenar los resultados
accessibility_data <- data.frame(Tissue = character(), Accessible_sites = integer())

# Calcular el número de sitios accesibles en cada tejido
for (tissue in tissue_list) {
  tissue_regions <- accessible_gr[accessible_gr$tissue == tissue]
  overlap <- findOverlaps(peaks_gr, tissue_regions)
  
  # Añadir los resultados al data.frame
  accessibility_data <- rbind(accessibility_data, data.frame(Tissue = tissue, Accessible_sites = length(overlap)))
}

# Ver los resultados que esten bien
print(accessibility_data)

# Visualizar la especificidad tisular de la accesibilidad cromática
library(ggplot2)

accessibility_data <- accessibility_data %>%
  arrange(Accessible_sites)  # Ordenar de menor a mayor

# En la siguiente linea donde pone reorder tmbn es pa ordenar si no lo quieres quitalo
ggplot(accessibility_data, aes(x = reorder(Tissue, Accessible_sites), y = Accessible_sites, fill = Tissue)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "Tissue specificity of chromatin accessibility in binding sites",
    x = "Tissue", y = "Number of chromatin accessible sites"
  ) +
  scale_fill_brewer(palette = "Set3")