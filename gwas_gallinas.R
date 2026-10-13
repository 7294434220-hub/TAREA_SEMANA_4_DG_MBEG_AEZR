# GWAS de peso corporal en gallinas indigenas de Etiopia
# Adaptacion del flujo de Flavio Lozano-Isla para la tarea de Alexei Zelada.
# Articulo fuente: Gheyas et al. (2022), doi:10.1038/s41597-022-01129-4.
# Datos: chickens_small.vcf.gz + Online-only Table 1 del mismo articulo.
# Abra su proyecto .Rproj y coloque este script junto al VCF. Pulse Source.
# El VCF original no se modifica. No se necesitan los archivos de arroz.

# ==== Configuracion y paquetes ====
options(stringsAsFactors = FALSE, timeout = 600, bigmemory.typecast.warning = FALSE)
set.seed(20261013)
CONFIG <- list(
  archivo_vcf = "chickens_small.vcf.gz",
  carpeta_salida = "resultados_gwas_gallinas",
  maf_min = 0.05,
  faltantes_max = 0.20,
  n_pc = 3L,
  alpha = 0.05,
  bloque_vcf = 2000L,
  reutilizar = TRUE,
  max_loop_farmcpu = 10L
)
nucleos <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(nucleos) || nucleos < 1L) nucleos <- 1L
cores <- max(1L, min(4L, as.integer(floor(nucleos * 0.8))))
# Tope de 4 hilos para un equipo personal; puede ajustar esta linea.

paquetes <- c("rMVP", "bigmemory")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) {
  install.packages(faltan, repos = "https://cloud.r-project.org")
}
if (!all(vapply(paquetes, requireNamespace, logical(1), quietly = TRUE))) {
  stop("No se pudieron instalar rMVP/bigmemory. Revise los mensajes de instalacion; ",
       "en Windows puede necesitar Rtools si CRAN no ofrece un binario para su R.")
}
suppressPackageStartupMessages({ library(rMVP); library(bigmemory) })

# Si Source se ejecuta desde otra carpeta, busque el VCF junto al script.
localizar_vcf <- function(nombre) {
  candidatos <- c(file.path(getwd(), nombre), file.path(getwd(), "datasets", nombre))
  marco <- sys.frames()
  for (f in marco) {
    if (!is.null(f$ofile) && length(f$ofile) == 1L) {
      candidatos <- c(candidatos, file.path(dirname(f$ofile), nombre))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  ruta_script <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(ruta_script)) {
    candidatos <- c(candidatos, file.path(dirname(ruta_script[1]), nombre))
  }
  candidatos <- unique(candidatos[file.exists(candidatos)])
  if (!length(candidatos)) {
    stop("No se encuentra ", nombre, ". Abra el .Rproj de la tarea y coloque ",
         "el script en la misma carpeta que ese VCF. Carpeta actual: ", getwd())
  }
  normalizePath(candidatos[1], winslash = "/", mustWork = TRUE)
}
vcf_original <- localizar_vcf(CONFIG$archivo_vcf)
proyecto <- dirname(vcf_original)
salida <- file.path(proyecto, CONFIG$carpeta_salida)
carpetas <- c("procesados", "tablas", "graficos/01_pca", "graficos/02_densidad",
              "graficos/03_fenotipo", "graficos/04_manhattan_clase",
              "graficos/05_manhattan_bonferroni", "graficos/06_qq")
for (x in carpetas) dir.create(file.path(salida, x), recursive = TRUE, showWarnings = FALSE)
procesados <- file.path(salida, "procesados")
tablas <- file.path(salida, "tablas")
graficos <- file.path(salida, "graficos")
cat("Proyecto:", proyecto, "\nHilos:", cores,
    "\nrMVP:", as.character(packageVersion("rMVP")), "\n")
guardar_csv <- function(x, archivo) {
  write.csv(x, file.path(tablas, archivo), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
# Compatibilidad con nombres cpu/cpus/ncpus entre versiones de rMVP.
llamar_rmvp <- function(fun, argumentos, hilos = cores) {
  formales <- names(formals(fun))
  hilo <- intersect(c("cpu", "cpus", "ncpus", "threads"), formales)
  if (length(hilo)) argumentos[[hilo[1]]] <- hilos
  do.call(fun, argumentos[names(argumentos) %in% formales])
}

# ==== Metadatos reales del articulo ====
# Datos transcritos de Online-only Table 1. Los campos vacios siguen siendo NA.
# URLs fuente: https://www.nature.com/articles/s41597-022-01129-4/tables/2
#             https://pmc.ncbi.nlm.nih.gov/articles/PMC8844291/
metadatos_texto <- paste(c(
  "ID_articulo,ENA_estudio,ENA_muestra,Region,Distrito,Localidad,Sexo,Edad_meses,Peso_kg",
  "ABB-2H_Pink,PRJEB39275,SAMEA7050617,Amhara,Fagita Lekoma,Batambie,F,12,1.5",
  "ABB-3C,PRJEB39275,SAMEA7050618,Amhara,Fagita Lekoma,Batambie,M,8,1.75",
  "ABB-4H1,PRJEB39275,SAMEA7050619,Amhara,Fagita Lekoma,Batambie,F,7,1",
  "ABB-4H2,PRJEB39275,SAMEA7050620,Amhara,Fagita Lekoma,Batambie,F,7,1",
  "ABB-5H1,PRJEB39275,SAMEA7050621,Amhara,Fagita Lekoma,Batambie,F,12,1.8",
  "ABB-6Cock1,PRJEB39275,SAMEA7050622,Amhara,Fagita Lekoma,Batambie,M,6,1",
  "ABB-6H2,PRJEB39275,SAMEA7050623,Amhara,Fagita Lekoma,Batambie,F,12,1.8",
  "ABB-6H3,PRJEB39275,SAMEA7050624,Amhara,Fagita Lekoma,Batambie,F,10,1",
  "ABS-1C,PRJEB39275,SAMEA7050625,Amhara,Banja,Surta,M,8,1.5",
  "ABS-1H,PRJEB39275,SAMEA7050626,Amhara,Banja,Surta,F,10,1.5",
  "ABS-2H,PRJEB39275,SAMEA7050627,Amhara,Banja,Surta,F,7,1.5",
  "ABS-3C,PRJEB39275,SAMEA7050628,Amhara,Banja,Surta,M,5,1.5",
  "ABS-4C,PRJEB39275,SAMEA7050629,Amhara,Banja,Surta,M,6,1.8",
  "ABS-5H1,PRJEB39275,SAMEA7050630,Amhara,Banja,Surta,F,12,1",
  "ABS-5H2,PRJEB39275,SAMEA7050631,Amhara,Banja,Surta,F,14,1.5",
  "ABS-6H,PRJEB39275,SAMEA7050632,Amhara,Banja,Surta,F,9,2",
  "ABS-7H,PRJEB39275,SAMEA7050633,Amhara,Banja,Surta,F,8,2",
  "AFA-1C,PRJEB39275,SAMEA7050634,Amhara,Fagita Lekoma,Amesha Shinkuri,M,12,1",
  "AFA-1H,PRJEB39275,SAMEA7050635,Amhara,Fagita Lekoma,Amesha Shinkuri,F,13,1.5",
  "AFA-2H,PRJEB39275,SAMEA7050636,Amhara,Fagita Lekoma,Amesha Shinkuri,F,12,1.5",
  "AFA-3C,PRJEB39275,SAMEA7050637,Amhara,Fagita Lekoma,Amesha Shinkuri,M,5,1",
  "AFA-4C,PRJEB39275,SAMEA7050638,Amhara,Fagita Lekoma,Amesha Shinkuri,M,6,1.5",
  "AFA-4H,PRJEB39275,SAMEA7050639,Amhara,Fagita Lekoma,Amesha Shinkuri,F,12,1.5",
  "AFA-5C,PRJEB39275,SAMEA7050640,Amhara,Fagita Lekoma,Amesha Shinkuri,M,7,1.5",
  "AFA-6H,PRJEB39275,SAMEA7050641,Amhara,Fagita Lekoma,Amesha Shinkuri,F,12,1",
  "AFA-7H,PRJEB39275,SAMEA7050642,Amhara,Fagita Lekoma,Amesha Shinkuri,F,8,1",
  "AFA-8H,PRJEB39275,SAMEA7050643,Amhara,Fagita Lekoma,Amesha Shinkuri,F,7,1",
  "AFDH-C1,PRJEB39275,SAMEA7050644,Afar,Dulecha,Hugub,M,24,1.4",
  "AFDH-C3,PRJEB39275,SAMEA7050645,Afar,Dulecha,Hugub,M,,1",
  "AFDH-C5,PRJEB39275,SAMEA7050646,Afar,Dulecha,Hugub,F,8,0.9",
  "AFDH-C7,PRJEB39275,SAMEA7050647,Afar,Dulecha,Hugub,M,,",
  "AFDH-H1,PRJEB39275,SAMEA7050648,Afar,Dulecha,Hugub,M,,",
  "AFDH-H2,PRJEB39275,SAMEA7050649,Afar,Dulecha,Hugub,F,12,1.6",
  "AFDH-H4,PRJEB39275,SAMEA7050650,Afar,Dulecha,Hugub,F,,",
  "AFDH-H6,PRJEB39275,SAMEA7050651,Afar,Dulecha,Hugub,F,5,0.8",
  "AFDH-H8,PRJEB39275,SAMEA7050652,Afar,Dulecha,Hugub,F,7,0.9",
  "AFDH-H9,PRJEB39275,SAMEA7050653,Afar,Dulecha,Hugub,F,24,1.3",
  "AFDK-C1,PRJEB39275,SAMEA7050654,Afar,Dulecha,Kefis,M,,",
  "AFDK-C2,PRJEB39275,SAMEA7050655,Afar,Dulecha,Kefis,M,,1.4",
  "AFDK-C3,PRJEB39275,SAMEA7050656,Afar,Dulecha,Kefis,M,,1.05",
  "AFDK-H10,PRJEB39275,SAMEA7050657,Afar,Dulecha,Kefis,F,12,0.97",
  "AFDK-H4,PRJEB39275,SAMEA7050658,Afar,Dulecha,Kefis,F,7,1.7",
  "AFDK-H5,PRJEB39275,SAMEA7050659,Afar,Dulecha,Kefis,F,7,1",
  "AFDK-H6,PRJEB39275,SAMEA7050660,Afar,Dulecha,Kefis,F,8,0.6",
  "AFDK-H7,PRJEB39275,SAMEA7050661,Afar,Dulecha,Kefis,F,8,1",
  "AFDK-H8,PRJEB39275,SAMEA7050662,Afar,Dulecha,Kefis,F,12,1.4",
  "AFDK-H9,PRJEB39275,SAMEA7050663,Afar,Dulecha,Kefis,F,,1.2",
  "AFG-1H,PRJEB39275,SAMEA7050664,Amhara,Fagita Lekoma,Gafera,F,6,1.5",
  "AFG-2C1,PRJEB39275,SAMEA7050665,Amhara,Fagita Lekoma,Gafera,M,6,1.5",
  "AFG-2H,PRJEB39275,SAMEA7050666,Amhara,Fagita Lekoma,Gafera,F,6,1.5",
  "AFG-3C,PRJEB39275,SAMEA7050667,Amhara,Fagita Lekoma,Gafera,M,5,1.5",
  "AFG-3H,PRJEB39275,SAMEA7050668,Amhara,Fagita Lekoma,Gafera,F,8,1.5",
  "AFG-4H_P_Y,PRJEB39275,SAMEA7050669,Amhara,Fagita Lekoma,Gafera,F,7,1.5",
  "AFG-5C,PRJEB39275,SAMEA7050670,Amhara,Fagita Lekoma,Gafera,M,6,1.5",
  "AFG-5H,PRJEB39275,SAMEA7050671,Amhara,Fagita Lekoma,Gafera,F,4,1",
  "AFG-6H,PRJEB39275,SAMEA7050672,Amhara,Fagita Lekoma,Gafera,F,4,1",
  "AFG-8C,PRJEB39275,SAMEA7050673,Amhara,Fagita Lekoma,Gafera,M,5,1",
  "AGT-10C,PRJEB39275,SAMEA7050674,Amhara,Gondar Zuria,Tsion Teguaz,M,8,1",
  "AGT-1C,PRJEB39275,SAMEA7050675,Amhara,Gondar Zuria,Tsion Teguaz,M,19,1.5",
  "AGT-2C,PRJEB39275,SAMEA7050676,Amhara,Gondar Zuria,Tsion Teguaz,M,8,1.5",
  "AGT-3H,PRJEB39275,SAMEA7050677,Amhara,Gondar Zuria,Tsion Teguaz,F,12,1.5",
  "AGT-4H,PRJEB39275,SAMEA7050678,Amhara,Gondar Zuria,Tsion Teguaz,F,12,1",
  "AGT-5H,PRJEB39275,SAMEA7050679,Amhara,Gondar Zuria,Tsion Teguaz,F,12,1",
  "AGT-6H,PRJEB39275,SAMEA7050680,Amhara,Gondar Zuria,Tsion Teguaz,F,15,1",
  "AGT-7H,PRJEB39275,SAMEA7050681,Amhara,Gondar Zuria,Tsion Teguaz,F,18,1",
  "AGT-8H,PRJEB39275,SAMEA7050682,Amhara,Gondar Zuria,Tsion Teguaz,F,9,1",
  "AGT-9C,PRJEB39275,SAMEA7050683,Amhara,Gondar Zuria,Tsion Teguaz,M,9,1.5",
  "AK025A-10C-059,PRJEB39275,SAMEA7050684,Amhara,Kalu,0-25Adane,M,6,1.1",
  "AK025A-1HC-55,PRJEB39275,SAMEA7050685,Amhara,Kalu,0-25Adane,F,12,1.1",
  "AK025A-2H-149,PRJEB39275,SAMEA7050686,Amhara,Kalu,0-25Adane,F,8,0.88",
  "AK025A-3H-85,PRJEB39275,SAMEA7050687,Amhara,Kalu,0-25Adane,F,7,1.4",
  "AK025A-4H-070,PRJEB39275,SAMEA7050688,Amhara,Kalu,0-25Adane,F,,1.06",
  "AK025A-5H-71,PRJEB39275,SAMEA7050689,Amhara,Kalu,0-25Adane,F,6,0.86",
  "AK025A-6H-97,PRJEB39275,SAMEA7050690,Amhara,Kalu,0-25Adane,F,11,1.14",
  "AK025A-7C1,PRJEB39275,SAMEA7050691,Amhara,Kalu,0-25Adane,M,8,1.2",
  "AK025A-8C,PRJEB39275,SAMEA7050692,Amhara,Kalu,0-25Adane,M,,",
  "AK025A-8C-50,PRJEB39275,SAMEA7050693,Amhara,Kalu,0-25Adane,M,8,1.34",
  "AKA-10H-11,PRJEB39275,SAMEA7050694,Amhara,Kalu,Arabo,F,,1.1",
  "AKA-1C,PRJEB39275,SAMEA7050695,Amhara,Kalu,Arabo,M,,0.9",
  "AKA-2H,PRJEB39275,SAMEA7050696,Amhara,Kalu,Arabo,F,5,0.82",
  "AKA-3C,PRJEB39275,SAMEA7050697,Amhara,Kalu,Arabo,M,5,1.14",
  "AKA-4H,PRJEB39275,SAMEA7050698,Amhara,Kalu,Arabo,F,12,1.2",
  "AKA-5H,PRJEB39275,SAMEA7050699,Amhara,Kalu,Arabo,F,6,0.85",
  "AKA-6H,PRJEB39275,SAMEA7050700,Amhara,Kalu,Arabo,F,8,1.2",
  "AKA-7C,PRJEB39275,SAMEA7050701,Amhara,Kalu,Arabo,M,8,1.2",
  "AKA-8H,PRJEB39275,SAMEA7050702,Amhara,Kalu,Arabo,F,14,1.3",
  "AKA-9C,PRJEB39275,SAMEA7050703,Amhara,Kalu,Arabo,M,8,1.5",
  "AMAM-10H-108,PRJEB39275,SAMEA7050704,Amhara,Menz Gera Midir,Alfa Midir/05/,F,7,1",
  "AMAM-1H-031,PRJEB39275,SAMEA7050705,Amhara,Menz Gera Midir,Alfa Midir/05/,F,24,1.12",
  "AMAM-2C-151,PRJEB39275,SAMEA7050706,Amhara,Menz Gera Midir,Alfa Midir/05/,M,12,1.01",
  "AMAM-3H-124,PRJEB39275,SAMEA7050707,Amhara,Menz Gera Midir,Alfa Midir/05/,F,14,0.7",
  "AMAM-4H-115,PRJEB39275,SAMEA7050708,Amhara,Menz Gera Midir,Alfa Midir/05/,F,8,1.26",
  "AMAM-5C-147,PRJEB39275,SAMEA7050709,Amhara,Menz Gera Midir,Alfa Midir/05/,M,9,1.38",
  "AMAM-6H-120,PRJEB39275,SAMEA7050710,Amhara,Menz Gera Midir,Alfa Midir/05/,F,9,1.08",
  "AMAM-7C-116,PRJEB39275,SAMEA7050711,Amhara,Menz Gera Midir,Alfa Midir/05/,M,6,1.16",
  "AMAM-8H-101,PRJEB39275,SAMEA7050712,Amhara,Menz Gera Midir,Alfa Midir/05/,F,12,0.78",
  "AMAM-9C-107,PRJEB39275,SAMEA7050713,Amhara,Menz Gera Midir,Alfa Midir/05/,M,9,1.52",
  "AMNA-10C,PRJEB39275,SAMEA7050714,Amhara,Menz Gera Midir,Negasi Amba/07/,M,12,1.32",
  "AMNA-1H-102,PRJEB39275,SAMEA7050715,Amhara,Menz Gera Midir,Negasi Amba/07/,F,12,1.18",
  "AMNA-2C-104,PRJEB39275,SAMEA7050716,Amhara,Menz Gera Midir,Negasi Amba/07/,M,12,1.2",
  "AMNA-3H-170,PRJEB39275,SAMEA7050717,Amhara,Menz Gera Midir,Negasi Amba/07/,F,12,1.04",
  "AMNA-4H-121,PRJEB39275,SAMEA7050718,Amhara,Menz Gera Midir,Negasi Amba/07/,F,12,1.1",
  "AMNA-5C-138,PRJEB39275,SAMEA7050719,Amhara,Menz Gera Midir,Negasi Amba/07/,M,15,0.98",
  "AMNA-6H-165,PRJEB39275,SAMEA7050720,Amhara,Menz Gera Midir,Negasi Amba/07/,F,11.5,1",
  "AMNA-7H-144,PRJEB39275,SAMEA7050721,Amhara,Menz Gera Midir,Negasi Amba/07/,F,9,1.2",
  "AMNA-8H,PRJEB39275,SAMEA7050722,Amhara,Menz Gera Midir,Negasi Amba/07/,M,,",
  "AMNA-9C-171,PRJEB39275,SAMEA7050723,Amhara,Menz Gera Midir,Negasi Amba/07/,M,12,1.7",
  "ASA-10H,PRJEB39275,SAMEA7050724,Amhara,South Achefer,Ashuda,F,8,1.5",
  "ASA-1C,PRJEB39275,SAMEA7050725,Amhara,South Achefer,Ashuda,M,6,1.5",
  "ASA-2H,PRJEB39275,SAMEA7050726,Amhara,South Achefer,Ashuda,F,8,1",
  "ASA-3H,PRJEB39275,SAMEA7050727,Amhara,South Achefer,Ashuda,F,8,1",
  "ASA-4C,PRJEB39275,SAMEA7050728,Amhara,South Achefer,Ashuda,M,5,1.5",
  "ASA-5H,PRJEB39275,SAMEA7050729,Amhara,South Achefer,Ashuda,F,8,1",
  "ASA-6H,PRJEB39275,SAMEA7050730,Amhara,South Achefer,Ashuda,F,9,1",
  "ASA-7H,PRJEB39275,SAMEA7050731,Amhara,South Achefer,Ashuda,F,11,1.5",
  "ASA-8C,PRJEB39275,SAMEA7050732,Amhara,South Achefer,Ashuda,M,7,1.5",
  "ASA-9C,PRJEB39275,SAMEA7050733,Amhara,South Achefer,Ashuda,M,6,1.5",
  "ASD-1C,PRJEB39275,SAMEA7050734,Amhara,South Achefer,Dikuli,M,9,1.5",
  "ASD-1H,PRJEB39275,SAMEA7050735,Amhara,South Achefer,Dikuli,F,9,1",
  "ASD-2C,PRJEB39275,SAMEA7050736,Amhara,South Achefer,Dikuli,M,18,2",
  "ASD-2H,PRJEB39275,SAMEA7050737,Amhara,South Achefer,Dikuli,F,9,1",
  "ASD-3H,PRJEB39275,SAMEA7050738,Amhara,South Achefer,Dikuli,F,9,1.5",
  "ASD-4H,PRJEB39275,SAMEA7050739,Amhara,South Achefer,Dikuli,F,10,1",
  "ASD-5C,PRJEB39275,SAMEA7050740,Amhara,South Achefer,Dikuli,M,12,1",
  "ASD-6H,PRJEB39275,SAMEA7050741,Amhara,South Achefer,Dikuli,F,12,1.5",
  "ASD-7H,PRJEB39275,SAMEA7050742,Amhara,South Achefer,Dikuli,F,12,1.5",
  "ASD-8C,PRJEB39275,SAMEA7050743,Amhara,South Achefer,Dikuli,M,9,1.5",
  "BGDG-10H,PRJEB39275,SAMEA7050744,Gumuz,Dibate,Gesses,F,,",
  "BGDG-1C-080,PRJEB39275,SAMEA7050745,Gumuz,Dibate,Gesses,M,,1.48",
  "BGDG-2C-95,PRJEB39275,SAMEA7050746,Gumuz,Dibate,Gesses,M,,0.92",
  "BGDG-3H-12,PRJEB39275,SAMEA7050747,Gumuz,Dibate,Gesses,F,8,1",
  "BGDG-4C-51,PRJEB39275,SAMEA7050748,Gumuz,Dibate,Gesses,M,,1.9",
  "BGDG-5H-2,PRJEB39275,SAMEA7050749,Gumuz,Dibate,Gesses,F,12,1",
  "BGDG-6C,PRJEB39275,SAMEA7050750,Gumuz,Dibate,Gesses,M,,",
  "BGDG-7H-23,PRJEB39275,SAMEA7050751,Gumuz,Dibate,Gesses,F,,1",
  "BGDG-8H,PRJEB39275,SAMEA7050752,Gumuz,Dibate,Gesses,F,7,1.25",
  "BGDG-9H,PRJEB39275,SAMEA7050753,Gumuz,Dibate,Gesses,F,9,1.4",
  "BGDK-10C-14,PRJEB39275,SAMEA7050754,Gumuz,Dibate,Kido,M,6,1.2",
  "BGDK-11H-56,PRJEB39275,SAMEA7050755,Gumuz,Dibate,Kido,F,,0.82",
  "BGDK-12C,PRJEB39275,SAMEA7050756,Gumuz,Dibate,Kido,M,,0.84",
  "BGDK-3C,PRJEB39275,SAMEA7050757,Gumuz,Dibate,Kido,M,,",
  "BGDK-5H-82,PRJEB39275,SAMEA7050758,Gumuz,Dibate,Kido,F,6,1.2",
  "BGDK-6H,PRJEB39275,SAMEA7050759,Gumuz,Dibate,Kido,F,,",
  "BGDK-7H,PRJEB39275,SAMEA7050760,Gumuz,Dibate,Kido,F,,",
  "BGDK-8H,PRJEB39275,SAMEA7050761,Gumuz,Dibate,Kido,F,,",
  "BGDK-9H,PRJEB39275,SAMEA7050762,Gumuz,Dibate,Kido,F,,0.6",
  "ODB-10H-003,PRJEB39275,SAMEA7050763,Oromia,Dugda,Bekele Girissa,F,6,1",
  "ODB-1H-076,PRJEB39275,SAMEA7050764,Oromia,Dugda,Bekele Girissa,F,5,1",
  "ODB-2C-083,PRJEB39275,SAMEA7050765,Oromia,Dugda,Bekele Girissa,M,6,1",
  "ODB-3H,PRJEB39275,SAMEA7050766,Oromia,Dugda,Bekele Girissa,F,9,1",
  "ODB-4H-002,PRJEB39275,SAMEA7050767,Oromia,Dugda,Bekele Girissa,F,12,1.5",
  "ODB-5H,PRJEB39275,SAMEA7050768,Oromia,Dugda,Bekele Girissa,F,5,1",
  "ODB-6C-093,PRJEB39275,SAMEA7050769,Oromia,Dugda,Bekele Girissa,M,12,1",
  "ODB-7C,PRJEB39275,SAMEA7050770,Oromia,Dugda,Bekele Girissa,M,5,1",
  "ODB-8H-013,PRJEB39275,SAMEA7050771,Oromia,Dugda,Bekele Girissa,F,5,1",
  "ODB-9H,PRJEB39275,SAMEA7050772,Oromia,Dugda,Bekele Girissa,F,6,1",
  "ODS-10H-068,PRJEB39275,SAMEA7050773,Oromia,Dugda,Shubi Gemo,F,7,1",
  "ODS-1H,PRJEB39275,SAMEA7050774,Oromia,Dugda,Shubi Gemo,F,,",
  "ODS-2H-88,PRJEB39275,SAMEA7050775,Oromia,Dugda,Shubi Gemo,F,,1",
  "ODS-3C,PRJEB39275,SAMEA7050776,Oromia,Dugda,Shubi Gemo,M,,",
  "ODS-4H-010,PRJEB39275,SAMEA7050777,Oromia,Dugda,Shubi Gemo,F,5,1",
  "ODS-5C-086,PRJEB39275,SAMEA7050778,Oromia,Dugda,Shubi Gemo,M,5,1",
  "ODS-6C-098,PRJEB39275,SAMEA7050779,Oromia,Dugda,Shubi Gemo,M,7,1",
  "ODS-7H,PRJEB39275,SAMEA7050780,Oromia,Dugda,Shubi Gemo,F,,",
  "ODS-8H-042,PRJEB39275,SAMEA7050781,Oromia,Dugda,Shubi Gemo,F,7,1.5",
  "ODS-9H-016,PRJEB39275,SAMEA7050782,Oromia,Dugda,Shubi Gemo,F,12,1",
  "SDK-10C-004,PRJEB39275,SAMEA7050783,SNNPR,Dara,Kumato,M,10,1",
  "SDK-1H-080,PRJEB39275,SAMEA7050784,SNNPR,Dara,Kumato,F,12,",
  "SDK-2C,PRJEB39275,SAMEA7050785,SNNPR,Dara,Kumato,M,18,1.7",
  "SDK-3H-054,PRJEB39275,SAMEA7050786,SNNPR,Dara,Kumato,F,18,1.7",
  "SDK-4C-062,PRJEB39275,SAMEA7050787,SNNPR,Dara,Kumato,M,12,1.5",
  "SDK-5H,PRJEB39275,SAMEA7050788,SNNPR,Dara,Kumato,F,,",
  "SDK-6H-052,PRJEB39275,SAMEA7050789,SNNPR,Dara,Kumato,F,12,1.5",
  "SDK-7H,PRJEB39275,SAMEA7050790,SNNPR,Dara,Kumato,F,12,1.6",
  "SDK-8C-021,PRJEB39275,SAMEA7050791,SNNPR,Dara,Kumato,M,12,1.6",
  "SDK-9H-026,PRJEB39275,SAMEA7050792,SNNPR,Dara,Kumato,F,10,1",
  "SDL-10H-017,PRJEB39275,SAMEA7050793,SNNPR,Dara,Loya,F,7,1.5",
  "SDL-1C-091,PRJEB39275,SAMEA7050794,SNNPR,Dara,Loya,M,12,1.5",
  "SDL-2H-064,PRJEB39275,SAMEA7050795,SNNPR,Dara,Loya,F,7,1.5",
  "SDL-3C-019,PRJEB39275,SAMEA7050796,SNNPR,Dara,Loya,M,8,2",
  "SDL-4H-066,PRJEB39275,SAMEA7050797,SNNPR,Dara,Loya,F,12,1",
  "SDL-5H-008,PRJEB39275,SAMEA7050798,SNNPR,Dara,Loya,F,7,1",
  "SDL-6H-057,PRJEB39275,SAMEA7050799,SNNPR,Dara,Loya,F,8,1",
  "SDL-7C-035,PRJEB39275,SAMEA7050800,SNNPR,Dara,Loya,F,9,1.5",
  "SDL-8C-090,PRJEB39275,SAMEA7050801,SNNPR,Dara,Loya,M,7,1",
  "SDL-9H-013,PRJEB39275,SAMEA7050802,SNNPR,Dara,Loya,F,8,1.5",
  "TGENF_110,PRJEB39275,SAMEA7050803,Tigray,Enderta,Meseret,F,25,1.25",
  "TGENF_128,PRJEB39275,SAMEA7050804,Tigray,Enderta,Meseret,F,26,0.8",
  "TGENF_134,PRJEB39275,SAMEA7050805,Tigray,Enderta,Meseret,F,25,1",
  "TGENF_145,PRJEB39275,SAMEA7050806,Tigray,Enderta,Meseret,F,26,1.3",
  "TGENF_177,PRJEB39275,SAMEA7050807,Tigray,Enderta,Meseret,F,28,0.85",
  "TGENF_178,PRJEB39275,SAMEA7050808,Tigray,Enderta,Meseret,F,30,1.7",
  "TGENM_135,PRJEB39275,SAMEA7050809,Tigray,Enderta,Meseret,M,28,1.55",
  "TGENM_146,PRJEB39275,SAMEA7050810,Tigray,Enderta,Meseret,M,30,1.8",
  "TGENM_150,PRJEB39275,SAMEA7050811,Tigray,Enderta,Meseret,M,26,1.45",
  "TGENM_175,PRJEB39275,SAMEA7050812,Tigray,Enderta,Meseret,M,26,1.2",
  "TMLHA_13H,PRJEB39275,SAMEA7050813,Tigray,Merebleke,Hadush Adi,F,8,1.2",
  "TMLHA-04C,PRJEB39275,SAMEA7050814,Tigray,Merebleke,Hadush Adi,M,,2.6",
  "TMLHA-19H,PRJEB39275,SAMEA7050815,Tigray,Merebleke,Hadush Adi,F,,1.2",
  "TMLHA-54H,PRJEB39275,SAMEA7050816,Tigray,Merebleke,Hadush Adi,F,,",
  "TMLHA-66C,PRJEB39275,SAMEA7050817,Tigray,Merebleke,Hadush Adi,M,,1.83",
  "TMLHA-86H,PRJEB39275,SAMEA7050818,Tigray,Merebleke,Hadush Adi,F,12,1.3",
  "TMLHA-90C,PRJEB39275,SAMEA7050819,Tigray,Merebleke,Hadush Adi,M,,",
  "TMLHA-91C,PRJEB39275,SAMEA7050820,Tigray,Merebleke,Hadush Adi,F,,1",
  "TMLHA-96H,PRJEB39275,SAMEA7050821,Tigray,Merebleke,Hadush Adi,F,10,1.43",
  "TMLM-_C94,PRJEB39275,SAMEA7050822,Tigray,Merebleke,Mihquan,M,12,1.45",
  "TMLM-26C,PRJEB39275,SAMEA7050823,Tigray,Merebleke,Mihquan,M,,",
  "TMLM-6H,PRJEB39275,SAMEA7050824,Tigray,Merebleke,Mihquan,F,12,1.1",
  "TMLM-C30,PRJEB39275,SAMEA7050825,Tigray,Merebleke,Mihquan,M,12,2.1",
  "TMLM-C89,PRJEB39275,SAMEA7050826,Tigray,Merebleke,Mihquan,M,10,1.1",
  "TMLM-H18,PRJEB39275,SAMEA7050827,Tigray,Merebleke,Mihquan,F,12,1.65",
  "TMLM-H53,PRJEB39275,SAMEA7050828,Tigray,Merebleke,Mihquan,F,10,1.6",
  "TMLM-H57,PRJEB39275,SAMEA7050829,Tigray,Merebleke,Mihquan,F,10,1.1",
  "TMLM-H72,PRJEB39275,SAMEA7050830,Tigray,Merebleke,Mihquan,F,12,1.3",
  "TMLM-H84,PRJEB39275,SAMEA7050831,Tigray,Merebleke,Mihquan,F,12,1.57",
  "TSSG-07H,PRJEB39275,SAMEA7050832,Tigray,Sahareti Samire,Gijet,F,,",
  "TSSG-29H,PRJEB39275,SAMEA7050833,Tigray,Sahareti Samire,Gijet,F,15,1.7",
  "TSSG-37C,PRJEB39275,SAMEA7050834,Tigray,Sahareti Samire,Gijet,M,9,1.5",
  "TSSG-40H,PRJEB39275,SAMEA7050835,Tigray,Sahareti Samire,Gijet,F,7,0.9",
  "TSSG-43C,PRJEB39275,SAMEA7050836,Tigray,Sahareti Samire,Gijet,M,15,1.2",
  "TSSG-44C,PRJEB39275,SAMEA7050837,Tigray,Sahareti Samire,Gijet,M,6,1.4",
  "TSSG-47H,PRJEB39275,SAMEA7050838,Tigray,Sahareti Samire,Gijet,F,12,1.4",
  "TSSG-48H,PRJEB39275,SAMEA7050839,Tigray,Sahareti Samire,Gijet,F,,",
  "TSSG-49C,PRJEB39275,SAMEA7050840,Tigray,Sahareti Samire,Gijet,M,5,1.7",
  "TSSM-17C,PRJEB39275,SAMEA7050841,Tigray,Sahareti Samire,Metkilimat,M,6,1.2",
  "TSSM-21H,PRJEB39275,SAMEA7050842,Tigray,Sahareti Samire,Metkilimat,F,6,0.9",
  "TSSM-35H,PRJEB39275,SAMEA7050843,Tigray,Sahareti Samire,Metkilimat,F,7,1.6",
  "TSSM-52C,PRJEB39275,SAMEA7050844,Tigray,Sahareti Samire,Metkilimat,M,8,1.5",
  "TSSM-62H,PRJEB39275,SAMEA7050845,Tigray,Sahareti Samire,Metkilimat,F,6,1.3",
  "TSSM-64H,PRJEB39275,SAMEA7050846,Tigray,Sahareti Samire,Metkilimat,F,12,1.2",
  "TSSM-68C,PRJEB39275,SAMEA7050847,Tigray,Sahareti Samire,Metkilimat,M,6,1.7",
  "TSSM-81H,PRJEB39275,SAMEA7050848,Tigray,Sahareti Samire,Metkilimat,F,8,1.6",
  "TSSM-92H,PRJEB39275,SAMEA7050849,Tigray,Sahareti Samire,Metkilimat,F,10,1.2",
  "TSSM-99C,PRJEB39275,SAMEA7050850,Tigray,Sahareti Samire,Metkilimat,M,8,1.3"
), collapse = "\n")
meta <- read.csv(text = metadatos_texto,
                 na.strings = c("", "NA"), check.names = FALSE)

meta$Taxa <- meta$ID_articulo
# Una sola discrepancia de escritura entre la tabla y el VCF proporcionado.
meta$Taxa[meta$ID_articulo == "TMLM-_C94"] <- "TMLM_-C94"
meta$Sexo <- toupper(trimws(meta$Sexo))
meta$Edad_meses <- as.numeric(meta$Edad_meses)
meta$Peso_kg <- as.numeric(meta$Peso_kg)
meta$Sexo_M <- ifelse(meta$Sexo == "M", 1, ifelse(meta$Sexo == "F", 0, NA_real_))
meta$Poblacion <- paste(meta$Region, meta$Distrito, meta$Localidad, sep = " / ")
if (anyDuplicated(meta$Taxa)) stop("Hay identificadores repetidos en los metadatos.")

leer_ids_vcf <- function(archivo) {
  con <- gzfile(archivo, open = "rt")
  on.exit(close(con))
  repeat {
    lineas <- readLines(con, n = 500L, warn = FALSE)
    if (!length(lineas)) stop("El VCF no contiene la cabecera #CHROM.")
    cab <- lineas[startsWith(lineas, "#CHROM\t")]
    if (length(cab)) return(strsplit(cab[1], "\t", fixed = TRUE)[[1]][-(1:9)])
  }
}
ids_vcf <- leer_ids_vcf(vcf_original)
if (anyDuplicated(ids_vcf)) stop("El VCF tiene identificadores duplicados.")
if (anyNA(match(ids_vcf, meta$Taxa))) {
  stop("Muestras sin metadatos: ", paste(setdiff(ids_vcf, meta$Taxa), collapse = ", "))
}
meta <- meta[match(ids_vcf, meta$Taxa), , drop = FALSE]
rownames(meta) <- NULL
stopifnot(identical(meta$Taxa, ids_vcf))
casos_completos <- complete.cases(meta[, c("Peso_kg", "Edad_meses", "Sexo_M")])
casos_completos <- casos_completos & meta$Peso_kg > 0 & meta$Edad_meses > 0
indices_gwas <- which(casos_completos)
if (length(indices_gwas) < 30L) stop("Muy pocos individuos con fenotipo y covariables completos.")
meta$Usado_en_GWAS <- casos_completos
meta$Motivo_exclusion <- vapply(seq_len(nrow(meta)), function(i) {
  faltantes <- c(if (is.na(meta$Peso_kg[i])) "peso faltante",
                if (is.na(meta$Edad_meses[i])) "edad faltante",
                if (is.na(meta$Sexo_M[i])) "sexo faltante")
  if (casos_completos[i]) "" else if (length(faltantes)) paste(faltantes, collapse = "; ") else "valor no valido"
}, character(1))
guardar_csv(meta, "01_metadatos_alineados.csv")
guardar_csv(meta[!casos_completos, ], "02_muestras_excluidas.csv")
guardar_csv(data.frame(ID_articulo = "TMLM-_C94", ID_VCF = "TMLM_-C94",
                      Motivo = "Diferencia de posicion del guion y guion bajo"),
            "03_correccion_identificador.csv")
cat("Muestras del VCF:", length(ids_vcf),
    "\nPeso registrado:", sum(!is.na(meta$Peso_kg)),
    "\nCasos completos para GWAS:", length(indices_gwas),
    "\nExcluidos por datos incompletos:", sum(!casos_completos), "\n")

# ==== Preparacion del VCF y control de calidad ====
# Tabla 2 del articulo: relacion exacta entre accesiones GRCg6a y cromosomas.
# En particular CM000124.5 es chr25; no se asigna por orden de aparicion.
diccionario <- data.frame(
  Accesion = c(sprintf("CM%06d.5", 93:116), "CM000124.5",
               sprintf("CM%06d.5", 117:119), "CM003637.2", "CM003638.2",
               "CM000120.4", "CM000123.5", "CM000121.5", "CM000122.5"),
  Cromosoma = c(as.character(1:28), "30", "31", "32", "33", "W", "Z")
)
guardar_csv(diccionario, "04_equivalencia_cromosomas.csv")
autosomas <- diccionario[!diccionario$Cromosoma %in% c("W", "Z"), ]
vcf_limpio <- file.path(procesados, "gallinas_autosomas_qc.vcf")
cache_qc <- file.path(procesados, "control_calidad.rds")
# El hash impide reutilizar una conversion de un VCF diferente.
firma <- list(vcf_md5 = unname(tools::md5sum(vcf_original)),
              muestras = ids_vcf, indices_gwas = indices_gwas,
              maf = CONFIG$maf_min, faltantes = CONFIG$faltantes_max,
              diccionario = diccionario, version_script = "1.0")

preparar_vcf <- function(entrada, destino) {
  con_in <- gzfile(entrada, "rt")
  on.exit(close(con_in), add = TRUE)
  temporal <- paste0(destino, ".tmp")
  con_out <- file(temporal, "wt")
  on.exit(try(close(con_out), silent = TRUE), add = TRUE)
  writeLines(c("##fileformat=VCFv4.2",
               "##source=GWAS_gallinas_Gheyas2022_autosomas_QC",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids_vcf), collapse = "\t")), con_out)
  contador <- c(Total_VCF = 0, Autosomas = 0, SNP_bialelicos_PASS = 0,
                Faltantes_aceptables = 0, MAF_aceptable = 0)
  qc_bloques <- list()
  faltantes_muestra <- numeric(length(ids_vcf))
  total_muestra <- 0
  codigos_gt <- c("0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1")
  dosis_gt <- c(0L, 1L, 1L, 2L, 0L, 1L, 1L, 2L)
  bloque <- 0L
  repeat {
    lineas <- readLines(con_in, n = CONFIG$bloque_vcf, warn = FALSE)
    if (!length(lineas)) break
    lineas <- lineas[!startsWith(lineas, "#")]
    if (!length(lineas)) next
    campos <- strsplit(lineas, "\t", fixed = TRUE)
    if (any(lengths(campos) != 9L + length(ids_vcf))) {
      stop("Hay una fila del VCF con una cantidad incorrecta de columnas.")
    }
    x <- matrix(unlist(campos, use.names = FALSE), nrow = length(campos), byrow = TRUE)
    contador["Total_VCF"] <- contador["Total_VCF"] + nrow(x)
    chr <- autosomas$Cromosoma[match(x[, 1], autosomas$Accesion)]
    contador["Autosomas"] <- contador["Autosomas"] + sum(!is.na(chr))
    mantener <- !is.na(chr) & x[, 7] %in% c("PASS", ".") &
      x[, 4] %in% c("A", "C", "G", "T") & x[, 5] %in% c("A", "C", "G", "T")
    x <- x[mantener, , drop = FALSE]
    chr <- chr[mantener]
    if (!nrow(x)) next
    contador["SNP_bialelicos_PASS"] <- contador["SNP_bialelicos_PASS"] + nrow(x)
    if (any(x[, 9] != "GT")) {
      stop("Este lector espera FORMAT=GT, como en chickens_small.vcf.gz. ",
           "Un VCF con otros campos requiere extraer primero GT.")
    }
    gt <- x[, -(1:9), drop = FALSE]
    dosis <- matrix(dosis_gt[match(gt, codigos_gt)], nrow = nrow(gt), ncol = ncol(gt))
    faltantes_muestra <- faltantes_muestra + colSums(is.na(dosis))
    total_muestra <- total_muestra + nrow(dosis)
    dosis_analisis <- dosis[, indices_gwas, drop = FALSE]
    ausentes <- rowMeans(is.na(dosis_analisis))
    n_llamados <- rowSums(!is.na(dosis_analisis))
    af <- rowSums(dosis_analisis, na.rm = TRUE) / (2 * pmax(1, n_llamados))
    maf <- pmin(af, 1 - af)
    pasa_faltantes <- ausentes <= CONFIG$faltantes_max & n_llamados > 0L
    contador["Faltantes_aceptables"] <- contador["Faltantes_aceptables"] + sum(pasa_faltantes)
    mantener <- pasa_faltantes & maf >= CONFIG$maf_min & maf > 0
    contador["MAF_aceptable"] <- contador["MAF_aceptable"] + sum(mantener)
    if (any(mantener)) {
      x <- x[mantener, , drop = FALSE]
      chr <- chr[mantener]
      gt <- gt[mantener, , drop = FALSE]
      # Identificadores unicos por posicion y alelos; los ID originales son '.'.
      snp_id <- paste(chr, x[, 2], x[, 4], x[, 5], sep = ":")
      gt[is.na(dosis[mantener, , drop = FALSE])] <- "./."
      gt <- gsub("|", "/", gt, fixed = TRUE)
      dim(gt) <- c(nrow(x), length(ids_vcf))
      nueva <- cbind(chr, x[, 2], snp_id, x[, 4], x[, 5], ".", "PASS", ".", "GT", gt)
      writeLines(apply(nueva, 1L, paste, collapse = "\t"), con_out)
      bloque <- bloque + 1L
      qc_bloques[[bloque]] <- data.frame(SNP = snp_id, Chromosome = as.integer(chr),
         Position = as.numeric(x[, 2]), REF = x[, 4], ALT = x[, 5],
         MAF_GWAS = maf[mantener], Fraccion_faltante_GWAS = ausentes[mantener])
    }
    if (contador["Total_VCF"] %% 50000L < CONFIG$bloque_vcf) {
      cat("Leidos:", contador["Total_VCF"], "; retenidos:", contador["MAF_aceptable"], "SNP\n")
    }
  }
  close(con_out)
  if (!length(qc_bloques)) stop("Ningun SNP supera los filtros.")
  qc_snps <- do.call(rbind, qc_bloques)
  if (anyDuplicated(qc_snps$SNP)) stop("Hay SNP repetidos por posicion y alelos.")
  if (!file.copy(temporal, destino, overwrite = TRUE)) stop("No se pudo guardar el VCF filtrado.")
  unlink(temporal)
  list(firma = firma, conteos = contador, snps = qc_snps,
       muestras = data.frame(Taxa = ids_vcf, Fraccion_faltante = faltantes_muestra / total_muestra))
}
qc <- NULL
if (CONFIG$reutilizar && file.exists(cache_qc) && file.exists(vcf_limpio)) {
  anterior <- readRDS(cache_qc)
  if (identical(anterior$firma, firma)) qc <- anterior
}
if (is.null(qc)) {
  qc <- preparar_vcf(vcf_original, vcf_limpio)
  saveRDS(qc, cache_qc)
} else cat("Se reutiliza el VCF filtrado con la misma firma.\n")
guardar_csv(data.frame(Etapa = names(qc$conteos), Numero_SNP = as.numeric(qc$conteos)),
            "05_control_calidad_conteos.csv")
guardar_csv(qc$muestras, "06_faltantes_por_muestra.csv")
guardar_csv(qc$snps, "07_snp_retenidos.csv")
print(data.frame(Etapa = names(qc$conteos), SNP = as.numeric(qc$conteos)))

# ==== Importacion a rMVP ====
prefijo <- file.path(procesados, "gallinas")
fenotipo_archivo <- file.path(procesados, "fenotipo_peso.csv")
write.csv(meta[, c("Taxa", "Peso_kg")], fenotipo_archivo, row.names = FALSE, na = "NA")
marca_conversion <- file.path(procesados, "firma_conversion.rds")
necesarios <- paste0(prefijo, c(".geno.desc", ".geno.bin", ".geno.map", ".geno.ind", ".phe"))
convertir <- TRUE
if (CONFIG$reutilizar && all(file.exists(necesarios)) && file.exists(marca_conversion)) {
  convertir <- !identical(readRDS(marca_conversion), firma)
}
if (convertir) {
  MVP.Data(fileVCF = vcf_limpio, filePhe = fenotipo_archivo, sep.phe = ",",
           SNP.impute = "Major", fileKin = FALSE, filePC = FALSE,
           priority = "memory", ncpus = cores, out = prefijo)
  saveRDS(firma, marca_conversion)
} else cat("Se reutiliza la conversion rMVP.\n")
genotype <- bigmemory::attach.big.matrix(paste0(prefijo, ".geno.desc"))
mapping <- read.table(paste0(prefijo, ".geno.map"), header = TRUE, check.names = FALSE)
mapping <- mapping[, 1:3, drop = FALSE]
names(mapping) <- c("SNP", "Chromosome", "Position")
mapping$Chromosome <- as.integer(mapping$Chromosome)
mapping$Position <- as.numeric(mapping$Position)
ids_rmvp <- read.table(paste0(prefijo, ".geno.ind"), header = FALSE, colClasses = "character")[[1]]
if (!identical(ids_rmvp, ids_vcf)) stop("La conversion altero el orden de los individuos.")
if (!identical(as.character(mapping$SNP), as.character(qc$snps$SNP))) {
  stop("El mapa y el control de calidad no estan alineados.")
}
marcadores_en_columnas <- nrow(genotype) == length(ids_vcf) && ncol(genotype) == nrow(mapping)
if (!marcadores_en_columnas && !(ncol(genotype) == length(ids_vcf) && nrow(genotype) == nrow(mapping))) {
  stop("Las dimensiones de genotipos, mapa y metadatos no coinciden.")
}
cat("Genotipos:", paste(dim(genotype), collapse = " x "),
    "\nSNP analizados:", nrow(mapping), "\n")

# ==== Kinship y PCA genotipico ====
# VanRaden sobre SNP autosomicos filtrados e imputados; todas las 234 muestras.
marca_kin <- file.path(procesados, "firma_kinship.rds")
kin_archivos <- paste0(prefijo, c(".kin.desc", ".kin.bin"))
calcular_kin <- TRUE
if (CONFIG$reutilizar && all(file.exists(kin_archivos)) && file.exists(marca_kin)) {
  calcular_kin <- !identical(readRDS(marca_kin), firma)
}
if (calcular_kin) {
  llamar_rmvp(MVP.Data.Kin, list(fileKin = TRUE, mvp_prefix = prefijo, out = prefijo,
                               mrk_bycol = marcadores_en_columnas, maxLine = 5000L))
  saveRDS(firma, marca_kin)
}
kinship_big <- bigmemory::attach.big.matrix(paste0(prefijo, ".kin.desc"))
kinship <- kinship_big[, ]
if (any(!is.finite(kinship))) stop("La matriz de parentesco tiene valores no finitos.")
kinship <- (kinship + t(kinship)) / 2
eigen_kin <- eigen(kinship, symmetric = TRUE)
varianza <- pmax(eigen_kin$values, 0)
varianza <- 100 * varianza / sum(varianza)
pcs <- eigen_kin$vectors[, seq_len(CONFIG$n_pc), drop = FALSE]
colnames(pcs) <- paste0("PC", seq_len(CONFIG$n_pc))
pca <- cbind(meta[, c("Taxa", "Region", "Poblacion", "Usado_en_GWAS")], as.data.frame(pcs))
guardar_csv(pca, "08_pca_coordenadas.csv")
guardar_csv(data.frame(PC = seq_along(varianza), Porcentaje_autovalor = varianza),
            "09_pca_varianza.csv")
write.csv(kinship, file.path(tablas, "10_kinship.csv"), row.names = FALSE)

dibujar_pca <- function() {
  regiones <- sort(unique(pca$Region))
  paleta <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")[seq_along(regiones)], regiones)
  par(mar = c(4.5, 4.5, 3.5, 1), las = 1)
  plot(pca$PC1, pca$PC2, pch = 19, cex = 1.1, col = paleta[pca$Region],
       xlab = sprintf("PC1 (%.2f%%)", varianza[1]),
       ylab = sprintf("PC2 (%.2f%%)", varianza[2]),
       main = paste("PCA genotipico -", nrow(pca), "gallinas"))
  legend("topright", legend = regiones, col = paleta, pch = 19, bty = "n", cex = 0.8)
}
grDevices::jpeg(file.path(graficos, "01_pca", "PCA_regiones.jpg"), width = 2400, height = 1800, res = 300, quality = 100)
dibujar_pca()
invisible(dev.off())
# Equivalente al facetado por subpoblacion del ejemplo de clase.
grDevices::jpeg(file.path(graficos, "01_pca", "PCA_poblaciones.jpg"), width = 3600, height = 3000, res = 300, quality = 100)
par(mfrow = c(4, 6), mar = c(2.5, 2.5, 2, 0.5), mgp = c(1.3, 0.4, 0), las = 1)
for (pop in unique(pca$Poblacion)) {
  seleccion <- pca$Poblacion == pop
  plot(pca$PC1, pca$PC2, pch = 19, col = "grey90", cex = 0.5,
       xlab = "PC1", ylab = "PC2", main = meta$Localidad[match(pop, meta$Poblacion)], cex.main = 0.75)
  points(pca$PC1[seleccion], pca$PC2[seleccion], pch = 19, col = "#0072B2", cex = 0.7)
}
invisible(dev.off())

# ==== Datos de entrada para GWAS ====
phenotype <- read.table(paste0(prefijo, ".phe"), header = TRUE, check.names = FALSE)
if (!identical(as.character(phenotype[[1]]), meta$Taxa)) stop("Fenotipos desordenados.")
trait <- data.frame(Taxa = meta$Taxa[indices_gwas], Peso_kg = meta$Peso_kg[indices_gwas])
covariables <- cbind(Edad_meses = meta$Edad_meses[indices_gwas],
                    Sexo_M = meta$Sexo_M[indices_gwas], pcs[indices_gwas, , drop = FALSE])
covariables <- scale(covariables, center = TRUE, scale = TRUE)
if (qr(cbind(Intercepto = 1, covariables))$rank < 1L + ncol(covariables)) {
  stop("Las covariables presentan dependencia lineal.")
}
K_gwas <- kinship[indices_gwas, indices_gwas, drop = FALSE]
# Copia a disco por bloques: las mismas muestras completas en cada matriz.
prefijo_gwas <- file.path(procesados, "gallinas_casos_completos")
marca_gwas <- paste0(prefijo_gwas, ".firma.rds")
archivos_gwas <- paste0(prefijo_gwas, c(".geno.desc", ".geno.bin"))
reusar_gwas <- CONFIG$reutilizar && all(file.exists(archivos_gwas)) &&
  file.exists(marca_gwas) && identical(readRDS(marca_gwas), firma)
if (reusar_gwas) {
  genotype_gwas <- bigmemory::attach.big.matrix(archivos_gwas[1])
} else {
  if (exists("genotype_gwas", inherits = FALSE)) { rm(genotype_gwas); invisible(gc()) }
  unlink(archivos_gwas)
  genotype_gwas <- bigmemory::filebacked.big.matrix(
    nrow = length(indices_gwas), ncol = nrow(mapping), type = "char",
    backingfile = basename(archivos_gwas[2]), descriptorfile = basename(archivos_gwas[1]),
    backingpath = procesados)
  for (inicio in seq.int(1L, nrow(mapping), by = 5000L)) {
    ind <- inicio:min(inicio + 4999L, nrow(mapping))
    if (marcadores_en_columnas) {
      genotype_gwas[, ind] <- genotype[indices_gwas, ind, drop = FALSE]
    } else {
      genotype_gwas[, ind] <- t(genotype[ind, indices_gwas, drop = FALSE])
    }
  }
  bigmemory::flush(genotype_gwas)
  saveRDS(firma, marca_gwas)
}
stopifnot(nrow(trait) == nrow(genotype_gwas), nrow(K_gwas) == nrow(trait),
          nrow(covariables) == nrow(trait), ncol(genotype_gwas) == nrow(mapping))
guardar_csv(trait, "11_fenotipo_usado.csv")
guardar_csv(cbind(Taxa = trait$Taxa, as.data.frame(covariables)), "12_covariables_usadas.csv")

# ==== GWAS: GLM, MLM y FarmCPU ====
# Peso en kg, ajustado por edad, sexo y 3 PCs. MLM ademas incorpora kinship.
# Los PCs ya estan en CV; nPC=0 evita agregarlos dos veces.
# file.type en MVP es un formato de imagen, no "csv".
# Las tablas se exportan explicitamente mas abajo.
imMVP <- MVP(phe = trait, geno = genotype_gwas, map = mapping, K = K_gwas,
             CV.GLM = covariables, CV.MLM = covariables, CV.FarmCPU = covariables,
             nPC.GLM = 0, nPC.MLM = 0, nPC.FarmCPU = 0,
             method = c("GLM", "MLM", "FarmCPU"), ncpus = cores,
             maxLoop = CONFIG$max_loop_farmcpu,
             file.output = FALSE, outpath = procesados)
saveRDS(imMVP, file.path(procesados, "resultados_GLM_MLM_FarmCPU.rds"))
invisible(gc())

modelos <- c(GLM = "glm.results", MLM = "mlm.results", FarmCPU = "farmcpu.results")
map_gwas <- as.data.frame(imMVP$map[, 1:3, drop = FALSE])
names(map_gwas) <- c("SNP", "Chromosome", "Position")
if (nrow(map_gwas) != nrow(mapping)) stop("El GWAS devolvio una cantidad inesperada de SNP.")
Pmap <- map_gwas
tabla_resumen <- list()
for (modelo in names(modelos)) {
  r <- imMVP[[modelos[[modelo]]]]
  if (is.null(r) || nrow(r) != nrow(map_gwas)) stop("El modelo ", modelo, " no devolvio resultados completos.")
  # En diferentes versiones el p-valor es la ultima columna (2 o 3 columnas).
  p <- as.numeric(r[, ncol(r)])
  p[!is.finite(p) | p < 0 | p > 1] <- NA_real_
  p[p == 0] <- .Machine$double.xmin
  if (sum(!is.na(p)) < 10L) stop("Muy pocos p-valores validos en ", modelo)
  Pmap[[paste0("Peso_kg.", modelo)]] <- p
  umbral_bonf <- CONFIG$alpha / nrow(map_gwas)
  tabla <- cbind(map_gwas, as.data.frame(r))
  tabla$P_valor <- p
  tabla$P_ajustada_Bonferroni <- p.adjust(p, method = "bonferroni", n = nrow(map_gwas))
  tabla$Significativo_Bonferroni <- !is.na(p) & p <= umbral_bonf
  tabla$Modelo <- modelo
  guardar_csv(tabla, paste0("13_GWAS_", modelo, ".csv"))
  top <- tabla[order(tabla$P_valor, na.last = NA), , drop = FALSE]
  guardar_csv(head(top, 20L), paste0("14_top20_", modelo, ".csv"))
  guardar_csv(tabla[tabla$Significativo_Bonferroni, , drop = FALSE],
              paste0("15_significativos_", modelo, ".csv"))
  p_validos <- p[!is.na(p)]
  lambda_gc <- median(qchisq(1 - p_validos, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  tabla_resumen[[modelo]] <- data.frame(Modelo = modelo, Muestras = nrow(trait),
      SNP_evaluados = nrow(map_gwas), SNP_p_valido = length(p_validos),
      P_min = min(p_validos), Lambda_GC = lambda_gc,
      Umbral_Bonferroni = umbral_bonf, SNP_significativos = sum(tabla$Significativo_Bonferroni))
}
resumen <- do.call(rbind, tabla_resumen)
rownames(resumen) <- NULL
guardar_csv(resumen, "16_resumen_modelos.csv")
guardar_csv(Pmap, "17_Pmap_tres_modelos.csv")

# ==== Grafico de densidad de SNP ====
# Densidad del subconjunto retenido, no los aproximadamente 15 M del articulo.
MVP.Report(MVP = map_gwas, plot.type = "d", bin.size = 1e6,
           col = c("darkgreen", "yellow", "red"),
           file.type = "jpg", dpi = 300, memo = "Gallinas_autosomas_QC",
           outpath = file.path(graficos, "02_densidad"))

# ==== Distribucion del fenotipo ====
MVP.Hist(phe = trait, breakNum = 18, file.type = "jpg", dpi = 300,
         memo = "Peso_corporal_kg_casos_completos",
         outpath = file.path(graficos, "03_fenotipo"))

# ==== Manhattan con los umbrales del docente ====
# Estos umbrales se conservan para comparar con la practica (1e-6 y 1e-4).
# La significancia corregida se decide con Bonferroni, en el grafico siguiente.
MVP.Report(MVP = Pmap, plot.type = "m", multracks = TRUE,
           threshold = c(1e-6, 1e-4), threshold.lty = c(1, 2),
           threshold.lwd = c(1, 1), threshold.col = c("black", "grey"),
           amplify = TRUE, bin.size = 1e6,
           chr.den.col = c("darkgreen", "yellow", "red"),
           signal.col = c("red", "green"), signal.cex = c(1, 1),
           file.type = "jpg", dpi = 300, memo = "Umbrales_clase",
           outpath = file.path(graficos, "04_manhattan_clase"))

# ==== Manhattan con correccion por pruebas multiples ====
umbral_bonf <- CONFIG$alpha / nrow(map_gwas)
MVP.Report(MVP = Pmap, plot.type = "m", multracks = TRUE,
           threshold = c(umbral_bonf, 1e-4), threshold.lty = c(1, 2),
           threshold.lwd = c(1, 1), threshold.col = c("black", "grey50"),
           amplify = TRUE, bin.size = 1e6,
           chr.den.col = c("darkgreen", "yellow", "red"),
           signal.col = c("red", "green"), signal.cex = c(1, 1),
           file.type = "jpg", dpi = 300, memo = "Bonferroni",
           outpath = file.path(graficos, "05_manhattan_bonferroni"))

# ==== Graficos Q-Q ====
# Se corrige threshold=1e6 del ejemplo: un p-valor debe estar entre 0 y 1.
MVP.Report(MVP = Pmap, plot.type = "q", multracks = TRUE,
           col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),
           threshold = umbral_bonf, signal.pch = 19, signal.cex = 1.5,
           signal.col = "red", box = FALSE, conf.int = TRUE,
           file.type = "jpg", dpi = 300, memo = "GLM_MLM_FarmCPU",
           outpath = file.path(graficos, "06_qq"))

# ==== Verificacion y registro final ====
archivos_graficos <- list.files(graficos, pattern = "\\.jpg$", full.names = TRUE, recursive = TRUE)
for (carpeta in carpetas[startsWith(carpetas, "graficos/")]) {
  if (!length(list.files(file.path(salida, carpeta), pattern = "\\.jpg$"))) {
    stop("Falta un grafico en ", carpeta, ". Revise los mensajes de esa seccion.")
  }
}
guardar_csv(data.frame(Archivo = basename(archivos_graficos),
                      Carpeta = basename(dirname(archivos_graficos))), "18_indice_graficos.csv")
writeLines(capture.output(sessionInfo()), file.path(salida, "sessionInfo.txt"))
writeLines(c(
  "GWAS exploratorio de peso corporal en gallinas indigenas de Etiopia.",
  "Fuente de genotipos: chickens_small.vcf.gz proporcionado para la tarea.",
  "Fuente de fenotipos: Gheyas et al. (2022), Online-only Table 1.",
  "DOI: https://doi.org/10.1038/s41597-022-01129-4",
  sprintf("VCF: %d muestras; %d variantes originales.", length(ids_vcf), qc$conteos["Total_VCF"]),
  sprintf("GWAS: %d individuos completos; %d SNP autosomicos filtrados.", nrow(trait), nrow(map_gwas)),
  "Fenotipo: peso corporal (kg). Covariables: edad (meses), sexo y 3 PCs.",
  "No se imputaron pesos ni edades. Se imputaron genotipos faltantes con Major (rMVP).",
  "Filtro SNP: bialelicos, PASS o sin filtro, MAF >=0.05 y faltantes <=0.20 en casos completos.",
  "Se excluyeron cromosomas W/Z y secuencias sin ubicacion autosomica.",
  "Los PCs proceden de la matriz VanRaden de los SNP filtrados; no se realizo poda por LD.",
  "Las tablas conservan la diferencia de ID TMLM-_C94 -> TMLM_-C94.",
  sprintf("Umbral de Bonferroni: 0.05/%d = %.8g.", nrow(map_gwas), umbral_bonf),
  "Los umbrales 1e-6 y 1e-4 del docente tambien se dibujan, como referencia de clase.",
  "El articulo describe datos genomicos; este GWAS de peso es una aplicacion de la tarea, no un resultado publicado alli.",
  "El VCF es un subconjunto. Los resultados no representan el conjunto de aproximadamente 15 M SNP del articulo.",
  "La muestra es pequena y el peso depende del manejo/ambiente. Una asociacion no demuestra causalidad."
), file.path(salida, "LEEME_resultados.txt"), useBytes = TRUE)
print(resumen)
cat("\nProceso finalizado.\nResultados:", normalizePath(salida, winslash = "/"),
    "\nGraficos JPG:", length(archivos_graficos),
    "\nAbra la carpeta graficos para ver PCA, densidad, histograma, Manhattan y Q-Q.\n")
