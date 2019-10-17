# Acceso-areas-verdes
accesibilidad áreas verdes 
# Script para optimizar el indicador de accesibilidad a Ã¡reas verdes
library(sf); library(tidyverse); library(glue); library(RCurl); library(spdep)
# Se definen carpetas de trabajo
rootdir <- "D:/ownCloud/Indicadores_de_Sustentabilidad_Expansion/Areas_Verdes/"
workdir <- glue("{rootdir}input/")
rawdir <- glue("{workdir}aavv/Pz_Pq_SIEDU_28062019.shp")
out_dir <- glue("{rootdir}Output/")
#Agregrar codigo de ciudades en estudio
codes <- read.csv(glue("{workdir}Ciudades_expansion.csv")) %>%
  mutate(Comuna = toupper(Comuna)) %>%
  rename(COMUNA = Comuna)
# Se lee shapefile y se crea un dataframe para almacenar walksheds 
shp <- read_sf(rawdir) %>%
  st_transform(5361)
#Filtrar areas verdes de comuna en estudio
shp <- shp[shp$CUT %in% codes$Codigo,]

#### Identificar Ã¡reas verdes compuestas ####

# Crear una red con las Ã¡reas verdes que estÃ©n a 12 metros por lo bajo.
# ElegÃ­ 12 porque sentÃ­ que fue el valor que dio el mejor resultado
nbs <- poly2nb(shp, snap = 12) %>%
  nb2listw(zero.policy = TRUE)

# La idea, dentro de todo, es generar una lista en que cada elemento sea un
# vector con los indices de las Ã¡reas verdes que lo componen. 
# Para eso tenemos que modificar un poco la red, primero, aÃ±adiendo indice
# del poligono a la lista de vecinos de dicho poligono. Segundo, ordenandolo
# De mayor a menor
lst <- nbs$neighbours

for (i in seq(lst)) {
  if (0 %in% lst[[i]]) {next()}
  lst[[i]] <- append(lst[[i]], i)
  lst[[i]] <- sort(lst[[i]])
}
# Luego, hacemos dos iteraciones de un sapply que combina las listas segÃºn
# los elementos que tienen en comÃºn
for (i in 1:2){
  lst <- unique(sapply(lst, function(x) 
    unique(unlist(lst[sapply(lst, function(y) 
      any(x %in% y))]))))
}

# Luego de esto creamos un nuevo vector con los indices de las Ã¡reas verdes
# compuestas que superen los 5000 metros 2 de superficie.
composite_big <- c()
for (i in 2:length(lst)){
  subset <- shp[lst[[i]],] %>%
    summarise(area = as.numeric(sum(st_area(.))))  %>%
    `st_geometry<-`(NULL) 
  if (subset$area >= 5000) {
    composite_big <- c(composite_big, lst[[i]])
  }
}

coords <- data.frame()
aavv_filtered <- st_sf(st_sfc(), crs = 5361)
# Loop que filtra Ã¡reas verdes segÃºn sus requisitos mÃ­nimos (10 metros ancho y 
# 0.5 ha superficie), y luego las secciona 
# con una grilla de 1ha. Las secciones que son igual a 1ha son eliminadas 
# (pues se considera que estÃ¡n dentro del Ã¡rea verde) Y sÃ³lo se dejan las 
# menores a 1ha y mayores a 0.05ha, que serÃ­an los bordes 
# del Ã¡rea verde. Luego se obtienen los centroides de estas secciones.
for(i in 1:nrow(shp)){
  aavv <- shp[i,] %>%
    mutate(Area = as.numeric(st_area(.)))
  # Evaluar si el polÃ­gono tiene el tamaÃ±o minimo.
  bbox <- st_bbox(aavv)
  ancholargo <- c(bbox[3] - bbox[1], bbox[4] - bbox[2])
  ancho <- min(ancholargo)
  largo <- max(ancholargo)
  # Si no tiene el tamaÃ±o mÃ­nimo, se pasa al siguiente
  if ((aavv$Area <= 5000| !(ancho > 10 & aavv$Area > largo*10)) & 
      !i %in%  composite_big) { print(i);next()}
  aavv_filtered <- rbind(aavv_filtered, aavv)
  # Generar grilla que segmenta las Ã¡reas verdes. Los pedazos de una hectarea
  # Son centros de un Ã¡rea verde. Los segmentos menores a 1ha son periferia del
  # area verde. Se filtran pedazos muy pequeÃ±os de 500 m2
  grilla <- st_make_grid(aavv, cellsize = c(100,100))
  aavv_sliced <- st_intersection(aavv, grilla) %>%
    select(Area) %>%
    mutate(slice_area = as.numeric(st_area(.))) %>%
    filter(slice_area < 10000 & slice_area > 500) %>%
    st_transform(4326) %>%
    st_centroid(.) %>%
    mutate(X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) %>%
    `st_geometry<-`(NULL) 
  coords <- rbind(coords, aavv_sliced)
}


# FunciÃ³n para obtener isocronas a partir de centroides
isocrona <- function(coordinates, mode, viaje.tiempo, viaje.velocidad = 1.138){
  consulta.url <- paste0('http://otpv2.cedeus.cl/otp/routers/enero2019/isochroneOld?fromPlace=', 
                         coordinates, 
                         '&walkTime=', viaje.tiempo, 
                         '&walkSpeed=', viaje.velocidad,
                         '&mode=', mode, 
                         "&time=", "2019-05-21T08:00:00",
                         '&toPlace=-33.5846161,-70.5410151&output=', 
                         "SHED") 
  
  print(consulta.url)
  walkshed.texto <- getURL(consulta.url)
  polygon <- read_sf(walkshed.texto) %>%
    st_transform(4326) %>%
    st_zm("ZM")
  return(polygon)
}
# Crear un sf vacÃ­o donde ir almacenando los datos
walksheds <- st_sf(st_sfc(), crs = 4326)
# Loop que calcula las isocronas. Para Ã¡reas de 2ha o mÃ¡s, 10 minutos caminando
# Para Ã¡reas menores a 2ha, 5 minutos caminando.
for (o in seq(nrow(coords))) {
  print(o)
  coords_subset <- coords[o,]
  coordenadas <- paste(coords_subset$Y, coords_subset$X, sep = ",")
  tiempo <- if_else(coords_subset$Area >= 20000, 10, 5)
  iso <- isocrona(coordenadas, mode = "WALK", viaje.tiempo = tiempo) %>%
    st_cast("POLYGON")
  if (st_geometry_type(iso) == "LINESTRING") {stop()}
  st_geometry(coords_subset) <- iso$geometry
  walksheds <- rbind(walksheds, coords_subset)
}
# Transformar a sirgas
walksheds <- st_transform(walksheds, 5361)

#### IntersecciÃ³n con manzanas
# Carpetas y shapefile
mzn_dir <- glue("{workdir}manzanas/")
mzn_shp_dir <- glue("{mzn_dir}Manzanas_expansion.shp")
# Leer shp. Orignalmente estÃ¡ en latlon, pasar a sirgas, calcular Ã¡rea
mzn_shp <- st_read(mzn_shp_dir) %>%
  st_transform(4326) %>%
  st_transform(5361) %>%
  mutate(area_mzn = as.numeric(st_area(.)))
# Intersectar isocronas con manzanas. Calcular el Ã¡rea de cobertura de cada 
# Isocrona relativo a su manzanas. Si una isocrona cubre menos del 75% de una 
# Manzana, esa manzana no se considera con acceso.


pop <-  st_intersection(mzn_shp, walksheds) %>% 
        mutate(area_mznsub = as.numeric(st_area(.)),
         completitud = round(area_mznsub/area_mzn*100)) %>%
  filter(completitud >= 75, TOTAL_P != 0) %>%
  group_by(MANZENT) %>%
  # Calcular la superficie total de Ã¡reas verdes a las cual una manzana tiene 
  # Acceso y dividirlo por el nÃºmero de habitantes de esa manzana.
  summarise(Habs = unique(TOTAL_P),
            area_av = sum(unique(Area))) %>%
  mutate(m2hab = area_av/Habs) %>%
  `st_geometry<-`(NULL) %>%
  # Unir con shp de manzanas total, y definir las manzanas que tienen acceso
  # Y las que no
  full_join(mzn_shp) %>%
  mutate(access = if_else(is.na(m2hab), 0, 1))

porcentaje <- pop %>%
  ungroup() %>%
  summarise(per = round(sum(TOTAL_PERS[access == 1])/sum(TOTAL_PERS)*100,2))
print(porcentaje)  
# Exportar resultados e isocronas
if_else(dir.exists(out_dir), dir.create(output_dir), NA)
write_sf(pop, glue("{out_dir}acceso_binario"), driver = "ESRI Shapefile")
write_sf(walksheds, glue("{out_dir}walksheds"), driver = "ESRI Shapefile")
write_sf(aavv_filtradas, glue("{out_dir}aavv_filtradas"), driver = "ESRI Shapefile")
