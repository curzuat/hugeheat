// Make heatmap for arko
params.outdir="."
process compile {
	storeDir "$params.outdir/bin/" 

	//publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_reduced_$filename" }
    
    input:
    val repo from "https://github.com/CGUTA/csvmipmap.git"

    output:
    file 'csvmipmap/target/release/csvmipmap' into binary


    """
    git clone $repo
    cd csvmipmap/
    cargo build --release
    #cp target/release/csvmipmap ../csvmipmap_b
    """

}

//log1 = Channel.create().subscribe { println "Log 1: $it" }

raw_file = Channel.fromPath("$params.input_file")
	.map { file -> tuple(file.baseName, file) }

params.size = '1,1'
params.pixel = '1,1'
process box_sampling_into_long {

	//publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_reduced_$filename" }
    
    input:
    set datasetID, file(dataset_to_convert) from raw_file
    file csvmipmap from binary

    output:
    set datasetID, file("data.csv") into to_normalize


    """
    ./$csvmipmap $params.size $dataset_to_convert > data.csv
    """

}

params.default_gap_size = 'NOT_PROVIDED'
params.column_gap_size = 'NOT_PROVIDED'
params.row_gap_size = 'NOT_PROVIDED'
params.intensity = 'NOT_PROVIDED'
params.column_gaps = 'NO_FILE_COLUMNS'
params.row_gaps = 'NO_FILE_ROWS'
params.grid = 'NO_FILE_GRID'

params.truncate_positive_at = 'NOT_PROVIDED'
params.truncate_negative_at = 'NOT_PROVIDED'
params.truncate_extremes_at = 'NOT_PROVIDED'
params.automatic_threshold = 'NOT_PROVIDED'

params.bicolor = 'NOT_PROVIDED'


process normalization_colorization {
	//echo true

	//publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_$filename" }

    label 'rscript'
    
    input:
    set datasetID, file(long_file) from to_normalize
    file tmatrix  from file("${baseDir}/data/T")
    file coldcurves from file("${baseDir}/data/blue_no_darkening.csv")
    file hotcurves from file("${baseDir}/data/fire_no_darkening.csv")
    file column_gaps from file(params.column_gaps)
    file row_gaps from file(params.row_gaps)
    file grid from file(params.grid) 

    output:
    set datasetID, file("data_output.csv") into to_heat


    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(magrittr)
    library(glue)



    melted <- fread("$long_file")
    blue <- fread("$coldcurves")
	fire <- fread("$hotcurves")
    Tmat <- read.table("$tmatrix", header = F) %>% as.matrix

    gamma_correction <- function(n){
	  if (n < 0.0031308){
	    n * 12.92
	  } else if (n < 1){
	    1.055 * n^(1/2.4) - 0.055
	  } else {
	    1
	  }
	}

	scale_to_u8 <- function(n){
	  round(n * 255)
	}

	into_rgb <- function(curve){
	  Tmat %*% curve %>%
	    apply(.,1,gamma_correction) %>%
	    scale_to_u8
	}

	into_hex <- function(rgb){
	    sprintf("%x",rgb) %>% # hex encoding
	    paste0(., collapse="") %>% # formatting
	    sprintf("#%s",.)
	}

	into_string <- function(rgb){
	    paste0(rgb, collapse="_") # formatting
	}

	colormix <- function(curve1, curve2, curve1_proportion){
	  
	  curve_of_mix <- curve1^curve1_proportion * curve2^(1-curve1_proportion)
	  
	  curve_of_mix %>% 
	    into_rgb 
	}

	choose_color <- function(x, output) {
	 color <- x[1]
	 color_negative <- x[2]
	 proportion <- x[3]
	 intensity <- x[4]

	if("$params.intensity" == "false"){
		intensity = 1.0
	}
	 
	 sRGB_to_linear <- function(rgb){
	    out <- numeric(length = length(rgb))
	    for(i in seq_along(rgb)) {
	        if (rgb[i] < 0.04045){
	          out[i] <- rgb[i]/12.92
	        } else {
	          out[i] <- ((rgb[i] + 0.055)/1.055)^2.4
	        }
	    }
	    out
	 }
	 
	 linear_to_sRGB <- function(rgb){
	   out <- numeric(length = length(rgb))
	    for(i in seq_along(rgb)) {
	        if (rgb[i] < 0.0031308){
	          out[i] <- rgb[i] * 12.92
	        } else {
	          out[i] <- 1.055 * rgb[i]^(1/2.4) - 0.055
	        } 
	    }
	    round(out)
	}
	 
	 mask <- function(rgb, intensity){
	   round(rgb * intensity)
	 }

	 gamma_correct_mask <- function(rgb, intensity){
	   linear_to_sRGB(sRGB_to_linear(rgb) * intensity)
	 }

	gamma_correct_rgb_mix <- function(rgb1, rgb2, rgb1_proportion){
  
		linear_to_sRGB(sRGB_to_linear(rgb1) * rgb1_proportion + sRGB_to_linear(rgb2) * (1 - rgb1_proportion))
		  
	}

	rgb_mix <- function(rgb1, rgb2, rgb1_proportion){
		  
		round(rgb1 * rgb1_proportion + rgb2 * (1 - rgb1_proportion))
		  
	}

	hex_to_rgb <- function(hex){
    	col2rgb(hex)[,1]
	}

	background_color <- hex_to_rgb("$params.background_color")
	 
	 if(!is.na(proportion)){
	   
	   if("$params.bicolor" == "NOT_PROVIDED"){ 
			color <- colormix(fire[[color]], blue[[color_negative]], proportion) %>% rgb_mix(background_color, intensity) %>% into_string
	   } 
	   else {
	   		color_pair <- "$params.bicolor" %>% strsplit(",") %>% unlist
	   		color_positive <- color_pair[1] %>% hex_to_rgb
	   		color_negative <- ifelse(length(color_pair) == 2, color_pair[2] %>% hex_to_rgb, color_positive)

	   		color <- rgb_mix(color_positive, color_negative, proportion) %>% rgb_mix(background_color, intensity) %>% into_string
	   }
	   
	 } else {
	   
	   color <- background_color %>% into_string
	 }
	 
	}


	head(melted)

	default_truncation <- ifelse("$params.truncate_extremes_at" == "NOT_PROVIDED", 1.0, $params.truncate_extremes_at) # hardcoded default
	truncate_positive_at <- ifelse("$params.truncate_positive_at" == "NOT_PROVIDED", default_truncation, $params.truncate_positive_at)
	truncate_negative_at <- ifelse("$params.truncate_negative_at" == "NOT_PROVIDED", default_truncation, $params.truncate_negative_at)

	threshold <- ifelse("$params.automatic_threshold" == "NOT_PROVIDED", 1, $params.automatic_threshold) # hardcoded default

	max <- ifelse("$params.automatic_threshold" == "NOT_PROVIDED", truncate_positive_at, quantile(melted[,value], probs = c(threshold)))
	min <- ifelse("$params.automatic_threshold" == "NOT_PROVIDED", -truncate_negative_at, quantile(melted[,value_negative], probs = c(1 - threshold)))

	melted[value > max, value := max]
	melted[value_negative < min, value_negative := min]

	if (max > 0){
		melted[, value := value/max]
	}
	if (min < 0){
		melted[, value_negative := value_negative/min]
	}


	melted[, color := 9]
	melted[, color_negative := 9]

	melted[, intensity := (value + value_negative)/2]
	melted[value == 0, intensity := value_negative]
	melted[value_negative == 0, intensity := value]
	melted[, intensity := round(intensity, 2)]

	for (i in 8:1){
	  melted[ value < i/9, color := i]
	  melted[ value_negative < i/9, color_negative := i]
	}

	melted[, proprotion_of_positive := round(value/(value + abs(value_negative)), 2)]

	conversion_table <- melted[,.N,by=.(color, color_negative, proprotion_of_positive, intensity)]
	hex_colors <- apply(conversion_table, 1, choose_color)

	conversion_table[[5]] <- hex_colors
	setkey(conversion_table, color, color_negative, proprotion_of_positive, intensity)
	setkey(melted, color, color_negative, proprotion_of_positive, intensity)

	merged <- melted[conversion_table, ]
	setnames(merged, "N", "hex_color")

	### Adding gaps in heatmap #####

	calcualte_gap_mapping <- function(gap_dir, compression_size){ 
		# Depending on kernel size gaps must be adjusted

		gaps <- fread(gap_dir)
		gaps[,.(mapped_position = floor(position/compression_size))][, mapped_position]
	}

	carve_gap <- function(what, gap_size, gap_locations){


		i_operation <- "{what} >= i" %>% glue %>% parse(text = .)
		j_operation <- "{what} := {what} + {gap_size}" %>% glue %>% parse(text = .)

		for (i in sort(gap_locations, decreasing = T)){
			merged[ eval(i_operation), eval(j_operation)] #  pixel gap
		}
		
	}


	# default_gap_size will be used if no gap specified
	gap_size <- ifelse("$params.default_gap_size" == "NOT_PROVIDED", 5, $params.default_gap_size) # hardcoded default gap size
	vertical_gap_size <- ifelse("$params.column_gap_size" == "NOT_PROVIDED", gap_size, $params.column_gap_size)
	horizontal_gap_size <- ifelse("$params.row_gap_size" == "NOT_PROVIDED", gap_size, $params.row_gap_size)

	horizontal_compression <- strsplit("${params.size}", split=",")[[1]][1] %>% as.numeric()
	vertical_compression <- strsplit("${params.size}", split=",")[[1]][2] %>% as.numeric()


	if("$params.column_gaps" != "NO_FILE_COLUMNS") {
		carve_gap("col", vertical_gap_size, calcualte_gap_mapping("$column_gaps", horizontal_compression))
	}
	if("$params.row_gaps" != "NO_FILE_ROWS") {
		carve_gap("id", horizontal_gap_size, calcualte_gap_mapping("$row_gaps", vertical_compression))
	}
	if("$params.grid" != "NO_FILE_GRID") {
		carve_gap("id", horizontal_gap_size, calcualte_gap_mapping("$grid", vertical_compression))
		carve_gap("col", vertical_gap_size, calcualte_gap_mapping("$grid", horizontal_compression))
	}



	###

	### pixel to parallelogram in render


	exageration <-c($params.pixel)
	max_row <- merged[,max(id)]
	#row_exageration <- data.table(id = rep(0:max_row, each=exageration[1]), exagerated_row = 0:(max_row*exageration[1]))
	row_exageration <- data.table(id = rep(0:max_row, each=exageration[1]), exagerated_row = 0:((max_row+1)*exageration[1] -1))

	max_col <- merged[,max(col)]
	col_exageration <- data.table(col = rep(0:max_col, each=exageration[2]), exagerated_col = 0:((max_col+1)*exageration[2] -1))

	exagerated_merged <- merged[row_exageration, on="id", allow.cartesian=T][col_exageration, on="col", allow.cartesian=T][order(-exagerated_row,-exagerated_col),.(exagerated_row,exagerated_col,hex_color)][!is.na(hex_color)]

	##


	fwrite(exagerated_merged, "data_output.csv", col.names = FALSE)
	#fwrite(merged[order(-id,-col),.(id, col, hex_color)], "data_output.csv", col.names = FALSE)
    """

}

process compile_display {
	storeDir "$params.outdir/bin/" 

	//publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_reduced_$filename" }
    
    input:
    val repo from "https://github.com/CGUTA/display_heatmap.git"

    output:
    file 'display_heatmap/target/release/display_heatmap' into binary_display


    """
    git clone $repo
    cd display_heatmap/
    cargo build --release
    """

}



process render_to_png {

	publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_t${params.threshold}_b${params.size}_p${params.pixel}_$filename" }

    
    input:
    set datasetID, file(long_file) from to_heat
    file display from binary_display

    output:
    file "heatmap.png" into out


	"""
	./$display "$long_file"
	"""

}


