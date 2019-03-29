// Make heatmap for arko

process compile {
	conda "rust git"
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

raw_file = Channel.fromPath("${params.input_file}")
	.map { file -> tuple(file.baseName, file) }

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



process normalization_colorization {
	conda "r-base r-magrittr r-data.table"

	//publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_$filename" }

    label 'rscript'
    
    input:
    set datasetID, file(long_file) from to_normalize
    file tmatrix  from file("${baseDir}/data/T")
    file coldcurves from file("${baseDir}/data/blue_no_darkening.csv")
    file hotcurves from file("${baseDir}/data/fire_no_darkening.csv")

    output:
    set datasetID, file("data_output.csv") into to_heat


    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(magrittr)



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

	#background_color <- hex_to_rgb("#e5e5e5")
	background_color <- hex_to_rgb("#000000")
	#background_color <- hex_to_rgb("#ffffff")
	 
	 if(!is.na(proportion)){
	   
	   color <- colormix(fire[[color]], blue[[color_negative]], proportion) %>% rgb_mix(background_color, intensity) %>% into_string
	   
	 } else {
	   
	   color <- background_color %>% into_string
	 }
	 
	}


	head(melted)



	threshold <- $params.threshold

	max <- quantile(melted[,value], probs = c(threshold))

	min <- quantile(melted[,value_negative], probs = c(1 - threshold))
	#max <- min * -1

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

	fwrite(merged[order(-id,-col),.(id,col,hex_color)], "data_output.csv", col.names = FALSE)
    """

}



process render_to_png {
	conda "scipy numpy pillow"

	publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${datasetID}_t${params.threshold}_b${params.size}_$filename" }

    label 'pythonscript'
    
    input:
    set datasetID, file(long_file) from to_heat

    output:
    file "heatmap.png" into out


"""
#!/usr/bin/env python3


import numpy as np
import scipy.misc as smp

def render(x, y, color):
    data[x, y][0], data[x, y][1], data[x, y][2]= [int(x) for x in color.split("_")]


with open("$long_file") as f:
    created = False
    for line in f:
        record = line.strip().split(",")
        x, y, = [int(x) for x in record[:2]]
        color = record[2]
        #print(x,y,color)
        if created:
            render(x,y,color)
        else:
            data = np.full( (x+1,y+1,3), $params.background_color_u8, dtype=np.uint8)
            render(x,y,color)
            created = True
        
img = smp.toimage( data )       # Create a PIL image
smp.imsave('heatmap.png', img)  
"""

}



