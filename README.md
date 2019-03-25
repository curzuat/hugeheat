A pipeline to make a visualization heatmap of large matrices with float values stored as csv

Procedure:

the pipeline reduces the size of the matrix by making an average of a reticule of --size (box sampling) treating negative and positive values separately


the values in the reduced matrix are scaled and truncated separately to 0-1 for negative and positive values where 1 is the --threshold percentile

a shade of blue is assigned to each cell coresponding to the  negative and positive values and a final color is created by linear interpolation

cell final colors are tinted with black if their absolute values were relatively low in the original matrix


finally the heatmap is plotted using --background_color_u8 as the color for cells with 0 positive and 0 negative values

Usage
nextflow run CGUTA/hugeheat \ 
--input_file \ #[path or glob] csv of numeric matrix
--outdir \ #[path] where you want the heatmap to be printed
--size \ #[u integer]
--background_color_u8 \ #[0-255] from black to white
--threshold=0.99 #[0-1] truncating threshold
 

