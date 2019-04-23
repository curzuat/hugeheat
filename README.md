A pipeline to make a visualization heatmap of large matrices with float values stored as csv

Procedure:

the pipeline reduces the size of the matrix by tiling it with bins of --size width,height (box sampling) cells covered by each bin get summarized into an average treating negative and positive values separately


the positive and negative values in the reduced matrix are separately scaled to fit a 0-255 range for displaying color intensities. The --threshold percentile value becomes 255 and higher values are truncated to 255. This facilitates contrast to be appreciated in the heatmap.

A shade of blue or red is assigned to each cell coresponding to the magnitude of the negative or positive values and a final color is created for the heatmap cell by linear interpolation of these two colors
Finally, the resulting colors are tinted with black inversely proportional to the absolute sum of the encoded negative and positive values


The heatmap is plotted using --background_color_u8 as the color for cells with 0 positive and 0 negative values should they exist. each cell is plotted over n,n --pixels, sometimes is desirable to exagerate one dimension of the cell for visualization purposes.

Usage
nextflow run CGUTA/hugeheat \ 
--input_file \ #[path or glob] csv of numeric matrix
--outdir \ #[path] where you want the heatmap to be printed 
--size \ #[Int,Int] # width and height of the binning reticule
--background_color_u8 \ #[0-255] from black to white
--threshold=0.99 #[0-1] truncating threshold

REQUIRED

--input_file \ #[path or glob] csv of numeric matrix

OPTIONAL

--outdir \ #[path] where you want the heatmap to be printed (default .)
--size \ #[Int,Int] # width and height of the binning tile (default 1,1)
--background_color_u8 \ #[0-255] from black to white (default 155)
--threshold #[0-1] truncating threshold (default 0.99)

EXAMPLE USAGE

nextflow run CGUTA/hugeheat \
--input_file my_big_matrix_with_header.csv \
--outdir . \
--size 2,1 \ # bin the rows in sets of two leave cols alone 
--threshold 1 \ # values are scaled to 0-255 but no truncation ocurrs

 

