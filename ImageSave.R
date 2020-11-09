#To export the plot results of skin a function is build for repeat calling
Image.Save <- function(x, counter, width.modifier1, width.modifier2, file.name, test.name, ...){
  tiff(paste0(file.name, counter,".tif"), width = ImageSize[paste0("RAG", counter), 1]/ConvSca1 * width.modifier1, height = ImageSize[paste0("RAG", counter), 2]/ConvSca1 * width.modifier2, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
  plot(x, xlim = NULL, ylim = c(ImageSize[paste0("RAG", counter), 2]/ConvSca1, 0), main = paste0("RAG", counter, " - ", test.name), ...);
  dev.off();
}