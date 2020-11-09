#To export the resulting graphs a function is created for recall
Graph.Save <- function(x, counter, width.modifier, file.name, test.name, ...) {
  tiff(paste0(file.name, counter, ".tif"), width = 150 * width.modifier, height = 150, units = "mm", res = GraphRes, pointsize = 12, compression = 'lzw');
  plot(x, main = paste0("RAG", counter, " - ", test.name), ...);
  dev.off();
}