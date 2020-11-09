# Model simulations (top half will generate single graphs, bottom will collate all sample graphs into a single image)

#This will simulate the fitted models to see if the outcome resembles the input (individual image version)
for (ww in 1:csvCount) {
  zzz <- simulate(XYCauList1[[paste0("RAG", ww)]], nsim = 4);
  tiff(paste0("CauSim", ww,".tif"), width = 100, height = 30, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
  plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 1, ncols = 4, pch = 1, cols = alpha("red", 0.4), maxsize = 1, main = paste0("RAG", ww, " - Cauchy Process Simulations"), main.panel = "");
  dev.off();
  
  zzz <- simulate(XYCoxList1[[paste0("RAG", ww)]], nsim = 4);
  tiff(paste0("CoxSim", ww,".tif"), width = 100, height = 30, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
  plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 1, ncols = 4, pch = 1, cols = alpha("red", 0.4), maxsize = 1, main = paste0("RAG", ww, " - Log-Gaussian Cox Process Simulations"), main.panel = "");
  dev.off();
  
  zzz <- simulate(XYMatList1[[paste0("RAG", ww)]], nsim = 4);
  tiff(paste0("MatSim(1)", ww,".tif"), width = 100, height = 100, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
  plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 2, ncols = 2, pch = 1, cols = alpha("red", 0.4), maxsize = 1, main = paste0("RAG", ww, " - Matern Cluster Process Simulations"), main.panel = "");
  dev.off();
  
  zzz <- simulate(XYTomList1[[paste0("RAG", ww)]], nsim = 4);
  tiff(paste0("TomSim", ww,".tif"), width = 100, height = 30, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
  plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 1, ncols = 4, pch = 1, cols = alpha("red", 0.4), maxsize = 1, main = paste0("RAG", ww, " - Thomas Process Simulations"), main.panel = "");
  dev.off();
  
  zzz <- simulate(XYVarList1[[paste0("RAG", ww)]], nsim = 4);
  tiff(paste0("VarSim", ww,".tif"), width = 100, height = 30, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
  plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 1, ncols = 4, pch = 1, cols = alpha("red", 0.4), maxsize = 1, main = paste0("RAG", ww, " - Varaince Gamma Cluster Process Simulations"), main.panel = "");
  dev.off();
}

rm(zzz);

#This will simulate the fitted models to see if the outcome resembles the input (single image version)
zzz <- anylapply(XYCauList1, simulate, nsim = 4);
tiff(paste0("CauSim.tif"), width = 280, height = 420, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 12, ncols = 8, pch = "20", cols = alpha("red", 0.4), maxsize = 1, main = "Cauchy Cluster Process Simulations", main.panel = "");
dev.off();

zzz <- anylapply(XYCoxList1, simulate, nsim = 4);
tiff(paste0("CoxSim.tif"), width = 280, height = 420, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 12, ncols = 8, pch = "20", cols = alpha("red", 0.4), maxsize = 1, main = "Log-Gaussian Cox Process Simulations", main.panel = "");
dev.off();

zzz <- anylapply(XYMatList1, simulate, nsim = 4);
tiff(paste0("MatSim.tif"), width = 280, height = 420, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 12, ncols = 8, pch = "20", cols = alpha("red", 0.4), maxsize = 1, main = "Matern Cluster Process Simulations", main.panel = "");
dev.off();

zzz <- anylapply(XYTomList1, simulate, nsim = 4);
tiff(paste0("TomSim.tif"), width = 280, height = 420, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 12, ncols = 8, pch = "20", cols = alpha("red", 0.4), maxsize = 1, main = "Thomas Process Simulations", main.panel = "");
dev.off();

zzz <- anylapply(XYVarList1, simulate, nsim = 4);
tiff(paste0("VarSim.tif"), width = 280, height = 420, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
plot(zzz, ylim = c(ImageSize[paste0("RAG", ww), 2]/ConvSca1, 0), arrange = TRUE, nrows = 12, ncols = 8, pch = "20", cols = alpha("red", 0.4), maxsize = 1, main = "Variance Gamma Cluster Process Simulations", main.panel = "");
dev.off();

rm(zzz);