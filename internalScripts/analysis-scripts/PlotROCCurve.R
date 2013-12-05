plot_roc <- function(results.file, title) {
  # Load the ROCR library.
  library(ROCR)
  
  # Get the data file produced by GeneratePhenotypeROCCurve script.
  results <- read.csv(results.file)
  
  # Calculate number of known reaction probabilities.
  match <- grep("unknown", results$reaction, invert=T)
  results.known <- results[match,]
  subtitle <- sprintf("%d genes with known probabilities", nrow(results.known))
  
  # Calculate the performance of the algorithm.
  pred <- prediction(results.known$prob, results.known$true)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, lwd=2, main=title, sub=subtitle) # type="b", lwd=1.5 
}

plot_rocx2 <- function(results.file1, results.file2, title) {
  # Load the ROCR library.
  library(ROCR)
  
  # Get the data file produced by GeneratePhenotypeROCCurve script.
  results1 <- read.csv(results.file1)
  results2 <- read.csv(results.file2)
  
  # Calculate number of known reaction probabilities.
  match1 <- grep("unknown", results1$reaction, invert=T)
  results1.known <- results1[match1,]
  match2 <- grep("unknown", results2$reaction, invert=T)
  results2.known <- results2[match2,]  
  subtitle <- sprintf("%d and %d genes with known probabilities", nrow(results1.known), nrow(results2.known))
  
  # Calculate the performance of the algorithm.
  pred1 <- prediction(results1.known$prob, results1.known$true)
  perf1 <- performance(pred1, "tpr", "fpr")
  pred2 <- prediction(results2.known$prob, results2.known$true)
  perf2 <- performance(pred2, "tpr", "fpr")
  
  # Show both on the same plot.
  plot(perf1, lwd=2, main=title, sub=subtitle, col="blue")
  plot(perf2, lwd=2, add=T, col="orange")
}
  
pdf("171101.1.pdf")
par(mfrow=c(2,2))
plot_roc("171101.1.model.std.int.minimal.int.SP4.knockoutsim.results.csv", "Standard gap fill")
plot_roc("171101.1.model.pa.int.minimal.int.SP4.knockoutsim.results.csv", "Probabilistic gap fill")
plot_roc("171101.1.model.std.iterative.int.minimal.int.SP4.knockoutsim.results.csv", "Standard iterative gap fill")
plot_roc("171101.1.model.pa.iterative.int.minimal.int.SP4.knockoutsim.results.csv", "Probabilistic iterative gap fill")
dev.off()

pdf("171101.1.x2.pdf", width=4)
par(mfrow=c(2,1))
plot_rocx2("171101.1.model.std.int.minimal.int.SP4.knockoutsim.results.csv", "171101.1.model.pa.int.minimal.int.SP4.knockoutsim.results.csv", "Standard vs. Probabilistic gap fill")
plot_rocx2("171101.1.model.std.iterative.int.minimal.int.SP4.knockoutsim.results.csv", "171101.1.model.pa.iterative.int.minimal.int.SP4.knockoutsim.results.csv", "Standard vs. Probabilistic iterative gap fill")
dev.off()
