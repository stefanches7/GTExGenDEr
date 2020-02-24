recoverTissueName <- function(filesystemname) {
  recoveredDashes <- gsub("___", " - ", filesystemname) %>% gsub(pattern = "_1", replacement = "-1") %>% gsub(pattern = "EBV_transformed", replacement = "EBV-transformed")
  recoveredBrackets <- recoveredDashes %>% gsub(pattern = "8_", replacement = "(") %>% gsub(pattern = "_9", replacement = ")")
  currentTissue <- gsub("_", " ", recoveredBrackets)
  currentTissue
}