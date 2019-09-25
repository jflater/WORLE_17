# Add OTU data to OTU_table
OTU_to_column <- function(physeq) {
	tax_table(physeq) <- cbind(tax_table(physeq), 
    		rownames(tax_table(physeq)))

	colnames(tax_table(physeq)) <- 
  		c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTUID")
	return(physeq)
}	
# Convert phyloseq OTU to vegan OTU table orientation (transpose)
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}