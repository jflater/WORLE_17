OTU_to_column <- function(physeq) {
	tax_table(physeq) <- cbind(tax_table(physeq), 
    		rownames(tax_table(physeq)))

	colnames(tax_table(physeq)) <- 
  		c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTUID")
	return(physeq)
}	
