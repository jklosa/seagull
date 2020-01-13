#' @title Simulated genotypes
#' 
#' @description A genotype matrix that contains information from single
#' nucleotide polymorphisms (SNPs). Dimensions are 1000 rows, 466 columns. The
#' data was simulated using the software
#' \href{https://alphagenes.roslin.ed.ac.uk}{AlphaSim} of which an R-package
#' is available. Each row corresponds to a single individual. 1000 individuals
#' were simulated, where 10 half sib families were created. Each family consists
#' of 100 half sibs. The half sibs share a common Sire. 466 SNPs are available,
#' distributed over 2 chromosomes to an equal amount, i.e., the first 233 SNPs
#' are located on chromosome 1, the remaining SNPs are on the second chromosome.
#' The complementary homozygote genotypes are coded as 0 and 2, respectively.
#' The heterozygote genotype as 1.
#' 
#' @name genotypes
#' 
#' @docType data
#' 
#' @author Jan Klosa \email{klosa@fbn-dummerstorf.de}
#' 
#' @format A matrix with 1000 rows and 466 columns. One row per individual.
#' One column per SNP.
"genotypes"
