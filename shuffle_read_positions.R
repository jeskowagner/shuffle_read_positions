### Author : Jesko Wagner
### Date: 2020.01.10
### Aim: shuffle positions of any genomic range across the chromosome
### to estimate if the read is positioned at its location by chance
### Note: currently this does not support strandedness.

shuffle_read_positions = 
  function(dat.df,
           genome = "hg19",
           chr.col = "chr",
           start.col = "start",
           end.col = "end", 
           chromosomes = paste0("chr", c(1:22, "X","Y")),
           chrom_sizes = NULL,
           keep.extra.columns = T, verbose = F) {
    
    # Dependencies
    library(GenomicRanges)
    library(data.table)
    
    # Sanity checks
    if("data.frame" %in% class(dat.df) || "GRanges" == class(dat.df)) {
      stop("Data must be data.frame or GRanges.")
    }
    
    if(dat.df[,any(is.na(get(chr.col)))]) {
      stop("Data must not have NA chromosome names.")
    }
    
    if(dat.df[,any(is.na(get(start.col)))]) {
      stop("Data must not have NA start positions.")
    }
    
    if(dat.df[,any(is.na(get(end.col)))]) {
      stop("Data must not have NA end positions.")
    }
    
    
    if(is.null(chrom_sizes)) {
      # Prepare, download chromosome sizes and obtain only those of required chromosomes
      chrom_sizes = fetchExtendedChromInfoFromUCSC(genome, quiet = T)
      chrom_sizes = as.data.table(chrom_sizes)
      chrom_sizes = chrom_sizes[UCSC_seqlevel %in% chromosomes,
                            .(UCSC_seqlevel, UCSC_seqlength)]
      setkey(chrom_sizes, UCSC_seqlevel)
    } else {
      chrom_sizes = as.data.table()  # if providing chromosomes sizes, column with sizes has to be either chr or UCSC_seqlevel
      # this needs further improvement to check also for chromosome names and data integrity
      if("chr" %in% names(chrom_sizes)) setnames(chrom_sizes, chr, UCSC_seqlevel)
      setkey(chrom_sizes, UCSC_seqlevel)
    }
    
    
    
    # Convert to dt for memory efficiency
    dat.df = as.data.table(dat.df)
    
   
    
    # Add temporary column containing the length of the read
    dat.df[,processing_width_tmp := (get(end.col) - get(start.col)), by = names(dat.df)]
    
    # Loop over chromosomes
    for(cur_chrom in dat.df[,as.character(unique(get(chr.col)))]) {
      
      dat.df[get(chr.col) == cur_chrom, # select only reads from the same chromosome
             
             # Generate random number between 1 and the size of the chromosome minus the width of the read
             eval(start.col) := round(digits = 0,
                                      runif(n = 1,
                                            min = 1,
                                            max = chrom_sizes[cur_chrom, UCSC_seqlength] - processing_width_tmp - 1 # honor 1-based coordinates
                                      )), 
             by = names(dat.df)
             ]
      
      # Add the new end postion of the read based on the previously calculated width
      dat.df[get(chr.col) == cur_chrom, eval(end.col) := get(start.col) + processing_width_tmp, ]
      
    }
    
    # Remove temporary column
    dat.df[, processing_width_tmp := NULL]
    
    # Return
    dat.df[]
  }
