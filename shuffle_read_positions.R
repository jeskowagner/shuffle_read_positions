### Author : Jesko Wagner
### Date: 2020.01.13
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
           chrom.sizes = NULL,
           chrom.sizes.name.col = "UCSC_seqlevel",
           chrom.sizes.length.col = "UCSC_seqlength",
           verbose = F) {
    
    # Dependencies
    library(GenomicRanges)
    library(data.table)
    
    ### Sanity checks
    if(verbose) message("Checking data integrity.")
    
    if(!any(is.data.frame(dat.df) | "GRanges" %in% class(dat.df))) {
      stop("Data must be data.frame or GRanges.")
    }
    
    if(dat.df[,
              any(is.na(get(chr.col))) |
              any(is.na(get(start.col))) |
              any(is.na(get(end.col)))
              ]) {
      stop("Data must not have NA values in chr, start, or end columnw.")
    }
    
    # Have chromosome sizes been provided? If no, download. If yes, check data
    if(is.null(chrom.sizes)) {
      # Prepare, download chromosome sizes and obtain only those of required chromosomes
      if(verbose) message("Downloading chromosome information from UCSC.")
      
      chrom.sizes = fetchExtendedChromInfoFromUCSC(genome, quiet = T)
      chrom.sizes = as.data.table(chrom.sizes)
      chrom.sizes = chrom.sizes[UCSC_seqlevel %in% chromosomes,
                                .(UCSC_seqlevel, UCSC_seqlength)]
      
    } else {
      chrom.sizes = as.data.table(chrom.sizes)  
      
      if(!is.null(chrom.sizes.name.col)) 
        setnames(chrom.sizes, old = chrom.sizes.name.col, new = "UCSC_seqlevel")
      
      if(!is.null(chrom.sizes.length.col)) 
        setnames(chrom.sizes, old = chrom.sizes.length.col, new = "UCSC_seqlength")
    }
    
    # Check that all entered chromosome information is correct
    if (chrom.sizes[, 
                    !is.numeric(UCSC_seqlength) | 
                    !(length(unique(UCSC_seqlevel)) == length(UCSC_seqlevel))
                    ]) {
      stop("chrom.sizes contains invalid or duplicate values.")
    }
    
    
    ### End sanity checks
    
    # Keying allows keyed join for additional speed
    setkey(chrom.sizes, UCSC_seqlevel)
    
    # Convert to dt for memory efficiency
    dat.df = as.data.table(dat.df)
    
    # Add temporary column containing the length of the read
    if(verbose) message("Calculating read length.")
    dat.df[,processing_width_tmp := (get(end.col) - get(start.col)), by = names(dat.df)]
    
    if(verbose) message("Assigning new read positions.")
    # Loop over chromosomes
    for(cur_chrom in dat.df[,as.character(unique(get(chr.col)))]) {
      
      dat.df[get(chr.col) == cur_chrom, # select only reads from the same chromosome
             
             # Generate random number between 1 and the size of the chromosome minus the width of the read
             eval(start.col) := round(digits = 0,
                                      runif(n = 1,
                                            min = 1,
                                            max = chrom.sizes[cur_chrom, UCSC_seqlength] - processing_width_tmp - 1 # honor 1-based coordinates
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
