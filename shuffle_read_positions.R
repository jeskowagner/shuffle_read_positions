### Author : Jesko Wagner
### Date: 2020.01.13

### Aim: shuffle positions of any genomic range across the chromosome
### to estimate if the read is positioned at its location by chance
### Note: currently this does not support strandedness.

shuffle_sv_positions = 
  function(sv.df,
           genome = "hg19",
           key.cols = key(sv.df),
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           chrom.sizes = NULL,
           chrom.sizes.name.col = "UCSC_seqlevel",
           chrom.sizes.length.col = "UCSC_seqlength",
           verbose = F) {
    
    # Dependencies
    library(GenomicRanges) # for getting chromosome sizes
    library(data.table)
    
    ### Sanity checks
    if(verbose) message("Checking data integrity.")
    
    if(!any(is.data.frame(sv.df) | "GRanges" %in% class(sv.df))) {
      stop("Data must be data.frame or GRanges.")
    }
    
    if(sv.df[,any(is.na(get(key.cols)))]) {
      stop("Data must not have NA values in chr, start, or end columns.")
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
    
    # Convert to dt for memory efficiency and joining
    sv.df = as.data.table(sv.df)
    
    # Keying allows keyed join for additional speed
    setkey(chrom.sizes, UCSC_seqlevel)
    setkeyv(sv.df, key.cols)

    # Add temporary column containing the length of the read
    if(verbose) message("Calculating read length.")
    
    ### For SVs other than TRA
    ### Note that SVs do not have an SV length, pos_alt is the new position
    # Get indices
    no_tra = sv.df[get(key.cols[1]) == get(key.cols[3]), which = TRUE]
    
    # calculate SV length
    sv.df[no_tra, processing_width_tmp := (get(key.cols[4]) - get(key.cols[2])), by = names(sv.df)]
    
    
    # Left join the chromosome sizes for chrom_ref
    sv.df[chrom.sizes, processing_width_chr1 := UCSC_seqlength]
    
    # Left join the chromosome sizes for chrom_alt
    setkeyv(sv.df, c(key.cols[3], key.cols[4], key.cols[1], key.cols[2]))
    sv.df[chrom.sizes, processing_width_chr2 := UCSC_seqlength]
    
    setkeyv(sv.df, key.cols) # return to original sorting
    
    # Calculate new start positions
    if(verbose) message("Assigning new read positions.")
    
    # For svs other than TRA
    sv.df[no_tra,
             # Generate random number between 1 and the size of the chromosome minus the width of the read
             eval(key.cols[2]) := round(digits = 0,
                                      runif(n = 1,
                                            min = 1,
                                            max = processing_width_chr1 - processing_width_tmp - 1 # honor 1-based coordinates
                                      )), 
             by = names(sv.df)
             ]
    
    # Add the new end postion of the read based on the previously calculated width
    sv.df[no_tra, eval(key.cols[4]) := get(key.cols[2]) + processing_width_tmp, ]
    
    # For TRA, just assign random positions
    sv.df[!no_tra, eval(key.cols[2]) := round(digits = 0, runif(n = 1, min = 1, max = processing_width_chr1)), names(sv.df)]
    sv.df[!no_tra, eval(key.cols[4]) := round(digits = 0, runif(n = 1, min = 1, max = processing_width_chr2)), names(sv.df)]
    
    # Remove temporary columns
    sv.df[, processing_width_tmp := NULL]
    sv.df[, processing_width_chr1 := NULL]
    sv.df[, processing_width_chr2 := NULL]
    
    # Return
    sv.df[]
  }



### Aim: shuffle positions of any genomic range across the chromosome
### to estimate if the read is positioned at its location by chance
### Note: currently this does not support strandedness.

shuffle_read_positions = 
  function(reads.df,
           genome = "hg19",
           chr.col = "chr",
           start.col = "start",
           end.col = "end", 
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           chrom.sizes = NULL,
           chrom.sizes.name.col = "UCSC_seqlevel",
           chrom.sizes.length.col = "UCSC_seqlength",
           verbose = F) {
    
    # Dependencies
    library(GenomicRanges)
    library(data.table)
    
    ### Sanity checks
    if(verbose) message("Checking data integrity.")
    
    if(!any(is.data.frame(reads.df) | "GRanges" %in% class(reads.df))) {
      stop("Data must be data.frame or GRanges.")
    }
    
    if(reads.df[,
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
    reads.df = as.data.table(reads.df)
    
    # Add temporary column containing the length of the read
    if(verbose) message("Calculating read length.")
    
    reads.df[,processing_width_tmp := (get(end.col) - get(start.col)), by = names(reads.df)]
    
    
    # Left join the chromosome sizes
    reads.df[chrom.sizes, processing_width_chr := UCSC_seqlength]
    
    # Calculate new start positions
    if(verbose) message("Assigning new read positions.")
    
    reads.df[,
           # Generate random number between 1 and the size of the chromosome minus the width of the read
           eval(start.col) := round(digits = 0,
                                    runif(n = 1,
                                          min = 1,
                                          max = processing_width_chr - processing_width_tmp - 1 # honor 1-based coordinates
                                    )), 
           by = names(reads.df)
           ]
    
    # Add the new end postion of the read based on the previously calculated width
    reads.df[, eval(end.col) := get(start.col) + processing_width_tmp, ]
    
    
    
    # Remove temporary column
    reads.df[, processing_width_tmp := NULL]
    reads.df[, processing_width_chr := NULL]
    
    # Return
    reads.df[]
  }

### Aim: update binned genome to count read ends per bin and add them to the possibly present count column
### Warning: this will  update your object in the calling environment to stay efficient
add_reads_per_bin = function(reads.df, bin.df, reads.key = key(reads.df), bin.key = key(bin.df)) {
  
  library(data.table)
  
  # Sanity check
  if(!is.data.table(reads.df)) reads.df = as.data.table(reads.df)
  if(!is.data.table(bin.df))   bin.df = as.data.table(bin.df)
  
  # Assign keys for fast overlapping
  if(!haskey(reads.df)) setkeyv(reads.df, reads.key)
  if(!haskey(bin.df))   setkeyv(bin.df, bin.key)
  
  # Copy start and end column to make ranges of: start to start, end to end
  # This allows to find only the read ends which overlap
  reads.df[,eval(paste0(reads.key[2], "_copy")) := get(key(reads.df)[2])]
  reads.df[,eval(paste0(reads.key[3], "_copy")) := get(key(reads.df)[3])]
  
  
  ## Find overlaps
  # Start position in bins
  start_in_bin = foverlaps(which = T, reads.df, bin.df, by.x=c(reads.key[1], reads.key[2], paste0(reads.key[2],"_copy")))
  
  # Prepare for fast join
  setkey(start_in_bin, yid)
  start_in_bin[,xid := NULL]

  # End position in bins
  end_in_bin   = foverlaps(which = T, reads.df, bin.df, by.x=c(reads.key[1], reads.key[3], paste0(reads.key[3],"_copy")))
  
  # Prepare for fast join
  setkey(end_in_bin, yid)
  end_in_bin[,xid := NULL]
  
  ## Merge
  read_in_bin = rbind(start_in_bin, end_in_bin)
  read_in_bin[,count := .N, by = yid]
  read_in_bin = unique(read_in_bin)
  setkey(read_in_bin, yid)
  
  ### Sanity check
  if(nrow(reads.df) * 2 != read_in_bin[,sum(count)])
    stop("Data failed sanity check. More of fewer read ends placed than there are reads (reads * 2 = possible ends to be placed).")
  ### End sanity check
  
  ## Add counts back onto bin.df
  # If the bins do not yet contain counts, then first add the column
  if(!"count" %in% names(bin.df)) {
    bin.df[,count:=0]
  }
  
  # Merge counts back onto bin df
  bin.df[read_in_bin, count := count+i.count, on = "yid"]
  
  # Remove temporary columns
  reads.df[,eval(paste0(reads.key[2], "_copy")) := NULL]
  reads.df[,eval(paste0(reads.key[3], "_copy")) := NULL]
}
