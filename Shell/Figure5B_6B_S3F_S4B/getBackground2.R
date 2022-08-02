library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
registerDoParallel(cores = 2)

ulimit::memory_limit(50000)
print(gc())

args=commandArgs(T)
file=args[1]
times=as.numeric(args[2])
targetRegionsFile=args[3]
  
targetRegions=read.table(targetRegionsFile, sep="\t", stringsAsFactors = F)
targetRegions$V2=targetRegions$V2+1
targetRegions <- GRanges(seqnames = targetRegions[,1], IRanges(start = targetRegions[,2], end = targetRegions[,3]), 
                seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))

data=read.table(file, sep="\t", stringsAsFactors = F)
dmrs <- GRanges(seqnames = data[,1], IRanges(start = data[,2], end = data[,3]), seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
dmrs <- sort(dmrs)
dmrs <- keepStandardChromosomes(dmrs, pruning.mode = "coarse")
start.time <- Sys.time() # start timing
background.set <- GRanges() # initialize
chr_matches <- vmatchPattern("CG", BSgenome.Hsapiens.UCSC.hg38) # find CpGs
chr_matches <- chr_matches[strand(chr_matches) == "+", ] # keep only plus strand
dmrs$nCpG <- countOverlaps(dmrs, chr_matches, ignore.strand = TRUE) # count overlaps between dmrs and cpgs
chr_tiles <- tileGenome(seqinfo(dmrs), tilewidth = 1000L, cut.last.tile.in.chrom = TRUE) # set up genome tiling
set.seed(4242, kind = "L'Ecuyer-CMRG") # Set random seed that stays consistent in parallel
chr_tiles <- trim(suppressWarnings(shift(chr_tiles, shift = sample.int(1000L, size = 1)))) # shift the start points randomly (all the same distance)
for (cpgCount in unique(sort(dmrs$nCpG))) { # loop for each distinct count value of CpGs per dmr
  message("cpg count is ", cpgCount,
          " . This is ", which(unique(sort(dmrs$nCpG)) %in% cpgCount),
          " out of ", length(unique(sort(dmrs$nCpG))))
  timenow <- Sys.time() - start.time # report time
  message(timenow, " ", attr(timenow, "units"), " has elapsed")
  dmr_sub <- dmrs[dmrs$nCpG == cpgCount, ] # operate only on dmrs with the current count
  sub.background.set <- foreach(dmrWidth = unique(sort(width(dmr_sub))), # loop in parallel over dmr widths
                                .combine = c,
                                .packages = "BSgenome.Hsapiens.UCSC.hg38" # dont forget to include your genome build here
  ) %dopar% {
    tiles <- trim(suppressWarnings(resize(chr_tiles, dmrWidth))) # resize initial tiles to the correct dmr width
    tiles$nCpG <- countOverlaps(tiles, chr_matches, ignore.strand = TRUE) # count all the overlaps
    tiles <- tiles[tiles$nCpG == cpgCount,] # keep only those tiles with the correct cpg count
    tiles=tiles[!is.na(findOverlaps(tiles, targetRegions,type="any", select="first",ignore.strand=T,minoverlap=dmrWidth)),]
    
    fg.count <- sum(width(dmr_sub) == dmrWidth)
    multiplier <- times
    tiles <- tryCatch(sort(sample(tiles, fg.count * multiplier)),
                      error = function(err) {
                        multiplier <- sum(tiles$nCpG == cpgCount) / fg.count
                        sort(sample(tiles, fg.count * multiplier))
                      }) # sample from these tiles to generate a 10X background, or barring that, as big a background as it can
  }
  background.set <- suppressWarnings(c(background.set, sub.background.set)); rm(sub.background.set); gc() # append and clean up
}
background.set <- sort(background.set) # done!
save(background.set, file=gsub(".bed",paste0(".BackgroundMotif.", times, "X.RData"),file))
timenow <- Sys.time() - start.time
message("ran in ", timenow, " ", attr(timenow, "units"))
