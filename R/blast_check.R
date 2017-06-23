#######################################################################
# Tue Jun 20 10:27:28 2017 ------------------------------
# BLAST to check mapping uniqueness of sequence.
#######################################################################
blast_check <- function(location, blast_dir, blast_db, blast_flank_bp, threads){
  # location: a data.frame to store position of variants, with colnames.
  #           column 1: chr, column 2: pos
  # balst_dir: path to blastn
  # blast_db: makeblastdb prepared blast databases.
  # blast_flank_bp: how many basepairs around the variant should be extracted
  # threads: threads used to do blast
  # output: a nonduplicated data.frame with columns
  #         1. chr
  #         2. pos
  #         3. uniquely mapped or not(TRUE or FALSE)?
  suppressMessages(library(Biostrings))
  tmpfile = tempfile(pattern = "file", tmpdir = tempdir(), fileext = "") # use temp file

  location = unique(location) # make query uniqueness.
  # set ranges to extract sequences
  location$start = location$pos - blast_flank_bp
  location$end = location$pos + blast_flank_bp
  location$index = with(location, str_c(chr,":",pos))
  location$pos = NULL # remove pos column
  # build GRanges object
  gr = makeGRangesFromDataFrame(location, ignore.strand=TRUE, keep.extra.columns=TRUE)

  ## get fasta and write to file
  fa = FaFile(ref_fasta) # load indexed fasta
  seq = getSeq(fa, gr)
  names(seq) = gr$index # rename seqs
  # write tmp fasta file
  writeXStringSet(seq, tmpfile, format = "fasta")

  ## do blast
  cat("Doing blastn, please wait for a long? time....\n")
  include_flags = "qseqid qlen sacc slen qstart qend sstart send bitscore score pident nident mismatch gaps evalue qcovs"
  # set command line string
  cmd = sprintf("%s/blastn -query %s -db %s -outfmt '6 %s' -num_threads %i -max_target_seqs 5 -qcov_hsp_perc 50",
                blast_dir, tmpfile, blast_db, include_flags, nt)
  blast_result <- read.table(pipe(cmd)) # pipe
  # add column name
  colnames( blast_result ) <- as.vector(str_split(include_flags," ",simplify = T))
  cat("Blastn done! Now filtering...\n")

  ## check multiple alignments
  output = location
  output$map.uniq = FALSE
  output$suppl.align = "" # supplemental alignment
  # foreach location
  for(i in 1:nrow(output)){
    idx = output[i,"index"]
    sub.out = blast_result[blast_result$qseqid == idx, ] # get blast results subset
    nm = nrow(sub.out) # how many match records?
    if(nm == 0){ # nothing found
      stop(str_c(idx," have nothing match!!!!"))
    } else if (nm == 1){ # uniquely mapped
      if(sub.out$qend - sub.out$qstart == 2*blast_flank_bp) {
        output[i,"map.uniq"] = TRUE
      } else {
        print(sub.out)
        stop(str_c(idx," have unknown problem !!!"))
      }
    } else if (nm > 1) { # multiple alignment found
      # remove itself match
      sub.out = subset(sub.out, !(sstart == output[i,"start"] & send == output[i,"end"] & sacc == output[i,"chr"]))
      output[i, c("map.uniq","suppl.align")] = with(sub.out, {
        is.sa = ((qstart <= blast_flank_bp / 2 & qend >= blast_flank_bp * 1.5)  & # read length
                   pident > 90) # similarity
        sa = ""
        if(any(is.sa)){
          sa = with(sub.out[is.sa,], str_c(
            sprintf("%s:%i-%i(%.1f)", sacc, sstart, send, pident),
            collapse = ";"
          ))
        }
        c(!any(is.sa), sa)
      })
    }
  }
  rm(tmpfile)
  ## return
  cat("Filtering done, return...\n")
  print(head(output))
  return(output)
}
