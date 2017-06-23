#######################################################################
#  check if a read carrying queried variants.
# Thu Jun 15 15:31:55 2017 ------------------------------
#   for this version, I have not yet check the exact sequence of insertion
# Fri Jun 16 13:09:08 2017 ------------------------------
#   add base quality (qual) filtering. If the mutated base is of low quality,
#   view it as not passed reads(i.e. let isMutated = FALSE)
#######################################################################
isMutatedRead <- function(start, mpos, ref, alt, MD, cigar ){
  # one read per input:
  # start: read start;
  # mpos: mutation position; ref: ref base; alt: alt base
  # Note: pos, ref and alt should be 1-based coordination
  # MD: mismatch tag; cigar: cigar;
  # return a list() with elements
  #   1. isMutated: TRUE: read carrying mutation, FALSE: read not carrying mutation
  #   2. mpos: the position of mutation on this read. 0 for notfound(non-mutated reads)
  type="SNP"
  if(ref == "-"){
    type = "INS"
  }else if(alt == "-"){
    type = "DEL"
  }
  if(is.na(MD)) {
    # Fri Jun 16 11:20:44 2017 ------------------------------
    # Some reads have no MD field (NA). See isHyperMutatedRead() for detail
    return(list(isMutated = FALSE, mpos = 0))
  }
  mp = mpos - start +1 # mutation point in a read
  len = end - start + 1 # length of a read
  clen = 0 # pos of mutation according for reference ( counts DEL)
  rlen = 0 # pos of mutation on a read (dont count DEL)
  isMutated = FALSE
  if(type %in% c("SNP","DEL")){
    mdz = str_match_all(MD,"[0-9]+\\^*[A-Z]*")[[1]]
    for(j in 1:nrow(mdz)){ # iterater for each mismatch
      x = mdz[j,]
      a = as.numeric(str_extract(x, "[0-9]+"))
      b = !is.na(str_extract(x, "\\^")) # is a DEL?
      c = str_extract(x,"[A-Z]+")
      if(is.na(c)){break} # the end of a read
      clen = clen + a
      rlen = rlen + a
      if(mp - clen == 1){ # meet the mutation
        if(type == "SNP" & !b & c == ref){
          isMutated = TRUE
          rlen = rlen + 1
          break
        }else if(type == "DEL" & b & c == ref){
          isMutated = TRUE
          rlen = rlen + 1
          break
        }
      }
      # if not meet the mutation, go on searching
      clen = clen + nchar(c)
      if(!b){# dont add rlen if DEL
        rlen = rlen + nchar(c)
      }
    }
  }else{
    # for INS, only count cigar
    # to update: need check the insertion bases equals to the altbases
    mdz = str_match_all(cigar,"[0-9]+[MID]")[[1]]
    for(j in 1:nrow(mdz)){ # iterater for each mismatch
      x = mdz[j,]
      a = as.numeric(str_extract(x, "[0-9]+"))
      c = str_extract(x,"[A-Z]+")
      if(c == "M"){
        clen = clen + a
        rlen = rlen + a
      }else if(c == "D"){
        clen = clen + a
      }else if(c == "I"){
        if(mp - clen == 1){
          if(a == nchar(alt)){
            isMutated = TRUE
            rlen = rlen + 1
            break
          }
        }
        rlen = rlen + a
      }
    }
  }
  return(list(isMutated = isMutated, mpos = ifelse(isMutated, rlen, 0)))
}
