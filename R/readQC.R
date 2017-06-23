##############################
# Functions
##############################
#######################################################################
# Thu Jun 15 14:53:20 2017 ------------------------------
# check if a read is clipping by examing the cigar (S or H)
#######################################################################
isClipped <- function(cigar, tol=3){
  # if clipping occurs in both side, return TRUE for not confident
  # if clipping occurs in one side, but  > tol bases, return TRUE
  # if no clipping or only < tol bases clipping, return FALSE
  flag = FALSE
  if (str_detect(cigar,"[SH]")) {
    clip_str = str_extract_all(cigar,"[0-9]+[SH]",simplify = TRUE)
    if (ncol(clip_str) > 1) {# has both sided clipping
      flag = TRUE
    } else {
      clip_base = as.numeric(str_extract_all(clip_str,"[0-9]+", simplify=TRUE))
      if (clip_base > tol) { # one sided clipping but bases > tol
        flag = TRUE
      }
    }
  }
  return(flag)
}

isClipped_mateClipped <- function(cigar, rcigar, tol=3){
  # check for reads and their mates
  n1 = length(cigar)
  n2 = length(rcigar)
  if(n1 != n2) {
    stop("cigar and rcigar has unequal length! ")
  }

  isFail = sapply(1:n1, function(i){
    isClipped(cigar[i], tol) | isClipped(rcigar[i], tol)
  })

  return(isFail)
}

#######################################################################
#  check if a read is unimapped read, that is, no multiple alignment or
#  supplementary alignment found for this read. By examining the 'SA' or
#  'XA' field
#######################################################################
isUniMapped <- function(SA){
  n = length(SA)
  f = !is.na(SA)
  if(length(f) == 1){
    return(rep(TRUE, n))
  }else{
    return(!f)
  }
}

#######################################################################
#  check if a read has sufficient mapping quality (MAPQ)
#######################################################################
isPass_MAPQ <- function(MAPQ, q=30){
  MAPQ >= q
}
isPass_MAPQ_mateMAPQ <- function(mq,rmq, q=30){
  isPass_MAPQ(mq, q) & isPass_MAPQ(rmq, q)
}

#######################################################################
# check if a read has correct insert size. (isize)
#######################################################################
isPass_isize <- function(isize, q=500){
  abs(isize) < q
}

#######################################################################
# check if a read carrying too many variants. (MD)
#######################################################################
isHyperMutatedRead <- function(MD, cigar, q=4){
  # Fri Jun 16 11:20:44 2017 ------------------------------
  # Some reads are weird as they have no "MD:Z" field in BAM.
  # Perhaps after indel realignmnet, something has changed?
  # for this case, I have to remove the reads with MD equals NA
  # and set the reads as non-pass.
  if(is.na(MD)){
    return(FALSE)
  }
  len.md = str_count(MD, "\\^*[A-Z]+")
  len.cigar = str_count(cigar,"I") # times of insertions
  len.nm = len.md + len.cigar
  return(len.nm >= q)
}

#######################################################################
# convert ASCII character to its Phred quality
#######################################################################
convertPhredCharacter <- function(qual, phred = 33) {
  if(!phred %in% c(33, 64)) {
    stop("phred should be only 33 or 64 !")
  }
  as.numeric(charToRaw(qual)) - phred
}

#######################################################################
#  check if a read has (q<20) bases > n%, e.g. 50% (qual)
#######################################################################
isTooManyLowQualBases <- function(qual, phred = 33, q = 20, tol = 0.5) {
  # qual: bam QUAL field, the same as FASTQ QUAL field. ASCII character
  # phred: 33 or 64, for illumina Phred+33 or +64, respectively
  # q: threshold for base quality.
  # tol: maximun allowed threshold
  # return: TRUE for Yes, FALSE for No
  qual = convertPhredCharacter(qual, phred)
  isFail = (sum(qual < q)/length(qual)) >= tol
  return(isFail)
}


#######################################################################
# Fri Jun 23 09:46:02 2017 ------------------------------
# calculate overlap bases between two regions
#######################################################################
computeOverlapBases <- function(reg1, reg2) {
  # compute how many bases in reg1 overlapping reg2
  # return overlapped bases
  # reg1 and reg2 is a data.frame, with column 1: chr, 2:start, 3:end
  width = reg1$end - reg1$start + 1
  if(class(reg2) != "data.frame") {
    reg2 = as.data.frame(reg2)
  }
  colnames(reg1)[1:3] = c("chr","start","end")
  colnames(reg2)[1:3] = c("chr","start","end")

  n.r2 = nrow(reg2) # howmany regions in reg2
  OB = rep(0, n.r2)

  for(i in 1:n.r2){
    f.chr = reg1[,"chr"] == reg2[i, "chr"]
    s.start = reg2[i, "start"] - reg1[,"start"] # dropped bases from start
    s.start[s.start < 0] = 0
    s.end = reg1[, "end"] - reg2[i, "end"] # dropped bases from end
    s.end[s.end < 0] = 0
    OB[i] = sum(width[f.chr] - s.start[f.chr] - s.end[f.chr])
  }

  return(OB)
}


#######################################################################
# Thu Jun 22 14:58:07 2017 ------------------------------
#  Wrapper filtering low-quality reads.
#######################################################################
wrapper_filter_reads <- function(which, bamName, what, flag, tag){
  TR = 0
  PR = 0
  TB = 0
  PB = 0

  # initialize the output list
  output = list(
    bam.pass = NULL, # return filter bam as a data.frame, as a format of 'samtools view *.bam chr*:*-*'.
    # if no read left, bam.pass will return NULL
    total.reads = TR, # total reads spanning the queried region
    pass.reads = PR, # total reads spanning the queried region that pass the read filters
    total.bases = TB, # total bases of a region be sequenced, note that TB = sum(read.length - outside_bases)
    pass.bases = PB  # total bases of a region be sequenced.
  )

  if(is.null(what)){
    flag = scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F,
                       isSecondaryAlignment = F, isNotPassingQualityControls = F,
                       isDuplicate = F)
  }
  if(is.null(flag)){
    what = c("qname","flag", "mapq", "isize", "seq", "qual")
  }
  if (is.null(tag)) {
    tag = c("MQ", "MC", "SA", "MD", "NM","XA")
  }

  # scan bam
  param <- ScanBamParam(flag=flag, which=which, what=what, tag=tag)
  bamCon = BamFile(bamName, asMates = FALSE) # see $mate_status for mated pairs (?BamFile)
  bam = as.data.frame(readGAlignments(bamCon,param=param,use.names = FALSE))

  TR = nrow(bam)
  if(TR < 1){ # no read cover this region
    return(output)
  }
  # compute overlap bases
  TB = computeOverlapBases(bam[,c("seqnames","start","end")], which)

  # filtering low confidence reads. TRUE for pass, FALSE for fail
  qc_flag = sapply(1:TR, function(i){
    flag = (!isClipped(bam$cigar[i], tol=5)) &
      (isPass_MAPQ(bam$mapq[i])) &
      isPass_isize(bam$isize[i]) &
      isUniMapped(bam$SA[i]) &
      isUniMapped(bam$XA[i]) &
      (!isHyperMutatedRead(bam$MD[i], bam$cigar[i], q = 4)) &
      bam$width[i] > 60
    if(!is.na(bam$MC[i])) {
      base_flag = flag & (!isClipped(bam$MC[i], tol=8))
    }
    base_flag
  })

  bam.pass = bam[qc_flag & !is.na(qc_flag),]
  PR = nrow(bam.pass)


  if(PR > 0){ # has pass reads
    output[["bam.pass"]] = bam.pass
  }
  # compute overlap bases
  PB = computeOverlapBases(bam.pass[,c("seqnames","start","end")], which)

  # set output
  output[["total.reads"]] = TR
  output[["pass.reads"]] = PR
  output[["total.bases"]] = TB
  output[["pass.bases"]] = PB

  return(output)
}

