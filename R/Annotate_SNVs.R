#!/usr/bin/env Rscript

library("optparse")

################ INPUT PARAMETERS
#define parameters
option_list = list(
  make_option(c("-p", "--panel.file"), type="character", default=NULL,help="Path to the smMIP design file [MANDATORY if -i was not supplied]", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,help="Path to a tab delimited table with the following columns: chr, pos, ref, alt. [MANDATORY if -p was not supplied]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,help="Path for the output (annotated panel). If not supplied, the output will be saved within the folder that contains the smMIP design file/input file", metavar="character"),
  make_option(c("-c", "--code"), type="character", default=getwd(),help="Path to smMIP tools source functions file. If not supplied, it is assumed that the file (smMIPs_Function.R) is located in your working directory", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=1,help="The number of cores to use for parallel processing", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="GRCh37",help="The genome assembly to be queried. For all the possible options, please refer to the cellbaseR package documentation", metavar="character"),
  make_option(c("-s", "--species"), type="character", default="hsapiens",help="The species to be queried. For all the possible options, please refer to the cellbaseR package documentation", metavar="character"),
  make_option(c("--skip-cellbase"), action="store_true", default=FALSE,help="Skip CellBase API annotation entirely. All annotation columns will be set to NA. Useful when CellBase is offline or unreachable."),
  make_option(c("--cellbase-timeout"), type="integer", default=30,help="Timeout in seconds for each CellBase API call (default: 30)"))

#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if ( is.null(opt$panel.file) & is.null(opt$input) ) {
  print_help(opt_parser)
  stop("Either the -p or -i parameter must be supplied", call.=FALSE)
}

if ( !is.null(opt$panel.file) & !is.null(opt$input) ) {
  print_help(opt_parser)
  stop("Either the -p or -i parameter must be supplied. Not both", call.=FALSE)
}

if ( is.null(opt$panel.file) & !is.null(opt$input)) {
  opt$output=dirname(opt$input)
  file.name=basename(opt$input)
}

if ( !is.null(opt$panel.file) & is.null(opt$input)) {
  opt$output=dirname(opt$panel.file)
  file.name=basename(opt$panel.file) 
}


if (is.null(opt$code)){ #load source functions
  source(paste0(getwd(),"/smMIPs_Functions.R"))
} else {
  source(paste0(opt$code,"/smMIPs_Functions.R"))
}


################  LOAD LIBRARIES
library("data.table")
library("cellbaseR")
library("IRanges")


#Print all the parameters
cat("###############################\n")
cat("        Run Parameters\n")
cat("###############################\n")
print(opt)
cat("###############################\n")
cat("            Running...\n")
cat("###############################\n")


################  CREATING A TABLE WITH ALL THE POSSIBLE OBSERVED ALLELES USING THE TARGETED GENOMIC LOCI. ALTERNATIVELY, LOADING THE USER INPUT FILE
if(!is.null(opt$panel.file)){
  panel<-load.panel(opt$panel.file)
  panel$target_seq[which(panel$probe_strand=="-")]=reverse(chartr("ATGC","TACG",panel$target_seq[which(panel$probe_strand=="-")])) #determine the reference alleles on the positive strand

  tmp1=data.table("smMIP"=rep(panel$id,unlist(lapply(strsplit(panel$target_seq,""),function(x) length(x)))),
                  "chr"=rep(panel$chr,unlist(lapply(strsplit(panel$target_seq,""),function(x) length(x)))),
                  "pos"=unlist(apply(cbind(panel$target_start,panel$target_stop),1,function(x) seq(x[1],x[2],1))), #we are not interested to call mutations from the smMIPs arms
                  "ref"=unlist(strsplit(panel$target_seq, "")),
                  "alt"="A")
  tmp2=tmp1
  tmp2$alt="T"
  tmp3=tmp1
  tmp3$alt="C"
  tmp4=tmp1
  tmp4$alt="G"

  genomicBases=rbind(tmp1,tmp2,tmp3,tmp4)
  rm(tmp1,tmp2,tmp3,tmp4)
  genomicBases=genomicBases[alt!=ref]
  colnames(genomicBases)[3]="pos"
  idx=1:nrow(genomicBases)
} else if (!is.null(opt$input)){
  genomicBases=fread(opt$input,header=T,sep="\t")
  if(length(grep("annotated",names(genomicBases)))>0){
    idx=which(genomicBases$annotated=="FALSE")
  } else {
    idx=1:nrow(genomicBases)
  }
}

################  POPULATE THE TABLE

# Helper: create a single NA annotation row (used as fallback)
.make_na_annotation <- function() {
  data.frame(gene=NA, protein=NA, cosmic=NA, maf=NA,
             variant_type=NA, cadd_scaled=NA, annotated=FALSE,
             stringsAsFactors=FALSE)
}

annotation_cols <- c("gene","protein","cosmic","maf","variant_type","cadd_scaled","annotated")

if (isTRUE(opt$`skip-cellbase`)) {
  #--- Offline mode: skip CellBase entirely ---
  cat("WARNING: --skip-cellbase flag set. All annotation columns will be NA.\n")
  a <- lapply(idx, function(i) .make_na_annotation())

} else {
  #--- Online mode: attempt CellBase API with timeout and fallback ---

  # Set connection timeout for CellBase HTTP calls
  cellbase_timeout <- opt$`cellbase-timeout`
  old_timeout <- getOption("timeout")
  options(timeout = cellbase_timeout)

  # Try to initialize CellBase client; fall back to NA if it fails
  cb_init_ok <- tryCatch({
    cb <<- CellBaseR(species = opt$species)
    cbparam <<- CellBaseParam(assembly = opt$genome)
    # Quick connectivity check: annotate a dummy variant to verify API is reachable
    invisible(capture.output(
      getVariant(object=cb, ids="1:1:A:T", resource="annotation", param=cbparam)
    ))
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!cb_init_ok) {
    cat("WARNING: CellBase API is unreachable or timed out. Falling back to NA annotations.\n")
    cat("  (To suppress this warning, use --skip-cellbase)\n")
    a <- lapply(idx, function(i) .make_na_annotation())
  } else {
    print("Annotating. Please wait...", quote = F)

    a <- mclapply(idx, mc.cores = opt$threads, mc.cleanup=T, mc.silent=F, function(i) {
      tryCatch({
        variant_annotation_using_cellbaseR(i)
      }, error = function(e) {
        warning(sprintf("CellBase annotation failed for variant %d: %s", i, conditionMessage(e)))
        .make_na_annotation()
      })
    })

    # Check if all results are NA (bulk connection failure mid-run)
    n_failed <- sum(sapply(a, function(x) all(is.na(x$gene))))
    if (n_failed == length(a)) {
      cat("WARNING: All CellBase API calls failed. Output will contain NA annotations.\n")
    } else if (n_failed > 0) {
      cat(sprintf("WARNING: %d of %d CellBase API calls failed. Those variants have NA annotations.\n",
                  n_failed, length(a)))
    }
  }

  # Restore original timeout
  options(timeout = old_timeout)
}

x=sapply(rbindlist(a),class)
idx2=which(!(names(x) %in% names(genomicBases)))
if(length(idx2)>0){
  for(j in names(x)[idx2]){
    genomicBases[,j]=""
  }
}
genomicBases[] <- lapply(genomicBases, as.character)
genomicBases[idx,annotation_cols]=rbindlist(a)


###############  RENAME AMINO ACIDS WITH THEIR SINGLE LETTER ABBREVIATIONS
f=c("STOP","ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
r=c("*","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

for(i in 1:length(f)){
  genomicBases$protein=gsub(f[i],r[i],genomicBases$protein)
}

system(paste0("printf '\\rDownloading SNV annotations :  100%%     '"))
cat("\n")
cat("Writing the annotated table to disk\n")

u=unique(genomicBases$annotated)
if(length(u)==1 & "TRUE" %in% u){genomicBases=genomicBases[,-"annotated"]}

write.table(genomicBases,file=paste0(opt$output,"/annotated_",file.name),col.names = T,row.names = F,quote = F,sep = '\t')

cat("\n")
cat("###############################\n")
cat("             DONE\n")
cat("###############################\n")

