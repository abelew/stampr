## ----options, include=FALSE---------------------------------------------------
library("stampr")
basedir <- system.file("extdata", package="stampr")
knitr::opts_knit$set(width=120,
                     root.dir=basedir,
                     progress=TRUE,
                     verbose=TRUE,
                     echo=TRUE)
knitr::opts_chunk$set(error=TRUE,
                      dpi=96)
old_options <- options(digits=4,
                       stringsAsFactors=FALSE,
                       knitr.duplicate.label="allow")
ggplot2::theme_set(ggplot2::theme_bw(base_size=10))
rundate <- format(Sys.Date(), format="%Y%m%d")
##tmp <- sm(loadme(filename=paste0(gsub(pattern="\\.Rmd", replace="", x=previous_file), "-v", ver, ".rda.xz")))
##rmd_file <- "03_expression_infection_20180822.Rmd"
library(ggplot2)

## ----installation, eval=FALSE-------------------------------------------------
#  devtools::install_github("abelew/stampr")

## ----sample_sheet, results='asis'---------------------------------------------
basedir <- system.file("extdata", package="stampr")
metadata_file <- file.path(basedir, "all_samples.xlsx")
meta <- hpgltools::extract_metadata(metadata_file)
knitr::kable(meta)

## ----extract_reads------------------------------------------------------------
untarred <- utils::untar(tarfile=system.file("extdata/trimmed_reads.tar", package="stampr"))

## ----extract_otus-------------------------------------------------------------
untarred <- utils::untar(tarfile=system.file("extdata/otus.tar", package="stampr"))

## ----extract_counts-----------------------------------------------------------
untarred <- utils::untar(tarfile=system.file("extdata/tag_counts.tar", package="stampr"))

## ----tally, eval=FALSE--------------------------------------------------------
#  ## This R script will take the output files containing representative STAMP barcode sequences and
#  ## their individual abundance in the sequencer output from all samples  that
#  ## you sequenced and copies them in a single .csv file in the format that is
#  ## expected from the Nb estimation script
#  
#  output_file <- "qiime_tally.csv"
#  ## files <- list.files(txt_dir)
#  files <- meta[["qiimeotus"]]
#  names <- files
#  v.length <- rep(NA, length(names))
#  for (f in 1:length(files)) {
#    file <- files[f]
#    ##mtrx <- read.table(file, sep="\t", header=FALSE, fill=TRUE)
#    ##colnames(mtrx) <- c(paste0(files[f], "_seq"), files[f])
#    eval(parse(text=paste("Data=as.matrix(read.table('",file,"', fill=TRUE, sep='\t'))",sep="")))
#    colnames(Data) <- c(paste(names[f],"_seq", sep=""), names[f])
#    name <- paste("Matrix", names[f], sep="")
#    assign(name, Data)
#    name <- paste("Matrix", names[f], sep="")
#    assign(name, Data)
#    v.length[f] <- nrow(Data)
#  }
#  max <- max(v.length) + 1
#  fill.up <- max - v.length
#  
#  FillUpMatrix <- matrix(data=NA, nrow=fill.up[1], ncol=2)
#  eval(parse(text=paste("OutputMatrix=rbind(Matrix", names[1],", FillUpMatrix)", sep="")))
#  for(ff in 2:length(files)) {
#    eval(parse(text=paste("MatrixX=Matrix", names[ff], sep="")))
#    FillUpMatrix <- matrix(data=NA, nrow=fill.up[ff], ncol=2)
#    Matrix <- rbind(MatrixX, FillUpMatrix)
#    OutputMarix <- cbind(OutputMatrix, Matrix)
#  }
#  write.csv(OutputMatrix, file=outputfilename)

## ----read_qiime_otus----------------------------------------------------------
qiime_tags <- read_qiime_otus(meta, otu_column="qiimeotus", trimmed_column="trimmedfastq")
head(qiime_tags[["modified_metadata"]])
meta <- qiime_tags[["modified_metadata"]]

## Number of reads / tag
plot_read_density(qiime_tags)

## ----read_perl_tags-----------------------------------------------------------
perl_tags <- read_tags(meta, index_column="indextable")
summary(perl_tags)
## read_idx returns a list with 3 elements:
##  the input metadata with some new columns.
##  the counts in a long format for easier plotting.
##  the counts in a wide format, which is easier to work with.
head(perl_tags[["modified_metadata"]])
## The new elements in the metadata include:
##  reads:  how many reads were used in this counting.
##  observed: how many tags were observed.
##  passed_cutoff: how many tags passed our simple cutoff for 'real'.
##  rejected_cutoff: passed + rejected = observed
##  max: What is the maximum number of reads for 1 tag.  Thus for a1, we see that one tag sucked up
##       almost half of all the reads.

## Given the above, we can trivially plot how many reads/tag there are for each sample:
plot_read_density(perl_tags)

## ----original_qiime_tally-----------------------------------------------------
original_tag_file <- system.file("extdata/original_qiime_tally.csv", package="stampr")
original_tags <- read.csv(original_tag_file)
head(original_tags)

## ----filter_tags--------------------------------------------------------------
new_meta <- perl_tags[["modified_metadata"]]
current_tags <- perl_tags[["sample_counts"]]
filtered_tags <- filter_tags(perl_tags)
summary(filtered_tags)
head(filtered_tags[["modified_metadata"]])
## Thus the new column 'passed_aggressive_cpm'

## ----quantify_nb_original, eval=FALSE-----------------------------------------
#  inputdata <- "../../original/BigMatAll_STAMP_CAUTI_calcurve1_working.csv"
#  ### read input sequence data file
#  ## Data=as.matrix(read.csv(inputdata,row.names=1))
#  Data <- read.csv(inputdata)
#  rownames(Data) <- make.names(Data[["X"]], unique=TRUE)
#  Data[["X"]] <- NULL
#  
#  #BigMatAll=DataII[,2:ncol(DataII)]
#  BigMatAll <- Data
#  
#  ### USER ### define output filename
#  name <- "New_EstimateMatrix_CAUTI_Calcurve1"
#  
#  ### shows colnames and col # to make it easier to identify the reference samples
#  ### and the data samples
#  cbind(1:ncol(BigMatAll), colnames(BigMatAll))
#  
#  ### should a calibration curve be used to normalize the data? TRUE/FALSE
#  usecalibrationcurve <- FALSE
#  
#  ### USER ### If calibration curve is used, give path to calibration data
#  ### ("CalibrationLookup")
#  if (usecalibrationcurve) {
#    load("/Users/kerbachta/Desktop/STAMP/calcurve3/CalibrationLookup.Rdata")
#  }
#  
#  ### USER ### define which columns of BigMatAll contain real read numbers and not
#  ### meta-data like binaries, average count number or absolute presence in
#  ### samples. Usually column 1 to 24 (1:24) contain real read numbers.
#  v.reads <- 1:23
#  
#  ### USER ### define which columns of BigMatAll contain reference samples
#  ### (usually the inoculum)? Define at least two columns e.g. 1:2 or c(1,2) would
#  ### mean column 1 and column 2 of BigMatAll contain reference samples.
#  ## v.reference <- 21:23
#  v.reference <- c(8, 15, 23)
#  
#  ### USER ### define which columns of BigMatAll contain data samples (that should
#  ### be the columns of BigMatAll that contain real read numbers minus the
#  ### reference columns)
#  ## v.bottlenecksamples <- 1:20
#  v.bottlenecksamples <- c(1:7, 9:14, 16:22)
#  
#  ### define threshold for inocula (e.g. sequence must be present at least once)
#  v.threshold <- which(rowSums(BigMatAll[, v.reference], na.rm=TRUE) > 0)
#  
#  ### choose real reads
#  ThresholdMatrix <- BigMatAll[v.threshold, v.reads]
#  
#  ## Why was this not done first?
#  ### NA= 0 (not present)
#  ThresholdMatrix[is.na(ThresholdMatrix) == TRUE] <- 0
#  
#  ### Create matrix with references only
#  ReferenceMatrix <- ThresholdMatrix[, v.reference]
#  head(ReferenceMatrix)
#  
#  ### Define FrequenzMatrix
#  fReferenceMatrix <- matrix(nrow=nrow(ReferenceMatrix),
#                             ncol=ncol(ReferenceMatrix),
#                             dimnames=list(rownames(ReferenceMatrix),
#                                           colnames(ReferenceMatrix)))
#  
#  ### Fill FrequenzMatrix with frequencies
#  for(x in 1:ncol(fReferenceMatrix)) {
#  	fReferenceMatrix[, x] <- ReferenceMatrix[, x] / sum(ReferenceMatrix[, x])
#  }
#  
#  ### Average frequencies
#  fTotalInoculum <- rowMeans(fReferenceMatrix, na.rm=TRUE)
#  
#  ### make matrix that connects bottleneckestimate to columnname
#  EstimateMatrix <- matrix(nrow=5, ncol=length(v.bottlenecksamples))
#  colnames(EstimateMatrix) <- colnames(BigMatAll[, v.bottlenecksamples])
#  rownames(EstimateMatrix) <- c("Nb","lower CI","Nb'_median","upper CI","# sequences")
#  
#  ### define expected counts and frequencies
#  expectedCounts <- rowMeans(ReferenceMatrix, na.rm=TRUE)
#  expvec <- fTotalInoculum[expectedCounts > 0]
#  inoculumSize <- sum(ReferenceMatrix, na.rm=TRUE)
#  
#  ### loop for estimating Ne for each sample
#  for (i in 1:length(v.bottlenecksamples)) {
#    bnvec <- ThresholdMatrix[, v.bottlenecksamples[i]]
#    bnvec <- bnvec[expectedCounts > 0]
#    sampleSize <- sum(bnvec)
#    binomialNenner <- expvec * (1 - expvec)
#    measured <- bnvec / sum(bnvec)
#    var <- (measured - expvec) ^ 2
#    binomial <- as.numeric(var) / as.numeric(binomialNenner)
#    BSize <- 1 / (mean(binomial) - 1 / sampleSize - 1 / inoculumSize)
#    EstimateMatrix[1, i] <- BSize
#    Log10BSize <- round(log10(BSize), 2)
#    EstimateMatrix[5, i] <- sampleSize
#  ### this part of the code compares to the lookupmatrix
#    if (usecalibrationcurve) {
#      if (Log10BSize < max(as.numeric(rownames(LookUpMatrix)))) {
#        if (usecalibrationcurve) {
#          EstimateMatrix[2:4, i] <- 10 ^ LookUpMatrix[as.character(Log10BSize), ]
#        }
#      }
#    }
#  }
#  
#  head(EstimateMatrix)
#  tmp <- t(EstimateMatrix)
#  meta[["Nb_original"]] <- 0
#  for (l in 1:nrow(tmp)) {
#    row <- tmp[l, ]
#    name <- rownames(tmp)[l]
#    name <- tolower(gsub(x=name, pattern="\\.txt", replacement=""))
#    meta[name, "Nb_original"] <- row[["Nb"]]
#  }

## ----my_nb_idx----------------------------------------------------------------
my_nb <- calculate_nb(perl_tags, metadata_column="MyNb")

## ----nb_filtered--------------------------------------------------------------
topn_filtered <- filter_topn_tags(my_nb, topn=600)
topn600_nb <- calculate_nb(topn_filtered, metadata_column="Nb600")
top500_filtered <- filter_topn_tags(topn600_nb, topn=500)
topn500_nb <- calculate_nb(top500_filtered, metadata_column="Nb500", write="modified_metadata-20210105.xlsx")

## ----plot_calibration---------------------------------------------------------
cali_curve <- plot_calibration(topn500_nb)
cali_curve$lm
hpgltools::pp(file="conservative_calibration_curve.png", image=cali_curve$plot)

