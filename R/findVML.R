# Map revmap column to probe names after reducing a GenomicRanges object

# Given a revmap row (e.g. 1 5 6), we map those positions to their
# corresponding probe names (and end up with something like "cg00000029",
# "cg00000158", "cg00000165".This is a helper function of findVML()).

# @param positions A revmap row in the form of a vector
# @param manifest_hvp the manifest of the highly variable probes used in the
# findVML() function with the probes as row names

# @return a vector with the names of the probes that conform one reduced region

# @examples
# target <- data.frame(row.names = c("a", "b", "c", "d"),
#                      values = c(1, 1, 1, 1))
# query <- c(2, 1)
#
# map_revmap_names(positions = query, manifest_hvp = target)
# ## Expected output: c("b", "a")
# }
map_revmap_names <- function(positions, manifest_hvp) {
  # We start with 1 5 6
  # We want to end with cg00000029,  cg00000158 cg00000165
  row.names(manifest_hvp)[positions]
}

#' Resolve the array manifest
#'
#' Turns the name of a supported Illumina array into its manifest, taken from
#' the corresponding Bioconductor annotation package. A manifest supplied
#' directly as a data frame is returned unchanged.
#'
#' @param array_manifest Either the name of a supported array or a data frame.
#'
#' @return A data frame with the probe locations, with the probes as row names.
#' @noRd
.resolveManifest <- function(array_manifest) {
  if (!is.character(array_manifest)) return(array_manifest)
  if (array_manifest == "IlluminaHumanMethylation450k") {
    if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                          quietly = TRUE)) {
      stop("Package \"IlluminaHumanMethylation450kanno.ilmn12.hg19\"  ",
           "must be installed to use this function.",
           call. = FALSE)
    }
    data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)
  } else if (array_manifest == "IlluminaHumanMethylationEPICv1") {
    if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                          quietly = TRUE)) {
      stop("Package 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19' ",
           "must be installed to use this function.",
           call. = FALSE)
    }
    data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations)
  } else if (array_manifest == "IlluminaHumanMethylationEPICv2") {
    if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                          quietly = TRUE)) {
      stop("Package 'IlluminaHumanMethylationEPICv2anno.20a1.hg38' ",
           "must be installed to use this function.",
           call. = FALSE)
    }
    data.frame(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations)
  }
}

#' Compute the variability threshold used to call Highly Variable Probes
#'
#' The threshold is a percentile of the variability scores. Depending on
#' `var_distribution` the percentile is taken over every probe in the data set,
#' or only over the ultrastable probes, which are the ones expected to carry
#' little biological variability.
#'
#' @param var_scores Data frame of variability scores, probes as row names.
#' @param var_distribution Either "all" or "ultrastable".
#' @param var_threshold_percentile Percentile used as the threshold.
#' @param array_manifest The manifest, or the name of a supported array, used
#' to resolve the ultrastable probe names for the array at hand.
#'
#' @return The variability score threshold.
#' @noRd
.varThreshold <- function(var_scores, var_distribution,
                          var_threshold_percentile, array_manifest) {
  Methyl450_Loci <- NULL # Binding for the tidyverse framework
  if (var_distribution == "all") {
    return(stats::quantile(var_scores$var_score, var_threshold_percentile))
  }
  ## var_distribution == "ultrastable": resolve the ultrastable probe names in
  ## the naming scheme of the array at hand
  if (is.data.frame(array_manifest)) {
    # A custom manifest is expected to use the 450k/EPICv1 probe names. The
    # number of ultrastable probes present was already checked further up.
    return(stats::quantile(
      var_scores[(row.names(var_scores) %in% RAMEN::ultrastable_cpgs), ],
      var_threshold_percentile
    ))
  }
  if (array_manifest == "IlluminaHumanMethylationEPICv2") {
    # The EPICv2 manifest renames the probes, so the ultrastable ones have to
    # be looked up through their 450k identifiers
    if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38",
                          quietly = TRUE)) {
      stop("Package 'IlluminaHumanMethylationEPICv2anno.20a1.hg38' must ",
           "be installed to use this function.",
           call. = FALSE)
    }
    ultrastable <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other |>
      data.frame() |>
      dplyr::filter(Methyl450_Loci %in% RAMEN::ultrastable_cpgs) |>
      rownames()
  } else {
    # EPICv1 and 450k use the same probe names as the ultrastable probes
    ultrastable <- RAMEN::ultrastable_cpgs
  }
  if (sum(row.names(var_scores) %in% ultrastable) < 25) {
    stop("Your data set should contain more than 25 ultrastable probes ",
         "(RAMEN::ultrastable_cpgs) to compute a variance threshold. If ",
         "not, please get the variability threshold based on all the ",
         "probes in your data set (var_distribution = 'all', ",
         "var_threshold_percentile = 0.9).")
  }
  var_threshold <- stats::quantile(
    var_scores[(row.names(var_scores) %in% ultrastable), ],
    var_threshold_percentile
  )
  var_threshold
}

#' Group probes into regions
#'
#' Builds a GenomicRanges object out of a probe manifest and reduces it, so
#' that probes closer to each other than `max_distance` end up in the same
#' region. Each region carries the number of probes it contains and their
#' names.
#'
#' @param probe_df Manifest with the columns "TargetID", "chr", "pos" and
#' "strand", sorted by chromosome.
#' @param max_distance Maximum distance between two probes of a region.
#' @param count_name Name given to the metadata column holding the number of
#' probes per region.
#'
#' @return A GRanges object with one range per region.
#' @noRd
.probeToRegions <- function(probe_df, max_distance, count_name) {
  seqnames_table <- table(probe_df$chr)
  gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(names(seqnames_table),
                              # Number of chromosome; as.numeric to convert
                              # from table to numeric vector
                              as.numeric(seqnames_table)),
    ranges = IRanges::IRanges(probe_df$pos,
                              end = probe_df$pos,
                              names = probe_df$TargetID
    ),
    # Rle() run length encodes the vector itself, so the strands only have to
    # be walked once rather than once for the values and again for the lengths
    strand = S4Vectors::Rle(as.character(probe_df$strand))
  )
  regions <- GenomicRanges::reduce(gr,
                                   with.revmap = TRUE,
                                   min.gapwidth = max_distance)
  # Add the number of probes in each region
  S4Vectors::mcols(regions)[[count_name]] <- vapply(
    S4Vectors::mcols(regions)$revmap,
    length,
    FUN.VALUE = numeric(1)
  )
  # Substitute revmap with the name of the probes in each region
  S4Vectors::mcols(regions)$probes <- lapply(S4Vectors::mcols(regions)$revmap,
                                             map_revmap_names,
                                             probe_df)
  # Remove revmap mcol
  S4Vectors::mcols(regions)$revmap <- NULL
  regions
}

#' Capture the canonical Variable Methylated Regions
#'
#' Keeps the candidate regions that hold more than one probe and, of those, the
#' ones whose probes are correlated with each other.
#'
#' @param candidate_VMRs Candidate regions, as returned by `.probeToRegions()`.
#' @param methylation_data The methylation data, used to compute correlations.
#' @param cor_threshold Minimum median pairwise correlation of a region.
#'
#' @return A GRanges object with the VMRs.
#' @noRd
.findVMRs <- function(candidate_VMRs, methylation_data, cor_threshold) {
  VMRs <- candidate_VMRs[
    (GenomicRanges::elementMetadata(candidate_VMRs)[, "n_VMPs"] > 1)]
  ### Check for correlation between probes only if we have VMRs
  if (length(VMRs) > 0) {
    # Compute the median correlation of each region
    VMRs <- medCorVMR(VML = VMRs, methylation_data = methylation_data)
    # Remove VMRs whose CpGs are not correlated
    VMRs <- VMRs[
      (GenomicRanges::elementMetadata(VMRs)[, "median_correlation"] >
         cor_threshold)]
  } else {
    warning("No VMRs were found in this data set")
  }
  VMRs
}

#' Capture the sparse Variable Methylated Probes
#'
#' A candidate region is a sparse Variable Methylated Probe when it holds a
#' single highly variable probe and that probe has no neighbour within
#' `max_distance` on the array, i.e. it is one of the lonely probes.
#'
#' @param candidate_VMRs Candidate regions, as returned by `.probeToRegions()`.
#' @param lonely_probes Names of the probes that sit alone in their region.
#'
#' @return A GRanges object with the sVMPs.
#' @noRd
.findSVMPs <- function(candidate_VMRs, lonely_probes) {
  probes_mcol <- GenomicRanges::elementMetadata(candidate_VMRs)[, "probes"]
  single_probe_region <- lengths(probes_mcol) == 1
  is_sVMP <- logical(length(candidate_VMRs))
  is_sVMP[single_probe_region] <-
    unlist(probes_mcol[single_probe_region]) %in% lonely_probes
  sVMPs <- candidate_VMRs[is_sVMP]
  # Add a column of NAs under the name of median_correlation to match the VMRs
  S4Vectors::mcols(sVMPs)$median_correlation <-
    rep(NA, nrow(S4Vectors::mcols(sVMPs)))
  sVMPs
}




#' Identify Variable Methylated Loci in microarrays
#'
#' Identifies Highly Variable Probes (HVP) and groups them into Variable
#' Methylated Loci (VML) given an Illumina manifest.The output of this function
#' provides the HVPs, and the identified VML, which are made of Variable
#' Methylated Regions and sparse Variable Methylated Probes. See Details below
#' for more information.
#'
#' This function identifies HVPs based on variance or MAD scores, and groups
#' them into VML, which are defined as genomic locations with high DNA methylation
#' variability.To best capture methylome variability patterns in microarrays,
#' we identify two types of VML: Variably Methylated Regions (VMRs) and sparse
#' Variably Methylated Probes (sVMPs) .
#'
#' In one hand, we defined VMRs as two or more proximal highly variable probes
#' (default: < 1kb apart) with correlated DNAme level (default: r > 0.15).
#' Modelling DNAme variability through regions rather than individual CpGs
#' provides several methodological advantages in association studies, since
#' CpGs display a significant correlation for co-methylation when they are
#' close (less than or equal to 1 kilobase). Modelling DNAme variability
#' through regions rather than individual CpGs provides several methodological
#' advantages in association studies, since CpGs display a significant
#' correlation for co-methylation when they are close (less than or equal to 1
#' kilobase). Modelling DNAme variability through regions rather than individual
#' CpGs provides several methodological advantages in association studies, since
#' CpGs display a significant correlation for co-methylation when they are close
#' (less than 1 kilobase). Some of these advantages include increasing statistical power
#' by testing redundant probes only once, reducing false-positives driven by one
#' problematic probe in a region, and improving comparability between studies
#' that analyze the same genomic region but measure distinct CpGs due to
#' microarray design differences. In contrast with Differentially Methylated
#' Regions (DMRs), a different class of regional construct used in the field,
#' VMRs denote regions with high inter-individual variability in methylation
#' levels within a single population, while DMRs represent regions where DNA
#' methylation differs significantly across a variable of interest.
#'
#' In addition to traditional VMRs, we also identified sparse Variably
#' Methylated Probes (sVMPs), a second type of VML that takes into account the
#' sparse and non-uniformly distributed coverage of CpGs in microarrays to
#' tailor our analysis to this DNAme platform. sVMPs aimed to retain genomic
#' regions with high DNAme variability measured by single probes, where probe
#' grouping based on proximity and correlation is therefore not applicable. This
#' is particularly relevant in the Illumina EPIC v1 array, where most covered
#' regulatory regions (up to 93%) are represented by just one probe. Notably,
#' based on empirical comparisons with whole-genome bisulfite sequencing data,
#' these single probes are mostly representative of local regional DNAme
#' levels due to their positioning (98.5-99.5%)
#'
#' This function uses GenomicRanges::reduce() to group the regions, which is
#' strand-sensitive. In the Illumina microarrays, the MAPINFO for all the probes
#' is usually provided for the + strand. If you are using this array, we
#' recommend to first convert the strand of all the probes to "+".
#'
#' This function supports parallel computing for increased speed. To do so, you
#' have to set the parallel backend
#' in your R session BEFORE running the function (e.g.,
#' *doParallel::registerDoParallel(4)*). After that, the function can be run as
#' usual. When working with big datasets, the parallel backend might throw an
#' error if you exceed the maximum allowed size of globals exported for future
#' expression. This can be fixed by increasing the allowed size (e.g. running
#' *options(future.globals.maxSize= +Inf)*)
#'
#' Note: this function does not exclude sex chromosomes. If you want to exclude
#' them, you can do so in the methylation_data object before running the
#' function.
#'
#' @param array_manifest Information about the probes on the array in a format
#' compatible with the Bioconductor annotation packages. The user can specify
#' one of the supported human microarrays ("IlluminaHumanMethylation450k" with
#' the hg19 genome build, "IlluminaHumanMethylationEPICv1" with the hg19 genome
#' build, or "IlluminaHumanMethylationEPICv2" with the hg38 genome build), or
#' provide a manifest. The manifest requires the probe names as row names, and
#' the following columns: "chr" (chromosome); "pos" (genomic location of the
#' probe in the genome); and "strand" (this is very important to set up, since
#' the VMRs will only be created based on CpGs on the same strand; if the
#' positions are reported based on a single DNA strand, this should contain
#' either a vector of only "+", "-" or "*" for all of the probes).
#' @param methylation_data A data frame containing M or B values, with samples
#' as columns and probes as rows. Data is expected to have already passed
#' through quality control and cleaning steps. Rows must be the CpG probe IDs.
#' @param cor_threshold Numeric value (0-1) to be used as the median pearson
#' correlation threshold for identifying VMRs (i.e. all VMRs will have a median
#' pairwise probe correlation higher than this threshold).
#' @param var_method A string indicating the metric to use to represent
#' variability in the data set. The options are "mad" (median absolute
#' deviation) or "variance".
#' @param var_distribution A string indicating which probes in the data set
#' should be used to create a variability distribution; the threshold to
#' identify Highly Variable Probes (determined also with the
#' var_threshold_percentile argument) is established based on this distribution.
#' The options 1 is "ultrastable" (a subset of CpGs that are stably
#' methylated/unmethylated across human tissues and developmental states
#' described by [Edgar R., et al.](https://doi.org/10.1186/1756-8935-7-28)
#' in 2014). This option is recommended, especially if you want to compare
#' different populations or tissues, as the threshold value should be
#' comparable. On the other hand, the user can use option 2: "all" (all probes
#' in the data set). The "ultrastable" option is only compatible with Illumina
#' human microarrays. The default is "ultrastable".
#' @param var_threshold_percentile The percentile (0-1) to be used as cutoff to
#' define Highly Variable Probes (which are then grouped into VML). If using
#' the variability of the "ultrastable" probes, we recommend a high threshold
#' (default is 0.99), since these probes are expected to display a very low
#' variation in human tissues. If using the variability of "all" probes, we
#' recommend using a percentile of 0.9 since it captures the top 10% most
#' variable probes, which has been traditionally used in studies. It is
#' important to note that the top 10% most variable probes will capture the
#' same amount of probes in a data set regardless of their overall variability
#' levels, which might differ between tissues or populations.
#' @param max_distance Maximum distance in base pairs allowed for two probes to
#' be grouped into a region. The default is 1000.
#'
#' @return A list with the following elements:
#'  - $var_score_threshold: threshold used to define Highly Variable Probes (mad
#'   or variance, depending on the specified choice).
#'  - $highly_variable_probes: a data frame with the probes that passed the
#'  variability score threshold imposed by the user, and their variability score
#'   (MAD score or variance).
#'  - $VML: a GRanges object listing VMRs (regions composed of two or more
#'   contiguous, correlated and proximal Highly Variable Probes), and sVMPs
#'   (highly variable probes without neighboring CpGs measured in *max_distance*
#'   on the array).
#'
#' @export
#' @examples
#' # Evaluate sequentially
#' foreach::registerDoSEQ()
#' VML <- RAMEN::findVML(
#'   methylation_data = RAMEN::test_methylation_data,
#'   array_manifest = "IlluminaHumanMethylationEPICv1",
#'   cor_threshold = 0.15,
#'   var_method = "variance",
#'   var_distribution = "ultrastable",
#'   var_threshold_percentile = 0.99,
#'   max_distance = 1000
#' )
findVML <- function(methylation_data,
                    array_manifest,
                    cor_threshold = 0.15,
                    var_method = "variance",
                    var_distribution = "ultrastable",
                    var_threshold_percentile = 0.99,
                    max_distance = 1000) {
  #### Binding of variables used within the tidyverse framework ####
  Methyl450_Loci <- TargetID <- chr <- pos <- strand <- NULL
  var_score <- NULL
  #### Argument checks ####
  # methylation_data
  argument_check(methylation_data, "data.frame")
  # array_manifest
  if (is.data.frame(array_manifest)) {
    columns_exist(array_manifest,
                  c("chr", "pos", "strand"),
                  paste("Otherwise, provide a string with one of the supported",
                        "human microarrays ('IlluminaHumanMethylation450k',",
                        "'IlluminaHumanMethylationEPICv1', or",
                        "'IlluminaHumanMethylationEPICv2')."))
    # Check that the array strand is in the format expected by the user
    if (base::length(base::unique(array_manifest$strand)) > 1) {
      warning("The manifest currently has more than one type of strands. ",
              "Please note that this function is strand sensitive. So, ",
              "probes in proximal coordinates but different strands on the ",
              "manifest will not be grouped together. Many array manifests ",
              "such as the Illumina EPIC v1 include the PROBE strand, but ",
              "the position of the actual CpGs (pos) is reported in the same ",
              "strand; in those cases we recommend setting all of the probes ",
              "to the same strand.")
    }
    if (var_distribution == "ultrastable") {
      # If the user provides their own manifest and is choosing to use the
      # ultrastable probes, make sure that a good number of them is present
      # in the data set. If not, throw an error
      if (sum(row.names(array_manifest) %in% RAMEN::ultrastable_cpgs) < 25) {
        stop("The var_distribution = 'ultrastable' option is only ",
             "compatible with Illumina human microarrays at the moment. If ",
             "you are using a human Illumina microarray please indicate it ",
             "with their corresponding string, or make sure that it contains ",
             "a more than 25 ultrastable probes (RAMEN::ultrastable_cpgs). ",
             "If not, please get the variability threshold based on all the ",
             "probes in your data set (var_distribution = 'all', ",
             "var_threshold_percentile = 0.9).")
      }
    }
  } else argument_char_options(array_manifest,
                               options = c("IlluminaHumanMethylation450k",
                                           "IlluminaHumanMethylationEPICv1",
                                           "IlluminaHumanMethylationEPICv2"),
                               extra_msg = paste("Otherwise, please provide a",
                                                 "custom manifest as a data frame."))
  # cor_threshold
  argument_check(cor_threshold, "numeric")
  if (!(cor_threshold >= 0 && cor_threshold <= 1)) {
    stop("'cor_threshold' must be a value between 0 and 1 (inclusive)")
  }
  # var_method
  argument_char_options(var_method, c("variance", "mad"))
  # var_distribution
  argument_char_options(var_distribution, c("ultrastable", "all"))
  # var_threshold_percentile
  argument_check(var_threshold_percentile, "numeric")
  if (!(var_threshold_percentile >= 0 && var_threshold_percentile <= 1)) {
    stop("'var_threshold_percentile' must be a value between 0 and 1 (inclusive)")
  }
  # max_distance
  argument_check(max_distance, "numeric")

  #### Compute variation metric ####
  if (var_method == "mad") {
    var_scores <- matrixStats::rowMads(as.matrix(methylation_data)) |>
      as.data.frame() |>
      dplyr::rename("var_score" = 1)
  } else if (var_method == "variance") {
    var_scores <- matrixStats::rowVars(as.matrix(methylation_data)) |>
      as.data.frame() |>
      dplyr::rename("var_score" = 1)
  }

  #### Get HVP threshold ####
  message("Identifying Highly Variable Probes...")
  # Get the variability threshold
  var_threshold <- .varThreshold(var_scores, var_distribution,
                                 var_threshold_percentile, array_manifest)
  # Replace the array manifest if the user provided a string with the name of
  # the array
  manifest <- .resolveManifest(array_manifest)
  #### Get HVPs ####
  full_manifest <- cbind(TargetID = rownames(manifest), manifest) |>
    # Select necessary columns in case there is more
    dplyr::select(c(TargetID, chr, pos, strand)) |>
    dplyr::filter(
      # Remove probes with no map info
      !is.na(pos),
      # keep only the probes where we have methylation information
      TargetID %in% row.names(var_scores)
    ) |>
    # important step for using Rle later when constructing the GenomicRanges
    # object
    dplyr::arrange(chr)
  if (is.factor(full_manifest$chr)) {
    full_manifest <- full_manifest |>
      dplyr::mutate(chr = droplevels(chr))
  }

  # Keep the probes that are highly variable
  manifest_hvp <- full_manifest |>
    dplyr::filter(
      TargetID %in% row.names(var_scores |>
        dplyr::filter(var_score >= var_threshold))
    ) |>
    dplyr::left_join(
      # Add variability information
      cbind(TargetID = rownames(var_scores), var_scores),
      by = "TargetID"
    )
  rownames(manifest_hvp) <- manifest_hvp$TargetID
  if (is.factor(manifest_hvp$chr)) {
    manifest_hvp <- manifest_hvp |>
      dplyr::mutate(chr = droplevels(chr))
  }

  #### Identify sparse Variable Methylated Probes####
  message("Identifying sparse Variable Methylated Probes")
  # Convert the full manifest to a GenomicRanges object
  # Group every probe into regions, so that the probes which sit alone can
  # be told apart from those that have neighbours
  regions_full_manifest <- .probeToRegions(full_manifest, max_distance,
                                         "n_probes")
  # Keep elements with only one probe
  lonely_probes <- regions_full_manifest[
    (GenomicRanges::elementMetadata(regions_full_manifest)[, "n_probes"] <= 1)]
  lonely_probes <- S4Vectors::mcols(lonely_probes)$probes |>
    unlist()

  #### Identify VMRs####
  message("Identifying Variable Methylated Regions...")
  # convert the highly variable probes data frame to a GenomicRanges object
  # Group the highly variable probes into candidate regions
  candidate_VMRs <- .probeToRegions(manifest_hvp, max_distance, "n_VMPs")

  ### Capture canonical VMRs ###
  message("Applying correlation filter to Variable Methylated Regions...")
  VMRs <- .findVMRs(candidate_VMRs, methylation_data, cor_threshold)

  ### Capture sVMPs ###
  sVMPs <- .findSVMPs(candidate_VMRs, lonely_probes)

  ### Make final GRanges object ####
  VML <- c(VMRs, sVMPs)
  S4Vectors::mcols(VML)$type <- c(rep("VMR", length(VMRs)),
                       rep("sVMP", length(sVMPs)))
  S4Vectors::mcols(VML)$VML_index <- paste0("VML", seq_along(VML))

  return(list(
    var_score_threshold = var_threshold,
    highly_variable_probes = cbind(TargetID = rownames(var_scores), var_scores) |>
      dplyr::filter(TargetID %in% manifest_hvp$TargetID),
    VML = VML
  ))
}
