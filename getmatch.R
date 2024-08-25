source("/home/rdc143/hapla-proj/sim/funcs.R")

basedir <- commandArgs(trailingOnly = TRUE)[3]

### numbers of haplotypes
n_ref <- nrow(read.table(paste0(basedir, "/ref.lst"), header = FALSE)) * 2
n_query <- nrow(read.table(paste0(basedir, "/query.lst"), header = FALSE)) * 2

### read true path
true_anc <- read.table(paste0(basedir, "/simulated_true_anc.csv"), header = TRUE)
true_anc["sample"] <- true_anc["sample"] + 1 # change indexing from 0 to 1

names(true_anc)[2] <- "pos"

seqlen <- max(true_anc$right)
npop <- length(unique(true_anc$population))

### read hapla path res
hapla <- read.table(commandArgs(trailingOnly=TRUE)[1], header = FALSE)
names(hapla) <- read.table(commandArgs(trailingOnly=TRUE)[2])[, 2]

### remap hapla population numbering based on order in ref inds
### as they are assigned arbitrarily at admixture, should now work for unequal ref no. and admixed refs
remap <- as.character(0:(npop - 1))
n_ref_per <- rle(read.table(paste0(basedir, "/samplemap.tsv"), header = FALSE)[, 2])$lengths * 2
n <- 1
for (i in 1:npop){ # findes mode of ancestry in first 100 wd's of each ref panel
    names(remap)[i] <- names(sort(-table(unlist(hapla[n:(n + n_ref_per[i] - 1), 1:100]))))[1]
    n <- n + n_ref_per[i]
}
hapla[] <- lapply(hapla, function(x) {
    as.numeric(unname(remap[as.character(x)])) })


### read and reformat rfmix path res
rf <- read.table(paste0(basedir, "/rfmix.fb.tsv"), header = TRUE)
rf_wide <- rf[5:ncol(rf)]
rf_max <- matrix(NA, nrow(rf), ncol(rf_wide) / npop)
for (i in 1:(ncol(rf_wide) / npop)) {
    chunk <- rf_wide[, ((i - 1) * npop + 1):(i * npop)]
    rf_max[, i] <- max.col(chunk, ties.method = "first") - 1
}

rf <- cbind(rf[1:4], rf_max)
names(rf)[5:(n_query + 4)] <- (n_ref + 1):(n_ref + n_query)


# get match proportions
matchdf <- data.frame()
for (haplo_index in (1 + n_ref):(n_ref + n_query)){
    match_props <- get_match_props(haplo_index)
    for (method in c("hapla", "rfmix")) {
        for (pop in names(match_props[[method]])){
            matchdf <- rbind(matchdf, list("haplotype" = haplo_index, "method" = method, "pop" = pop, "match" = match_props[[method]][[pop]][1], "true" = match_props[[method]][[pop]][2]))
        }
    }
}

write.table(sum(matchdf[matchdf$method == "hapla", "match"]) / sum(matchdf[matchdf$method == "hapla", "true"]),
            "meanhaplamatch.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)
