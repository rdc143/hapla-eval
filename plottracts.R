library(ggplot2)
source("/home/rdc143/hapla-proj/sim/funcs.R")

### numbers of haplotypes
n_ref <- nrow(read.table("ref.lst", header = FALSE)) * 2
n_query <- nrow(read.table("query.lst", header = FALSE)) * 2


### read true path
true_anc <- read.table("simulated_true_anc.csv", header = TRUE)
true_anc["sample"] <- true_anc["sample"] + 1 # change indexing from 0 to 1

names(true_anc)[2] <- "pos"

seqlen <- max(true_anc$right)
npop <- length(unique(true_anc$population))


### read hapla path res
hapla <- read.table("hapla.path", header = FALSE)
names(hapla) <- read.table("hapla.simulated_data.w.info")[, 2]


### remap hapla population numbering based on order in ref inds
### as they are assigned arbitrarily at admixture, should now work for unequal ref no. and admixed refs
remap <- as.character(0:(npop - 1))
n_ref_per <- rle(read.table("samplemap.tsv", header = FALSE)[, 2])$lengths * 2
n <- 1
for (i in 1:npop){ # finds mode of ancestry in first 100 wd's of each ref panel
    names(remap)[i] <- names(sort(-table(unlist(hapla[n:(n + n_ref_per[i] - 1), 1:100]))))[1]
    n <- n + n_ref_per[i]
}

hapla[] <- lapply(hapla, function(x) {
    as.numeric(unname(remap[as.character(x)])) })


### read and reformat rfmix path res, TODO: read .msp.tsv instead of taking max of fb.tsv, should've read the manual :))
rf <- read.table("rfmix.fb.tsv", header = TRUE)
rf_wide <- rf[5:ncol(rf)]
rf_max <- matrix(NA, nrow(rf), ncol(rf_wide) / npop)
for (i in 1:(ncol(rf_wide) / npop)) {
    chunk <- rf_wide[, ((i - 1) * npop + 1):(i * npop)]
    rf_max[, i] <- max.col(chunk, ties.method = "first") - 1
}

rf <- cbind(rf[1:4], rf_max)
names(rf)[5:(n_query + 4)] <- (n_ref + 1):(n_ref + n_query)


### plots
png("tract.png", width = 4000, height = 400)
haplo_index <- 1 + n_ref
subset_hap(haplo_index)
plot_hap(haplo_index, remap)
dev.off()

pdf("tracts.pdf", width = 60, height = 6)
for (haplo_index in (1 + n_ref):(n_ref + 20)) {
    subset_hap(haplo_index)
    plot_hap(haplo_index, remap)
}
dev.off()


### boxplot of match for all query haplotypes
matchdf <- data.frame()
for (haplo_index in (1 + n_ref):(n_ref + n_query)){
    for (method in c("hapla", "rfmix")) {
        match_props <- get_match_props(haplo_index)
        for (pop in names(match_props[[method]])){
            matchdf <- rbind(matchdf, list("haplotype" = haplo_index, "method" = method, "pop" = pop, "match" = match_props[[method]][[pop]][1] / match_props[[method]][[pop]][2]))
        }
    }
}

matchplot <- ggplot(matchdf, aes(pop, match, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    stat_summary(aes(group = method), fun=mean, geom="point", shape=23, size=3, color="#000000", fill="#4dff64", position = position_dodge(width = 0.75)) +
    ylim(.7, 1.0)

ggsave("match.png", matchplot, width = 5, height = 5)

### old plots
# png("tract.png", width=4000, height=400)
# haplo_index <- 1 + n_ref
# plot(true_anc[true_anc["sample"] == haplo_index, ]$pos,
#         true_anc[true_anc["sample"] == haplo_index, ]$population,
#         type = "l", lwd = 3, lty = 1, col = "#181818")
# lines(rf$physical_position, unlist(rf[as.character(haplo_index)]), type = "l", lwd = 3, lty = 1, col="#e82626bd")
# lines(names(hapla), hapla[haplo_index,], type = "l", lwd = 3, lty = 1, col = "#25bf28bd")
# dev.off()

# pdf("tracts.pdf", width=60, height=6)

# for (haplo_index in (1+n_ref):(n_ref + n_query)) {

#     plot(true_anc[true_anc["sample"] == haplo_index, ]$pos,
#         true_anc[true_anc["sample"] == haplo_index, ]$population,
#         type = "l", lwd = 3, lty = 1, col = "#181818")
#     lines(rf$physical_position, unlist(rf[as.character(haplo_index)]), type = "l", lwd = 3, lty = 1, col="#e82626bd")
#     lines(names(hapla), hapla[haplo_index,], type = "l", lwd = 3, lty = 1, col = "#25bf28bd")
#     }
# dev.off()
