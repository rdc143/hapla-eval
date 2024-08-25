### colors for plotting
popcols <- c("0" = "#ee4d25", "1" = "#2ee28e", "2" = "#2e8be2", "3" = "#a543cf")


### extract single haplotype from results
subset_hap <- function(haplo_index) {
    true_anc_sub <<- data.frame(pos = true_anc[true_anc["sample"] == haplo_index, ]$pos,
                               pop = true_anc[true_anc["sample"] == haplo_index, ]$population)
    true_anc_sub$len <<- c(diff(true_anc_sub$pos), seqlen - true_anc_sub[nrow(true_anc_sub), 1])

    hapla_sub <<- data.frame(pos = as.numeric(names(hapla)),
                            pop = unlist(hapla[haplo_index, ]))
    hapla_sub[1, 1] <<- 0
    hapla_sub <<- hapla_sub[c(1, 1 + which(diff(hapla_sub$pop) != 0)), ] # only keep first instances of each pop in tract to remove redundancy and make plots look nicer
    hapla_sub$len <<- c(diff(hapla_sub$pos), seqlen - hapla_sub[nrow(hapla_sub), 1])

    rf_sub <<- data.frame(pos = rf$physical_position,
                         pop = round(rf[, as.character(haplo_index)])) # rounding perhaps removes a little nuance
    rf_sub[1, 1] <<- 0
    rf_sub <<- rf_sub[c(1, 1 + which(diff(rf_sub$pop) != 0)), ] # only keep first instances of each pop in tract
    rf_sub$len <<- c(diff(rf_sub$pos), seqlen - rf_sub[nrow(rf_sub), 1])
}


get_match_props <- function(haplo_index) {
    subset_hap(haplo_index)

    pops <- as.character(unique(true_anc_sub$pop))
    hapla_sums <- integer(length(pops))
    names(hapla_sums) <- pops
    rf_sums <- integer(length(pops))
    names(rf_sums) <- pops
    for (tract in seq(2, nrow(true_anc_sub))){ # split inferred tracts where true tracts change to make sure no inferred tract extends beyond any true tract
        split_index_hapla <- sum(hapla_sub$pos < true_anc_sub[tract, "pos"])
        newlen <- true_anc_sub[tract, "pos"] - hapla_sub[split_index_hapla, "pos"]
        diff <- hapla_sub[split_index_hapla, "len"] - newlen
        hapla_sub[split_index_hapla, "len"] <- newlen
        hapla_sub <- rbind(hapla_sub, c(true_anc_sub[tract, "pos"], hapla_sub[split_index_hapla, "pop"], diff))
        hapla_sub <- hapla_sub[order(hapla_sub$pos), ]

        split_index_rf <- sum(rf_sub$pos < true_anc_sub[tract, "pos"])
        newlen <- true_anc_sub[tract, "pos"] - rf_sub[split_index_rf, "pos"]
        diff <- rf_sub[split_index_rf, "len"] - newlen
        rf_sub[split_index_rf, "len"] <- newlen
        rf_sub <- rbind(rf_sub, c(true_anc_sub[tract, "pos"], rf_sub[split_index_rf, "pop"], diff))
        rf_sub <- rf_sub[order(rf_sub$pos), ]
    }
    for (tract in seq(1, nrow(true_anc_sub))){ # go through all true ancestry tracts
        startpos <- true_anc_sub[tract, 1]
        if (tract != nrow(true_anc_sub)) {
            stoppos <- true_anc_sub[tract + 1, 1]
        } else {
            stoppos <- seqlen
        }
        true_pop <- as.character(true_anc_sub[tract, 2])

        # find tracts wholly inside true tract
        overlap_hapla <- hapla_sub[hapla_sub$pos >= startpos & hapla_sub$pos < stoppos, ]
        overlap_rf <- rf_sub[rf_sub$pos >= startpos & rf_sub$pos < stoppos, ]

        # get length of tracts that match pop and add sum to vector
        hapla_sums[true_pop] <- hapla_sums[true_pop] + sum(overlap_hapla[true_pop == overlap_hapla[, 2], 3])
        rf_sums[true_pop] <- rf_sums[true_pop] + sum(overlap_rf[true_pop == overlap_rf[, 2], 3])
    }
    hapla_res <- list()
    rf_res <- list()
    for (pop in pops){ ### report matching proportion per pop
        hapla_res[[pop]] <- c(hapla_sums[pop], sum(true_anc_sub[true_anc_sub$pop == pop, 3]))
        rf_res[[pop]] <- c(rf_sums[pop], sum(true_anc_sub[true_anc_sub$pop == pop, 3]))
    }
    hapla_res <- hapla_res[order(names(hapla_res))]
    rf_res <- rf_res[order(names(rf_res))]
    list("hapla" = hapla_res, "rfmix" = rf_res)
}


plot_hap <- function(haplo_index, remap) {
    ind_index <- (haplo_index + 1) %/% 2

    q <- read.table(list.files(pattern = "hapla.*.Q")[1])
    admix_str <- "Admixture proportions inferred from hapla: "
    for (i in 1:dim(q)[2]){
        admix_str <- paste0(admix_str, i - 1, ": ", signif(q[ind_index, as.numeric(names(remap[i])) + 1], 4), "   ")
    }

    true_admix <- read.table("trueQ.csv", header = FALSE)
    true_admix_str <- "True admix proportions: "
    for (i in 1:dim(q)[2]){
        true_admix_str <- paste0(true_admix_str, i - 1, ": ", true_admix[ind_index, i], "   ")
    }

    match_props <- get_match_props(haplo_index)
    for (pop in names(match_props[["hapla"]])){
        match_props[["hapla"]][pop] <- signif(match_props[["hapla"]][[pop]][1] / match_props[["hapla"]][[pop]][2], 3)
        match_props[["rfmix"]][pop] <- signif(match_props[["rfmix"]][[pop]][1] / match_props[["rfmix"]][[pop]][2], 3)
    }

    titlestr <- paste("Hapla match per population:", paste(names(match_props[["hapla"]]), match_props[["hapla"]], sep = ": ", collapse = ", "), "  |  ",
                      "RFmix match per population:", paste(names(match_props[["rfmix"]]), match_props[["rfmix"]], sep = ": ", collapse = ", "))
    
    plot(c(0, seqlen), c(0, 1), type = "n", main = titlestr,
        sub = paste(admix_str,  "  |  ", true_admix_str),
        ylab = "", yaxt = "n", xlab = "Position\n", cex.main = 3, cex.lab = 2, cex.sub = 3)
    text(-seqlen * 0.02, .2, "True Ancestry", cex = 1.8)
    text(-seqlen * 0.015, .4, "RFmix", cex = 1.8)
    text(-seqlen * 0.015, .6, "Hapla", cex = 1.8)

    for (start in seq(1, nrow(true_anc_sub))){
        segments(true_anc_sub[start, 1], 0.2, true_anc_sub[start + 1, 1], 0.2, col = popcols[as.character(true_anc_sub[start, 2])], lend = 1, lwd = 64)
    }
    segments(true_anc_sub[nrow(true_anc_sub), 1], 0.2, seqlen, 0.2, col = popcols[as.character(true_anc_sub[nrow(true_anc_sub), 2])], lend = 1, lwd = 64)

    for (start in seq(1, nrow(hapla_sub))){
        segments(hapla_sub[start, 1], .6, hapla_sub[start + 1, 1], .6, col = popcols[as.character(hapla_sub[start, 2])], lend = 1, lwd = 64)
    }
    segments(hapla_sub[nrow(hapla_sub), 1], .6, seqlen, .6, col = popcols[as.character(hapla_sub[nrow(hapla_sub), 2])], lend = 1, lwd = 64)

    for (start in seq(1, nrow(rf_sub))){
        segments(rf_sub[start, 1], .4, rf_sub[start + 1, 1], .4, col = popcols[as.character(rf_sub[start, 2])], lend = 1, lwd = 64)
    }
    segments(rf_sub[nrow(rf_sub), 1], .4, seqlen, .4, col = popcols[as.character(rf_sub[nrow(rf_sub), 2])], lend = 1, lwd = 64)

    legend("right", legend = names(popcols)[1:length(unique(true_anc_sub$pop))],
                   col = popcols[1:length(unique(true_anc_sub$pop))],
                   lty = 1, cex = 2.5, lwd = 16)

     if (file.exists("hapla.prob")) {
        text(-seqlen * 0.02, .8, "Hapla Posterior\nProbability", cex = 1.7)
        prob <- unlist(read.table("hapla.prob")[haplo_index, ])
        wdpos <- unlist(read.table("hapla.simulated_data.w.info")[, 2])
        lines(c(0, seqlen), c(0.7, 0.7))
        lines(c(0, seqlen), c(0.9, 0.9))
        lines(wdpos, prob * 0.2 + 0.7, lwd = 1.5, col = "red")
    }
}
