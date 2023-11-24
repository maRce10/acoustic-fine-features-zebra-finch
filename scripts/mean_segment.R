mean_segment <-
    function(wave,
             cores = 1,
             plot = TRUE,
             pb = TRUE,
             thinning = 1,
             col = wave_col,
             mean = TRUE,
             type = "ac",
             npeak = 20) {
        # thin
        if (thinning < 1) {
            if (length(wave@left) * thinning < 10) {
                stop2("thinning is too high, no enough samples left for at least 1 sound file")
            }
            
            # reduce size of envelope
            wavefrm <-
                stats::approx(
                    x = seq(0, duration(wave), length.out = length(wave@left)),
                    y = wave@left,
                    n = round(length(wave@left) * thinning),
                    method = "linear"
                )$y
        } else {
            wavefrm <- wave@left
        }
        
        # get empirical mode decomposition
        if (type == "EMD") {
            emds <-
                EMD::emd(wavefrm, seq_len(length(wavefrm)), boundary = "wave")
            
            perd <- emds$imf[, 4] / max(emds$imf[, 4])
            # plot(x = seq_len(length(wavefrm)), y = perd, type = "l")
            # lines(y = wavefrm / max(wavefrm), x = seq_len(length(wavefrm)), col = "gray", lty = 2)
        }
        
        if (type == "ac") {
            ac <-
                acf(
                    x = wavefrm,
                    lag.max = length(wavefrm),
                    type = "covariance",
                    demean = FALSE,
                    plot = FALSE
                )
            perd <- ac$acf / max(ac$acf)
        }
        
        tpks <-
            seewave::fpeaks(cbind(seq_len(length(perd)), perd),
                            plot = FALSE,
                            threshold = 0.5
            )
        
        if (nrow(tpks) > npeak) {
            tpks <- tpks[1:npeak, ]
        }
        
        segment_df <-
            data.frame(
                selec = seq_len(nrow(tpks)),
                pos = tpks[, 1],
                peak = tpks[, 2]
            )
        
        # get mean number of sample between peaks
        mean_dist_peak <- round(mean(diff(segment_df$pos)))
        
        segment_df$start <- segment_df$pos - mean_dist_peak / 2
        segment_df$end <- segment_df$pos + mean_dist_peak / 2
        
        
        # fix if values are out of wavefrm size
        if (segment_df$start[1] > 0) {
            segment_df$start[1] <- 0
        }
        if (segment_df$end[nrow(segment_df)] > length(wavefrm)) {
            segment_df$end[nrow(segment_df)] <- length(wavefrm)
        }
        
        # extract segments into a list
        segments <- lapply(seq_len(nrow(segment_df)), function(x) {
            wavefrm[segment_df$start[x]:segment_df$end[x]]
        })
        
        
        # make all the same number of samples
        segments <-
            lapply(segments, function(x) {
                approx(x, n = max(sapply(
                    segments, length
                )))$y
            })
        
        # normalize between 1, -1
        segments <- lapply(segments, function(x) {
            x / max(x)
        })
        
        # put all segments in a data frame
        segments <-
            as.data.frame(segments, col.names = seq_len(length(segments)))
        
        # compute mean segment
        mean_segment <- rowMeans(segments)
        
        if (plot) {
            mean_segment_df <-
                data.frame(
                    time = seq(0, 1, length.out = nrow(segments)),
                    mean.amp = rowMeans(segments),
                    sd.amp = apply(segments, 1, sd)
                )
            
            gg <- ggplot(
                data = mean_segment_df,
                mapping = aes(x = time, y = mean.amp)
            ) +
                geom_line(color = wave_col) +
                geom_ribbon(aes(ymin = mean.amp - sd.amp, ymax = mean.amp + sd.amp),
                            alpha = 0.2
                ) +
                theme_classic(base_size = 25)
            
            print(gg)
        }
        if (mean) {
            return(list(mean_segment = mean_segment, mean_dist_peak = mean_dist_peak))
        } else {
            return(segments)
        }
    }