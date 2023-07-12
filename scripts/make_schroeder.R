qwindpc <- function(rftime, srate, sig) {
  ramppts <- round((rftime / 1000) * srate)
  hold1 <- seq(0, 1, length.out = ramppts + 1)
  onramp <- sin(0.5 * pi * hold1)^2
  offramp <- cos(0.5 * pi * hold1)^2
  steady <- rep(1, length(sig) - (2 * ramppts) - 2)
  wind <- c(onramp, steady, offramp)

  if (length(wind) > length(sig)) {
    print("Window not applied, block too short for window")
  } else {
    sig <- wind * sig
  }

  return(sig)
}


make_schroeder <-
  function(samp.rate = 44100,
           f0, # fundamental freq
           dur = 40, # in ms
           save.wave = FALSE,
           plot = TRUE,
           n.components,
           scalar = 1,
           color = "#3E4A89FF",
           path = ".",
           file.name = NULL) {
    t <- seq(1 / samp.rate, 1, by = 1 / samp.rate)
    sumwavepos <- rep(0, samp.rate)


    # number of samples
    durpts <- round((dur / 1000) * samp.rate)
    N <- n.components

    phase <- rep(0, N)
    amplin <- rep(1, N)
    f <- rep(0, N)

    Ns <- 1:N

    for (n in Ns) {
      compnum <- n
      phase[compnum] <- scalar * pi * compnum * (compnum - 1) / N
      f[n] <- f0 + ((compnum - 1) * f0)
    }

    wave <- matrix(0, nrow = N, ncol = length(t))

    for (i in Ns) {
      wave[i, ] <- amplin[i] * cos((2 * pi * f[i] * t) + phase[i])
      sumwavepos <- sumwavepos + wave[i, ]
    }

    posschr <- sumwavepos[1:durpts] / max(abs(sumwavepos[1:durpts]))

    # plot(t[1:length(posschr)], posschr, type = "l")

    spectr <- fft(sumwavepos)
    power <- abs(spectr)^2
    sample <- length(t) / max(t)
    freq <- (1:length(t)) * sample / length(t)
    powerdb <- 10 * log10(power[1:(length(t) / 2)])

    rftime <- 10
    posschr <- qwindpc(rftime, samp.rate, posschr)
    posschr <- posschr / max(abs(posschr))

    # if (plot)
    #     plot(t[1:length(posschr)], posschr, type = "l", xlab = "Time (s)", ylab = "Amplitude")
    # max(abs(posschr))

    # Sound playback
    # Note: Sound playback in R can be system-dependent and may require additional setup.
    # Adjust the code accordingly or save the audio to a file for playback in a suitable player.
    # sound(posschr, samp.rate)
    # Sys.sleep(2)
    # sound(posschr, samp.rate)

    wave_obj <- tuneR::normalize(tuneR::Wave(posschr, samp.rate = samp.rate, bit = 16), unit = "16")

    if (plot) {
      opar <- par()
      opar$cin <- opar$cra <- opar$sci <- opar$cxy <- opar$din <- opar$csi <- opar$page <- NULL

      # on.exit(par(opar))
      par(mfrow = c(2, 1), mar = c(5, 4, 0, 0) + 0.1)
      plot(
        freq[1:(length(t) / 2)],
        powerdb,
        type = "l",
        xlim = c(0, 5200),
        xlab = "Frequency (Hz)",
        ylab = "Power(dB)",
        col =  color
      )

      seewave::oscillo(wave = wave_obj, colwave = color)
    }

    # Save audio to a file
    if (save.wave) {
      if (is.null(file.name)) {
        file.name <- paste0("f0-", f0, "_ncomp-", N, "_", if (scalar == 1) "pos" else "neg", ".wav")
      }

      writeWave(object = wave_obj, filename = file.path(path, file.name), extensible = FALSE)
    } else {
      return(wave_obj)
    }
  }
