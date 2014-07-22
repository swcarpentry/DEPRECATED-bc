# These settings control the behavior of all chunks in the novice R materials.
# For example, to generate the lessons with all the output hidden, simply change
# `results` from "markup" to "hide".
# For more information on available chunk options, see
# http://yihui.name/knitr/options#chunk_options

library("knitr")
opts_chunk$set(tidy = FALSE, results = "markup")
