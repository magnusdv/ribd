alc = function(x, id, gen) {
  anc = ancestors(x, id, maxGen = gen, inclusive = TRUE)
  length(anc)/sum(2^seq_len(gen))
}

avk = alc

