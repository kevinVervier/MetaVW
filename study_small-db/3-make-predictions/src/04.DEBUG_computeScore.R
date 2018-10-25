

# read parameters
args = commandArgs(trailingOnly = TRUE)
input.kmers  = args[1]
input.weights = args[2]


# read kmers
kmers = readLines(input.kmers)
kmers = strsplit(kmers, split = " " )[[1]]
kmers = kmers[-c(1,2)]
kmers = table(kmers)

# read weights
w = as.character(read.table(input.weights)$V1)
w = strsplit(w, split = ":")
weights = as.numeric(sapply(w, function(x){x[[3]]}))
k = sapply(w, function(x){x[[1]]})
k = strsplit(k, split = "\\[")
k  = sapply(k, function(x){x[[1]]})
names(weights) = k

# compute score
ind.match = match(names(kmers), names(weights))
score = sum(kmers * weights[ind.match])

# add intercept
score = score + weights[grep("Constant", names(weights))]
cat("score found =", score, "\n")

