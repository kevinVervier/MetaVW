#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.

# script converting genome ids into (species level) taxids based on a meta-data file

# read input parameters
args = commandArgs(trailingOnly = TRUE)
gis.file = args[1]
meta.file = args[2]
output.file = args[3]

# read data
meta = read.table(meta.file, header = T)
gis = read.table(gis.file)$V1

# extract taxids
ind.match = match(gis, meta$genome.id)
if(sum(is.na(ind.match)) > 0){
	stop("gis not found")
}
taxids = meta$taxid.species[ind.match]

# create output file
write.table(taxids,  file = output.file, row.names = F, col.names = F, quote = F)



