#Copyright (c) by respective owners including bioMerieux, and
#individual contributors Kevin Vervier and Pierre Mahe. All #rights reserved.  Released under a BSD (revised)
#license as described in the file LICENSE.


# specify model parameters #
#--------------------------#
DB = "small"
NBATCHES = "10"


#-----util function------------------
computePerf = function(ref,pred){
	pred.correct = ref == pred
	acc.micro = mean(pred.correct)
	tmp = split(pred.correct, ref)
        acc.micro = round( 100*mean(pred.correct), digits = 2)
	acc.species = round( 100*sapply(tmp, mean), digits = 2)
	acc.macro = round(mean(acc.species), digits = 2)
	acc.median = round(median(acc.species), digits = 2)
	return(list("acc.micro"=acc.micro,"acc.macro"=acc.macro,"acc.median"=acc.median,"acc.species"=acc.species))
}
#------------------------------------


perf = list()

# read results fragments #
#------------------------#
	# read reference labels
ref = read.table("../../1-generate-test-datasets/output/test.fragments.taxid")$V1
	# read predictions
pred = read.table(paste("../output/test.fragments.", DB, "-db.preds.taxid", sep = ""))$V1
	# compute perf
perf.tmp = computePerf(ref, pred)
perf[["fragments"]] = perf.tmp


# read results homo reads #
#-------------------------#
model.list = c("Balzer","Richter","Margulies")
for(model in model.list){
	# read reference labels
	ref = read.table(paste("../../1-generate-test-datasets/output/test.homo-",model,"-reads.taxid", sep = ""))$V1
	# read predictions
	pred = read.table(paste("../output/test.homo-",model,"-reads.",DB,"-db.preds.taxid", sep = ""))$V1
	# compute perf
	perf.tmp = computePerf(ref, pred)
	perf[[model]] = perf.tmp
}


# read results mutation reads #
#-----------------------------#
freq.list = c(seq(1,10),10.7)
for(freq in freq.list){
	# read reference labels
	ref = read.table(paste("../../1-generate-test-datasets/output/test.mutation-",freq,"-reads.taxid", sep = ""))$V1
	# read predictions
	pred = read.table(paste("../output/test.mutation-",freq,"-reads.",DB,"-db.preds.taxid", sep = ""))$V1
	# compute perf
	perf.tmp = computePerf(ref, pred)
	perf[[paste("mutation-",freq,sep="")]] = perf.tmp
}




# plot #
#------#
acc.median = sapply(perf, function(x){x$acc.median})
pdf(paste("../output/median-species-level-accuracy_", DB, "-db.pdf", sep = ""), width = 12)
par(mar = c(8,4,4,2))
bargraph  = barplot(acc.median, ylim = c(50,100), xpd = F, las = 2, xlab = "", ylab = "median species-level accuracy", main = "performance per dataset")
text(bargraph, acc.median , acc.median, pos = 1)
dev.off()


