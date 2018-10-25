
tab1 = read.table("../output/preds.fragments_L-200.from-vw.preds.vw")$V1
tab2 = read.table("../output/preds.fragments_L-200.from-spectrum")$V1


cat("*** fraction of predictions in agreement =", mean(tab1==tab2), "***\n")

# focus on "reachable and stringent" reads
labels = read.csv2("../../../data/small/test-dataset/fragments_L-200.labels.csv")
ind = which(labels$to.keep == 1)
tab1.reach = tab1[ind]
tab2.reach = tab2[ind]
cat("*** fraction of predictions in agreement among reads to keep =", mean(tab1.reach==tab2.reach), "***\n")
ind.dis = which(tab1.reach != tab2.reach)


