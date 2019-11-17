readname="e7e72eea-c040-4d34-9c2b-bb3b0e7f39d2"
#deletion
test1 <- subset(df, query=="5S_unit2_del2" &subject==readname & qstart <10 & qend > 45)
test1 <-test1[order(test1$sstart),]
View(test1)
#unit 1
test2<-subset(df, query=="cel-5s-unit1" &subject==readname & sstart >16000 &sstart <17000)
test2
test3<-subset(df, query=="cel-5s-unit1" &subject==readname & sstart >39000 &sstart <40500)
test3
