source("cormad.r")

dat = read.table("data.txt",sep="\t")
st = system.time({
		for(i in  1:20000){
			cc =cormad(dat[,1],dat[,2])
		}
})

cat(st[[1]],"\n")
cat(cc,"\n")
