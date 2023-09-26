# /public/home/zhaojinqian/anaconda3/envs/R3.6.3/bin/R
library(SPARK)
sample <- c("E135A", "E155A", "E165A", "E175A1",
           "P0A1","E135B", "E155B", "E165B", "E175A2", "P0A2")
for (i in sample){
    rawcount <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                              i, "_rawcount.csv"), header=T, row.names = 1, check.names = F)
    coor <- read.csv(paste0("/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/",
                              i, "-coor.csv"), header=T, row.names = 1)
    rownames(coor) <- colnames(rawcount)
    spark <- CreateSPARKObject(counts=rawcount, 
                             location=coor[,1:2],
                             percentage = 0.1, 
                             min_total_counts = 10)
    spark@lib_size <- apply(spark@counts, 2, sum)
    spark <- spark.vc(spark, 
                   covariates = NULL, 
                   lib_size = spark@lib_size, 
                   num_core = 5,
                   verbose = F)
    spark <- spark.test(spark, 
                     check_positive = T, 
                     verbose = F)
    saveRDS(spark, paste0("/public/home/zhangzb/test/work/spark/", i, ".rds"))
    print(paste0(i, " done"))
    # a <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
    # a <- a[which(a$adjusted_pvalue <= 0.01),]
    # b <- data.frame()
}
for (i in sample){
    spark <- readRDS(paste0("/public/home/zhangzb/test/work/spark/", i, ".rds"))
    a <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
    b <- row.names(a)
    b <- as.data.frame(b)
    colnames(b) <- paste0(i, "_gene")
    c <- as.data.frame(a[,2])
    colnames(c) <- paste0(i, "_adjusted_pvalue")
    d <- cbind(b, c)
    if (i=="E135A"){
        d2 <- d
    }
    else {
       d2 <- merge(d2, d, by = "row.names", all = T)
       d2 <- d2[,-1]
    }     
}
write.csv(d2, "/public/home/zhangzb/test/work/spark/all_spark.csv")