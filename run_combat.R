library(ggfortify)
library(sva)
library(bladderbatch)
a <- c("E175A2", "E175A2")
b <- c("E175A1", "E175B")
all <- read.csv("/public/home/zhangzb/test/work/raw_data/raw_count-2/all_rawcount.csv",
                header = T, row.names = 1)
for (i in c(1:2)){
    # mh070 <- t(all[, c(grep(a[i], colnames(all)))])
    mh125 <- t(all[, c(grep(b[i], colnames(all)))])

    mh070 <- t(read.csv(paste0("/public/home/zhangzb/test/work/combat/combat_data2/",
                        a[i], "-combat.csv"),
                    header = T,
                    row.names = 1))
    # mh125 <- t(read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count-2/",
    #                     b[i], ".csv"),
    #                 header = T,
    #                 row.names = 1))
    # mh070 <- log2(mh070 + 1)
    # mh125 <- log2(mh125 + 1)
    gene <- intersect(colnames(mh070), colnames(mh125))
    rst1 <- as.data.frame(rep(a[i], length(rownames(mh070))))
    colnames(rst1) <- "label"
    mh070 <- cbind.data.frame(rst1, mh070)


    rst2 <- as.data.frame(rep(b[i], length(rownames(mh125))))
    colnames(rst2) <- "label"
    mh125 <- cbind.data.frame(laabel = rst2, mh125)

    combins <- rbind(mh070, mh125)
    #combins <- combins[which(colSums(combins) > 0),]

    out_pca <- prcomp(combins[, -1]) #删除第一列
    #plot(out_pca, type="l")
    setwd("/public/home/zhangzb/test/work/combat/pca")
    f1 <- paste0(a[i], "-", b[i], "-rawPCA.jpeg")
    jpeg(filename = f1)
    f2 <- autoplot(out_pca, data=combins, colour ='label', size = 0.1,)
            #  label=TRUE, label.size=3)
    # print(f2)
    # dev.off()
    edata <- t(combins[, -1])
    pheno <- combins
    combat_edata <- ComBat(dat = edata, batch = pheno$label)
    mh070_combat <- combat_edata[, c(grep(a[i], colnames(combat_edata)))]
    mh125_combat <- combat_edata[, c(grep(b[i], colnames(combat_edata)))]
    # write.csv(mh070_combat,
    #         paste0("/public/home/zhangzb/test/work/combat/combat_data2/",
    #         a[i], "-combatlog.csv"))
    write.csv(mh125_combat,
        paste0("/public/home/zhangzb/test/work/combat/combat_data2/",
        b[i], "-combat.csv"))
    # write.csv(combat_edata, "/public/home/zhangzb/test/work/combat/combat_data2/E175A2-combatlog.csv")
    dat <- as.data.frame(t(combat_edata))
    dat <- as.data.frame(cbind(pheno$label, dat), stringsAsFactors = FALSE)
    colnames(dat)[1] <- "label"

    combat_pca <- prcomp(dat[, -1])
    f3 <- paste0(a[i], "-", b[i], "-combatPCA.jpeg")
    jpeg(filename = f3)
    f4 <- autoplot(combat_pca, data = dat, colour ='label', size = 0.1)
    # print(f4)
    # dev.off()
}

library(sva)

sample <- c("E135", "E155", "E165", "E175", "P0")
for (i in sample){
    if (i == "E175" | i == "P0"){
        a <- paste0(i, "A1")
        b <- paste0(i, "A2")
        count1 <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                  a, "_rawcount.csv"), header = T, row.names = 1, check.names =F)
        count2 <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                  b, "_rawcount.csv"), header = T, row.names = 1, check.names =F)
        count <- cbind(count1, count2)
        batch1 <- as.data.frame(rep(a, length(colnames(count1))), row.names = colnames(count1))
        colnames(batch1) <- "label"
        batch2 <- as.data.frame(rep(b, length(colnames(count2))), row.names = colnames(count2))
        colnames(batch2) <- "label"
        batch <- rbind(batch1, batch2)
        combat_count <- ComBat(dat = count, batch = batch$label)       
    }
    else{
        a <- paste0(i, "A")
        b <- paste0(i, "B")
        count1 <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                  a, "_rawcount.csv"), header = T, row.names = 1, check.names =F)
        count2 <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                  b, "_rawcount.csv"), header = T, row.names = 1, check.names =F)
        count <- cbind(count1, count2)
        batch1 <- as.data.frame(rep(a, length(colnames(count1))), row.names = colnames(count1))
        colnames(batch1) <- "label"
        batch2 <- as.data.frame(rep(b, length(colnames(count2))), row.names = colnames(count2))
        colnames(batch2) <- "label"
        batch <- rbind(batch1, batch2)
        combat_count <- ComBat(dat = count, batch = batch$label)
    }
    if(i == "E135"){
        combat_count2 <- combat_count
    }
    else{
        combat_count2 <- cbind(combat_count2, combat_count)
    }
}
write.csv(combat_count2, "/public/home/zhangzb/test/work/combat/all_combat.csv")