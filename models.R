library(MCMCglmm)
library(plyr)
library(ggplot2)
set.seed(1)
ntime <- 10		#模擬總時間
jnt <- 14  		#初始稚體數
ant <- 56  		#初始成體數
abund <-c(jnt,ant) 	#初始族群數量
j11 <- 0.04 		#稚體於苦茶粕影響下存活且停留於稚體的機率平均值
sd_j11 <- 0.05 	#稚體於苦茶粕影響下存活且停留於稚體的機率標準差
j12 <- 0.13		#稚體於苦茶粕影響下存活且成長為成體的機率平均值
sd_j12<- 0.05		#稚體於苦茶粕影響下存活且成長為成體的機率標準差
jsurv <- j11+j12	#稚體於苦茶粕影響下存活機率平均值
fer <-1.64		#成體每單位時間繁殖率平均值
sd_fer <- 0.59		#成體每單位時間繁殖率標準差
a22 <- 0.88		#成體每單位時間繁殖率平均值
sd_a22 <- 0.03		#成體每單位時間繁殖率標準差
nrep <- 1000		#重複模擬次數

ant_arr <-matrix(0,2,ntime+1)			#儲存族群動態與結構
colnames(ant_arr) <- c(seq(0,ntime/2,0.5))
rownames(ant_arr) <- c("juvenile","adult")
ant_arr[,1] <- abund
h30_list <- list()
h30_arr <- array(0,dim=c(2,ntime+1,nrep),
                 dimnames = list(c("jvenile","adult"),
                                 c(seq(0,ntime/2,0.5)),
                                 c(seq(1,nrep))))
#模擬各生命力為截尾常態分布
j11_trn <- round(rtnorm(nrep,mean=j11,sd=sd_j11,lower = 0),3)
j12_trn <- round(rtnorm(nrep,mean=j12,sd=sd_j12,lower = 0),3)
#假若稚體存活率>1則進行線性調整
for(i in 1:nrep){
  if((j11_trn[i]+j12_trn[i])>1){
    c <- j11_trn[i]+j12_trn[i]
    j11_trn[i] <- j11_trn[i]/c
    j12_trn[i] <- j12_trn[i]/c
  }
}
fer_trn <- round(rtnorm(nrep,mean=fer,sd=sd_fer,lower = 0),3)
a22_trn <- round(rtnorm(nrep,mean=a22,sd=sd_a22,lower = 0,
                        upper = 1),3)
h30_total_arr <- array(0,dim=c(1,ntime+1,nrep))
h30_arr <- array(0,dim=c(2,ntime+1,nrep))
for(ii in 1:nrep){   
  for(i in 1:ntime){
    #mpm為隨機抽樣的矩陣族群模式
    mpm <- matrix(c(j11_trn[ii],j12_trn[ii],fer_trn[ii],
                    a22_trn[ii]),2,2)
    ant_arr[,i+1] <- mpm %*% ant_arr[,i]
    ant_arr[,i+1] <- ant_arr[,i+1]
  } 
  h30_list[[ii]] <- ant_arr[1,]+ant_arr[2,]
  h30_total_arr[,,ii] <- ant_arr[1,]+ant_arr[2,]
  h30_arr[,,ii] <- ant_arr
}
#經過模擬後族群動態的平均值
sto_pop_mean = aaply(laply(h30_list, as.matrix),c(2), mean)
#經過模擬後族群動態的標準差
sto_pop_sd = aaply(laply(h30_list, as.matrix), 2, sd)
sto_pop_mean =as.data.frame(sto_pop_mean)
sto_pop_sd =as.data.frame(sto_pop_sd)
ymin= sto_pop_mean - sto_pop_sd
ymin[ymin<0]=0
colnames(ymin)=c("ymin")
ymax= sto_pop_mean + sto_pop_sd
colnames(ymax)=c("ymax")
#儲存族群模擬的平均值與標準差
ABUND_SD <- data.frame(sto_pop_mean = sto_pop_mean,ymin=ymin,ymax=ymax,
                       Days=seq(seq(0,ntime/2,0.5)))
# expo_day <- ggplot(ABUND_SD,aes(x=Days,y=ans2,
# 原代碼
expo_day <- ggplot(ABUND_SD,aes(x=Days,y=sto_pop_mean,
                                group=1))+
  geom_line(lwd=1)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= ymin,ymax= ymax),width=.3,
                position=position_dodge(.9),size=1)+
  xlab("Days")+ylab("Abundance")+
  geom_hline(yintercept=1, linetype="dashed", color = "red",lwd=1)+
  theme(text=element_text(family = "D"),
        axis.text.x = element_text(face="bold", size=I(20)),
        axis.text.y = element_text(face="bold", size=I(20)),
        axis.title.x = element_text(size = 15,face="bold"),
        axis.title.y = element_text(size = 15,face="bold"))
#微型裸腹蚤受苦茶粕、溫度與食物量影響的族群動態模擬
expo_day
#計算滅絕率，族群量小於1時滅絕
ep <- sum(h30_total_arr[,11,]<1)/nrep*100
ep
