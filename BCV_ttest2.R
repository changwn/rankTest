BCV_ttest2<-function(data0,rounds=20,slice0=2,maxrank0=30,msep_cut=0.001)
{
          x<-data0
        fff_cc<-c()
        for(kk in 1:rounds)
        {
                cv_result <- cv.svd.gabriel(x, slice0, slice0, maxrank = maxrank0)
                fff_cc<-rbind(fff_cc,cv_result$msep)
        }
        pp<-c()
        ddd<-apply(fff_cc,2,mean)
          ddd<-ddd/sum(ddd)
        for(kk in 1:(ncol(fff_cc)-1))
        {
                pp_c<-1
                if(mean(fff_cc[,kk],na.rm=T)>mean(fff_cc[,kk+1],na.rm=T))
                {
                                if(ddd[kk]>msep_cut)
                                {
                          pp_c<-t.test(fff_cc[,kk],fff_cc[,kk+1])$p.value
                                }
                }
                pp<-c(pp,pp_c)
        }
          #boxplot(fff_cc)
        return(pp)       
}



