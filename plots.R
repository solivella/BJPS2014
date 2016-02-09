library(maps)
library(lattice)

##Map of the world
map(fill=TRUE, col="grey")
map(regions=c("Colombia"
              , "Croatia"
              , "Italy"
              , "Israel"
              , "Japan"
              , "Macedonia"
              , "Moldova"
              , "New Zealand"
              , "Romania"
              , "South Africa"
              , "Taiwan"
              , "Ukraine"
              , "Canada"
              ,"UK"),fill=TRUE,col="red",add=TRUE)

##Plots
data.for.plots1 <- array(NA,c(1330*2,3))
#data.for.plots1[,1] <- round(c(post.samples.agg[[1]][,
#                                        grep("theta",
#                                             colnames(post.samples.agg[[1]]))]),3)
data.for.plots1[,1] <- round(c(rnorm(1330,0.308,0.295)
                         ,rnorm(1330,0.716,0.373)),3)
data.for.plots1[,2] <- rep(c("SMD to PR","PR to SMD")
                           ,each=1330)
data.for.plots1[,3] <- "Aggregated"


data.for.plots2 <- array(NA,c(1330*4,3))
#data.for.plots2[,1] <- round(c(post.samples.dis[[1]][,
#                                        c(grep("theta1",
#                                             colnames(post.samples.dis[[1]])),
#                                          grep("theta2\\[1\\]",
#                                             colnames(post.samples.dis[[1]])))
#                                          ]),3)
data.for.plots2[,1] <- round(c(rnorm(1330,0.324, 0.287)
                         ,rnorm(1330,0.317, 51.692)
                         ,rnorm(1330,0.688, 0.379)
                         ,rnorm(1330,0.690, 0.458)),3)
data.for.plots2[,2] <- rep(c("SMD to RPR","SMD to NPR","RPR to SMD","NPR to SMD")
                           ,each=1330)
data.for.plots2[,3] <- "Disaggregated"

data.for.plots3 <- array(NA,c(1330*4,3))
#data.for.plots3[,1] <- round(c(post.samples.inter[[1]][,
#                                        grep("theta",
#                                             colnames(post.samples.inter[[1]]))]),3)
data.for.plots3[,1] <- round(c(rnorm(1330,0.234, 0.311)
                         ,rnorm(1330,0.712, 0.538)
                         ,rnorm(1330,0.645, 0.404)
                         ,rnorm(1330,0.680, 0.399)),3)
data.for.plots3[,2] <- rep(c("SMD to PR, No Nomination Change"
                             ,"SMD to PR, Nomination Change"
                             ,"PR to SMD, No Nomination Change"
                             ,"PR to SMD, Nomination Change")
                           ,each=1330)
data.for.plots3[,3] <- "Conditional"


data.for.plots <- rbind(data.for.plots1
                        ,data.for.plots2
                        ,data.for.plots3)
(data.for.plots[,3]) <- factor(data.for.plots[,3],levels=c("Conditional","Disaggregated","Aggregated"))

bwplot(data.for.plots[,2]~data.for.plots[,1]|data.for.plots[,3],
       scales=list(y=list(relation="free"))
     ,ylab="Treatment"
       , xlab="Effect", 
        main="Effects of Mandate Change By Model", 
        layout=(c(1,3)),
       panel=function(x,y){
         panel.abline(v=0,lty=3)
         panel.bwplot(x,y
                      ,fill="black"
                      ,do.out=FALSE
                      ,stats=function(x,...){
                        return(list(stats=c(quantile(x,probs=0.1)
                                            ,rep(median(x),3)
                                            ,quantile(x,probs=0.9)),
                                    n=length(x)
                                    ,conf=NA
                                    ,out=NA))
                      }
                      ,drop.unused.levels=TRUE
                      )
       }
               ,par.settings = list(box.umbrella = list(lty=1,col = 'black'))
       ,xlim=c(-1,2)
       ,index.cond=list(c(2,3,1))
       ,strip=strip.custom(bg="black",par.strip.text=list(col="white"))
       )

                           



