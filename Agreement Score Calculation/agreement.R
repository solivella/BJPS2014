
weights.rice <- function(votevec){
  votevec <- votevec[complete.cases(votevec)]
  Y <- length(votevec[votevec==1])
  N <- length(votevec[votevec==2])
	return(1 - abs((Y-N)/(Y+N)))
}

agreement <- function(votevec){
  tab.res <- sort(-table(c(votevec)))
  #cat("votevec is:", votevec,"\n")
  mod.choice <- names(tab.res)[1]
  multimod <- any(tab.res[-1]==tab.res[1])
  #cat("multimod is:", multimod,"\n")
  #cat("mod.choice is:", mod.choice,"\n")
  if(multimod|(mod.choice<1|mod.choice>3)){
    #cat("forced comparison is: ",rep(0,length(votevec)),"\n")
    return(rep(0,length(votevec)))
  }else{
    votevec[which(votevec<1|votevec>3)] <- mod.choice
    #cat("comparison is: ", as.numeric(votevec!=mod.choice),"\n")
    return(votevec!=mod.choice)    
  }  
}

tot.vot <- function(votemat){
  #print(ifelse(votemat>0&votemat<4,1,0))
  rowSums(ifelse(votemat>0&votemat<4,1,0))
}

get.scores <- function(){
  score.by.party <- list()
  seq.of.years <- seq(1994,2006,by=4)
  for(i in seq.of.years){
    cat("I'm working on",i,"\n")
    temp.votes <- read.csv(paste(i,".csv",sep="")
                           ,encoding=ifelse(i==1994,"latin1","utf8"))
    names.parties <- by(temp.votes[,c(2:3)],temp.votes[,3],function(x)x)
    weights <- apply(temp.votes[,-c(1:3)],2,weights.rice)
    dis.mat.list <- by(temp.votes[,-c(1:3)],temp.votes[,3],FUN=apply,2,agreement)
    tot.vec.list <- by(temp.votes[,-c(1:3)],temp.votes[,3],FUN=tot.vot)
    by.party.ind.scores.w <- lapply(dis.mat.list,function(x){x%*%weights})
    by.party.ind.scores.uw <- lapply(dis.mat.list,function(x){x%*%rep(1,length(weights))})
    for(p in 1:length(dis.mat.list)){
      score.by.party[[p*i]] <- names.parties[[p]]
      score.by.party[[p*i]][,3] <- i
      score.by.party[[p*i]][,4] <- by.party.ind.scores.w[[p]]/tot.vec.list[[p]]
      score.by.party[[p*i]][,5] <- by.party.ind.scores.uw[[p]]/tot.vec.list[[p]]
      colnames(score.by.party[[p*i]]) <- c("MP","party","year","disagree.score.w","disagree.score.uw")
    }
  }
  return(do.call(rbind,score.by.party))
}

agreement.scores <- get.scores()






  
