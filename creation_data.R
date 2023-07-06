creation_data <- function(comm1,comm2,mode,points_number,tol.abs = 0.5,tol.rel = NULL, plot=FALSE) {
  
 
  #3 different modes
  if  (mode == "Different") {
    if(is.null(tol.rel)){
      first <- data.frame(X = c(comm1[,1],min(comm1[,1]) + runif(points_number)*tol.abs), 
                          Y = c(comm1[,2],max(comm1[,2]) - runif(points_number)*tol.abs))
      
      second <- data.frame(X = c(comm2[,1],max(comm2[,1]) - runif(points_number)*tol.abs), 
                           Y = c(comm2[,2],min(comm2[,2]) + runif(points_number)*tol.abs))
    }else{
      first <- data.frame(X = c(comm1[,1],min(comm1[,1]) + runif(points_number)*tol.rel*(max(comm1[,1])-min(comm1[,1]))), 
                          Y = c(comm1[,2],max(comm1[,2]) - runif(points_number)*tol.rel*(max(comm1[,2])-min(comm1[,2]))))
      
      second <- data.frame(X = c(comm2[,1],max(comm2[,1]) - runif(points_number)*tol.rel*(max(comm2[,1])-min(comm2[,1]))), 
                           Y = c(comm2[,2],min(comm2[,2]) + runif(points_number)*tol.rel*(max(comm2[,2])-min(comm2[,2]))))
    }
  }else if (mode == "Uniform") {
    
    first <- data.frame(X = c(comm1[,1],min(comm1[,1]) + runif(points_number)*(max(comm1[,1])-min(comm1[,1]))), 
                        Y = c(comm1[,2],min(comm1[,2]) + runif(points_number)*(max(comm1[,2])-min(comm1[,2]))))
    
    second <- data.frame(X = c(comm2[,1],min(comm2[,1]) + runif(points_number)*(max(comm2[,1])-min(comm2[,1]))), 
                         Y = c(comm2[,2],min(comm2[,2]) + runif(points_number)*(max(comm2[,2])-min(comm2[,2]))))
  }else if (mode =="Similar") {
    
    if (min(comm2[,2]) >= min(comm1[,2])) {
      
      if (min(comm2[,1]) >= min(comm1[,1])) {
        first <- data.frame(X = c(comm1[,1],min(comm2[,1]) + runif(points_number)*tol.abs), 
                            Y = c(comm1[,2],max(comm2[,2]) - runif(points_number)*tol.abs))
        
        second <- data.frame(X = c(comm2[,1],min(comm2[,1]) + runif(points_number)*tol.abs), 
                             Y = c(comm2[,2],max(comm2[,2]) - runif(points_number)*tol.abs))
      }else{
      first <- data.frame(X = c(comm1[,1],min(comm2[,1]) + runif(points_number)*(max(comm1[,1])-min(comm2[,1]))), 
                          Y = c(comm1[,2],min(comm2[,2]) + runif(points_number)*(max(comm2[,2])-min(comm2[,2]))))
      
      second <- data.frame(X = c(comm2[,1],min(comm2[,1]) + runif(points_number)*(max(comm1[,1])-min(comm2[,1]))), 
                           Y = c(comm2[,2],min(comm2[,2]) + runif(points_number)*(max(comm2[,2])-min(comm2[,2]))))
      }
    }
    else {
      
      first <- data.frame(X = c(comm1[,1],min(comm2[,1]) + runif(points_number)*(max(comm1[,1])-min(comm2[,1]))), 
                          Y = c(comm1[,2],min(comm1[,2]) + runif(points_number)*(max(comm2[,2])-min(comm1[,2]))))
      
      second <- data.frame(X = c(comm2[,1],min(comm2[,1]) + runif(points_number)*(max(comm1[,1])-min(comm2[,1]))), 
                           Y = c(comm2[,2],min(comm1[,2]) + runif(points_number)*(max(comm2[,2])-min(comm1[,2]))))
      
    }
    
  }else if  (mode == "Disjoint similar") {
    if(is.null(tol.rel)){
      first <- data.frame(X = c(comm1[,1],max(comm1[,1]) - runif(points_number)*tol.abs), 
                          Y = c(comm1[,2],min(comm1[,2]) + runif(points_number)*tol.abs))
      
      second <- data.frame(X = c(comm2[,1],min(comm2[,1]) + runif(points_number)*tol.abs), 
                           Y = c(comm2[,2],max(comm2[,2]) - runif(points_number)*tol.abs))
    }else{
      first <- data.frame(X = c(comm1[,1],max(comm1[,1]) - runif(points_number)*tol.rel*(max(comm1[,1])-min(comm1[,1]))), 
                          Y = c(comm1[,2],min(comm1[,2]) + runif(points_number)*tol.rel*(max(comm1[,2])-min(comm1[,2]))))
      
      second <- data.frame(X = c(comm2[,1],min(comm2[,1]) + runif(points_number)*tol.rel*(max(comm2[,1])-min(comm2[,1]))), 
                           Y = c(comm2[,2],max(comm2[,2]) - runif(points_number)*tol.rel*(max(comm2[,2])-min(comm2[,2]))))
    }
  }
  else {
    print("Wrong mode")
  }
  
  ###plot
  if(plot==TRUE){
    xlim= c(min(c(comm1[,1],comm2[,1]))-0.5,max(c(comm1[,1],comm2[,1]))+0.5)
    ylim= c(min(c(comm1[,2],comm2[,2]))-0.5,max(c(comm1[,2],comm2[,2]))+0.5)
    hull1 <- chull(first)
    hull2 <- chull(second)
    hull1 <- c(hull1,hull1[1])
    hull2 <- c(hull2,hull2[1]) 
    
    plot(first, col = "blue",ylim=ylim,xlim=xlim,xlab="Trait1", ylab="Trait2")
    points(second, col = "red")
    lines(first[hull1,], col = "blue")
    lines(second[hull2,],col = "red")
  }
  #put OTU names
  OTUs <- rep("OTU",2*nrow(first))
  for (i in 1:(2*nrow(first))) { OTUs[i] <- paste(OTUs[i],i,sep = "") }
  
  new_trait <- data.frame(Trait1 = c(first[,1],second[,1]),
                          Trait2 = c(first[,2],second[,2]))
  rownames(new_trait) <-  OTUs
  
  new_comm <- data.frame(Site1 = c(rep(1,nrow(first)),rep(0,nrow(second))),
                         Site2 = c(rep(0,nrow(first)),rep(1,nrow(second))))
  new_comm <- t(new_comm)
  new_comm <- as.matrix(new_comm)
  colnames(new_comm) <- rownames(new_trait)
  
  res <- list(new_comm,new_trait)
  return(res)
  
}
