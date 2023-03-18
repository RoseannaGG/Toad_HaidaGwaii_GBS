
remove_sibs_genlight <- function(genlight_object,
                      names){
  
  initial_length  <- length(genlight_object@ind.names)
  
  print(paste0("Initial object made up of ", initial_length, " individuals"))
  
  for(i in 1:length(names)){
    
    genlight_object <- genlight_object[indNames(genlight_object) != names[i]]
    
  }  
  
 
  
  
  final_length <- length(genlight_object@ind.names)
  
  print(paste0("Final object made up of ", final_length, " individuals"))
  

  return(genlight_object)
  
  
}

