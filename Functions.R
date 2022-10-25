# Transition table function ----

tt.function <- function(max.days.of.vigilance.EI = 10,
                        mean.time.to.infection.EI = 3,
                        max.days.of.vigilance.ES = 12,
                        mean.time.to.susceptible.ES = 7,
                        infection.risk.not.vaccinated = 0.375) {
  
  infection.risks <- infection.risk.not.vaccinated
                     
  #' Cumulative Probability of Infection
  cpi <- c(0, 
           ppois(1:(max.days.of.vigilance.EI-1), mean.time.to.infection.EI), 
           1)
  
  #' Transition probabilities from Exposed to Infection
  d <- length(cpi)
  transition.probabilities.EI <- (cpi[2:d]-cpi[1:(d-1)])/(1-cpi[1:(d-1)])
  
  
  #' Cumulative Probability of Susceptible
  cps <- c(0,
           ppois(1:(max.days.of.vigilance.ES-1), mean.time.to.susceptible.ES), 
           1)
  
  #' Transition probabilities from Exposed to Susceptible
  d <- length(cps)
  transition.probabilities.ES <- (cps[2:d]-cps[1:(d-1)])/(1-cps[1:(d-1)])
  
  
  #' Assume both transition arrays last the same time in days  
  ifelse(length(transition.probabilities.EI) <= length(transition.probabilities.ES), 
         transition.probabilities.EI[length(transition.probabilities.EI):length(transition.probabilities.ES)] <- 1, 
         transition.probabilities.ES[length(transition.probabilities.ES):length(transition.probabilities.EI)] <- 1)
  
  
  
  #' ---------------------------
  #' Transition table
  for (j in 1:length(infection.risks)) {

    intermediate.table <- data.frame(tpEI = transition.probabilities.EI,
                                     tpES = transition.probabilities.ES)
    
    for (i in 1:nrow(intermediate.table)) {
      if (i == 1){
        intermediate.table$IRisk[i] = infection.risks[j]
        intermediate.table$SRisk[i] = 1 - infection.risks[j]
      } else {
        intermediate.table$IRisk[i] = intermediate.table$Prob.EE_I[i-1] 
        intermediate.table$SRisk[i] = intermediate.table$Prob.EE_S[i-1]
      }
      
      intermediate.table$Prob.EE_I[i] = (1 - intermediate.table$tpEI[i]) * (intermediate.table$IRisk[i]) #Prob. of staying in Exposed, given it will be INFECTED in the end the follow up (intermediate value for calculation)
      intermediate.table$Prob.EE_S[i] = (1 - intermediate.table$tpES[i]) * (intermediate.table$SRisk[i]) #Prob. of staying in Exposed, given it will be SUSCEPTIBLE in the end the follow up (intermediate value for calculation)
      intermediate.table$Prob.EE[i] = intermediate.table$Prob.EE_I[i] + intermediate.table$Prob.EE_S[i]
      intermediate.table$Prob.EI[i] = (intermediate.table$tpEI[i]) * (intermediate.table$IRisk[i])
      intermediate.table$Prob.ES[i] = (intermediate.table$tpES[i]) * (intermediate.table$SRisk[i])
    }
    
    if (j == 1) {
      a <- intermediate.table[,7:9]
      colnames(a) <- paste((colnames(infection.risks)[j]),colnames(a), sep=".")
      transition.table <- a
    } else {
      a <- intermediate.table[,7:9]
      colnames(a) <- paste((colnames(infection.risks)[j]),colnames(a), sep=".")
      transition.table <- cbind(transition.table,a)
    }
  }
  
  return(transition.table)
}

# Expand nrows of TABLES 
expand.to.nrow <- function(probability.table, max.length){
  if (nrow(probability.table) < max.length){probability.table[(nrow(probability.table)+1):max.length,] <- 0}
  return(probability.table)
}



# Vaccine efficacy functions ----
prop.vaccinated <- function(dataset,
                            t0.efficacy = 0.31,
                            max.efficacy = 0.70,
                            waned.efficacy = 0.50,
                            time.to.max.efficacy = 15, #in days
                            time.to.waned.efficacy = 180, #in days
                            population = 10000000){
  
  rate.to.max <- -log(1-(max.efficacy - t0.efficacy)) / time.to.max.efficacy
  rate.to.waned <- -log(1-(waned.efficacy - max.efficacy)) / time.to.waned.efficacy
  
  efficacy <- NULL #vector of vaccine efficacy per day
  for (t in 0:time.to.max.efficacy){
    efficacy[length(efficacy)+1] = t0.efficacy + (1-exp(-rate.to.max*t))
  }
  for (t in 1:time.to.waned.efficacy){
    efficacy[length(efficacy)+1] =  max.efficacy + (1-exp(-rate.to.waned*t))
  }
  
  days = max(nrow(dataset), length(dataset)) #lenght of dataset in number of days
  tot.vac <- rep(0L, days)
  efficacy[(length(efficacy)+1):length(tot.vac)] <- efficacy[length(efficacy)]
  
  for (day in 1:length(tot.vac)){
    tot.vac[day:length(tot.vac)] = tot.vac[day:length(tot.vac)] + (dataset[day] * efficacy[1:(length(efficacy)-day+1)])
  }
  
  return(tot.vac/population)
} #-- end of prop.function



# Model functions ----
# Dataset fine tuning
fine.tuning <- function(data,
                        vaccines,
                        variants) {
  if (vaccines == FALSE) {data[,3] <- 0L}
  if (variants == FALSE) {data[,10:(ncol(data)-1)] <- 0; data[,ncol(data)] <- 1}
  return(data)
}


# Exposed->E/S/I Function
transitions <- function(input.output=model.daily.table,
                        list.transition.tables=tables.list){
  for (row in rownames(input.output)[1:(nrow(input.output)-1)]){
    if (row == "S(V)") {
      table <- eval(as.name(list.transition.tables[grep("other", list.transition.tables)]))
      table[,2] <- 0
    } else {
      table <- eval(as.name(list.transition.tables[grep(row, list.transition.tables)]))
    } #Choose the transition table to follow
    
    input.output[row,"EI"] <- 0L

    for (transition in nrow(table):1) {
      if (input.output[row,transition]==0) {next}
      if (sum(table[transition,])==0) {next}
      input.output[row,(transition+1)] <- sum(sample(x=0:1,
                                                     size=input.output[row,transition],
                                                     replace=T,
                                                     prob=c(table[transition,2]+table[transition,3],table[transition,1]))) #Calculate E->E

      remainder <- input.output[row,transition] - input.output[row,(transition+1)] #People who have not transitioned to the next E state.
      if (remainder==0) {next}
      if (sum(table[transition,2:3])==0) {next}
      ES <- sum(sample(x=0:1,
                       size=remainder,
                       replace=T,
                       prob=c(table[transition,2], table[transition,3]))) #Calculate E->S
      
      input.output[row,"ES"] <- input.output[row,"ES"] + ES #Summation off all E->S
      
      input.output[row,"EI"] <- input.output[row,"EI"] + remainder - ES #Calculate E->I
    }
  }
  input.output[nrow(input.output),] <- colSums(input.output[(1:(nrow(input.output)-1)),])

  return(input.output)
}


# EI-new cases proportion calculation
prop.cases.exposed <- function(daily.cases=today$confirmados_novos,
                               expected.daily.cases=model.daily.table["total","EI"],
                               rule=keep.proportions) {
  proportion <- expected.daily.cases / daily.cases
  
  if (proportion >=1 & rule==T) {
    return(1)
  } else {
    return(proportion)
  }
}

# E0 calculations
new.exposed <- function(data=today, 
                        input.output=model.daily.table) {
  
  input.output[,1] <- 0L
  already.exposed <- sum(input.output["total",(1:(ncol(input.output)-2))])
  exposed <- data[,"E"]
  
  if (exposed < already.exposed){ #If we have less people exposed than those inside model...
    for (col in (ncol(input.output)-2):1){ #... those on the later days of follow-up... 
      input.output[,col] <- 0 #... end their follow-up and return to susceptible...
      already.exposed <- sum(input.output["total",(1:(ncol(input.output)-2))])
      if(exposed >= already.exposed) {break} #... until we can introduce new exposed (E0).
    }
  }
  
  newly.exposed <- input.output["total","E0"] <- exposed - already.exposed
  group.exposed <- rownames(input.output)[1:(nrow(input.output)-1)]
  probs.group.exposed <- data[,group.exposed]
  probs.group.exposed[,(2:ncol(probs.group.exposed))] <- (1-probs.group.exposed[,1]) * probs.group.exposed[,(2:ncol(probs.group.exposed))] # Probability of having any variant AFTER accounting for vaccination.
  new.exposed.per.group <- data.frame(table(sample(colnames(probs.group.exposed),
                                                   size=newly.exposed,
                                                   replace=T,
                                                   prob=probs.group.exposed[1,])))
  
  for(row in as.character(new.exposed.per.group[,1])) {
    input.output[row, 1] <- new.exposed.per.group[which(row==new.exposed.per.group$Var1),2] 
  }

  return(input.output)
}
