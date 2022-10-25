#'==========================================================
#' title: COVID-19 cases from high-risk contacts simulation
#' 
#' file: Main
#' authors: marquescruz , jomigalves  
#' last_version: Oct. 2022
#'==========================================================

source("Functions.R")

# Parameterizations ####
## Population ####
population = 1e+07 # Number of susceptible individuals at the start. Usual approximation for the Portuguese population.
## Transition table for variants ####
# Each variant included must have a "transition.table.VARIANT" created through tt.function().
# The name for the variable must be kept as "transition.table.variant"
# One of the variants must be named "other" and must always be created.
# In this example code, we included 5 variants (and the "other" variant), in concordance to the data on variants used in this simulation.
transition.table.other <- tt.function(max.days.of.vigilance.EI = 10,
                                      mean.time.to.infection.EI = 3,
                                      max.days.of.vigilance.ES = 12,
                                      mean.time.to.susceptible.ES = 7,
                                      infection.risk.not.vaccinated = 0.375)
transition.table.alpha <- tt.function(7,4,14,9)
transition.table.beta <- tt.function(7,4,14,9)
transition.table.gamma <- tt.function(7,4,14,9)
transition.table.delta <- tt.function(6,3,14,9)
transition.table.omicron <- tt.function(14,2,14,9)
tables.list <- ls()[grep('transition.table.', ls())]

# All transition tables must have the same number of rows.
max.days <- NULL
for (table in tables.list){
  max.days <- max(max.days, nrow(eval(as.name(table))))
}

for (table in tables.list){
  assign(table, expand.to.nrow(eval(as.name(table)), max.days)) 
}; rm(max.days, table)

## Protection from vaccines ####
# This variables are used to determine the number of individuals protected from infection through vaccines.
# Refer to README.md for more information on this topic.

# Efficacy measures for first inoculation
t0.efficacy = 0.31 # Vaccine efficacy at inoculation.
max.efficacy = 0.70 # Maximal vaccine efficacy 
waned.efficacy = 0.50 # Vaccine efficacy after waning period
time.to.max.efficacy = 15 # Number of days from inoculation to maximal vaccine efficacy
time.to.waned.efficacy = 180 # Number of days from maximal vaccine efficacy to waned vaccine efficacy


# Efficacy measures for booster doses
booster.t0.efficacy = 0.70
booster.max.efficacy = 0.80 
booster.waned.efficacy = 0.60
booster.time.to.max.efficacy = 15
booster.time.to.waned.efficacy = 180


# Data ####
## Importation ####
# Data included in this example is from 2nd March 2020 to 20th January 2022.
# Datasets included ("data.csv", vaccines.csv", "variant.csv") are in the main branch of the GitHub repository (marquescruz/covid19-contact-tracing-cases-simulation).
# These datasets are adapted and include only the data needed for the simulation to run.
# Original data sources are described in README.md/Data sources and acknowledgements.
# Refer to README.md/How to use the simulation/Input data for information on the expected input data structure.
ds <- read.csv("datasets/data_preview.csv")
vaccines <- read.csv("datasets/vaccines_preview.csv") 
variant <- read.csv("datasets/variant_preview.csv")

## Preprocessing ####
### Vaccines ####
vaccines$vaccinated_protected <- prop.vaccinated(vaccines$vaccinated_daily,
                                                 population = population,
                                                 t0.efficacy = t0.efficacy,
                                                 max.efficacy = max.efficacy,
                                                 waned.efficacy = waned.efficacy,
                                                 time.to.max.efficacy = time.to.max.efficacy,
                                                 time.to.waned.efficacy = time.to.waned.efficacy)

vaccines$booster_protected <- prop.vaccinated(vaccines$vaccinated_daily,
                                              population = population,
                                              t0.efficacy = booster.t0.efficacy,
                                              max.efficacy = booster.max.efficacy,
                                              waned.efficacy = booster.waned.efficacy,
                                              time.to.max.efficacy = booster.time.to.max.efficacy,
                                              time.to.waned.efficacy = booster.time.to.waned.efficacy) * (1-vaccines$vaccinated_protected) #Only considering those who could have gained protection from booster, since "vaccines$vaccinated_protected" were already with protection from the vaccine.

vaccines$S.V <- vaccines$vaccinated_protected + vaccines$booster_protected


### Data merge ####
ds <- merge(ds,vaccines[,c(1,6)], by ="data", all.x=T)
ds$S.V[which(is.na(ds$S.V))] <- 0
ds$S.V <- ds$S.V * ds$S / population #Set protected individuals according to susceptible individuals at each day.
ds$S <- 1e+07 - rowSums(ds[,4:7]) #Susceptible as aprox. national population
ds <- ds[,c(1,9,8,6,7,4,5,2,3)] #Reorder columns

colnames(ds) <- c("Data", "S", "S(V)", "E", "I", "R", "D", "confirmados", "confirmados_novos") #rename columns
ds <- merge(ds, variant, by="Data")
rm(vaccines, variant)

## Templates and fine-tuning ####
### Daily transition table ####
# A table to keep daily simulated data inside the simulation. 
cols <- c(paste0("E", 0:(nrow(transition.table.other)-1)), "EI", "ES")
variant_names <- c(colnames(ds)[10:ncol(ds)], "total")
model.daily.table <- data.frame(matrix(data=0L,
                                       nrow = length(variant_names)+1,
                                       ncol = length(cols),
                                       dimnames=list(c("S(V)",variant_names),
                                                     cols)))
zero.model.daily.table <- model.daily.table
rm(cols, variant_names)

### Dataframe to keep data from repeated simulations ####
simulation.table <- data.frame(Data = ds$Data, 
                               sum = 0L,
                               sum.squares = 0L)

### Fine tuning of the simulation ####
# Fine tuning dataset
ds.simulation <- fine.tuning(data=ds,
                             vaccines=TRUE, #If TRUE, vaccine protection will be taken into account
                             variants=TRUE) #If TRUE, all variants created will be included in the simulation; If FALSE, only the variant "other" will be used.

# Simulation ####
## Model control ####
set.seed(26)
simulations = 10 # Number of simulations to run
keep.proportions = FALSE #If TRUE, maximum proportion of cases from exposed is 1. If FALSE, "proportions" can be bigger than 1. Refer to README.md for more information on this issue.

## Simulation iterations ####
for (simulation in 1:simulations) { ## SIMULATION ITERATION
  for (day in 1:nrow(ds.simulation)) { ## DAILY MODEL ITERATIONS
    # Import information from ds about each day
    today <- ds.simulation[day,]
    
    # En -> En+1 / EI / ES transitions with vaccination or per variant and total
    model.daily.table <- transitions()
    
    # EI/new_confirmed metrics
    simulation.table$sum[day] <-  simulation.table$sum[day] + prop.cases.exposed()
    simulation.table$sum.squares[day] <- simulation.table$sum.squares[day] + prop.cases.exposed()**2
    
    # S -> E0 transitions per vaccinated and per variant
    model.daily.table <- new.exposed()
    
  } #End of 'day' loop
  model.daily.table <- zero.model.daily.table
} #End of 'simulation' loop

# Results ####
# Daily simulated data was treated like a counter. Hence, it is necessary to convert the results into a mean and a standard deviation pertaining to each simulation day.
simulation.table$mean <- simulation.table$sum/simulations #Mean
simulation.table$sd <- sqrt((simulation.table$sum.squares - simulation.table$sum**2 / simulations) / (simulations - 1))  #Standard deviation