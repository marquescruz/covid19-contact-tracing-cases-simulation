***NOTE: Still under construction*** 🚧🧱👷

# Simulation of COVID-19 cases arising from high-risk contacts
## Summary
Simulation work 

## Data sources and acknowledgements
### COVID-19 cases and vaccines
Data for Portugal :portugal: on COVID-19 cases (exposed, infected, recovered and deaths) and vaccines (number of individuals completely vaccinated and number of booster doses) was collected through the DSSG data respository (https://github.com/dssg-pt/covid19pt-data).

*We would like to aknowledge that this simulation would have been impossible to build without DSSG’s generous contribution in compiling and curating the data reported in PDF format by DGS ( our :portugal: NHS), specially pertaining to exposed individuals identified.*

### COVID-19 variant prevalence
Data on weekly prevalence of COVID-19 variants in Portugal :portugal: was collected through the ECDC open repository (https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea). We used data collected through TESSy.

## Simulation description
### Epidemiological compartmental model
![plot](https://github.com/marquescruz/covid19-cases-contact-tracing-simulation/blob/a9e45ece85e3bc17e08fab7524eb5ab5f01219a8/Modelo%20VSEI.png)

Compartmental model used. 
Legend: S – Susceptible individuals; E – Exposed individuals; E0, E1, E2, …, En – Sub-compartments of the exposed compartment. Each number represents the number of days since exposure; I - Infected individuals; V - Subset of susceptible individuals immune to SARS-CoV-2 infection through vaccination.

### Transitions
The following equations show the calculations behind each compartment transition

    E-S = ∑_(i=0)^n▒〖Ei-S〗
    
    Ei-S = ∑_(w=0)^n▒〖P(Ei-S)w .Ei〗
    
    E-I = ∑_(i=0)^n▒〖Ei-I〗
    
    Ei-I = ∑_(w=0)^n▒〖P(Ei-I)w .Ei〗
    
    Proportion of new cases from exposed = E-I ÷ New cases
    
    S-I = New cases - E-I
    
    S-E = Efinal - Einitial = E0
    
    V = Pprotection . Vaccinated

E-S : Exposed -> Susceptible

Ei-S : Sub-compartment i of compartment E -> Susceptible

E-I : Exposed -> Infected

Ei-I : Sub-compartment i of compartment E -> Infected

S-I : Susceptible -> Infected

S-E : Susceptible -> Exposed

w : COVID-19 variant

### Assumptions of the simulation
The transitions between E sub-compartments were defined by a Markov chain Monte Carlo simulation: in each simulation day, an individual in a sub-compartment of the E compartment transitions probabilistically to the following sub-compartment, to compartment S, or to compartment I.

A constant rate was assumed between two events: 1) the efficacy after inoculation and maximum efficacy, and 2) maximum efficacy and waned efficacy. Vaccinated individuals were probabilistically placed in subset V each day. The same procedure was replicated to individuals with booster doses. 
Individuals in subset V could transition to compartment E, but their infection risk was set at 0. These individuals were included in transition S-E, following the assumption that the proportion of individuals immune due to vaccination in that transition was equal to the proportion of said individuals in compartment S.

## Files included and description
**functions.R**
  Includes all functions needed to run the simulation

**main_simulation.R**
  Builds an example simulation (fully functioning), based on collected data.

**datasets**
  Folder with data needed for the simulation, namely the number of individuals to keep in each compartment daily (*data_preview.csv*), data on daily vaccinated individuals (*vaccines_preview.csv*) and data on weekly prevalence of each variant (transformed in daily prevalence) (*variant_preview.csv*).
  *datasets* include only a preview of the first 100 rows of each dataset, to show the needed structure for the simulation to run.

## How to use the simulation
### Data inputs
To properly run, the simulation needs several inputs, namely:
- **ds.simulation**, a dataframe which includes a minimum of 10 columns
  1. *Date*
  2. *S* - susceptible
  3. *S(V)* - proportion of susceptible (S) individuals protected from infection through vaccination
  4. *E* - exposed/high risk contacts identified
  5. *I* - infected
  6. *R* - recovered
  7. *D* - deaths
  8. *confirmed* - total of confirmed COVID-19 cases (= I)
  9. *new_confirmed* - daily new confirmed cases
  10. and through: *variants* - must at least have one column named "**other**"

- **transition.table**, one dataframe for each variant included. Must at least have one "other"

- **model.daily.table**, a dataframe to keep daily simulated data. Should have as many columns as the number of 

### Function inputs
**tt.function()**
  - For each variant included, 5 parameters can be inputed
  

### Simulation Outputs
A single dataframe, named **simulation.table**, which includes 5 columns:
  1. *Date*
  2. *Sum* - summation of the estimated proportion each day
  3. *Sum of squares* - summation of the estimation squared proportion each day
  4. *Mean* - Mean proportion each day
  5. *SD* - Standard deviation for the proportion each day

### Special cases
Some parameterizations can lead to "proportions" higher than 1. In such cases, what the model is telling is that we would expect more cases from the exposed individuals alone, than those ultimately observed.
In order to correct for this result, a boolean variable was created (*keep.proportions*). If TRUE, maximum proportion of cases from exposed is 1. If FALSE, "proportions" can be bigger than 1.
In the **main_simulation.R** file (line 132), *keep.proportions* is set to **FALSE**.


## Works using this simulation
This simulation was built in order to implement the research protocol that won the 1st Edition of the Scholarship Award Amélia Leitão (https://www.anmsp.pt/post/anmsp-premeia-investigacao-na-covid-19-realizada-por-medicos-de-saude-publica).

An article based on the simulation shown in this repository has been submitted to a scientific journal with peer-review.

To allow a more personalized experience of simulation with different parameters, a web app (built using Shiny) is currently in production.
