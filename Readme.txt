##################################################

Paper: 
A Hidden Markov model addressing missing ordinal responses 
and replicated covariates to monitor Parkinson’s disease 
progression based on recorded speech


Authors:
Lizbeth Naranjo (1), Carlos J. Perez (2), Yolanda Campos-Roca (3).

(1) Departamento de Matemáticas, Facultad de Ciencias, 
Universidad Nacional Autonoma de Mexico (UNAM), Mexico

(2) Departamento de Matematicas, Facultad de Veterinaria, 
Universidad de Extremadura, Spain

(3) 3Departamento de Tecnologias de los Computadores y de las Comunicaciones, 
Escuela Politecnica, Universidad de Extremadura, Spain


Journal: 
Submitted. Under Revision. 

##################################################

Instructions to run the codes in R and JAGS are provided. 
The codes are applied to obtain a similar analysis as in Section 4 ‘Simulation based case’, but without cross-validation. 

##################################################

##################################################
FILES 

The file ‘HMMprogressionOrdinalMissing.R’ contains the R code. The JAGS code is run from this R file.

The file ‘HMMprogressionOrdinalMissing.bug' contains the JAGS model. 


##################################################

To run the files, do the following.
 
1.- Download JAGS from www.mcmc-jags.sourceforge.net/

2.- Install the packages necessary to run the R file. 
These are indicated in the R file. 

3.- Change the address indicated in ‘setwd()’. 
setwd("HERE"). 
This is the address where the file ‘HMMprogressionOrdinalMissing.bug’ is in.

##################################################


