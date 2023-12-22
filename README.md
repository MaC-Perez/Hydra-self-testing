# Hydra-self-testing
The objective of this project is to perform self-tests by fitting the model to data simulated from the operating model (simulation version of Hydra) for a subset of stocks for the Georges Bank ecosystem in the Northwest Atlantic. The operating model was conditioned based on simpler version for two predator and two prey stocks (Spiny dogfish, Atlantic cod, Atlantic mackerel, and Atlantic herring) using 5 length bins, 2 fleets and one survey.

Estimated parameters: 
•	Initial year abundance for each stock 
•	Average annual recruits for each stock 
•	Recruitment deviations for each stock
•	Annual F for each fleet 
•	Fishery catchability for fleet
•	Fishery logistic selectivity for each fleet
•	Survey catchability for each stock and survey 
•	Survey logistic selectivity for each stock and survey


START WITH THE FILE sim_code.R to read the GB real data (hydraDataList), read the outputs from the initial run or OM, create the simulated data sets with the OM and save the ts.dat files, then run the model nsim times and save the .rep and .par files for each simulation. 

Folder R: some functions to read data and create ts files
Folder inputs: original files for differents scenarios (files for 4 species created by Sarah)   
Folder inputs --> initial run: OM used to create the simulations 
Folder sims: nsim ts.dat files 
Folder sims--> rep nsim .rep files (outputs with the simulated data)
Folder sims--> par nsim .par files (outputs with the simulated data)

