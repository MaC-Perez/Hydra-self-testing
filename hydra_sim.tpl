// Copyright (c) 2008, 2009, 2010 Regents of the University of California.
//
// ADModelbuilder and associated libraries and documentations are
// provided under the general terms of the "BSD" license.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2.  Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3.  Neither the name of the  University of California, Otter Research,
// nor the ADMB Foundation nor the names of its contributors may be used
// to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
// MultiSpecies Size Structured Assessment Model MS3AM
// 		based on LeMANS as modified by M. Fogarty and S. Gaichas
//		coded by S. Gaichas, help from K. Curti, G. Fay, T. Miller
//		testing version, February 2012
//		initial working version, July 2012
//
// Renamed hydra_sim and modified for use as size structued simulation model
//		operating model for testing produciton, nonlinear time series models
//		May 2013
//
// Dec 29, 2016: Modified and extended by Andy Beet
//             : see word doc, documenting changes


//=======================================================================================
GLOBALS_SECTION
//=======================================================================================
  //Including C++ libraries
  #include "statsLib.h"
  #include <iostream>
  #include <cmath>
  #include <fstream>
  #include <string>
  #include <sstream>
  #include <time.h>
  time_t baseTime;
  clock_t startTime = clock();

  // ofstream test("test.csv"); // for debugging only
  ofstream gavjunk("gav.junk");
  ofstream pmse_predvals("pmse_predvals.out");


//=======================================================================================
DATA_SECTION
//=======================================================================================

//Debug statements from Kiersten Curti's model
  //  0's = Main body of program
      //  =  1:  Exits after data section
      //  =  2:  Exits after parameter section
      //  =  3:  Exits at end of procedure section
      //  =  4:  Prints checkpoints after each function in the procedure section, except those within the year loop
      //  =  5:  Prints checkpoints after each function within the year loop
      //  =  6:  Prints parameter estimates after each iteration
      //  =  7:  Prints pop dy variables after each iteration
      //  =  8:  Prints trophic matrices at end of year loop, then exits
      //  =  9:  Prints out predicted indices that go into the objective function after each iteration
  //10's = Initial states function
      //  = 10:  Outputs food-selection parameters at the end of the function and then exits
      //  = 11:  Outputs food selection paramters at the end of each iteration
      //  = 12:  Outputs fishery and survey selectivity matrices at the end of the function and then exits
      //  = 13:  Outputs abundance and biomass arrays and then exits
      //  = 14:  Outputs Yr1, Age1 and iFt matrices to ensure 'means + devt'ns' parameterized correctly and then exits
      //  = 15:  Outputs initial N and proportion mature arrays
  //20's = Suitability function
      //  = 20:  Output at end of function
      //  = 21:  Output at end of function, then exits
      //  = 22:  Outputs suitability and scaled suitability for each predator sp and then exits
      //  = 23:  Outputs Eta, Sigma, WtRatio, G for pred, prey combos and then exits
  //30's = Predation mortality function
      //  = 30:  Prints checkpoints throughout the function
      //  = 31:  Prints intermediate arrays (Avail, scAv, Consum) for each predator sp and then exits
      //  = 32:  Prints intermediate arrays for each prey sp (Consum, sumCon) and then exits
      //  = 33:  Prints sumCon and B right before M2 is calculated, and M2 after it is calculated
      //  = 34:  Debug 31 but does not exit
      //  = 35:  Prints nf and Other Food in every year
  //40's = Population dynamics function
      //  = 40:  Outputs N, C_hat and Cprop_hat at end of function
      //  = 41:  Outputs N, C_hat and Cprop_hat at end of function and then exits
      //  = 42:  Outputs mortality components for each species in each year and exits at end of year loop after trophic =1
      //  = 43:  Outputs mortality components at end of function and then exits
  //50's = Survey abundance function
      //  = 50:  Prints intermediate arrays for species where survey data is one contiguous time series and then exits
      //  = 51:  Prints intermediate arrays for species where the survey data is split into multiple segments and then exits
      //  = 52:  Prints predicted q, FICs, FIC_hat and N for each species and then exits
      //  = 53:  Prints estimated q matrix at the end of each iteration
  //60's = Log likelihood function
      //  = 60: Prints checkpoints after each likelihood component
      //  = 61: Prints checkpoints for multinomial components within the year loop
      //  = 62: Prints predicted and log predicted indices for TotC and TotFIC
      //  = 63: Prints predicted and log predicted indices for Cprop
      //  = 64: Prints predicted and log predicted indices for Sprop
      //  = 65: FHprop, when added
      //  = 66: Prints summary of objective function components at end of each iteration
  //70's = Food habits function
      //  = 70:  Prints Avpd (scAv) and FHtmp for each predator, year and then exits
      //  = 71:  Prints Avpd, Avpy for each prey species within Avpd, the colsum of Avpy, and FHtmp; exits after all predator species
      //  = 72:  Prints bin years, FHtmp, FHsum, and average FH, FH_hat, for each FH bin, and exits after all predator species
      //  = 73:  Prints total %W for each pred age, summed across prey sp, for each pred sp and FH bin.  This value should always == 1.
  //80's = Penalty functions
      // = 80: Biomass penalty function: Prints pre- and post- B for every species and year, regardless of whether a penalty is imposed
      // = 81: Biomass penalty function: Prints pre- and post- biomass and assorted arrays when biomass falls below threshold
      // = 82: Yr1 penalty function: Prints avgZ, thYr1, Yr1 and Ypen arrays and then exits at end of function
      // = 83: Recruitment penalty function: Prints Age1 parameter estimates, calculated CV and recruitment penalty

//take random number seeds from command line,
//code from https://groups.nceas.ucsb.edu/non-linear-modeling/projects/skate/ADMB/skatesim.tpl
//#1: spinydog
//#2: winterskate
//#3: Aherring
//#4: Acod
//#5: haddock
//#6: yellowtailfl
//#7: winterfl
//#8: Amackerel
//#9: silverhake
//#10: goosefish

  int sim;
  int rseed;
 LOC_CALCS


    int on,opt;
    sim=0;
     rseed=1;
//    rseed=123456;
    //the following line checks for the "-sim" command line option
    //if it exists the if statement retreives the random number seed
    //that is required for the simulation model
    if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
    {
      sim=1;
      rseed=atoi(ad_comm::argv[on+1]);
    }

 END_CALCS

  init_int debug;                          //Debugger switch, 0 off, other above, now in dat file

//read in bounding indices from .dat file
  init_int Nyrs  												//number of years
  init_int Nspecies												//number of species
  init_int Nsizebins											//number of size bins
  init_int Nareas												//number of areas
  init_int Nfleets												//number of fleets
  init_int Nsurveys                       //number of surveys (e.g. NEFSC spring, NEFSC fall, NEAMAP, etc)  
//  init_int Nages												//number of age classes
  int Totsizebins
  !!  Totsizebins = Nspecies*Nsizebins;
  int spp                         //loop counter species
  int size                         //loop counter lengthbins
  int area                        //loop counter areas
  int pred								 //loop counter predators
  int prey								 //loop counter prey
  int t									 //loop counter timestep
  int yr								 //loop counter year
  int fleet              //loop counter fleet
  int Nstepsyr           //model timesteps per year, set using max(growthprob_phi) below
  int Tottimesteps       //total timesteps for dimensioning state variable arrays
  int yrct               //year counter for dynamics
  int iguild             // loop counter in calc_survey_abundance
  int iassess            // loop counter in calc_assessment_strategy
  int ithreshold         // loop counter in calc_assessment_strategy
  ivector maxThreshold(1,Nareas) // stored most severe  threshold detected

  int Nprey
  !! Nprey = Nspecies + 1;

  number o
  !!  o = 0.001;         //small number to add before taking logs in objective function calcs
  // should add 1 since log of 1 = 0. objecive function not changed! Andy Beet

  init_number wtconv     //multiply by weight in grams to desired units for biomass, catch

 // time series data file
  init_adstring datfilename;

//read in length bin sizes, weight at length parameters from .dat file
  init_matrix binwidth(1,Nspecies,1,Nsizebins) 					//length bin width in cm
  //init_3darray binwidth(1,Nareas,1,Nspecies,1,Nsizebins) 		//length bin width in cm, area specific
  init_vector lenwt_a(1,Nspecies)				//weight = lenwt_a * length ^ lenwt_b UNITS cm to g
  init_vector lenwt_b(1,Nspecies)				//weight = lenwt_a * length ^ lenwt_b UNITS cm to g
  //init_matrix lenwt_a(1,Nareas,1,Nspecies) 					//weight = lenwt_a * length ^ lenwt_b, area specific
  //init_matrix lenwt_b(1,Nareas,1,Nspecies) 					//weight = lenwt_a * length ^ lenwt_b, area specific

//calculate length and weight attributes of bins
//dimension all by area and add outside area loop to make these area specific
  matrix lbinmax(1,Nspecies,1,Nsizebins)						//upper end of length bins
  matrix lbinmin(1,Nspecies,1,Nsizebins)						//lower end of length bins
  matrix lbinmidpt(1,Nspecies,1,Nsizebins)						//midpoint of length bins
  matrix wtbinmax(1,Nspecies,1,Nsizebins)						//max weight of length bins
  matrix wtbinmin(1,Nspecies,1,Nsizebins)						//min weight of length bins
  matrix wtatlbinmidpt(1,Nspecies,1,Nsizebins)					//wt at length midpoint of length bins
  matrix binavgwt(1,Nspecies,1,Nsizebins)						//average weight for length bins
  matrix powlbinmaxb(1,Nspecies,1,Nsizebins)
  !!	for (spp=1; spp<=Nspecies; spp++){
  !!		lbinmax(spp,1)   = binwidth(spp,1);
  !!		lbinmin(spp,1)   = 0.0;						//lowest bin assumed to start at 0!
  !!		lbinmidpt(spp,1) = binwidth(spp,1)/2.0;
  !!		for (size=2; size<=Nsizebins; size++){
  !!			lbinmax(spp, size)   = binwidth(spp, size) + lbinmax(spp, size-1);
  !!			lbinmin(spp, size)   = lbinmax(spp, size-1);
  !!			lbinmidpt(spp, size) = binwidth(spp, size)/2.0 + lbinmax(spp, size-1);
  !!		}
  !!    	wtbinmax(spp) = lenwt_a(spp)* pow(lbinmax(spp), lenwt_b(spp));
  !!    	wtbinmin(spp) = lenwt_a(spp)* pow(lbinmin(spp), lenwt_b(spp));
  !!    	wtatlbinmidpt(spp) = lenwt_a(spp)* pow(lbinmidpt(spp), lenwt_b(spp));
  !!	}
  !!	binavgwt = (wtbinmin + wtbinmax)/2.0;

//read in covariate information from .dat file
  init_int Nrecruitment_cov  									//number of recruitment covariates
  init_int Nmaturity_cov  										//number of maturity covariates
  init_int Ngrowth_cov  										//number of growth covariates
  init_matrix recruitment_cov(1,Nrecruitment_cov,1,Nyrs)		//time series of recruitment covariates
  init_matrix maturity_cov(1,Nmaturity_cov,1,Nyrs)				//time series of maturity covariates
  init_matrix growth_cov(1,Ngrowth_cov,1,Nyrs)
  
   //time series of growth covariates
  //init_3darray recruitment_cov(1,Nareas,1,Nrecruitment_cov,1,Nyrs)  //time series of recruitment covariates, area specific
  //init_3darray maturity_cov(1,Nareas,1,Nmaturity_cov,1,Nyrs)  //time series of maturity covariates, area specific
  //init_3darray growth_cov(1,Nareas,1,Ngrowth_cov,1,Nyrs)   	//time series of growth covariates, area specific

//read in survey and catch observations from .dat file
// SKG: for now, sub in Nsurveys for the unused Nareas dimension for input survey indices and comps
  //init_3darray obs_survey_biomass(1,Nareas,1,Nspecies,1,Nyrs)  	//spring or fall? units needed

  init_3darray obs_effort(1,Nareas,1,Nfleets,1,Nyrs)  	//standardized effort units needed
  //init_4darray for survey size comp by area, species, year?
  //init_4darray obs_survey_size(1,Nsurveys,1,Nspecies,1,Nyrs,1,Nsizebins)  //numbers, uncomment when in dat file
  //init_5darray for catch at size by area, species, fleet, year?

//read in mean stomach content weight time series from .dat file for intake calculation
  init_4darray mean_stomwt(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)

  //want variance for this in this for fitting?


//read in temperature time series from .dat file for intake calculation
  init_matrix obs_temp(1,Nareas,1,Nyrs)       //want variance measure? data source?

//read in estimation phases from .dat file
  init_int yr1Nphase            //year 1 N at size estimation phase
  init_int recphase				//recruitment parameter estimation phase
  init_int avg_rec_phase		//average recruitment estimation phase (could make species specific, currently global)
  init_int recsigmaphase   // std dev of recruit resids phase (could make species specific, currently global)
  init_int avg_F_phase			//average fishing mort estimation phase (could make species specific, currently global)
  init_int dev_rec_phase		//recruitment deviation estimation phase (could make species specific, currently global)
  init_int dev_F_phase			//fishing mort deviation estimation phase (could make species specific, currently global)
  init_int fqphase              //fishery q estimation phase
  init_int fsphase              //fishery selectivity estimation phase
  init_int sqphase              //survey q estimation phase
  init_int ssphase              //survey selectivity estimation phase  
  init_int ssig_phase           //survey sigma (obs error) phase
  init_int csig_phase           //catch sigma (obs error) phase
  init_int m1_phase            // M1 phase
  init_int oF1_phase           // amount of other food included in the M2 term for the base (predator 1) phase
  init_int oFdev_phase          //deviation from base other food for predators 2+ phase
  init_int vuln_phase           // phase for vulnerability parameters

//read in lists of species names, area names, fleet names, actual years from .dat file

//the following are parameters that will be fixed in initial runs, so read in as "data" from .dat
//to estimate any of them within the model, place in parameter section, initialize from .pin file
//recruitment parameters from .dat file
  init_matrix recGamma_alpha(1,Nareas,1,Nspecies)			//eggprod gamma Ricker model alpha
  init_matrix recGamma_shape(1,Nareas,1,Nspecies)			//eggprod gamma Ricker model shape parameter
  init_matrix recGamma_beta(1,Nareas,1,Nspecies)			//eggprod gamma Ricker model beta

  init_matrix recDS_alpha(1,Nareas,1,Nspecies)			//SSB Deriso-Schnute model alpha
  init_matrix recDS_shape(1,Nareas,1,Nspecies)			//SSB Deriso-Schnute model shape parameter
  init_matrix recDS_beta(1,Nareas,1,Nspecies)			//SSB Deriso-Schnute model beta

  init_matrix recGamSSB_alpha(1,Nareas,1,Nspecies)			//SSB gamma alpha
  init_matrix recGamSSB_shape(1,Nareas,1,Nspecies)			//SSB gamma shape parameter
  init_matrix recGamSSB_beta(1,Nareas,1,Nspecies)			//SSB gamma beta

  init_matrix recRicker_alpha(1,Nareas,1,Nspecies)			//SSB Ricker model alpha
  init_matrix recRicker_shape(1,Nareas,1,Nspecies)			//SSB Ricker model shape parameter=1.0 not used
  init_matrix recRicker_beta(1,Nareas,1,Nspecies)			//SSB Ricker model beta

  init_matrix recBH_alpha(1,Nareas,1,Nspecies)			//SSB Beverton Holt model alpha
  init_matrix recBH_shape(1,Nareas,1,Nspecies)			//SSB Beverton Holt model shape parameter=1.0 not used
  init_matrix recBH_beta(1,Nareas,1,Nspecies)			//SSB Beverton Holt model beta

  init_matrix recShepherd_alpha(1,Nareas,1,Nspecies)			//SSB Shepherd model alpha
  init_matrix recShepherd_shape(1,Nareas,1,Nspecies)			//SSB Shepherd model shape parameter=1.0 not used
  init_matrix recShepherd_beta(1,Nareas,1,Nspecies)			//SSB Shepherd model beta

  init_matrix recHockey_alpha(1,Nareas,1,Nspecies)			//SSB Hockey Stick model alpha
  init_matrix recHockey_shape(1,Nareas,1,Nspecies)			//SSB Hockey stick  model shape (S*) breakpoint
  init_matrix recHockey_beta(1,Nareas,1,Nspecies)			//SSB Hockey Stick  model beta parameter=1.0 not used

  init_matrix recSegmented_alpha(1,Nareas,1,Nspecies)			//SSB Segmented regression model alpha
  init_matrix recSegmented_shape(1,Nareas,1,Nspecies)			//SSB  Segmented regression model shape. use this for the breakpoint
  init_matrix recSegmented_beta(1,Nareas,1,Nspecies)			//SSB  Segmented regression model beta
  init_ivector rectype(1,Nspecies)  //switch for alternate recruitment functions
  
  init_ivector stochrec(1,Nspecies)  //switch for stochastic recruitment

  matrix rec_alpha(1,Nareas,1,Nspecies)//
  matrix rec_shape(1,Nareas,1,Nspecies)//
  matrix rec_beta(1,Nareas,1,Nspecies)//

  //initialize recruitment parameters for each type (case 9 no functional form, dummy pars)
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!	  switch (rectype (spp)){
  !!       case 1:	  				//egg production based recruitment, 3 par gamma (Ricker-ish)
  !!		  rec_alpha(area,spp) = recGamma_alpha(area,spp);
  !!		  rec_shape(area,spp) = recGamma_shape(area,spp);
  !!		  rec_beta(area,spp) = recGamma_beta(area,spp);
  !!	   break;
  !!	   case 2:                   //SSB based recruitment, 3 par Deriso-Schnute; see Quinn & Deriso 1999 p 95
  !!          rec_alpha(area,spp) = recDS_alpha(area,spp);
  !!          rec_shape(area,spp) = recDS_shape(area,spp);
  !!          rec_beta(area,spp) = recDS_beta(area,spp);
  !!       break;
  !!	   case 3:                   //SSB based recruitment, 3 par gamma
  !!          rec_alpha(area,spp) = recGamSSB_alpha(area,spp);
  !!          rec_shape(area,spp) = recGamSSB_shape(area,spp);
  !!          rec_beta(area,spp) = recGamSSB_beta(area,spp);
  !!       break;
  !!	   case 4:                   //SSB based recruitment, 2 par Ricker
  !!          rec_alpha(area,spp) = recRicker_alpha(area,spp);
  !!          rec_shape(area,spp) = recRicker_shape(area,spp);
  !!          rec_beta(area,spp) = recRicker_beta(area,spp);
  !!       break;
  !!	   case 5:                   //SSB based recruitment, 2 par BevHolt
  !!          rec_alpha(area,spp) = recBH_alpha(area,spp);
  !!          rec_shape(area,spp) = recBH_shape(area,spp);
  !!          rec_beta(area,spp) = recBH_beta(area,spp);
  !!       break;
  !!       case 6:                // SSB based recruitment, 3 parameters Shepherd
  !!          rec_alpha(area,spp) = recShepherd_alpha(area,spp);
  !!          rec_shape(area,spp) = recShepherd_shape(area,spp);
  !!          rec_beta(area,spp) = recShepherd_beta(area,spp);
  !!       break;
  !!       case 7:                // SSB based rcruitment. hockeyStick
  !!          rec_alpha(area,spp) = recHockey_alpha(area,spp);
  !!          rec_shape(area,spp) = recHockey_shape(area,spp);
  !!          rec_beta(area,spp) = recHockey_beta(area,spp);
  !!       break;
  !!       case 8:                // SSB based recruitment,segmented regresion with breakpoint
  !!          rec_alpha(area,spp) = recSegmented_alpha(area,spp);
  !!          rec_shape(area,spp) = recSegmented_shape(area,spp);
  !!          rec_beta(area,spp) = recSegmented_beta(area,spp);
  !!       break;
  !!	   case 9:                   //no functional form, uses average+devs in .pin file
  !!          rec_alpha(area,spp) = 0;
  !!          rec_shape(area,spp) = 0;
  !!          rec_beta(area,spp) = 0;
  !!       break;
  !!       default:
  !!          cout<<"undefined recruitment type, check .dat file"<<endl;
  !!          exit(1);
  !!		}
  !!    }
  !! }

  init_matrix sexratio(1,Nareas,1,Nspecies)  //this is proportion female
  init_matrix recruitment_covwt(1,Nspecies,1,Nrecruitment_cov)	//recruitment covariate weighting factor
  //init_3darray recruitment_covwt(1,Nareas,1,Nspecies,1,Nrecruitment_cov) //area specific weighting

//fecundity parameters from .dat file and calculate fecundity at length
  init_matrix fecund_d(1,Nareas,1,Nspecies)
  init_matrix fecund_h(1,Nareas,1,Nspecies)
  init_3darray fecund_theta(1,Nareas,1,Nspecies,1,Nsizebins)
  3darray fecundity(1,Nareas,1,Nspecies,1,Nsizebins)
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!    	fecundity(area, spp) = elem_prod(fecund_theta(area, spp),
  !!											(fecund_d(area, spp)
  !!                            	 			* pow(lbinmidpt(spp),fecund_h(area, spp))));
  !!	}
  !!  }

//maturity parameters from .dat file
  init_matrix maturity_nu(1,Nareas,1,Nspecies)
  init_matrix maturity_omega(1,Nareas,1,Nspecies)
  init_matrix maturity_covwt(1,Nspecies,1,Nmaturity_cov)
  matrix covariates_M(1,Nspecies,1,Nyrs) // intermediate calculation to obtain maturity covariates //AndyBeet
//growth parameters from .dat file and calculate simple (no cov) prob of growing through length interval
  init_matrix growth_psi(1,Nareas,1,Nspecies)    //power function growth length=psi*age^kappa
  init_matrix growth_kappa(1,Nareas,1,Nspecies)  //power function growth length=psi*age^kappa
  init_matrix growth_covwt(1,Nspecies,1,Ngrowth_cov)
  init_matrix vonB_Linf(1,Nareas,1,Nspecies)    //alternate parameterization, vonB growth
  init_matrix vonB_k(1,Nareas,1,Nspecies)       //alternate parameterization, vonB growth
  init_ivector growthtype(1,Nspecies)                          //switch for alternate growth types

  init_number phimax
  4darray growthprob_phi(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)
  vector delta_t(1,Nsizebins)
  vector lmax_test(1,2);
  number lmax_use
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!     for(yr=1; yr<=Nyrs; yr++){
  !!      switch (growthtype (spp)){
  !!        case 1:	  	 //exponential no covariates
  !!              // these yimes  of growing out of bin relative to size zero. SGaichas
  !!    	  delta_t = pow((lbinmax(spp)/growth_psi(area, spp)),(1.0/growth_kappa(area, spp)));
  !!              // these are "probabilities" of growing into size bin, r given in size bin r-1. ABeet
  !!              growthprob_phi(area, spp, yr,1) = 1/delta_t(1);
  !!              for (int isize=2;isize<=Nsizebins; isize++) {
  !!                growthprob_phi(area, spp, yr, isize) = 1/(delta_t(isize) - delta_t(isize-1));
  !!              }
  !!        break;
  !!        case 2:       //exponential with covariates
  !!              // these "probabilitlies of growing out of bin relative to size zero. SGaichas
  !!    	  delta_t = pow((lbinmax(spp)/growth_psi(area, spp)* mfexp(growth_covwt(spp)*trans(growth_cov)(yr))),
  !!	        										(1.0/growth_kappa(area, spp)));
  !!              // these are "probabilities" of growing into size bin, r given in size bin r-1. ABeet
  !!              growthprob_phi(area, spp, yr,1) = 1/delta_t(1);
  !!              for (int isize=2;isize<=Nsizebins; isize++) {
  !!                growthprob_phi(area, spp, yr, isize) = 1/(delta_t(isize) - delta_t(isize-1));
  !!              }
  !!        break;
  !!        case 3:       //VonB no covariates          
//  !!          lmax_test(1) = 1e-09; 
//  !!          for (size=1;size<=Nsizebins;size++) {
//  !!          lmax_test(2) = vonB_Linf(area, spp)-lbinmax(spp,size);
//  !!          lmax_use = max(lmax_test);
//  !!          cout << spp << " " << yr << " " << size << " " << lmax_use << endl;
//  !!          growthprob_phi(area, spp, yr, size) = vonB_k(area, spp)/log(
//  !!                                          ((vonB_Linf(area, spp)-lbinmin(spp,size))/
//  !!                                                   (lmax_use)));
//  !!          }
    !!          for (size=1;size<=Nsizebins;size++) {
    !!            if (lbinmin(spp,size)>=vonB_Linf(area, spp) || lbinmax(spp,size)>=vonB_Linf(area, spp)) {
    !!             growthprob_phi(area, spp, yr, size) = 0.0;
    !!            }
    !!            if (lbinmax(spp,size)<vonB_Linf(area, spp)) {
    !!               growthprob_phi(area, spp, yr, size) = vonB_k(area, spp)/log(
    !!                                                   (vonB_Linf(area, spp)-lbinmin(spp,size))/
    !!                                                   (vonB_Linf(area, spp)-lbinmax(spp,size)));
    !!            }
    !!          }
//  !!          growthprob_phi(area, spp, yr) = vonB_k(area, spp)/log(
//  !!                                          elem_div((vonB_Linf(area, spp)-lbinmin(spp)),
//  !!                                                   (vonB_Linf(area, spp)-lbinmax(spp))));
//  !! cout << "Growthprob " << spp << " " << growthprob_phi(area,spp,yr) << endl;
  !!        break;
  !!        case 4:       //VonB with covariates
  !!          growthprob_phi(area, spp, yr) = vonB_k(area, spp)/log(
  !!                                          elem_div((vonB_Linf(area, spp)*mfexp(growth_covwt(spp)*trans(growth_cov)(yr))-lbinmin(spp)),
  !!                                                   (vonB_Linf(area, spp)*mfexp(growth_covwt(spp)*trans(growth_cov)(yr))-lbinmax(spp))));
  !!        break;
  !!        default:
  !!          cout<<"undefined growth type, check .dat file"<<endl;
  !!          exit(1);
  !!        }
  !!      growthprob_phi(area, spp, yr)(Nsizebins) = 0.0; //set prob of outgrowing highest bin to 0
  !!      double tempmax =  max(growthprob_phi(area, spp, yr));
  !!      phimax = max(tempmax,phimax);
  !!	  }
  !!	}
  !!    growthprob_phi(area) /= phimax;  //rescale so no group has >1 prob growing out
  !!  }

  !!//  growthprob_phi /= phimax;   //rescale so no group has >1 prob growing out--not working on 4d array
  !!  Nstepsyr = round(phimax);            //set model timestep to phimax
  !!  Tottimesteps = Nstepsyr*Nyrs;        //for scaling state variable arrays

//intake parameters from .dat file
  init_matrix intake_alpha(1,Nareas,1,Nspecies)
  init_matrix intake_beta(1,Nareas,1,Nspecies)
  4darray intake(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!      for(yr=1; yr<=Nyrs; yr++){
  !!        intake(area, spp, yr) = 24.0 * (intake_alpha (area, spp) *
  !!                              mfexp(intake_beta (area, spp) * obs_temp (area,yr))) *
  !!                              mean_stomwt(area, spp,yr) * //daily intake in g
  !!                              365.0 / //annual intake
  !!                              Nstepsyr;  //intake per model timestep
  !!      }
  !!    }
  !!  }

////natural mortality parameters from .dat file and calculate weight ratio, size preference, suitability
//  init_3darray M1ann(1,Nareas,1,Nspecies,1,Nsizebins)
//  3darray M1(1,Nareas,1,Nspecies,1,Nsizebins) 
//  !!  for (area=1; area<=Nareas; area++){
//  !!	  for(spp=1; spp<=Nspecies; spp++){
//  !!          M1(area, spp)  = 1.0 - pow((1.0 - M1ann(area, spp)), (1.0 / Nstepsyr)) ; //scale for steps per year to equal annual input from dat  
////  !!!!            M1(area, spp) = M1ann(area, spp)/Nstepsyr;
//  !!    }
//  !!  }
    
  init_3darray isprey(1,Nareas,1,Nspecies,1,Nspecies)    //preds in columns, prey in rows
  int Npred 
  int Npreypar
  !! Npreypar = 0;
  !! Npred = 0;
  !! Npreypar = sum(isprey(1));
  !! for (spp=1;spp<=Nspecies;spp++) {
  !!  if(sum(extract_column(isprey(1),spp))>0) Npred += 1;
  !! }
  
  init_matrix preferred_wtratio(1,Nareas,1,Nspecies)     //pred specific, not size

  init_vector sd_sizepref(1,Nspecies)              //pred specific, not size
  4darray wtratio(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins)  //2nd dim pred spp, 3rd all spp as prey lengths, 4th pred lengths
  !!  for (area=1; area<=Nareas; area++){
  !!  	for (pred=1; pred<=Nspecies; pred++){
  !!    	for(prey=1; prey<=Nspecies; prey++){
  !!			dmatrix wttemp = outer_prod(wtatlbinmidpt(prey), 1.0/wtatlbinmidpt(pred));// ijth =  prey size i/pred size j
  !!			wttemp.rowshift(prey*Nsizebins-(Nsizebins-1));
  !!                    // since wttemp is inserted into a larger matrix you need to set the .rowmin() property to the row it will be inseted into
  !!                  wtratio(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins) = wttemp;
  !!		}
  !!    }
  !!  } //

  4darray sizepref(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins) //2nd dim pred spp, 3rd all spp as prey lengths, 4th pred lengths
  !!  // log normal distribution mu's and sigmas, read in. wtratioprey variable.
  !!  for (area=1; area<=Nareas; area++){
  !!  	for (pred=1; pred<=Nspecies; pred++){
  !!    	for(prey=1; prey<=Totsizebins; prey++){
  !!     	    for(int isize=1; isize<=Nsizebins; isize++){
  !!              double wtratioprey = wtratio(area, pred, prey, isize);
  !!			  sizepref(area, pred, prey, isize) =
  !!              1/(wtratioprey*sd_sizepref(pred)*sqrt(2*M_PI))*exp(-square(log(wtratioprey)-preferred_wtratio(area, pred))/(2*square(sd_sizepref(pred))));
  !!
  !!
  !!              }
  !!		}
  !!	}
  !!  } //ok
// GF - moved to procedure section
//  4darray suitability(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins)
//  !!  //suitability = sizepref * isprey
//  !!  for (area=1; area<=Nareas; area++){
//  !!  	for (pred=1; pred<=Nspecies; pred++){
//  !!    	for(prey=1; prey<=Nspecies; prey++){
//  !!                    double sumOfSizePrefs = sum(sizepref(area,pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins));
//  !!			suitability(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins) =
//  !!						sizepref(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins)*
//  !!					//	(sizepref(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins)/sumOfSizePrefs)* // normalize rates
//  !!						isprey(area, prey, pred);
//  !!             }
//  !!	}
//  !!  }

// GF - moved to parameter section
//  //fishery selectivity pars from dat file, for now not area specific
//  init_matrix fishsel_c(1,Nspecies,1,Nfleets)  //fishery selectivity c par
//  init_matrix fishsel_d(1,Nspecies,1,Nfleets)  //fishery selectivity d par

  // All inputs below have been added by andyBeet
  init_matrix B0(1,Nareas,1,Nspecies) // Equilibrium biomass. Obtained from a baseline run with zero fishing effort and no Errors added(recruitment, survey, catch)
  init_int Nguilds // number of guilds
  imatrix catchToDiscardsGuild(1,Nareas,1,Nguilds) // binary vector indicating which guilds have dropped below the most severe threshold
  imatrix catchToDiscardsSpecies(1,Nareas,1,Nspecies) // binary vector indicating which species have dropped below the most severe threshold

  imatrix maxGuildThreshold(1,Nareas,1,Nguilds) // most severe exceedence for each guild
  imatrix maxSpeciesThreshold(1,Nareas,1,Nspecies) // most severe exceedence for each species. currently a binary response
  init_ivector guildMembers(1,Nspecies) // assign each species to a guild. 1,2,3 etc.
  init_ivector fleetMembers(1,Nguilds) // assign each guild to a fleet (1,2,...,Nfleets). A fleet that predominantly catches the guild
  init_int AssessmentPeriod // time (yrs) when we assess guildlevel biomass levels

  matrix B0_guilds(1,Nareas,1,Nguilds) // equilibrium biomass for guild.
  // calculates the unfished equilibrium biomass of guild
  !! for (area=1; area<=Nareas; area++) {
  !!     for (iguild=1; iguild<=Nguilds; iguild++ ) {
  !!          for (spp=1; spp<=Nspecies; spp++) {
  !!               if (guildMembers(spp)== iguild) {
  !!                  // sum up equilibr biomass for each guild
  !!                  B0_guilds(area,iguild) += B0(area,spp);
  !!
  !!              }
  !!          }
  !!      }
  !! }
  init_int flagRamp // flag for type of response function. 0 = step function , 1 = linear ramp
  init_vector minExploitation(1,Nfleets) //min exploitation for each fleet
  init_vector maxExploitation(1,Nfleets) //max Exploitation for each fleet
  init_vector minMaxExploitation(1,2) //[MinExploitation,MaxExploitation
  init_vector minMaxThreshold(1,2) // MinThreshold,MaxThreshold]
  init_int Nthresholds // number of thresholds used for corrective fishing rules

  init_vector threshold_proportion(1,Nthresholds) //levels at which action is taken
  init_vector exploitation_levels(1,Nthresholds) //levels to drop exploitation to if threshold is exceeded
  init_vector threshold_species(1,Nspecies) // individual species thresholds (fraction)
  init_int AssessmentOn // binary, yes or no
  init_int speciesDetection // binary yes or no. Determins if species level should influence exploitation rate change during assessment

  init_int LFI_size // determins the size deemed a large fish. Used in LFI calc_health_indices module
  init_number scaleInitialN // scale the initial values of yr1N individuals
  //init_number otherFood // amount of other food included in the M2 term. previously hard coded  //GF moving to parameter section
  init_matrix effortScaled(1,Nareas,1,Nspecies) // species specific scaling of effort due to S-R function data. eg. herring/mackerel use NE effort, others GB effort
  init_4darray discard_Coef(1,Nareas,1,Nspecies,1,Nfleets,1,Nsizebins) // proportion of fleet catch that are discarded
  init_4darray discardSurvival_Coef(1,Nareas,1,Nspecies,1,Nfleets,1,Nsizebins) // of the discards what proportion survives

  init_vector predOrPrey(1,Nspecies) // binary vector indicating which species are considered predators. used in pred:prey indices
  init_int bandwidth_metric // moving window for catch variance
  init_number baseline_threshold // value of threshold that we stop landing catch. Typically 0.2
 //  !!  baseline_threshold = 0.2;
  // GF March 2022 - modifying fishing for estimation
  init_3darray indicator_fishery_q(1,Nareas,1,Nfleets,1,Nspecies) // used to determin which species used to calculate updated effort under assessment
  //!! cout << indicator_fishery_q << endl;
  int Nqpars
  !! Nqpars = sum(indicator_fishery_q)-(Nfleets*Nareas);
  imatrix f_map(1,Nareas,1,Nfleets)  //primary species for each fleet (i.e. which species the F refers to)
  int dim2
  !! dim2 = Nqpars;
  !! if (dim2 == 0) dim2 = 1;
  imatrix q_map(1,dim2,1,3)  //object that maps the catchability parameters to area, species, and fleet
  int dum
  !! f_map = 0;
  !! q_map = 0;
  !! dum = 0;
  !! for(int area=1;area<=Nareas;area++) {
  !!   for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
  !!     for (int species=1;species<=Nspecies;species++) {
  !!       if (Nqpars >0) {
  !!       if (f_map(area,ifleet)!=0 && indicator_fishery_q(area,ifleet,species) == 1) {
  !!        dum += 1;
  !!        q_map(dum,1) = area;
  !!        q_map(dum,2) = species;
  !!        q_map(dum,3) = ifleet;
  !!       }
  !!       }
  !!       if (f_map(area,ifleet)==0 && indicator_fishery_q(area,ifleet,species) == 1) f_map(area,ifleet) = species;
  !!     }
  !! 
  !!   } 
  !! }
  !! cout << "q par map" << endl;
  !! cout << Nqpars << endl;
  !! cout << f_map << endl;
  !! cout << q_map << endl;
  //!! exit(-1);

  init_vector AR_parameters(1,3) // rho (autoregressive parameters) for survey, recruitment, catch
  number rho_AR_Survey
  !! rho_AR_Survey =  AR_parameters(1);
  number rho_AR_Recruitment
  !! rho_AR_Recruitment = AR_parameters(2);
  number rho_AR_Catch
  !! rho_AR_Catch = AR_parameters(3);
  3darray sim_survey_error(1,Nareas,1,Nspecies,1,Nyrs) // used to calculate AR process for survey
  3darray sim_recruit_error(1,Nareas,1,Nspecies,1,Nyrs) // used to calculate AR process for recruitment
  4darray sim_catch_error(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs) // used to calculate AR process for catch
  3darray sim_extreme_recruit_error(1,Nareas,1,Nspecies,1,Nyrs) // used to calculate errors for extreme event

  init_int flagMSE //flag to determine output created. If MSE = 1, otherwise = 0
  init_matrix residentTime(1,Nareas,1,Nspecies) // proportion of time each species spent in management area
    // .10 1 .17 1 1 1 1 .14 1 1
  init_matrix areaMortality(1,Nareas,1,Nspecies) // total mortality of pop outside management area
  // 0 0 0 0 0 0 0 0 0 0  

    //survey obs error
  // init_matrix ln_surv_sigma(1,Nareas,1,Nspecies)
  // matrix surv_sigma(1,Nareas,1,Nspecies)
  // 3darray surv_obsError(1,Nareas,1,Nspecies,1,Nyrs)

  // //catch obs error
  // init_3darray ln_catch_sigma(1,Nareas,1,Nspecies,1,Nfleets)
  // 3darray catch_sigma(1,Nareas,1,Nspecies,1,Nfleets)
  // 4darray catch_obsError(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs)

  
  // time varying mortality (fixed multiplier for prey items)
  init_int m1_change_yr  // yr to implement change, if -99 no change applied
  init_number m1_multiplier // value to scale M1 by

  // time varying other food (fixed multiplier, same for all)
  init_int of_change_yr  // yr to implement change, if -99 no change applied
  init_number of_multiplier  // value to scale oF by

  // whether to have varying q or a fixed input
  init_int flag_tvar_q  // 0 = fixed value over time, 1 = annual values
  int Nqpar_vec 
  !! Nqpar_vec = 1;
  !! if(flag_tvar_q==1) Nqpar_vec = Nyrs;



 //flag marking end of file for data input
  init_int eof;

  //change to read in raw data
  !!ad_comm::change_datafile_name(datfilename);


  //init_3darray obs_survey_biomass(1,Nsurveys,1,Nspecies,1,Nyrs)   //input spring, fall separately, in tons
  //init_int Nsurvey_obs  //number of survey observations
  init_int Nsurvey_obs;
  !!cout << Nsurvey_obs << endl;
  init_matrix obs_survey_biomass(1,Nsurvey_obs,1,5)   //GF May 2021 revised structure

  //read in survey size comp data
  init_int Nsurvey_size_obs;
  !! int ncol = Nsizebins+5;
  init_matrix obs_survey_size(1,Nsurvey_size_obs,1,ncol)   //legnth comps for surveys

  // GF May 2021 - alternate data structure
  init_int Ncatch_obs  //number of catch observations
  init_matrix obs_catch_biomass(1,Ncatch_obs,1,6)   //total catch in tons
  //init_3darray obs_catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)    //total catch in tons
  init_int Ncatch_size_obs
  !! ncol = Nsizebins+6;
  init_matrix obs_catch_size(1,Ncatch_size_obs,1,ncol)

  //read in diet proportion data
  init_int Ndietprop_obs;
  !! ncol = Nspecies+6;
  init_matrix obs_dietprop(1,Ndietprop_obs,1,ncol)   //diet proportions by weight
  //!! cout << obs_dietprop << endl;



//debugging section, check inputs and initial calculations
	LOCAL_CALCS

  if (debug == 1)
    {
    cout<<"Nyrs\n"<<Nyrs<<endl;
    cout<<"Nspecies\n"<<Nspecies<<endl;
    cout<<"Nsizebins\n"<<Nsizebins<<endl;
    cout<<"Nareas\n"<<Nareas<<endl;
    cout<<"Nfleets\n"<<Nfleets<<endl;
    cout<<"wtconv\n"<<wtconv<<endl;
    cout<<"Totsizebins\n"<<Totsizebins<<endl;
    cout<<"binwidth\n"<<binwidth<<endl;
    cout<<"lenwt_a\n"<<lenwt_a<<endl;
    cout<<"lenwt_b\n"<<lenwt_b<<endl;
    cout<<"lbinmax\n"<<lbinmax<<endl;
    cout<<"lbinmin\n"<<lbinmin<<endl;
    cout<<"lbinmidpt\n"<<lbinmidpt<<endl;
    cout<<"wtbinmax\n"<<wtbinmax<<endl;
    cout<<"wtbinmin\n"<<wtbinmin<<endl;
    cout<<"wtatlbinmidpt\n"<<wtatlbinmidpt<<endl;
    cout<<"binavgwt\n"<<binavgwt<<endl;
    cout<<"Nrecruitment_cov\n"<<Nrecruitment_cov<<endl;
    cout<<"Nmaturity_cov\n"<<Nmaturity_cov<<endl;
    cout<<"Ngrowth_cov\n"<<Ngrowth_cov<<endl;
    cout<<"recruitment_cov\n"<<recruitment_cov<<endl;
    cout<<"maturity_cov\n"<<maturity_cov<<endl;
    cout<<"growth_cov\n"<<growth_cov<<endl;
    cout<<"obs_survey_biomass\n"<<obs_survey_biomass<<endl;
    cout<<"obs_survey_size\n"<<obs_survey_size<<endl;
    cout<<"obs_catch_biomass\n"<<obs_catch_biomass<<endl;
    cout<<"obs_catch_size\n"<<obs_catch_size<<endl;    
    cout<<"obs_dietprop\n"<<obs_dietprop<<endl;        
    cout<<"mean_stomwt\n"<<mean_stomwt<<endl;
    cout<<"obs_temp\n"<<obs_temp<<endl;
    cout<<"recruitment_covwt\n"<<recruitment_covwt<<endl;
    cout<<"rectype\n"<<rectype<<endl;
    cout<<"stochrec\n"<<stochrec<<endl;
    cout<<"rec_alpha\n"<<rec_alpha<<endl;
    cout<<"rec_shape\n"<<rec_shape<<endl;
    cout<<"rec_beta\n"<<rec_beta<<endl;
    cout<<"fecund_d\n"<<fecund_d<<endl;
    cout<<"fecund_h\n"<<fecund_h<<endl;
    cout<<"fecund_theta\n"<<fecund_theta<<endl;
    cout<<"fecundity\n"<<fecundity<<endl;
    cout<<"maturity_nu\n"<<maturity_nu<<endl;
    cout<<"maturity_omega\n"<<maturity_omega<<endl;
    cout<<"maturity_covwt\n"<<maturity_covwt<<endl;
    cout<<"growth_psi\n"<<growth_psi<<endl;
    cout<<"growth_kappa\n"<<growth_kappa<<endl;
    cout<<"growth_covwt\n"<<growth_covwt<<endl;
    cout<<"vonB_Linf\n"<<vonB_Linf<<endl;
    cout<<"vonB_k\n"<<vonB_k<<endl;
    cout<<"growthtype (1 power, 2 power/covariates, 3 vonB, 4 vonB covariates)\n"<<growthtype<<endl;
    cout<<"growthprob_phi\n"<<growthprob_phi<<endl;
    cout<<"phimax\n"<<phimax<<endl;
    cout<<"Nstepsyr\n"<<Nstepsyr<<endl;
    cout<<"Tottimesteps\n"<<Tottimesteps<<endl;
    cout<<"intake_alpha\n"<<intake_alpha<<endl;
    cout<<"intake_beta\n"<<intake_beta<<endl;
    cout<<"intake\n"<<intake<<endl;
   // cout<<"M1ann\n"<<M1ann<<endl;
   // cout<<"M1\n"<<M1<<endl;
    cout<<"isprey\n"<<isprey<<endl;
    cout<<"Npred\n"<<Npred<<endl;
    cout<<"preferred_wtratio\n"<<preferred_wtratio<<endl;
    cout<<"sd_sizepref\n"<<sd_sizepref<<endl;
    cout<<"wtratio\n"<<wtratio<<endl;
    cout<<"sizepref\n"<<sizepref<<endl;
    //cout<<"suitability\n"<<suitability<<endl;
    cout<<"B0\n"<<B0<<endl;
    cout<<"Nguilds\n"<<Nguilds<<endl;
    cout<<"guildMembers\n"<<guildMembers<<endl;
    cout<<"AssessmentPeriod\n"<<AssessmentPeriod<<endl;

//    cout<<setprecision(10)<<"FH\n"<<FH<<endl;
//    cout<<setprecision(10)<<"FHideal\n"<<FHideal<<endl;
    cout<<"eof\n"<<eof<<endl;
    }

  if(eof != 54321) {cout<<"Stop, data input failed"<<endl<<"eof: "<<eof<<endl; exit(1);}

  if (debug == 1) {cout<<"\nManually exiting at end of data section..."<<endl;  exit(-1);}

	END_CALCS

//=======================================================================================
INITIALIZATION_SECTION
//=======================================================================================

//=======================================================================================
PARAMETER_SECTION
//=======================================================================================
  
  matrix effort_updated(1,Nareas,1,Nfleets)     // calculated new Effort for each fleet following threshold exceedence calc_assessment_strategy
  3darray obs_effortAssess(1,Nareas,1,Nfleets,1,Nyrs)  	//standardized effort units needed
  !!obs_effortAssess = obs_effort;
  //Initial N year 1
  init_3darray ln_yr1N(1,Nareas,1,Nspecies,1,Nsizebins,yr1Nphase)       //initial year N at size, millions
  3darray yr1N(1,Nareas,1,Nspecies,1,Nsizebins);       //initial year N at size, millions


 // recruitment parameters all read in from Dat file. These are redundant. Andy Beet
 //recruitment parameters (alts in .dat file read in with switch for rec form by spp)
  init_matrix recruitment_alpha(1,Nareas,1,Nspecies,recphase)			//recruitment model alpha
  init_matrix recruitment_shape(1,Nareas,1,Nspecies,recphase)			//recruitment model shape parameter
  init_matrix recruitment_beta(1,Nareas,1,Nspecies,recphase)			//recruitment model beta

  //proportion mature
  4darray propmature(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) //from maturity pars and covs

  //egg production
  3darray eggprod(1,Nareas,1,Nspecies,1,Nyrs) //from fecundity, propmature, sexratio, N, in millions

  //recruitment: average annual, annual devs, actual (avg+dev)
  init_matrix ln_avg_recruitment(1,Nareas,1,Nspecies,avg_rec_phase)  //average annual recruitment by area, species
  matrix avg_recruitment(1,Nareas,1,Nspecies)  //average annual recruitment by area, species
  init_3darray recruitment_devs(1,Nareas,1,Nspecies,2,Nyrs,dev_rec_phase)  //recruitment deviations by area, species

  3darray recruitment(1,Nareas,1,Nspecies,2,Nyrs)  //by definition into first size bin for each species, millions

  //recruitment simulation
  init_matrix ln_recsigma(1,Nareas,1,Nspecies,recsigmaphase)    //sigma for stochastic recruitment from SR curve
  matrix recsigma(1,Nareas,1,Nspecies)    //sigma for stochastic recruitment from SR curve

  3darray rec_procError(1,Nareas,1,Nspecies,1,Nyrs)   //to generate deviations from SR curve

  //growth options
  //independent of bioenergetics--age based, use growthprob_phi above based on pars of age predicting length
  //4darray length(1,Nareas,1,Nspecies,1,Nages,1,Nyrs) //not needed in basic calculations so leave aside for now
  //bioenergetics based, depends on consumption--to be added

  //fishing mort: average annual, annual devs, actual (avg+dev)
  //**needs to be done by fleet, currently assuming each fleet has the same avg_F, F_devs and Ftot**
  //**then they sum to give F by species**
  //init_3darray avg_F(1,Nareas,1,Nspecies,1,Nfleets,avg_F_phase)  //logspace average annual fishing mort by area, species
  init_matrix avg_F(1,Nareas,1,Nfleets,avg_F_phase)  //logspace average annual fishing mort by area, species
  //init_3darray F_devs(1,Nspecies,1,Nfleets,1,Nyrs,dev_F_phase)  //logspace F deviations by area, species--NEEDS TO BE 4D, CANT DO?, FIX
  init_3darray F_devs(1,Nareas,1,Nfleets,1,Nyrs,dev_F_phase)  //logspace F deviations by area, species--NEEDS TO BE 4D, CANT DO?, FIX
  //
  //***********June 2014 replace with F = q*E formulation*****************************
  4darray Fyr(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs)  //array to get annual Fs by fleet from either avg/devs or q*effort, logspace

  4darray suitpreybio(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins);  //suitable prey for each predator size and year, see weight unit above for units
  
  !!cout << "Nprey " << Nprey << endl;
  init_vector logit_vuln(1,Npreypar,vuln_phase);
  3darray vulnerability(1,Nareas,1,Nspecies,1,Nspecies);
  4darray suitability(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins); //GF moved here 01/09/2023 because vulnerabilities shifted

  5darray biomass_prey_avail_no_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins,1,Nprey);

  //N, B, F, Z, M2, C, need total prey consumed per pred, suitability, available prey, food habits?
  4darray N(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //numbers by area of total pop, species, size,  timestep , in millions
  4darray Narea(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //numbers by area of pop in area, species, size,  timestep , in millions
  4darray Nnotarea(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //numbers by area of pop outside, species, size,  timestep , in millions
  4darray B(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //biomass by area, species, size,  timestep , see weight unit above for units
  4darray F(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //fishing mort by area, species, size,  timestep
  4darray C(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //catch numbers by area, species, size, timestep , in millions
  4darray Z(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //total mort by area, species, size,  timestep
  4darray M2(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //predation mort by area, species, size,  timestep
  4darray eatN(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //number eaten by predators by area, species, size,  timestep
  4darray discardN(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //number discarded by vessels by area, species, size,  timestep
  4darray otherDead(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //number unknown deaths by area, species, size,  timestep
  4darray D(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //discard mortality by area, species, size,  timestep

  //Fishery selectivities, fleet F and catch, fishery q
  init_matrix fishsel_pars(1,2,1,Nfleets,fsphase)
  matrix fishsel_c(1,Nspecies,1,Nfleets)  //fishery selectivity c par
  matrix fishsel_d(1,Nspecies,1,Nfleets)  //fishery selectivity d par
  4darray fishsel(1,Nareas,1,Nspecies,1,Nfleets,1,Nsizebins)  //fishery selectivity
  5darray Ffl(1,Nareas,1,Nspecies,1,Nfleets,1,Tottimesteps,1,Nsizebins) //fleet specific  Fs
  5darray Dfl(1,Nareas,1,Nspecies,1,Nfleets,1,Tottimesteps,1,Nsizebins) //fleet specific Discard mortaliy s
  5darray Cfl(1,Nareas,1,Nspecies,1,Nfleets,1,Tottimesteps,1,Nsizebins) //fleet specific Catch in numbers
  //init_matrix ln_fishery_q(1,Nyrs,1,Nqpars,fqphase) //Nareas,1,Nspecies,1,Nfleets,fqphase)
  init_matrix ln_fishery_q(1,Nqpar_vec,1,Nqpars,fqphase) //Nareas,1,Nspecies,1,Nfleets,fqphase)
  4darray fishery_q(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs)
  3darray  mean_guild_fishery_q(1,Nareas,1,Nguilds,1,Nfleets) // mean q for guild and fleet// andybeet
//  matrix  mean_fishery_q(1,Nareas,1,Nfleets) // mean q for fleet. ignore values of zero //andybeet

  // gavinfay March 2022 - moving this code to procedure section as function of estimated parameters
//  // calculates the mean fishery_q for each guild (over fleets)
//  !! for (area=1; area<=Nareas; area++) {
//  !!     for (iguild=1; iguild<=Nguilds; iguild++ ) {
//  !!           for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
//  !!             int icount = 0;
//  !!               for (spp=1; spp<=Nspecies; spp++) {
//  !!                 if (guildMembers(spp)== iguild) {
//  !!                    icount++;
//  !!                    // sum up q's
//  !!                    mean_guild_fishery_q(area,iguild,ifleet) += fishery_q(area,spp,ifleet);
//  !!                 }
//  !!              }
//  !!                    mean_guild_fishery_q(area,iguild,ifleet) =  mean_guild_fishery_q(area,iguild,ifleet)/icount;
//  !!          }
//  !!      }
//  !! }
//
//
//  // calculates the mean q for each fleet ignoring zero q's.
//  // this is used to update effort when an assessment dictates
//  !! for (area=1; area<=Nareas; area++) {
//  !!           for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
//  !!             int icount = 0;
//  !!             mean_fishery_q(area,ifleet) = 0;
//  !!               for (spp=1; spp<=Nspecies; spp++) {
//  !!                    if (fishery_q(area,spp,ifleet) < 1e-29) {
//  !!                      //ignore
//  !!                    } else {
//  !!                       icount = icount + indicator_fishery_q(area,spp,ifleet);
//  !!                       mean_fishery_q(area,ifleet) += indicator_fishery_q(area,spp,ifleet)*fishery_q(area,spp,ifleet);
//  !!                    }
//  !!              }
//  !!              if (icount == 0) { // then all q's are < 1-e29. This occurs during testing a new fleet with no information
//  !!                 mean_fishery_q(area,ifleet) = 0;
//  !!              } else {
//  !!                 mean_fishery_q(area,ifleet) =  mean_fishery_q(area,ifleet)/icount;
//  !!              }
//  !!          }
//  !! }


  //Survey qs 
  init_bounded_matrix ln_survey_q(1,Nsurveys,1,Nspecies,-30,2,sqphase)
  matrix survey_q(1,Nsurveys,1,Nspecies)

  //Survey selectivity (will want to be derived based on estimated parameters)
  init_matrix survey_selpars(1,2,1,Nsurveys,ssphase)
  3darray survey_sel(1,Nsurveys,1,Nspecies,1,Nsizebins)

  // //survey obs error
  // init_matrix ln_surv_sigma(1,Nareas,1,Nspecies,ssig_phase)
  // matrix surv_sigma(1,Nareas,1,Nspecies)
  // 3darray surv_obsError(1,Nareas,1,Nspecies,1,Nyrs)

  // //catch obs error
  // init_3darray ln_catch_sigma(1,Nareas,1,Nspecies,1,Nfleets,csig_phase)
  // 3darray catch_sigma(1,Nareas,1,Nspecies,1,Nfleets)
  // 4darray catch_obsError(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs)

  //annual total B and SSB
  3darray avByr(1,Nareas,1,Nspecies,1,Nyrs) //uses B units--"true" simulated biomass
  3darray SSB(1,Nareas,1,Nspecies,1,Nyrs) //uses B units--"true" SSB

  //annual eaten B
  3darray eaten_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units. M2
  3darray discard_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units. discards
  3darray otherDead_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units. M1
  3darray total_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units. N
  4darray eaten_biomass_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)  //uses B units. M2
  4darray discard_biomass_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)  //uses B units. discards
  4darray otherDead_biomass_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)  //uses B units. M1
  4darray total_biomass_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)  //uses B units. N
  4darray catch_biomass_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)  //uses B units--"true" simulated catch
  3darray predation_mortality(1,Nareas,1,Nspecies,1,Nyrs) // uses B units = eaten/total biomass
  3darray fishing_mortality(1,Nareas,1,Nspecies,1,Nyrs) // uses B units  = (landing+discards) / total biomass
  4darray predation_mortality_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) // uses B units = eaten/total biomass
  4darray fishing_mortality_size(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) // uses B units  = (landing+discards) / total biomass

  //estimated fishery catch and survey catch
  3darray catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units--"true" simulated catch
  4darray fleet_catch_biomass(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs) //B units--"true" catch by fleet
  //4darray est_fleet_catch_biomass(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs) //B units, can have obs error and q
 // 3darray est_catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units, sums over est_fleet_catch_biomass
  //4darray est_fleet_catch_guild_biomass(1,Nareas,1,Nguilds,1,Nfleets,1,Nyrs)// uses B units. sums over species. calc_catc_etx
  //4darray est_fleet_catch_guild_assessment(1,Nareas,1,Nguilds,1,Nfleets,1,Nyrs)// moving average of catch over  NumAssessment years
  //3darray est_catch_guild_biomass(1,Nareas,1,Nguilds,1,Nyrs) // sums over species and fllet for guild
  // 3darray est_survey_biomass_assessment(1,Nareas,1,Nspecies,1,Nyrs)// moving average for each species over NumAssessment years

  // 3darray est_survey_biomass(1,Nareas,1,Nspecies,1,Nyrs) //uses B units, can have q, obs error
  // 3darray est_survey_guild_biomass(1,Nareas,1,Nguilds,1,Nyrs) //uses B units, sums over species calc_survey_abundance
  // 3darray est_survey_guild_biomass_assessment(1,Nareas,1,Nguilds,1,Nyrs)   // moving average of guild biomass over NumAssessment years

  //objective function, penalties, components of objective function (placeholder, not yet developed)
  //3darray resid_catch(1,Nareas,1,Nspecies,1,Nyrs)   //log(obs)-log(est) catch
  3darray resid_bio(1,Nareas,1,Nspecies,1,Nyrs)     //log(obs)-log(est) survey bio
  matrix totcatch_fit(1,Nareas,1,Nspecies)  //fit to total catch in weight by area and species
  matrix catchcomp_fit(1,Nareas,1,Nspecies) //fit to catch at length composition
  matrix totbio_fit(1,Nareas,1,Nspecies)    //fit to total survey biomass by area and species
  matrix biocomp_fit(1,Nareas,1,Nspecies)   //fit to survey catch at length composition
  //matrix agelencomp_fit(1,Nareas,1,Nspecies) //fit to age at length composition, where available


  //GF objective function piece placeholders - check consistency with those above and remove redundancies
  // catch observations
  vector pred_catch_biomass(1,Ncatch_obs);
  vector resid_catch(1,Ncatch_obs);
  vector nll_catch(1,Ncatch_obs);


  //vector pred_catch_biomass(1,Ncatch_obs);
  //vector resid_catch(1,Ncatch_obs);
  !! int Nsize_obs = 0;
  !! for (int i=1;i<=Ncatch_size_obs;i++) {
  !!    for (int ilen=1;ilen<=Nsizebins;ilen++) 
  !!      if (obs_catch_size(i,6+ilen)>=0) Nsize_obs += 1;
  !! }
  vector pred_catch_size(1,Nsize_obs);
  vector nll_catch_size(1,Nsize_obs);

  vector pred_survey_index(1,Nsurvey_obs);
  vector resid_survey(1,Nsurvey_obs);
  vector nll_survey(1,Nsurvey_obs);

  !! Nsize_obs = 0;
  !! for (int i=1;i<=Nsurvey_size_obs;i++) {
  !!    for (int ilen=1;ilen<=Nsizebins;ilen++) 
  !!      if (obs_survey_size(i,5+ilen)>=0) Nsize_obs += 1;
  !! }
  vector pred_survey_size(1,Nsize_obs);
  vector nll_survey_size(1,Nsize_obs);

  //4darray est_survey_size(1,Nsurvey,1,Nyrs,1,Nspecies,1,Nsizebins);


  !! Nsize_obs = 0;
  !! for (int i=1;i<=Ndietprop_obs;i++) {
  !!    for (int ilen=1;ilen<=Nspecies;ilen++) 
  !!      if (obs_dietprop(i,5+ilen)>=0) Nsize_obs += 1;
  !!    if (obs_dietprop(i,6+Nspecies)>=0) Nsize_obs += 1;
  !! }
  vector pred_dietprop(1,Nsize_obs);
  vector nll_dietprop(1,Nsize_obs);

  !! Nsize_obs = Nspecies*Nareas*Nyrs;
  vector recdev(1,Nsize_obs);
  number nll_recruit;


// calc_health_indices variables AndyBeet
  matrix index_Simpsons_N(1,Nareas,1,Nyrs); // simpsons_N. index values
  matrix index_Simpsons_Nrecip(1,Nareas,1,Nyrs); //1/simpsons for N. index values
  matrix index_Simpsons_C(1,Nareas,1,Nyrs); // simpsons for catch. index values
  matrix index_Simpsons_Crecip(1,Nareas,1,Nyrs); // 1/simpsons for catch. index values
  3darray index_LFI_Biomass(1,Nareas,1,Nspecies,1,Nyrs) // large fish defined as # in top size class. For each species
  3darray index_LFI_Catch(1,Nareas,1,Nspecies,1,Nyrs) // large fish defined as # in top size class. For each species
  3darray index_LFI_N(1,Nareas,1,Nspecies,1,Nyrs) // large fish index defines as number of fish in largest size class
  matrix LFI_threshold(1,Nareas,1,Tottimesteps) // large fish index. large fish defined as exceeding x cm (parameter read in)
  vector prob_species(1,Nspecies) //function of N
  number LF_Biomass
  vector B_total(1,Nspecies) // total B by species at each time t
  vector B_largestClass(1,Nspecies) // biomass of largest size class
  4darray N_tot(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) // acumulative N over the year
  4darray B_tot(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) // total B over year
  4darray C_tot(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) // total C over year
  5darray Cfl_tot(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs,1,Nsizebins) // total C by fleet over year
  matrix index_predBio(1,Nareas,1,Nyrs) // annual predator biomass. used in health indices module
  matrix index_preyBio(1,Nareas,1,Nyrs) // annual prey biomass. used in health indices module
  matrix index_predToPreyRatio(1,Nareas,1,Nyrs) // annual predator:prey ratio in terms of biomass. used in health indices module
  matrix index_plankToPiscRatio(1,Nareas,1,Nyrs) // annual predator:prey ratio in terms of biomass. used in health indices module
  vector index_catch(1,bandwidth_metric) // temp vector used in indices module. stores catch
  vector index_biomass(1,bandwidth_metric) // temp vector used in indices module. stores catch
  3darray index_stdev_catch(1,Nareas,1,Nspecies,1,Nyrs) /// std of catch uing a window of size = bandwidth_metric
  3darray index_stdev_biomass(1,Nareas,1,Nspecies,1,Nyrs) // std of biomass using window of size = bandwidth_metric
  3darray index_ExploitationRate(1,Nareas,1,Nspecies,1,Nyrs) // species exploitation rate
  matrix index_SystemExploitationRate(1,Nareas,1,Nyrs) // system exploitation rate
  3darray  exploitation_update(1,Nareas,1,Nfleets,1,Nyrs) // changing exploitation rate
  3darray index_status_species(1,Nareas,1,Nspecies,1,Nyrs) // species biomass < .2*B0
  3darray index_status_guild(1,Nareas,1,Nguilds,1,Nyrs) // species biomass < .2*B0
  matrix exploitationLevelSpecies(1,Nareas,1,Nspecies) // stores adjusted exploitation levels
  matrix exploitationLevelGuild(1,Nareas,1,Nguilds) // stores adjusted exploitation levels
  matrix newExploitationLevel(1,Nareas,1,Nfleets)// the new exploitation rate after an assessment
  matrix objfun_areaspp(1,Nareas,1,Nspecies) //sum over components for area and species
  3darray rec_EventError(1,Nareas,1,Nspecies,1,Nyrs) // error for extreme recruitment event

  objective_function_value objfun

  	LOCAL_CALCS
    if (debug == 2)
      {
       cout<<"rectype\n"<<rectype<<endl;
       cout<<"recruitment_alpha\n"<<recruitment_alpha<<endl;
       cout<<"recruitment_shape\n"<<recruitment_shape<<endl;
       cout<<"recruitment_beta\n"<<recruitment_beta<<endl;
      //cout<<"aAge1\n"<<aAge1<<endl<<"aFt\n"<<aFt<<endl;
      //for (i=1; i<=nsp; i++)  {
      //  cout<<"species: "<<i<<endl;
      //  cout<<"idAge1\n"<<idAge1(i)<<endl<<"idFt\n"<<idFt(i)<<endl;
      //  cout<<"iagesel\n"<<iagesel(i)<<endl<<"iFICsel\n"<<iFICsel(i)<<endl;
      //  cout<<"iYr1\n"<<iYr1(i)<<endl<<endl;
      //  }
      //cout<<"iRho\n"<<iRho<<endl;
      cout<<"\nManually exiting at the end of the parameter section...\n"<<endl;
      exit(-1);
      }
	END_CALCS

  init_matrix ln_M1ann(1,Nareas,1,Nspecies,m1_phase)
  //3darray M1(1,Nareas,1,Nspecies,1,Nsizebins)
  4darray M1(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) 

  init_number ln_otherFood_base(oF1_phase) // amount of other food included in the M2 term for the base (predator 1)
  init_vector otherFood_dev(1,Npred-1,oFdev_phase)  //deviation from base other food for predators 2+  (same other food for all size classes of each predator)
  matrix otherFood(1,Nspecies,1,Nyrs)   // vector of other food included in the M2 term



//=======================================================================================
PRELIMINARY_CALCS_SECTION
//=======================================================================================
  recruitment_alpha = rec_alpha;
  recruitment_shape = rec_shape;
  recruitment_beta = rec_beta;

  // for (int i=1;i<=Nsurveys;i++)
  //  for (int k=1;k<=Nspecies;k++)
  //  for (int j=1;j<=Nsizebins;j++)
  // // need to add survey selectivity at length, but for now set at 1 for all.
  //     survey_sel(i,k,j) = 1.;

// GF moved below out of procedure section May 09 2022
  //propmature do once for whole cov time series, estimate maturity params, covwt, or both in later phases
  //propmature = 1/(1+exp(-(maturity_nu+maturity_omega*lbinmidpt)*sum(maturity_covwt*maturity_cov)))
  // first calculate the covarate part to add. andybeet
  for (spp=1; spp<=Nspecies; spp++) {
      for (yr=1; yr<=Nyrs; yr++) {
          for (int icov=1; icov<=Nmaturity_cov;icov++) {
              covariates_M(spp,yr) += maturity_covwt(spp,icov)*maturity_cov(icov,yr);
          }
      }
  }




//=======================================================================================
PROCEDURE_SECTION
//=======================================================================================

  transform_parameters(); if (debug == 4) {cout<<"completed parameter transform"<<endl;}

  for (t=1;t<=Tottimesteps;t++){
  for (area=1; area<=Nareas; area++){
   for(spp=1; spp<=Nspecies; spp++){
//    //      M1(area, spp)  = 1.0 - pow((1.0 - M1ann(area, spp)), (1.0 / Nstepsyr)) ; //scale for steps per year to equal annual input from dat  
          M1(area, spp, t)  = 1.0 - pow((1.0 - mfexp(ln_M1ann(area, spp))), (1.0 / Nstepsyr)) ; //scale for steps per year to equal annual input from dat  
//    //        M1(area, spp) = mfexp(ln_M1ann(area, spp))/Nstepsyr;
    }
   }
   }

  // Other Food - fill in the other food vector
   otherFood = 0.;
   int j=0;
   for(spp=1; spp<=Nspecies; spp++){
    if (sum(extract_column(isprey(1),spp))>0) {
      if (j==0) otherFood(spp) = mfexp(ln_otherFood_base);
      if (j>=1) otherFood(spp) = mfexp(ln_otherFood_base + otherFood_dev(j));
      j = j+1;
    }
    else
    {otherFood(spp) = mfexp(ln_otherFood_base);}   //GF 11/14/22 added trap to avoid 0 other food for non-predators, need to check full implications of why this is needed [nans in M2 calcs otherwise but want to check it's not using this info and this is just a product of looping over all species]
   }
   for(spp=1; spp<=Nspecies; spp++){
    for (yr=1;yr<=Nyrs;yr++) {
     if (of_change_yr !=-99 && yr >= of_change_yr) 
       otherFood(spp,yr) = otherFood(spp,yr)*of_multiplier;;
    }
   }
    
  //cout << "other food" << endl;
  //cout << otherFood << endl;


  //suitability calcs - previously in data section
    //suitability = sizepref * vulnerability
    for (area=1; area<=Nareas; area++){
      for (pred=1; pred<=Nspecies; pred++){
        for(prey=1; prey<=Nspecies; prey++){
                      double sumOfSizePrefs = sum(sizepref(area,pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins));
        suitability(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins) =
              sizepref(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins)*
            //  (sizepref(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins)/sumOfSizePrefs)* // normalize rates
              vulnerability(area, prey, pred);
               }
    }
    }


  //ofstream popout("popstructure.out");
  //ofstream recout("recstructure.out");

//  calc_fishery_qs(); if (debug == 4) {cout<<"completed fishery qs"<<endl;}  //gavinfay March 2022 - moved from PARAMETER_SECTION 
  
  calc_initial_states();  if (debug == 4) {cout<<"completed Initial States"<<endl;}

  yrct=1;


  for (t=1; t<=Tottimesteps; t++)
     {

		//if (debug == 3) {cout<<yrct<<" "<<t<<endl;}
                if (t>1 && (Nstepsyr == 1 || t % Nstepsyr == 1)) {yrct++;} // first time step in new year = > increment year

                if (t>1 && m1_change_yr !=-99 && yrct >= m1_change_yr) calc_tv_m1();

                if (t>1) calc_update_N(); // N(t) = N(t-1)

                  // add recruits at start of year.  update N to add recruits to bin 1
                if (t>1) calc_recruitment(); if (debug == 4) {cout<<"completed Recruitment"<<endl;}
                //if (t % Nstepsyr == 1) recout << yrct << " " << recruitment(1,2,yrct) << endl;
                calc_available_N();

                calc_pred_mortality(); if (debug == 4) {cout<<"completed Predation Mortality"<<endl;}

                calc_fishing_mortality(); if (debug == 4) {cout<<"completed Fishing Mortality"<<endl;}

                calc_total_mortality(); // We calculate Z(t) = M1 + M2 + F
  
    //cout << t << " " << Z(1,1,t) << endl;
    //cout << t << " " << Z(1,1,t) << " " << M1(1,1) << " " << M2(1,1,t) << " " << F(1,1,t) << " " << D(1,1,t) << endl;

		calc_catch_etc(); if (debug == 4) {cout<<"completed Catch"<<endl;} // split F among fleets
    //for (int spp=1;spp<=Nspecies;spp++)
		//popout << t << " " << spp << " " << N(1,spp,t) << " " << Z(1,spp,t) << " " << M1(1,spp) << " " << M2(1,spp,t) << " " << F(1,spp,t) << " " << D(1,spp,t) << endl;
    calc_pop_dynamics(); if (debug == 4) {cout<<"completed Pop Dynamics"<<endl;} // update N - death + growth
    
                calc_SSB();

		calc_movement(); if (debug == 4) {cout<<"completed Movement"<<endl;}

                calc_survey_abundance();  if (debug == 4) {cout<<"completed Survey Abundance"<<endl;}

//   GF commented out below code during debugging 04/28/2022
//                calc_health_indices();  if (debug == 4) {cout<<"completed Survey Abundance"<<endl;}

                // // enter assessment module if turned on in data file, if end of year, if curent year is a multiple of assessmentPeriod
                // if (AssessmentOn == 1) {
                //  // if end of year and enough years have passed to perform assessment.
                //  // every AssessmentPeriod we monitor stocks and adjust the effort for the future
                //   if ((t % Nstepsyr == 0) && (yrct <= (Nyrs-AssessmentPeriod))) {
                //    if( yrct % AssessmentPeriod == 0) {

                //     if (flagRamp == 0) {    // step function
                //      calc_assessment_strategy_Step(); if (debug == 4) {cout<<"completed calc_assessment_strategy"<<endl;}
                //     } else { // linear function
                //      calc_assessment_linear_independent_fleet(); if (debug == 4) {cout<<"completed calc_assessment_strategy"<<endl;}
                //     }
                //    }
                //   }
                // }

	 }

  if (debug == 4) {cout<<"completed timestep loop"<<endl;}

  calculate_predicted_values(); {cout<<"completed predicted values for surveys"<<endl;}

  evaluate_the_objective_function(); if (debug == 4) {cout<<"completed Log Likelihood"<<endl;}

//   GF commented out below code during debugging 04/28/2022
//  ///////////////////////////////////// OUTPUT //////////////////////////////////

//     // write out indices to a file

//    if (debug == 3 && flagMSE  == 1) // MSE runs only. Limited output
//     {
//       write_outIndices(); // MSE output
//       exit(0);
//     }
//    if (debug == 3 && flagMSE  == 2) // MSE darwin runs only. biomass and catch output only
//     {
//       write_outDarwin(); // MSE output
//       exit(0);
//     }
//    if (debug == 3 && flagMSE == 0) // All diagnostic plots
//     {
//       write_simout_KRAKEN(); // kraken output
//       write_outIndices(); // MSE output
//       write_outDiagnostics(); // diagnostic plot output
//       exit(0);
//     }
// //////////////////////////////////////////////////////////////////////////////

// //----------------------------------------------------------------------------------------
// FUNCTION calc_fishery_qs
// //----------------------------------------------------------------------------------------

//   // calculates the mean fishery_q for each guild (over fleets)
//    for (area=1; area<=Nareas; area++) {
//      for (iguild=1; iguild<=Nguilds; iguild++ ) {
//            for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
//              int icount = 0;
//                for (spp=1; spp<=Nspecies; spp++) {
//                  if (guildMembers(spp)== iguild) {
//                     icount++;
//                     // sum up q's
//                     mean_guild_fishery_q(area,iguild,ifleet) += fishery_q(area,spp,ifleet);
//                  }
//               }
//                     mean_guild_fishery_q(area,iguild,ifleet) =  mean_guild_fishery_q(area,iguild,ifleet)/icount;
//           }
//       }
//  }


 //  // calculates the mean q for each fleet ignoring zero q's.
 //  // this is used to update effort when an assessment dictates
 // for (area=1; area<=Nareas; area++) {
 //           for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
 //             int icount = 0;
 //             mean_fishery_q(area,ifleet) = 0;
 //               for (spp=1; spp<=Nspecies; spp++) {
 //                    if (fishery_q(area,spp,ifleet) < 1e-29) {
 //                      //ignore
 //                    } else {
 //                       icount = icount + indicator_fishery_q(area,ifleet,spp);
 //                       mean_fishery_q(area,ifleet) += indicator_fishery_q(area,ifleet,spp)*fishery_q(area,spp,ifleet);
 //                    }
 //              }
 //              if (icount == 0) { // then all q's are < 1-e29. This occurs during testing a new fleet with no information
 //                 mean_fishery_q(area,ifleet) = 0;
 //              } else {
 //                 mean_fishery_q(area,ifleet) =  mean_fishery_q(area,ifleet)/icount;
 //              }
 //          }
 // }

//----------------------------------------------------------------------------------------
FUNCTION calc_tv_m1
//----------------------------------------------------------------------------------------

  for (area=1; area<=Nareas; area++){
   //for(spp=1; spp<=Nspecies; spp++){
          M1(area, 1, t, 1) = M1(area, 1, t, 1)*m1_multiplier;
          M1(area, 2, t)  = M1(area, 2, t)*m1_multiplier;
          M1(area, 3, t)  = M1(area, 3, t)*m1_multiplier;
          M1(area, 5, t, 1) = M1(area, 5, t, 1)*m1_multiplier;
          M1(area, 5, t, 2) = M1(area, 5, t, 2)*m1_multiplier;
          M1(area, 6, t, 1) = M1(area, 6, t, 1)*m1_multiplier;
          M1(area, 6, t, 2) = M1(area, 6, t, 2)*m1_multiplier;
          M1(area, 10, t, 1) = M1(area, 10, t, 1)*m1_multiplier;
          M1(area, 10, t, 2) = M1(area, 10, t, 2)*m1_multiplier;          
   // }
   }

//----------------------------------------------------------------------------------------
FUNCTION transform_parameters
//----------------------------------------------------------------------------------------

  yr1N = mfexp(ln_yr1N);
  avg_recruitment = mfexp(ln_avg_recruitment);
  recsigma = mfexp(ln_recsigma);
  //fishery_q = mfexp(ln_fishery_q);
  fishery_q.initialize();
  // fishery catchabilities  //gavinfay March 2022
  for (yr=1;yr<=Nyrs;yr++) {
  for (area=1;area<=Nareas;area++) {
    for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
      for (int species=1;species<=Nspecies;species++) fishery_q(area,species,ifleet,yr) = 1e-15; //0.;
      fishery_q(area,f_map(area,ifleet),ifleet,yr) = 1.;
    }
   }
  if (Nqpars > 0) {
  for (int ipar=1;ipar<=Nqpars;ipar++) {
    if (flag_tvar_q==1 ) fishery_q(q_map(ipar,1),q_map(ipar,2),q_map(ipar,3),yr) = mfexp(ln_fishery_q(yr,ipar));
    if (flag_tvar_q==0 ) fishery_q(q_map(ipar,1),q_map(ipar,2),q_map(ipar,3),yr) = mfexp(ln_fishery_q(1,ipar));
  }
  }
  }
  //cout << "fishery q" << endl;
  //cout << fishery_q << endl;
  //exit(-1);

  //fishery selectivity
  for (int species=1;species<=Nspecies;species++) {
   for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
     fishsel_c(species,ifleet) = mfexp(fishsel_pars(1,ifleet));
     fishsel_d(species,ifleet) = mfexp(fishsel_pars(2,ifleet));
   }
  }
  //cout << "fishery selectivit pars" << endl;
  //cout << fishsel_c << endl;
  //cout << fishsel_d << endl;
  //exit(-1);
  
  dvariable neglog19 = -1.*log(19.);
  for (int i=1;i<=Nsurveys;i++)
   for (int k=1;k<=Nspecies;k++) {
   for (int j=1;j<=Nsizebins;j++) {
  // need to add survey selectivity at length, but for now set at 1 for all.
  //survey_sel(i,k,j) = 1/(1 + mfexp(-1.*mfexp(survey_selpars(1,i))*(lbinmidpt(k,j)-mfexp(survey_selpars(2,i)))));
  survey_sel(i,k,j) = 1./(1.+mfexp(neglog19*(lbinmidpt(k,j)-mfexp(survey_selpars(1,i)))/mfexp(survey_selpars(2,i))));
                                     }
    survey_sel(i,k) /= max(survey_sel(i,k));
    }

   survey_q = mfexp(ln_survey_q);
  // surv_sigma = mfexp(ln_surv_sigma);
  // catch_sigma = mfexp(ln_catch_sigma);

  //vulnerabilities
  int iprey = 0;
  vulnerability.initialize();
  for (area=1; area<=Nareas; area++){
    for (pred=1; pred<=Nspecies; pred++){
      for(prey=1; prey<=Nspecies; prey++){
        if (isprey(area,prey,pred)==1) {
         iprey+=1;
         vulnerability(area,prey,pred) = mfexp(logit_vuln(iprey))/(1.+mfexp(logit_vuln(iprey)));
        }
      }
    }
  }
  //cout << "vulnerability" << endl;
  //cout << logit_vuln << endl;
  //cout << vulnerability << endl;
  //exit(-1);



//----------------------------------------------------------------------------------------
FUNCTION calc_initial_states
//----------------------------------------------------------------------------------------

  propmature.initialize();
  eggprod.initialize();
  recruitment.initialize();
  rec_procError.initialize();
  suitpreybio.initialize();
  Fyr.initialize();
  fishsel.initialize(); Ffl.initialize(); Dfl.initialize(); Cfl.initialize();
  N.initialize(); B.initialize(); F.initialize(); D.initialize(); C.initialize(); Narea.initialize();Nnotarea.initialize();
  Z.initialize(); M2.initialize(); eatN.initialize(); otherDead.initialize();discardN.initialize();
  avByr.initialize(); SSB.initialize();
  Cfl_tot.initialize(); C_tot.initialize();
  N_tot.initialize(); B_tot.initialize();
  eaten_biomass.initialize();
  otherDead_biomass.initialize();
  discard_biomass.initialize();
  // surv_obsError.initialize();
  // catch_obsError.initialize();
  // est_survey_biomass.initialize(); est_catch_biomass.initialize();
  //  // andy beet
  // est_survey_guild_biomass.initialize();
  // est_fleet_catch_guild_biomass.initialize();
  // est_catch_guild_biomass.initialize();
//  est_survey_guild_biomass_assessment.initialize();
  // est_survey_biomass_assessment.initialize();
  // est_fleet_catch_guild_assessment.initialize();
  covariates_M.initialize();
  index_predBio.initialize();
  index_preyBio.initialize();
  index_predToPreyRatio.initialize();
  index_plankToPiscRatio.initialize();
  exploitation_update.initialize();
  catchToDiscardsSpecies.initialize(); catchToDiscardsGuild.initialize();
  eaten_biomass_size.initialize();
  otherDead_biomass_size.initialize();
  discard_biomass_size.initialize();
  total_biomass_size.initialize();
  catch_biomass_size.initialize();
  rec_EventError.initialize();
  //andybeet
  fleet_catch_biomass.initialize(); //est_fleet_catch_biomass.initialize();
  catch_biomass.initialize();
  index_stdev_catch.initialize();
  index_stdev_biomass.initialize();
  index_status_species.initialize();
  index_status_guild.initialize();


  //gavinfay
  biomass_prey_avail_no_size.initialize();


  totcatch_fit.initialize(); catchcomp_fit.initialize();
  totbio_fit.initialize(); biocomp_fit.initialize();

  //need year 1 N to do recruitment, pred mort, N for following years
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
         N(area, spp, 1) = yr1N(area, spp)*scaleInitialN;
         B(area, spp, 1) = wtconv*elem_prod(N(area, spp, 1),binavgwt(spp));
         avByr(area, spp, 1) = sum(B(area,spp,1))/Nstepsyr;
//         est_survey_biomass(area,spp,1) = avByr(area, spp, 1);  //perfect surveys as placeholder
      }
  }

// GF May 9 2022 moving out of procedure section as not dependent on parameters
  // //propmature do once for whole cov time series, estimate maturity params, covwt, or both in later phases
  // //propmature = 1/(1+exp(-(maturity_nu+maturity_omega*lbinmidpt)*sum(maturity_covwt*maturity_cov)))
  // // first calculate the covarate part to add. andybeet
  // for (spp=1; spp<=Nspecies; spp++) {
  //     for (yr=1; yr<=Nyrs; yr++) {
  //         for (int icov=1; icov<=Nmaturity_cov;icov++) {
  //             covariates_M(spp,yr) += maturity_covwt(spp,icov)*maturity_cov(icov,yr);
  //         }
  //     }
  // }


  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
		for(yr=1; yr<=Nyrs; yr++){
//			propmature(area, spp)(yr) = 1/(1+exp(-(maturity_nu(area, spp) +
// Andybeet                                          maturity_omega(area, spp)*lbinmidpt(spp)) +
//                                          maturity_covwt(spp)*trans(maturity_cov)(yr)));
                   for (int isizebin=1; isizebin<=Nsizebins; isizebin++) {
                    // make an exception for dogfish, (species 1). SR relationship used only females. therefore SSB should only use females.
                    // females considered to be only members of largest size class and only class that can reproduce
                    // this isn't smart coding. ideally we'd want a function that would create this and we'd just pass parameter values in dat file
                     if ((spp == 1) && (isizebin < Nsizebins)) {
                        propmature(area,spp,yr,isizebin) = 0;
                     } else if ((spp == 1) && (isizebin == Nsizebins)) {
                        propmature(area,spp,yr,isizebin) = 1;// eventually code for covariates on dogfish
                     } else {
			propmature(area,spp,yr,isizebin) = 1/(1+mfexp(-1.*(maturity_nu(area, spp) +
                                          maturity_omega(area, spp)*lbinmidpt(spp,isizebin)) +
                                          covariates_M(spp,yr)));
                    }
                   }
		}
    }
  }

  //as long as growth is fit outside the model and transition probs done once at the beginning, could move here
  
  //fill F arrays; either start with avg F and devs by fleet, or calculate from q and effort by fleet
  // effort is scaled by species depending on how S-R data was assembled. GB or region wide
  for (int iyr =1; iyr<=Nyrs; iyr++) {
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
	  for(fleet=1; fleet<=Nfleets; fleet++){
          //fill Fyr array with area, species, fleet, year specific Fs
          // Fyr(area,spp,fleet) = avg_F(area,spp,fleet) + F_devs(spp,fleet); //WARNING ONLY WORKS WITH 1 AREA:REDO
//            Fyr(area,spp,fleet) = log(fishery_q(area,spp,fleet)*obs_effort(area,fleet)); //Andy Beet
//            effordScaled redundant. we now use proportion of GB effort to shelf effort rather than assume constant scaling
//            Fyr(area,spp,fleet) = fishery_q(area,spp,fleet)*obs_effort(area,fleet)*effortScaled(area,spp); //Andy Beet
//                  Fyr(area,spp,fleet) = fishery_q(area,spp,fleet)*obs_effort(area,fleet); //Andy Beet
            Fyr(area,spp,fleet,iyr) = fishery_q(area,spp,fleet,iyr)*mfexp(avg_F(area,fleet)+F_devs(area,fleet,iyr));  //gavinfay March 2022 - modding for F
            //Fyr(area,spp,fleet,iyr) = fishery_q(area,spp,fleet)*mfexp(F_devs(area,fleet,iyr));  //gavinfay March 2022 - modding for F
            //cout << iyr << " " << area << " " << spp << " " << fleet << " " << Fyr(area,spp,fleet,iyr) << endl;
      }
    }
  }
  }
  //exit(-1);
//   cout << "HI GAVIN" << endl;
// ///////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////   Random Number Generators ////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////

//   // Add simulated observation errors for survey
//     random_number_generator rng (rseed);
//     dvector obsError(1,Nyrs);
//     for (area=1; area<=Nareas; area++){
//   	  for(spp=1; spp<=Nspecies; spp++){
//                obsError.fill_randn(rng);
//                // N(0,1) values
//                surv_obsError(area,spp) = obsError;
//       }
//     }
//     // create AR(1) error structure for survey
//     // note if AR parameter = 0 then we have white noise
//     for (area=1; area<=Nareas; area++){
//         for (spp=1; spp<=Nspecies; spp++){
//           sim_survey_error(area,spp,1) =  value(surv_sigma(area,spp)*surv_obsError(area,spp,1))*
//                               pow(1-pow(rho_AR_Survey,2),0.5)  - 0.5 *value( surv_sigma(area,spp) * surv_sigma(area,spp));

//           for (int iy=2; iy<=Nyrs; iy++) {
//               sim_survey_error(area,spp,iy) = rho_AR_Survey*sim_survey_error(area,spp,iy-1) +  value(surv_sigma(area,spp)*surv_obsError(area,spp,iy))*
//                               pow(1-pow(rho_AR_Survey,2),0.5)  - 0.5 *value( surv_sigma(area,spp) * surv_sigma(area,spp))*(1-pow(rho_AR_Survey,2))  ;
//           }
//         }
//      }

//     // create AR(1) error structure for catch (fleet specific)
//     // note if AR parameter = 0 then we have white noise
//     random_number_generator rng2 (rseed+10000);
//     dvector CobsError(1,Nyrs);
//     for (area=1; area<=Nareas; area++){
//   	  for(spp=1; spp<=Nspecies; spp++){
// 		 for(fleet=1; fleet<=Nfleets; fleet++){
//                         CobsError.fill_randn(rng2); //N(0,1) values
//                         catch_obsError(area,spp,fleet) = CobsError;
// 	         }
//           }
//     }
//     for (area=1; area<=Nareas; area++){
//      for(spp=1; spp<=Nspecies; spp++){
//       for(fleet=1; fleet<=Nfleets; fleet++){
//          sim_catch_error(area,spp,fleet,1) = value(catch_sigma(area,spp,fleet)*catch_obsError(area,spp,fleet,1))*
//                               pow(1-pow(rho_AR_Catch,2),0.5)  - 0.5 *value(catch_sigma(area,spp) * catch_sigma(area,spp));

//         for (int iy=2; iy<=Nyrs; iy++) {
//           sim_catch_error(area,spp,fleet,iy) = rho_AR_Catch*sim_catch_error(area,spp,fleet,iy-1) +  value(catch_sigma(area,spp,fleet)*catch_obsError(area,spp,fleet,iy))*
//               pow(1-pow(rho_AR_Catch,2),0.5)  - 0.5 *value( catch_sigma(area,spp) * catch_sigma(area,spp))*(1-pow(rho_AR_Catch,2))  ;
//         }
//        }
//      }
//    }



//     // Add simulated process errors for recruitment
//     // create AR(1) error structure for recruitment
//     // note if AR parameter = 0 then we have white noise
//     random_number_generator rng3 (rseed+20000);
//     dvector RprocError(1,Nyrs);
//     for (area=1; area<=Nareas; area++){
//      for(spp=1; spp<=Nspecies; spp++){
//          RprocError.fill_randn(rng3);
//          rec_procError(area,spp) = RprocError;
//      }
//     }


//     for (area=1; area<=Nareas; area++){
//         for (spp=1; spp<=Nspecies; spp++){
//           sim_recruit_error(area,spp,1) =  value(recsigma(area,spp)*rec_procError(area,spp,1))*
//                               pow(1-pow(rho_AR_Recruitment,2),0.5)  - 0.5 *value(recsigma(area,spp) * recsigma(area,spp));

//           for (int iy=2; iy<=Nyrs; iy++) {
//               sim_recruit_error(area,spp,iy) = rho_AR_Recruitment*sim_recruit_error(area,spp,iy-1) +  value(recsigma(area,spp)*rec_procError(area,spp,iy))*
//                               pow(1-pow(rho_AR_Recruitment,2),0.5)  - 0.5 *value( recsigma(area,spp) * recsigma(area,spp))*(1-pow(rho_AR_Recruitment,2))  ;
//           }
//         }
//      }

//   // p(extreme event) for recruitment
//     random_number_generator rngRecruits (rseed);
//     dvector recruitEventError(1,Nyrs);
//     for (area=1; area<=Nareas; area++){
//      for (spp=1; spp<=Nspecies; spp++){
//        recruitEventError.fill_randu(rngRecruits);
//        rec_EventError(area,spp) = recruitEventError;
//      }
//     }

    // now given the distribution of errors for extreme events we calculate the error
    //sim_extreme_recruit_error

   //////////// HERE /////////////////////

  if (debug == 15){
    cout<<"Ninit\n"<<N<<endl;
    cout<<"propmature\n"<<propmature<<endl;
    cout<<"Fyr\n"<<Fyr<<endl;
    // cout<<"surv_obsError\n"<<surv_obsError<<endl;
    // cout<<"catch_obsError\n"<<catch_obsError<<endl;
    cout<<"rec_procError\n"<<rec_procError<<endl;
    cout<<endl<<"manually exiting after calc_initial_states...\n"<<endl;
    exit(-1);
  }
  //other covariate sums here too? time series will be input


//----------------------------------------------------------------------------------------
FUNCTION calc_update_N
//----------------------------------------------------------------------------------------

// simply make N(t) = N(t-1)
 for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
            for(int isize=1; isize<=Nsizebins; isize++){
               N(area,spp,t,isize) = N(area,spp,t-1,isize);
            }
         }
  }



//----------------------------------------------------------------------------------------
FUNCTION calc_recruitment
//----------------------------------------------------------------------------------------

  //recruitment(t) =  recruitment_alpha  * pow (egg production(t-1),recruitment_shape) *
  //              exp(-recruitment_beta * egg production(t-1) +
  //              sumover_?(recruitment_covwt * recruitment_cov(t)))


  if ((Nstepsyr == 1 || t % Nstepsyr == 1) && (yrct <= Nyrs)) {  // recruits enter at start of year
    // simulate a vector of size Nspecies from uniform distribution between 0, 1 - probabilities
    // if prob < threshold then extreme event occurs and we sample from alternative distribution otherwise from ricker, beverton etc
//    for (spp =1; spp<=Nspecies;spp+) {



    for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){

	//	 switch (rectype(spp)){
  //         case 1:	  				//egg production based recruitment, 3 par gamma (Ricker-ish)
		// 	eggprod(area,spp)(yrct-1) /= Nstepsyr; //average egg production for a single "spawning" timestep
		// 	//eggprod(area,spp)(yrct) = recruitment_shape(area,spp)/recruitment_beta(area,spp);
		// 	recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * pow(eggprod(area,spp)(yrct-1), recruitment_shape(area,spp)) *
  //                                         mfexp(-recruitment_beta(area,spp) * eggprod(area,spp)(yrct-1) +
  //                                              recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));
  //     //cout << spp << " " << yrct << " " << recruitment(area,spp)(yrct) << " " <<
  //     //recruitment_alpha(area,spp) << " " << eggprod(area,spp)(yrct-1) << " " << recruitment_shape(area,spp) << " " <<
  //     //                                    recruitment_beta(area,spp) << " " << 
  //     //                                         recruitment_covwt(spp) << " " << recruitment_cov(yrct-1) << endl;
  //     //exit(-1);
		//   break;
	 //  case 2:                   //SSB based recruitment, 3 par Deriso-Schnute; see Quinn & Deriso 1999 p 95
		//     //SSB(area,spp)(yrct) /= Nstepsyr; //use? average spawning stock bio for a single "spawning" timestep, now SSB is at time t
  //                       recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1) *
  //                                          pow((1-recruitment_beta(area,spp)*recruitment_shape(area,spp)*SSB(area,spp)(yrct-1)),
  //                                           (1/recruitment_shape(area,spp)));
  //                                    //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92
  //                       recruitment(area,spp)(yrct) *= mfexp(-recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));
		//   break;
  //         case 3:	  				//SSB based recruitment, 3 par gamma (Ricker-ish)
		// 	//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
		// 	recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * pow(SSB(area,spp)(yrct-1), recruitment_shape(area,spp)) *
  //                                         mfexp(-recruitment_beta(area,spp) * SSB(area,spp)(yrct-1) +
  //                                              recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));
		//   break;
  //         case 4:	  				//SSB based recruitment, 2 par Ricker
		// 	//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
		// 	recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1) *
  //                                         mfexp(-recruitment_beta(area,spp) * SSB(area,spp)(yrct-1) +
  //                                              recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));
		//   break;
  //         case 5:	  				//SSB based recruitment, 2 par Beverton Holt
		// 	//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
		// 	recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1) /
  //                                        (1 + (recruitment_beta(area,spp) * SSB(area,spp)(yrct-1)));
  //                                    //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92
  //             /////////////////////////////////////////////// WHY -ve recruitment_covwt /////////////////////////////////////////////////////
  //                                     recruitment(area,spp)(yrct) *= mfexp(-recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));
		//   break;


  //          case 6:
		// 	//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
		// 	recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1) /
  //                                        (1 +  pow(value( SSB(area,spp)(yrct-1)/recruitment_beta(area,spp)),recruitment_shape(area,spp)) );
  //                                    //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92
  //                                     recruitment(area,spp)(yrct) *= mfexp(recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));

  //                break;

  //          case 7:  // Hockey Stick Stock recruitment
		// 	//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
  //                       if(SSB(area,spp)(yrct-1) <= recruitment_shape(area,spp)) { // S*
		// 	               recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1); // alpha.SSB
  //                       }  else {
		// 	               recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * recruitment_shape(area,spp); // alpha.SSB*
  //                       }
  //                                    //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92
  //                                     recruitment(area,spp)(yrct) *= mfexp(recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));

  //                break;


  //          case 8:  // Segmented Regression with breakpoint
		// 	//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
  //                       if(SSB(area,spp)(yrct-1) <= recruitment_shape(area,spp)) { // breakpoint
		// 	               recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1); // alpha.SSB
  //                       }  else {
		// 	               recruitment(area,spp)(yrct) = (recruitment_alpha(area,spp) * SSB(area,spp)(yrct-1)) +
  //                                            (recruitment_beta(area,spp) *( SSB(area,spp)(yrct-1)-recruitment_shape(area,spp)) ); // alpha.SSB + beta(ssB-breakpoint)
  //                       }

  //                       // check to see if recruitment goes below zero. which it is possible to do
  //                       if ( recruitment(area,spp)(yrct) < 0) {
  //                          recruitment(area,spp)(yrct) = 0.0;
  //                       }
  //                                    //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92
  //                                     recruitment(area,spp)(yrct) *= mfexp(recruitment_covwt(spp) * trans(recruitment_cov)(yrct-1));

  //                break;


  //          case 9:                   //Average recruitment plus devs--giving up on functional form
                       //recruitment(area,spp)(yrct) = mfexp(avg_recruitment(area,spp)+recruitment_devs(area,spp,yrct));
                       recruitment(area,spp)(yrct) = avg_recruitment(area,spp)*mfexp(recruitment_devs(area,spp,yrct)-0.5*recsigma(area, spp)*recsigma(area, spp));  //GF 2022/03/04, avg_recruitment is already in real space. This equation does not include lognormal bias correction (yet)
      //cout << spp << " " << yrct << " " << recruitment(area,spp)(yrct) << " " << avg_recruitment(area,spp) << " " << recruitment_devs(area,spp,yrct) << endl;
      //exit(-1);

		//   break;

  //          default:
  //           exit(1);
		// } //end switch


        // if(stochrec(spp)){                //simulate devs around recruitment curve
        //  // we allow the option of a "large" recruitment event (larger than under log normal) every so often as recommended by
        //  // CIE review team (Daniel Howell). We sample a random number from uniform distribution and based on frequency of large event
        //  // (from literature) determine if event should occur for species. We then sample from a distribution of event magnitudes.
        //  // NOT YET IMPLEMENTED

        //   if (rec_EventError(area,spp,yrct) >= 0) { // U~[0,1] to determine an extreme event x% of time

        //    // assumes log normal error. R = S.exp(Z) where Z = N(-sig2/2, sig2). In expectation R has mean = S. Slightly diff if Z = AR1
        //      recruitment(area,spp)(yrct) *=  mfexp(sim_recruit_error(area,spp,yrct));

        //    } else {   // extreme event
        //     //  recruitment(area,spp)(yrct)  = some other transformation depending on error structure (sim_extreme_recruit_error)
        //     ///////////  place holder for now ///////////////////
        //     recruitment(area,spp)(yrct) *=  mfexp(sim_recruit_error(area,spp,yrct));
        //    }
        // }  //end if stochastic
        // Now add recruitment to 1st size class
        N(area,spp,t,1) = N(area,spp,t,1) + recruitment(area,spp,yrct);


      }  //end spp
    }  //end area

  }  //end if last timestep in year



//----------------------------------------------------------------------------------------
FUNCTION calc_available_N
//----------------------------------------------------------------------------------------
// Since some of the species are only in certain management areas for a small percentage of the time
// fishing and predation mortality should be applied to a smaller proportion of the stock.

 
 // simply make Narea(t) = N(t)
 for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
            for(int isize=1; isize<=Nsizebins; isize++){
               Narea(area,spp,t,isize) = N(area,spp,t,isize);
            }
         }
  }

  //adjust Narea based on proportion of population in management area
  for (area=1; area<=Nareas; area++){
     for(spp=1; spp<=Nspecies; spp++){
         Narea(area,spp,t) = N(area,spp,t)*residentTime(area,spp);
     }
  }

//     cout << N(1,8,t,3)<<"-"<<N(1,6,t,5)<<endl;
//     cout << Narea(1,8,t,3)<<"-"<<Narea(1,6,t,5)<<endl;
//     cout << "___________"<<endl;
  



//----------------------------------------------------------------------------------------
FUNCTION calc_pred_mortality
//----------------------------------------------------------------------------------------
  
  //M2.initialize();
  //totalconsumedbypred = allmodeledprey(pred,predsize) + otherprey

  for (area=1; area<=Nareas; area++){
  	for(pred=1; pred<=Nspecies; pred++){
	    for(prey=1; prey<=Nspecies; prey++){
               // select the rows of suitability given predator (in blocks of Nsizebins )
               // suittemp is a Nsizebins x Nsizebins matrix.each value is suitability of pred size (col) on prey size(row)
		  dvar_matrix suittemp = suitability(area,pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins);
               // when using .sub on a higher order array, even if a matrix is the result the .rowmin() value is not set to 1
               // it uses the value of the row it occupied in the large array. This .rowmin() value determins if matrices can be multiplied
               // so we need to use .rowshif to designate the matrix to have rows starting from .rowmin()=1
		  suittemp.rowshift(1); //needed to match up array bounds
                  // vector * matrix = vector. result = [ sum(v*mat[,1]),sum(v*mat[,2]),sum(v*mat[,3]),sum(v*mat[,4]),sum(v*mat[,5])]
                  // standard  matrix  multiplication
                  //cout << pred << " " << prey << " " << suittemp << endl;
	      suitpreybio(area,pred,t) += wtconv*(elem_prod(binavgwt(prey),Narea(area,prey,t)) *  suittemp);
        
        for (int ipredsize=1;ipredsize<=Nsizebins;ipredsize++)     
          biomass_prey_avail_no_size(area,pred,yrct,ipredsize,prey) += sum(elem_prod(elem_prod(binavgwt(prey),Narea(area,prey,t)), column(suittemp,ipredsize)));
             }
        for (int ipredsize=1;ipredsize<=Nsizebins;ipredsize++)     
          biomass_prey_avail_no_size(area,pred,yrct,ipredsize,Nprey) += otherFood(pred,yrct); //suitability of other food needed?  //mod GF 11/14/22 to allow for predator-specific other food
        }
  } //ok

  


  /// NB: IN INITAL GAICHAS CODE ISSUE WITH THE MATRIX MULTIPLICATION OVER "INCORRECT" DIMENSION. THIS WAS CHANGED TO INCLUDE NESTED SIZE LOOPS
  //M2(area, prey, preysize,t) = sumover_preds_predsizes(intake*N(area,pred,predsize,t)*suitability(area,predpreysize)/
  //								sumover_preds_predsizes(totalconsumedbypred))

  //cout << "PM1 " << suitability(1,1) << endl;
  //cout  << "PM2 " << intake(1,1,yrct) << endl;
  //cout  << "PM3 " << suitpreybio(1,1,t) << endl;
  //cout  << "PM4 " << otherFood << endl;
  //cout  << "PM5 " << M2(1,1,t) << endl;

   for (area=1; area<=Nareas; area++){
  	for(prey=1; prey<=Nspecies; prey++){
               for(pred=1; pred<=Nspecies; pred++){
		  dvar_matrix suittemp2 = suitability(area,pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins);
                  // see above description of why row shif is needed
		  suittemp2.rowshift(1); //needed to match up array bounds
                  for (int ipreysize =1; ipreysize<=Nsizebins; ipreysize++) {
                   for (int ipredsize =1; ipredsize<=Nsizebins; ipredsize++) {                     
                     M2(area,prey,t,ipreysize) += (intake(area,pred,yrct,ipredsize)*Narea(area,pred,t,ipredsize) * suittemp2(ipreysize,ipredsize)) /
                           (suitpreybio(area,pred,t,ipredsize) + otherFood(pred,yrct));    //Hall et al 2006 other prey too high
                     //if (prey ==1) cout << "M2 " << pred << " " << ipreysize << " " << ipredsize << " " << M2(area,prey,t,ipreysize) << " " << intake(1,pred,yrct,ipredsize) << " " << suittemp2(ipreysize,ipredsize) << " " << suitpreybio(1,pred,t,ipredsize) << endl;
                    }
                  }

               } //pred

    } //prey
  } // ok

  //cout  << "PM6 " << M2(1,1,t) << endl;
  // Beet big dumb loop was written to explore issues with original M2. See 1_1_2 for details




//----------------------------------------------------------------------------------------
FUNCTION calc_fishing_mortality
//----------------------------------------------------------------------------------------

  //NOTE: Ftots are by area, species, and should be separated by fleet
  //not currently set up that way, assuming each fleet has same Ftot and they sum to F
  //selectivities are not currently by area, assuming fleet selectivity same in each area
  dvariable neglog19 = -1.*log(19.);
  for (area=1; area<=Nareas; area++){
      for(spp=1; spp<=Nspecies; spp++){


           for(fleet=1; fleet<=Nfleets; fleet++){

               for(int isizebin=1; isizebin<=Nsizebins; isizebin++) { //abeet added this loop to avoid compilation warnings.
               // could be created in calc_initial_sattes since it is not time dependent
                //fishsel(area,spp,fleet,isizebin) = 1/(1 + mfexp(-1.*fishsel_c(spp,fleet)*
                 //                    (lbinmidpt(spp,isizebin)-fishsel_d(spp,fleet))));
                fishsel(area,spp,fleet,isizebin) = 1./(1.+mfexp(neglog19*(lbinmidpt(spp,isizebin)-fishsel_c(1,fleet))/fishsel_d(2,fleet)));                                     
                                     }
                                     fishsel(area,spp,fleet) /= max(fishsel(area,spp,fleet));
                for(int isizebin=1; isizebin<=Nsizebins; isizebin++) { 
                 // Ffl(area,spp,fleet,t,isizebin) = fishsel(area,spp,fleet,isizebin)*(Fyr(area,spp,fleet,yrct))/Nstepsyr; //Andy Beet
                  // "catch mortality by fleet"  multiply by 1-p(discard) .i.e all landable catch:  p(no discard)
                  Ffl(area,spp,fleet,t,isizebin) = (1-discard_Coef(area,spp,fleet,isizebin))*fishsel(area,spp,fleet,isizebin)*(Fyr(area,spp,fleet,yrct))/Nstepsyr; //Andy Beet
                  // discard mortality by fleet" multiply by p(discard).(1-p(survive|discard))
                  Dfl(area,spp,fleet,t,isizebin) = (discard_Coef(area,spp,fleet,isizebin)*(1-discardSurvival_Coef(area,spp,fleet,isizebin))
                                                  *fishsel(area,spp,fleet,isizebin))*(Fyr(area,spp,fleet,yrct))/Nstepsyr; //Andy Beet

                  // partition fishing mortality,F  into "landings mortality", and "discard mortality", D
                  //  (1-p(discard)).F    +   p(discard).(1-p(discard|survive)).F
                  // == F - F.p(discard).p(survive|discard). See documentation
               }
               // sum mortalities over fleet
               D(area,spp,t) += Dfl(area,spp,fleet,t);
               F(area,spp,t) += Ffl(area,spp,fleet,t);
      }
    }
  }


//----------------------------------------------------------------------------------------
FUNCTION calc_total_mortality
//----------------------------------------------------------------------------------------
 for (area=1; area<=Nareas; area++){
     for(spp=1; spp<=Nspecies; spp++){

       //mort components, with all fleets already in F and D
       // F is mortality due to Fishing - landed species, D is discard mortality
       // Split F up in calc_catch_etc. Catch = F (landings) + D (discards)
       Z(area,spp,t) = M1(area,spp,t) +  M2(area,spp,t) +  F(area,spp,t) + D(area,spp,t);

     }
 }

//----------------------------------------------------------------------------------------
FUNCTION calc_catch_etc
//----------------------------------------------------------------------------------------

  //calculate Catch numbers at size (C), total catch_biomass, and N/biomass eaten and dead of other causes

  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
      //temp vectors for holding proportions
      dvar_vector Fprop = elem_div(F(area,spp,t),Z(area,spp,t)); //prop of death due to fishing of each size class
      dvar_vector M2prop = elem_div(M2(area,spp,t),Z(area,spp,t)); //prop of death due to predation of each size class
      dvar_vector M1prop = elem_div(M1(area,spp,t),Z(area,spp,t)); // prop of death due to other mortality. M1 read in from Data file
      dvar_vector Dprop = elem_div(D(area,spp,t),Z(area,spp,t)); // prop of death due to Discards of each size class
      dvar_vector Ndeadtmp = elem_prod((1-exp(-Z(area,spp,t))),Narea(area,spp,t));// total number dead in each size class
      // note: Z = total mortality

      //these are numbers at size dying each timestep from fishing, predation, and M1 (other)
      C(area,spp,t) = elem_prod(Fprop, Ndeadtmp); //fishing catch on GB
      eatN(area,spp,t) = elem_prod(M2prop, Ndeadtmp); // predation, M2
      discardN(area,spp,t) = elem_prod(Dprop,Ndeadtmp); // discards on vessel either not target species or not allowed to land
      otherDead(area,spp,t) = Ndeadtmp - C(area,spp,t) - eatN(area,spp,t)- discardN(area,spp,t); // M1

//      GF commented chunk out 04/28/2022 to remove discontinutity in objective function (associated with projections)
      // // all catch is considered discard since can not be landed if found to be so in assessment.
      // // catchTtoDiscards is a binary vector indicating threshold exceedance. Default all = 0
      // if (catchToDiscardsSpecies(area,spp) == 1) {
      //    discardN(area,spp,t) =  discardN(area,spp,t) + C(area,spp,t);
      //    C(area,spp,t) = 0.0;
      // }

//      GF commented chunk out 04/28/2022 to remove discontinutity in objective function (associated with projections)
//      // check to see if species part of a guild in trouble. if so set catch to discards and catch(landings) = 0
//      //      Default flag:  all = 0
//      if (catchToDiscardsGuild(area,guildMembers(spp)) == 1) {
//
//         // this species is a member of a guild whose guild biomass has exceeded threshold
//        discardN(area,spp,t) =  discardN(area,spp,t) + C(area,spp,t);
//        C(area,spp,t) = 0.0;
//      }

      // size class level
      eaten_biomass_size(area,spp,yrct) += wtconv*elem_prod(eatN(area,spp,t),binavgwt(spp));
      discard_biomass_size(area,spp,yrct) +=  wtconv*elem_prod(discardN(area,spp,t),binavgwt(spp));
      otherDead_biomass_size(area,spp,yrct) += wtconv*elem_prod(otherDead(area,spp,t),binavgwt(spp));
      total_biomass_size(area,spp,yrct) += wtconv*elem_prod(Narea(area,spp,t),binavgwt(spp));

      //these are annual total biomass losses for comparison with production models
      eaten_biomass(area,spp,yrct) += sum(wtconv*elem_prod(eatN(area,spp,t),binavgwt(spp)));
      discard_biomass(area,spp,yrct) +=  sum(wtconv*elem_prod(discardN(area,spp,t),binavgwt(spp)));
      otherDead_biomass(area,spp,yrct) += sum(wtconv*elem_prod(otherDead(area,spp,t),binavgwt(spp)));
      total_biomass(area,spp,yrct) += sum(wtconv*elem_prod(Narea(area,spp,t),binavgwt(spp)));

      //do fleet specific catch in numbers, biomass, sum for total catch
      for(fleet=1; fleet<=Nfleets; fleet++){
	  dvar_vector Fflprop = elem_div(Ffl(area,spp,fleet,t),F(area,spp,t));// proportion dead due to fleet in each size class. vec length= num classes
          Cfl(area,spp,fleet,t) = elem_prod(Fflprop, C(area,spp,t)); // numbers dying from fleet by sizeclass
          fleet_catch_biomass(area,spp,fleet,yrct) += sum(wtconv*elem_prod(Cfl(area,spp,fleet,t),binavgwt(spp)));
          catch_biomass_size(area,spp,yrct) += wtconv*elem_prod(Cfl(area,spp,fleet,t),binavgwt(spp));
          catch_biomass(area,spp,yrct) += sum(wtconv*elem_prod(Cfl(area,spp,fleet,t),binavgwt(spp)));

//          //add obs error for est fleet catch, sum for est total catch for simulations
//          est_fleet_catch_biomass(area,spp,fleet,yrct) = fleet_catch_biomass(area,spp,fleet,yrct) * exp(sim_catch_error(area,spp,fleet,yrct) ); //add obs error
         // est_fleet_catch_biomass(area,spp,fleet,yrct) = fleet_catch_biomass(area,spp,fleet,yrct) * exp(catch_sigma(area,spp,fleet)
         //                                * catch_obsError(area,spp,fleet,yrct)
         //                                - 0.5 * catch_sigma(area,spp,fleet) * catch_sigma(area,spp,fleet)  ); //add obs error
   //        if (t % Nstepsyr == 0){//if we are just about to end the year //andybeet
   //              est_catch_biomass(area,spp,yrct) += est_fleet_catch_biomass(area,spp,fleet,yrct);
	  // }//end if

          for (int isize=1; isize<=Nsizebins; isize++){
              // used for indices- LFI
              Cfl_tot(area,spp,fleet,yrct,isize) +=  Cfl(area,spp,fleet,t,isize);// fishing. total number by size/species/fleet each year. summed over t
              C_tot(area,spp,yrct,isize) += Cfl(area,spp,fleet,t,isize); // fishing. total number by size/species each yr summed over fleet and Nstepsyr timesteps
          }

      }//end fleet loop
    }//end species loop
  }//end area loop

 // // aggregate catch to guild level at end of year Andy Beet
 // // also calculate predation rate and fishingRate
 //   if ((t % Nstepsyr == 0) && (yrct <= Nyrs)){

 //   for(area=1; area<=Nareas; area++){
 //          for (spp=1; spp<=Nspecies; spp++) {
 //             for (int isize = 1; isize<=Nsizebins; isize++ ) {
 //             // look out for nans
 //                if (total_biomass_size(area,spp,yrct,isize) < .001) { // zero
 //                  predation_mortality_size(area,spp,yrct,isize) = 0;
 //                  fishing_mortality_size(area,spp,yrct,isize) = 0;
 //                } else {
 //                  predation_mortality_size(area,spp,yrct,isize) = eaten_biomass_size(area,spp,yrct,isize)/(total_biomass_size(area,spp,yrct,isize)/Nstepsyr);
 //                  fishing_mortality_size(area,spp,yrct,isize) = (catch_biomass_size(area,spp,yrct,isize)+discard_biomass_size(area,spp,yrct,isize))/(total_biomass_size(area,spp,yrct,isize)/Nstepsyr);
 //               }
 //             }
 //             if (total_biomass(area,spp,yrct) <.001) {
 //               predation_mortality(area,spp,yrct) = 0;
 //               fishing_mortality(area,spp,yrct) = 0;
 //             } else {
 //                predation_mortality(area,spp,yrct) = eaten_biomass(area,spp,yrct)/(total_biomass(area,spp,yrct)/Nstepsyr);
 //                fishing_mortality(area,spp,yrct) = (catch_biomass(area,spp,yrct)+discard_biomass(area,spp,yrct))/(total_biomass(area,spp,yrct)/Nstepsyr);
 //             }
 //          }
 //   }
 //   //  test<< "ttt = "<< t <<", yr =  "<<yrct <<endl;
 //   for(area=1; area<=Nareas; area++){
 //       for (int iguild=1; iguild<=Nguilds; iguild++) {
 //          for (spp=1; spp<=Nspecies; spp++) {
 //            for (fleet=1; fleet<=Nfleets; fleet++) {
 //                 if (guildMembers(spp) == iguild) {
 //                // cout<<iguild<<","<<spp<<","<<fleet<<endl;
 //                    est_fleet_catch_guild_biomass(area,iguild,fleet,yrct) += est_fleet_catch_biomass(area,spp,fleet,yrct);
 //                    est_catch_guild_biomass(area,iguild,yrct) += est_fleet_catch_biomass(area,spp,fleet,yrct);
 //                 }
 //            }
 //          }
 //       }
 //   }
 //  } // end if t%


//----------------------------------------------------------------------------------------
FUNCTION calc_pop_dynamics
//----------------------------------------------------------------------------------------

 // For species not resident in area.
 // Assume contant total mortality rate for population not in management area.
 // Adjust that proportion of population

  //ofstream popout("popstructure.out");

  for(area = 1; area <=Nareas; area++) {
       for (spp = 1; spp <=Nspecies; spp++) {
           for(int isize=1; isize <= Nsizebins; isize++){
              // For pop outside of management area => entire population * proportion of population out of area * mortality rate
              Nnotarea(area,spp,t,isize) = N(area,spp,t,isize)*(1.0-residentTime(area,spp))*(1.0-areaMortality(area,spp));
           }
       }
  }
  //cout << "A "<< t << " " << Nnotarea(1,1,t) << endl;

  // POP DYNAMICS for Area of interest
  //for all older than recruits,
  //pop is composed of survivors from previous size growing into current size and staying in area
  //plus survivors in current size not growing out of current size and staying in area
  //plus immigrants of current size from other areas
  //minus emigrants of current size to other areas
  //movement is not yet specified, so we leave out the immigration and emigration parts for now

  // Recruits added already in recuitment module

  //N(area, spp,t,bin) +=
  //                       N(area, spp,t,bin-1)*S(area,spp,t,bin-1)*growthprob_phi(area,spp,bin-1) +
  //                       N(area, spp,t,bin)*S(area,spp,t,bin)*(1-growthprob_phi(area,spp,bin)))

  //cout << "C "<< t << " " << Z(1,1,t) << endl;
  //cout << "D "<< t << " " << growthprob_phi(1,1,yrct) << endl;

  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){

        // For all bins except smallest. execute from largest to smallest. Remember N(t) = N(t-1) in first step
        //N = surviving and growing from smaller bin and surviving and staying in current bin
        for(int isize=Nsizebins; isize>=2; isize--){
       	 Narea(area,spp,t,isize) = Narea(area,spp,t,isize-1) * exp(-Z(area,spp,t,isize-1)) * growthprob_phi(area,spp,yrct,isize-1)
                            +  Narea(area,spp,t,isize)* exp(-Z(area,spp,t,isize))  * (1-growthprob_phi(area,spp,yrct,isize));
         //popout << spp << " " << t << " " << isize << " " << Narea(area,spp,t,isize) << endl;
         N_tot(area,spp,yrct,isize) += Narea(area,spp,t,isize);// cumulate sum. averaged in indices

        }//end size loop


        // smallest size class. Survivors that stay in same size class
        // we added recruits at start of year to current time. they were then fished in this time period
//	N(area,spp,t,1) = N(area,spp,t-1,1)* exp(-Z(area,spp,t,1))*(1-growthprob_phi(area,spp,yrct,1));
	Narea(area,spp,t,1) = Narea(area,spp,t,1)* exp(-Z(area,spp,t,1))*(1-growthprob_phi(area,spp,yrct,1));

        //popout << spp << " " << t << " " << 1 << " " << Narea(area,spp,t,1) << endl;

        N_tot(area,spp,yrct,1) += Narea(area,spp,t,1); // running total for the year. Used for indices

      }//end species loop

  }//end area loop


  // Add two populations (in management area, outside management area)
   for(area = 1; area <=Nareas; area++) {
       for (spp = 1; spp <=Nspecies; spp++) {
           N(area,spp,t) = Nnotarea(area,spp,t) + Narea(area,spp,t);
       }
   }

       //  cout<<endl;

//----------------------------------------------------------------------------------------
FUNCTION calc_SSB
//-------------------------------------------------------------------------------------
  //egg production(t) = sumover_length(fecundity*prop mature(t)*sexratio*N(t))
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
	  dvar_vector fecundmature = elem_prod(fecundity(area,spp), propmature(area, spp)(yrct));
          dvar_vector sexratioN = sexratio(area, spp) * N(area,spp)(t);
          eggprod(area,spp)(yrct) += sum(elem_prod(fecundmature, sexratioN));  //accumulates eggs all year--appropriate?

          dvar_vector Nmature = elem_prod(propmature(area, spp)(yrct), N(area,spp)(t));
                //SSB(area,spp)(yrct) += sum(wtconv*elem_prod(Nmature,binavgwt(spp)));  //accumulates SSB all year; not appropriate
          SSB(area,spp)(yrct) = sum(wtconv*elem_prod(Nmature,binavgwt(spp)));  //SSB in this timestep, overwrites previous

          // Final SSB(year) = SSB in time step 5, 1- etc
    }
  }

//----------------------------------------------------------------------------------------
FUNCTION calc_movement
//----------------------------------------------------------------------------------------

  //not yet specified will probably have migration in array and add random too
  int probmovein = 0.0;        //will be an array for area to area movement
  int probmoveout = 0.0;       //as will this
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
			  //N(area,spp,t) += N(area,spp,t) * probmovein(area,spp);
			 // N(area,spp,t) -= N(area,spp,t) * probmoveout(area,spp);
   // ******************************************************************************
   // weight/length is not linear. The following line  will underestimate B
      B(area,spp,t) = wtconv*elem_prod(Narea(area,spp,t),binavgwt(spp));  //do after movement

      for (int isize=1;isize<=Nsizebins;isize++) {
          // add up B over t for each year keeping size class structure. used in indices
          B_tot(area,spp,yrct,isize) += B(area,spp,t,isize);
      }
    }
  }

//----------------------------------------------------------------------------------------
FUNCTION calc_survey_abundance
//----------------------------------------------------------------------------------------


  for (area=1; area<=Nareas; area++){
      for(spp=1; spp<=Nspecies; spp++){
	   avByr(area,spp)(yrct) += sum(B(area,spp,t))/Nstepsyr;
//            est_survey_biomass(area,spp,yrct) =  avByr(area,spp,yrct)*survey_q(area,spp); //add surv q

//        GF commented next two lines out 04/28/2022 as relate to generating index values
//            // multiplicative LN error
//            est_survey_biomass(area,spp,yrct) *= exp(sim_survey_error(area,spp,yrct)); // error created in calc_initial_states

          //  if(yrct == 1) {
          //    est_survey_biomass(area,spp,yrct) *= exp(rho_AR_Survey*xts(area,spp,1) +
          //                               surv_sigma(area,spp)*surv_obsError(area,spp,yrct)*pow(1.0 -pow(rho_AR_Survey,2),0.5)
          //                               - 0.5 * surv_sigma(area,spp) * surv_sigma(area,spp)  ); //add obs error
          //  } else {
          //    est_survey_biomass(area,spp,yrct) *= exp(rho_AR_Survey*xts(area,spp,yrct-1) +
          //                               surv_sigma(area,spp)*surv_obsError(area,spp,yrct)*pow(1.0 -pow(rho_AR_Survey,2),0.5)
          //                               - 0.5 * surv_sigma(area,spp) * surv_sigma(area,spp)  ); //add obs error
          //  }
//       est_survey_biomass(area,spp,yrct) *= exp(surv_sigma(area,spp)
//                                         * surv_obsError(area,spp,yrct)
//                                         - 0.5 * surv_sigma(area,spp) * surv_sigma(area,spp)  ); //add obs error


      }
  }

// GF 04/28/2022 below code not needed for estimation but not affecting differentiability of likelihood

 // Added by Andy Beet
 // we need to sum up the biomass over each guild and check for excedences.
 // Do at end of year only. Used in assessment module and health indices module
  // if ((t % Nstepsyr == 0) && (yrct <= Nyrs)){
  // //  test<< "ttt = "<< t <<", yr =  "<<yrct <<endl;
  //  for(area = 1; area<=Nareas; area++){
  //      for (iguild=1; iguild<=Nguilds; iguild++) {
  //         for (spp=1; spp<=Nspecies; spp++) {
  //            if (guildMembers(spp) == iguild) {
  //              est_survey_guild_biomass(area,iguild,yrct) += est_survey_biomass(area,spp,yrct);
  //            }
  //         }
  //      }
  //  }
  // } // end if t%


// GF 04/28/2022 commenting out as part of debugging....
// //----------------------------------------------------------------------------------------
// FUNCTION calc_health_indices
// //----------------------------------------------------------------------------------------
// // Here we calculate several indices: measures of system health at the end of the year
// // These metrics could all be calculated in R after the run but we may want to use these in management so we need them real time

//  if ((t % Nstepsyr == 0) && (yrct <= Nyrs)){

// // note: we could combine these indices into the same loop, but chose not to ease readability
// // 1. Simpsons Diversity Index (Richness (number of species) and evenness(relative numbers)) = sum((N_i/N)^2) = sum(p_i^2)
// //  we use the  mean N over Nstepsyr as the annual value of N
// // simpsons for N

//  for(int iarea=1;iarea<=Nareas;iarea++){
//       prob_species.initialize();
//       dvariable N_total = 0;
//       for (int isp=1; isp<=Nspecies;isp++) {

//         if (sum(N_tot(iarea,isp,yrct)) < .0001) {
//             prob_species(isp) = 0;
//         } else {
//            prob_species(isp) = pow(sum(N_tot(iarea,isp,yrct))/Nstepsyr,2);
//         }
//         N_total += sum(N_tot(iarea,isp,yrct))/Nstepsyr;
//       }
//       index_Simpsons_N(iarea,yrct) = sum(prob_species)/pow(N_total,2);
//       index_Simpsons_Nrecip(iarea,yrct) =1/index_Simpsons_N(iarea,yrct);
//    }

// // simpsons for Catch - summed over fleet
//    for(int iarea=1;iarea<=Nareas;iarea++){
//       prob_species.initialize();
//       dvariable C_total = 0;
//       for (int isp=1; isp<=Nspecies;isp++) {
//           if (sum(C_tot(iarea,isp,yrct))< 0.0001) {
//              prob_species(isp) = 0;
//           } else {
//              prob_species(isp) = pow(sum(C_tot(iarea,isp,yrct))/Nstepsyr,2);
//           }
//         C_total += sum(C_tot(iarea,isp,yrct))/Nstepsyr;
//       }
//       if (C_total < .0001) {
//           index_Simpsons_C(iarea,yrct) = sum(prob_species);
//       } else {
//           index_Simpsons_C(iarea,yrct) = sum(prob_species)/pow(C_total,2);
//       }
//       index_Simpsons_Crecip(iarea,yrct) = 1/index_Simpsons_C(iarea,yrct);
//    }

//  //    Cfl_tot(area,spp,fleet,yrct,isize)
// // 2. Large Fish Indices
// //  i. LFI_Biomass = %biomass of largest sizeclass relative to total biomass for each species
// // ii. LFI_Catch = same for catch data
// //iii. LFI_N = number of large fish
//    for (int iarea=1;iarea<=Nareas;iarea++) {
// //       LF_Biomass = 0;
//        for (int isp=1; isp<=Nspecies;isp++) {
//           if (sum(B_tot(iarea,isp,yrct)) < .0001 ) {
//             index_LFI_Biomass(iarea,isp,yrct) = 0 ;
//            } else {
//              index_LFI_Biomass(iarea,isp,yrct) = B_tot(iarea,isp,yrct,Nsizebins)/sum(B_tot(iarea,isp,yrct)); // large fish in top size category for each fish. Biomass
//           }
//           if (sum(C_tot(iarea,isp,yrct)) < .0001) {
//             index_LFI_Catch(iarea,isp,yrct) = 0;
//           } else {
//             index_LFI_Catch(iarea,isp,yrct) = C_tot(iarea,isp,yrct,Nsizebins)/sum(C_tot(iarea,isp,yrct)); // large fish in top size category for each fish. Catch
//           }
//           index_LFI_N(iarea,isp,yrct) = N_tot(iarea,isp,yrct,Nsizebins)/Nstepsyr; // number of large fish per year
//        }
//    }

// // 3. Predator to prey biomass ratio. (Pred = dogfish, skate, goosefish, cod, silverhake ) / (Prey = herring, mackerel,haddock, yellowtail, winter flounder)
// //      uses average annual biomass for each species (over all sizeclasses)
// // 4. planktivore : piscivore ratio
//    for (int iarea=1;iarea<=Nareas;iarea++) {
//      for (int isp=1; isp<=Nspecies ; isp++){
//      // pred:prey
//       if (predOrPrey(isp) == 1) { // predator
//          index_predBio(iarea,yrct) += avByr(iarea,isp,yrct);
//        } else { // prey
//          index_preyBio(iarea,yrct) += avByr(iarea,isp,yrct);
//        }
//      }
//      index_predToPreyRatio(iarea,yrct) = index_predBio(iarea,yrct)/index_preyBio(iarea,yrct);
//      // planktivore:piscivore
//      if (est_survey_guild_biomass(iarea,1,yrct) < .0001) {
//         index_plankToPiscRatio(iarea,yrct) = 0;
//      } else {
//           index_plankToPiscRatio(iarea,yrct) = est_survey_guild_biomass(iarea,2,yrct)/est_survey_guild_biomass(iarea,1,yrct);
//      }
//    }

// // 5. variance of catch using a m-year moving window (bandwidth_metric)
//    if (yrct >= bandwidth_metric) { // take mean of last bandwidth_metric years
//       for (int iarea=1; iarea<=Nareas; iarea++){
//             for (int isp=1; isp<=Nspecies; isp++) {
//            index_catch.initialize();
//             index_biomass.initialize();
//                int ic = 0;
//                 for (int iyear = yrct-bandwidth_metric+1;  iyear<=yrct; iyear++) {
//                   // test << iyear << "," << est_catch_biomass(iarea,isp,iyear)  <<endl;
//                    ic++;
//                    index_catch(ic) = est_catch_biomass(iarea,isp,iyear);
//                    index_biomass(ic) = avByr(iarea,isp,iyear);
//                   //test << iyear << "," << catch_data(ic)  << endl;
//                 }
//                 // check to see if all elements are same abs(mean - geometric mean). if so std_dev() fails
//                 if ((sum(index_catch) < 1e-6) || (abs(value(mean(index_catch) - exp(sum(log(index_catch))/bandwidth_metric))) < 1e-6 ) ) {
//                    index_stdev_catch(iarea,isp,yrct) = 0;
//                 } else {
//                    index_stdev_catch(iarea,isp,yrct) = std_dev(index_catch);
//                 }
//                 if ((sum(index_biomass) < 1e-6) || (abs(value(mean(index_biomass) - exp(sum(log(index_biomass))/bandwidth_metric))) < 1e-6 ) ) {
//                    index_stdev_biomass(iarea,isp,yrct) = 0;
//                 } else {
//                    index_stdev_biomass(iarea,isp,yrct) = std_dev(index_biomass);
//                 }
//               //  test << yrct << "," << mean(catch_data)<<","<< std_dev(catch_data) << "\n"<<endl;
//             }
//       }
//    }

// // 6. Exploitation Rate
//       for (int iarea=1; iarea<=Nareas; iarea++){
//          dvariable total_catch = 0.0;
//          dvariable total_bio = 0.0;

//             for (int isp=1; isp<=Nspecies; isp++) {
//               if(total_biomass(iarea,isp,yrct) < .0001) {
//                 index_ExploitationRate(iarea,isp,yrct) = 0;
//               } else {
//                 index_ExploitationRate(iarea,isp,yrct) = est_catch_biomass(iarea,isp,yrct)/total_biomass(iarea,isp,yrct);
//               }
//                 total_catch += est_catch_biomass(iarea,isp,yrct);
//                 total_bio += total_biomass(iarea,isp,yrct);
//             }

//             index_SystemExploitationRate(iarea,yrct) = total_catch/total_bio;
//       }

// // 7. individual Species Biomass < 20% B0.
// // guild biomass < 20%
// // We indicate "yes" or "no". Does a species fall below threshold in each year yrct
//    for (int iarea=1; iarea<=Nareas; iarea++){
//        // species level
//         for (int isp=1; isp<=Nspecies; isp++) {
//             if (est_survey_biomass(iarea,isp,yrct) <= (B0(iarea,isp)* (baseline_threshold + threshold_species(isp)))) {
//                  index_status_species(iarea,isp,yrct) = 1;
//             }
//         }
//         // guild level
//         for (int iguild=1; iguild<=Nguilds; iguild++) {
//             if (est_survey_guild_biomass(iarea,iguild,yrct) <= (B0_guilds(iarea,iguild)* baseline_threshold)) {
//                  index_status_guild(iarea,iguild,yrct) = 1;
//             }
//         }

//     }



//   } // end of year if



// //----------------------------------------------------------------------------------------
// FUNCTION calc_assessment_linear_independent_fleet
// //----------------------------------------------------------------------------------------
// // We enter this loop every AssessmentPeriod years, during the last time period of the year.
// // We recalculate the new exploitation baased on the current biomass levels using a linear relationship (no Step)
// // Each species or complex is flagged to indicate which have breached the baseline threshold (indicating big trouble).
// // Consequently landings are not allowed and all catch is considered discards. This is dealt with in catch module.

// // ALL FLEETS ARE IMPACED INDEPENDENTLY OF EACH OTHER.
// // FLEET BASED EXPLOITATION

// // Move from system wide exploitation to fleet based exploitation
// // Defines the ramp properties based on min and max exploitation. The range of exploitations permitted post assessment.
// // This is independent of starting exploitation which could be outside the range

//    // ramp properties for each species
//    dmatrix slopeSpecies(1,Nareas,1,Nspecies);
//    dmatrix interceptSpecies(1,Nareas,1,Nspecies);

//    // ramp down for each fleet. Need to compare species in functional group against fleet they are predominantly caught in
//    dmatrix slopeGuild(1,Nareas,1,Nfleets);
//    dmatrix interceptGuild(1,Nareas,1,Nfleets);
//    for (area=1; area<=Nareas; area++) {
//        for (int ifleet = 1; ifleet<=Nfleets; ifleet++) {
//            slopeGuild(area,ifleet) = (maxExploitation(ifleet)-minExploitation(ifleet))/(minMaxThreshold(2) - minMaxThreshold(1));
//            interceptGuild(area,ifleet) =  maxExploitation(ifleet) - slopeGuild(area,ifleet)*minMaxThreshold(2);
//             //             cout<<slopeGuild(area,ifleet) <<" - "<<interceptGuild(area,ifleet)<<endl;
//        }
//    }

//    // Species ramp are assigned to the fleet which predominantly catches them. So if we need to protect
//    // the species we can ramp the fleet which impacts them the most
//    for (area=1; area<=Nareas; area++) {
//        for (spp = 1; spp<=Nspecies; spp++) {
//            int ifleet = fleetMembers(guildMembers(spp));
//            slopeSpecies(area,spp) = (maxExploitation(ifleet)-minExploitation(ifleet))/(minMaxThreshold(2) - minMaxThreshold(1));
//            interceptSpecies(area,spp) =  maxExploitation(ifleet) - slopeSpecies(area,spp)*(minMaxThreshold(2)+threshold_species(spp));
//        }
//    }

//  // Assess the guild/functional group biomass levels. What % of B0 are they at

//    // guild level calcs
//    // now average the guild biomass values over AssessmentPeriod yrs and then we adjust the exploitation rate
//      exploitationLevelGuild.initialize();
//        for (area=1 ; area<=Nareas; area++) {
//          for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//              newExploitationLevel(area,ifleet) = .000001;
//          }
//          for (iguild=1; iguild<=Nguilds; iguild++) {
//              catchToDiscardsGuild(area,iguild) = 0;//resets flag to indicate species has not exceeded min threshold
//              for (iassess=1; iassess<=AssessmentPeriod;iassess++){
//                   // calculate the mean biomass and catch over the Assessment period
//                   est_survey_guild_biomass_assessment(area,iguild,yrct) += est_survey_guild_biomass(area,iguild,yrct-iassess+1)/AssessmentPeriod;
//                   for (int ifleet=1;ifleet<=Nfleets;ifleet++) {// catch by fleet over last AssessmentPeriod Years
//                       est_fleet_catch_guild_assessment(area,iguild,ifleet,yrct) += est_fleet_catch_guild_biomass(area,iguild,ifleet,yrct-iassess+1)/AssessmentPeriod;
//                   }
//              }

//              // a fleetMember is the fleet most asscoiated with fishing the guild
//              int ifleet = fleetMembers(iguild);
//              // check to see if average < min threshold or > max threshold and assign new exploitation
//              // otherwise adjust exploitation linearly
//              dvariable biomassLevel = est_survey_guild_biomass_assessment(area,iguild,yrct)/B0_guilds(area,iguild);
//              if (biomassLevel <= minMaxThreshold(1) ) {// min threshold
//                 exploitationLevelGuild(area,iguild) = minExploitation(ifleet);
//              } else if (biomassLevel >= minMaxThreshold(2) ) {// max threshold
//                 exploitationLevelGuild(area,iguild) = maxExploitation(ifleet);
//              } else { // linear ramp
//                exploitationLevelGuild(area,iguild) = (biomassLevel * slopeGuild(area,ifleet)) + interceptGuild(area,ifleet);
//              }
//              if (biomassLevel <= baseline_threshold) { // this is bad. no longer allowed to land caught fish in this guild
//                catchToDiscardsGuild(area,iguild) = 1;
//              }
//             // cout<<yrct<<"  " <<iguild <<"   "<<biomassLevel<<"  "<<exploitationLevelGuild(area,iguild)<<endl;
//          } // guild loop

//          // For each fleet, take the min exploitation of all guilds caught by that fleet.
//          for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//              int icountf = 0;
//              dvariable exploitRate;
//              exploitRate.initialize();
//              for(iguild=1; iguild<=Nguilds; iguild++){
//                  if(fleetMembers(iguild) == ifleet) {
//                    icountf = icountf + 1;
//                    if (icountf == 1){
//                       exploitRate = exploitationLevelGuild(area,iguild);
//                    } else {
//                       exploitRate = min(value(exploitRate),value(exploitationLevelGuild(area,iguild)));
//                    }
//                  }
//               }
//              newExploitationLevel(area,ifleet) = exploitRate;
//          }

//      }  // area loop

   
//    // Check at species level also
//    // Include species detection level in determining rate change
//    if (speciesDetection == 1) {
//      // we check for exceedances at the species level. take the mean abundance over last AssessmentPeriod yrs for each species
//      exploitationLevelSpecies.initialize();
//      for (area=1; area<=Nareas; area++){
//          for (spp=1; spp<=Nspecies; spp++){
//              catchToDiscardsSpecies(area,spp) = 0; //resets flag to indicate species has not exceeded min threshold
//              for (iassess=1; iassess<=AssessmentPeriod; iassess++){
//                  // mean of last few years
//                  est_survey_biomass_assessment(area,spp,yrct) +=  est_survey_biomass(area,spp,yrct-iassess+1)/AssessmentPeriod;
//              }
//              int ifleet = fleetMembers(guildMembers(spp));
//               // now check for exceedances. if average < min threshold or > max threshold and assign new exploitation
//              // otherwise adjust exploitation linearly
//              dvariable biomassLevel =  est_survey_biomass_assessment(area,spp,yrct)/B0(area,spp);
//              if (biomassLevel <= minMaxThreshold(1) ) {// min threshold
//                 exploitationLevelSpecies(area,spp) = minExploitation(ifleet);
//              } else if (biomassLevel >= (minMaxThreshold(2)+threshold_species(spp)) ) {// max threshold
//                 exploitationLevelSpecies(area,spp) = maxExploitation(ifleet) ;
//              } else { // linear ramp
//                exploitationLevelSpecies(area,spp) = (biomassLevel * slopeSpecies(area,spp)) + interceptSpecies(area,spp);
//              }

//              if (biomassLevel <= baseline_threshold) { // this is bad. no longer allowed to land caught fish
//                catchToDiscardsSpecies(area,spp) = 1;
//              }
//              // Adjust the exploitation level found at functional group level.
//              // If any species are in trouble the exploitation will need to be reduced further
//              newExploitationLevel(area,ifleet) = min(value(exploitationLevelSpecies(area,spp)),value(newExploitationLevel(area,ifleet)));

//          } // spp  loop

//      } // area loop
//    } // species detection
 
//     // now we found new exploitation rates we need to act on them and adjust effort to correspond to rate
//     // If the scenario is a FixedRate scenario (all minExploitation == maxExploitation across fleets)
//     // we set exploitation rates all equal to fixed rate in data file, therefore when we
//     // encounter this phase the effort is unchanged

//     // store current exploitation level and set it for next few years until new assessment is due. used as output only
//     // if first time set exploitation_update for the first few years prior to assessment
//     // This section is purely for reporting out
//       for (area=1; area<=Nareas; area++) {
//          if (t==(Nstepsyr*AssessmentPeriod)) { //first assessment. assign 1st 3 yrs (not effected by assessment) to actual starting value of exploitation
//            for (int iassess=1; iassess <= AssessmentPeriod; iassess++) {
//               for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//                 exploitation_update(area,ifleet,iassess) =  maxExploitation(ifleet);// starting exploitation, maximum
//               }
//            }
//          }
//          // set all subsequent yrs to new exploitation otherwise last few years will revert to original rate/
//          for (int iy = yrct+1; iy<=Nyrs; iy++) {
//             for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//                exploitation_update(area,ifleet,iy) = newExploitationLevel(area,ifleet);
//             }
//          }
//      }


//      // Calculate the new effort for each fleet.
//      // Note that each  fleets is pacted based on the functional group/guild it fishes
//      for (area=1 ; area<=Nareas ; area++) {

//         for (int ifleet=1;ifleet<=Nfleets;ifleet++){
//          if ( mean_fishery_q(area,ifleet) < 1e-29){ // a fleet doesn't fish a particluar guild. keep effort same
//               // this will only happen at guild q not fleet q.
//               effort_updated(area,ifleet) = obs_effort(area,ifleet,yrct);
//          } else {
//               effort_updated(area,ifleet) = newExploitationLevel(area,ifleet)/mean_fishery_q(area,ifleet);
//          }
//         }
//      }
     
//      // Effort is used just once in initial.calcs() to obtain the Fyr terms for the whole simulation  so
//      // we need to use this new effort and create updated values for Fyr
//      for (area=1; area<=Nareas ; area++) {
//         for (spp=1; spp<=Nspecies; spp++) {
//             for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//                 for (int iy = yrct+1; iy<=Nyrs; iy++) {// set all subsequent yrs to new effort otherwise last few years will revert to original rate
//                  // this will all be updated during next assessment
// //               for (int iassess=1; iassess <= AssessmentPeriod; iassess++) {
//                   obs_effortAssess(area,ifleet,iy) = effort_updated(area,ifleet); // output only
//                   Fyr(area,spp,ifleet,iy) = fishery_q(area,spp,ifleet)*effort_updated(area,ifleet)*effortScaled(area,spp); //Andy Beet
//                 }
//             }
//         }
//      }



// //----------------------------------------------------------------------------------------
// FUNCTION calc_assessment_equal_fleet
// //----------------------------------------------------------------------------------------
// // We enter this loop every AssessmentPeriod years, during the last time period of the year.
// // We recalculate the new exploitation baased on the current biomass levels using a linear relationship (no Step)
// // Each species or complex is flagged to indicate which have breached the baseline threshold (indicating big trouble).
// // Consequently landings are not allowed and all catch is considered discards. This is dealt with in catch module

// // ALL FLEETS ARE IMPACED THE SAME WAY.
// // EG. IF EXPLOITATION IS REDUCED FROM 10% TO 5% ALL FLEETS ARE REDUCED TO 5%
// // regardless of which species/guild caused breach.
// // SYSTEM WIDE EXPLOITATION


//    // ramp properties
//    dmatrix slopeSpecies(1,Nareas,1,Nspecies);
//    dmatrix interceptSpecies(1,Nareas,1,Nspecies);
//    dvariable newExploitationLevel=0;
//    for (area=1; area<=Nareas; area++) {
//        for (spp = 1; spp<=Nspecies; spp++) {
//            slopeSpecies(area,spp) = (minMaxExploitation(2)-minMaxExploitation(1))/(minMaxThreshold(2) - minMaxThreshold(1));
//            interceptSpecies(area,spp) =  minMaxExploitation(2) - slopeSpecies(area,spp)*(minMaxThreshold(2)+threshold_species(spp));
//        }
//    }
//    // allows for the propects of extra protection for a guild
//    dmatrix slopeGuild(1,Nareas,1,Nguilds);
//    dmatrix interceptGuild(1,Nareas,1,Nguilds);
//    for (area=1; area<=Nareas; area++) {
//        for (iguild = 1; iguild<=Nguilds; iguild++) {
//            slopeGuild(area,iguild) = (minMaxExploitation(2)-minMaxExploitation(1))/(minMaxThreshold(2) - minMaxThreshold(1));
//            interceptGuild(area,iguild) =  minMaxExploitation(2) - slopeGuild(area,iguild)*minMaxThreshold(2);
//             //             cout<<slopeGuild(area,iguild) <<" - "<<interceptGuild(area,iguild)<<endl;
//        }
//    }

//    // guild level calcs
//    // now average the guild biomass values over AssessmentPeriod yrs and then we adjust the exploitation rate
//      exploitationLevelGuild.initialize();
//      for (area=1 ; area<=Nareas; area++) {
//          for (iguild=1; iguild<=Nguilds; iguild++) {
//              catchToDiscardsGuild(area,iguild) = 0;//resets flag to indicate species has not exceeded min threshold
//              for (iassess=1; iassess<=AssessmentPeriod;iassess++){
//                   // calculate the mean biomass and catch over the Assessment period
//                   est_survey_guild_biomass_assessment(area,iguild,yrct) += est_survey_guild_biomass(area,iguild,yrct-iassess+1)/AssessmentPeriod;
//                   for (int ifleet=1;ifleet<=Nfleets;ifleet++) {// catch by fleet over last AssessmentPeriod Years
//                       est_fleet_catch_guild_assessment(area,iguild,ifleet,yrct) += est_fleet_catch_guild_biomass(area,iguild,ifleet,yrct-iassess+1)/AssessmentPeriod;
//                   }
//              }


//              // check to see if average < min threshold or > max threshold and assign new exploitation
//              // otherwise adjust exploitation linearly
//              dvariable biomassLevel = est_survey_guild_biomass_assessment(area,iguild,yrct)/B0_guilds(area,iguild);
//              if (biomassLevel <= minMaxThreshold(1) ) {// min threshold
//                 exploitationLevelGuild(area,iguild) = minMaxExploitation(1);
//              } else if (biomassLevel >= minMaxThreshold(2) ) {// max threshold
//                 exploitationLevelGuild(area,iguild) = minMaxExploitation(2);
//              } else { // linear ramp
//                exploitationLevelGuild(area,iguild) = (biomassLevel * slopeGuild(area,iguild)) + interceptGuild(area,iguild);
//              }
//              if (biomassLevel <= baseline_threshold) { // this is bad. no longer allowed to land caught fish in this guild
//                catchToDiscardsGuild(area,iguild) = 1;
//              }
//             // cout<<yrct<<"  " <<iguild <<"   "<<biomassLevel<<"  "<<exploitationLevelGuild(area,iguild)<<endl;
//          } // guild loop
//          // take the smallest of all recalculated levels. this will be the new level
//         // THIS NEEDS TO CHANGE, NEED TO TARGET FLEET THAT FISH ON GUILD
//           newExploitationLevel = min(exploitationLevelGuild(area));
//         //  cout<<"nlevel = "<<newExploitationLevel<<endl;

//      }  // area loop



//    // check at species level also
//    if (speciesDetection == 1) { // include species detection level in determining rate change
//      // we check for exceedances at the species level. take the mean abundance over last AssessmentPeriod yrs for each species
//      exploitationLevelSpecies.initialize();
//      for (area=1; area<=Nareas; area++){
//          for (spp=1; spp<=Nspecies; spp++){
//              catchToDiscardsSpecies(area,spp) = 0; //resets flag to indicate species has not exceeded min threshold
//              for (iassess=1; iassess<=AssessmentPeriod; iassess++){
//                  // mean of last few years
//                  est_survey_biomass_assessment(area,spp,yrct) +=  est_survey_biomass(area,spp,yrct-iassess+1)/AssessmentPeriod;
//              }
//              // now check for exceedances. if average < min threshold or > max threshold and assign new exploitation
//              // otherwise adjust exploitation linearly
//              dvariable biomassLevel =  est_survey_biomass_assessment(area,spp,yrct)/B0(area,spp);
//              if (biomassLevel <= minMaxThreshold(1) ) {// min threshold
//                 exploitationLevelSpecies(area,spp) = minMaxExploitation(1);
//              } else if (biomassLevel >= (minMaxThreshold(2)+threshold_species(spp)) ) {// max threshold
//                 exploitationLevelSpecies(area,spp) = minMaxExploitation(2) ;
//              } else { // linear ramp
//                exploitationLevelSpecies(area,spp) = (biomassLevel * slopeSpecies(area,spp)) + interceptSpecies(area,spp);
//              }

//              if (biomassLevel <= baseline_threshold) { // this is bad. no longer allowed to land caught fish
//                catchToDiscardsSpecies(area,spp) = 1;
//              }

//          } // spp  loop
//          // if species detection is on then the new level will be determined by the species until we can target fleets based on guild exceedences
//          newExploitationLevel = min(exploitationLevelSpecies(area));
//      } // area loop
//    } // species detection


//     // now we found new exploitation rates we need to act on them and adjust effort to correspond to rate
//     // Also if the scenario is a FixedRate scenario we set exploitation rates all equal to fixed rate in data file, therefore when we
//     // encounter this phase the efort is unchanged

//     // store current exploitation level and set it for next few years until new assessment is due. used as output only
//      // if first time set exploitation_update for the first few years prior to assessment
//       for (area=1; area<=Nareas; area++) {
//          if (t==(Nstepsyr*AssessmentPeriod)) { //first assessment. assign 1st 3 yrs (not effected by assessment) to actual starting value of exploitation
//            for (int iassess=1; iassess <= AssessmentPeriod; iassess++) {
//              exploitation_update(area,iassess) =  minMaxExploitation(2);// starting exploitation, maximum
//            }
//          }
//          // set all subsequent yrs to new exploitation otherwise last few years will revert to original rate/
//          for (int iy = yrct+1; iy<=Nyrs; iy++) {
//           exploitation_update(area,iy) = newExploitationLevel;
//          }
//      }
//     //   cout << exploitation_update << endl;
//       // now we calculate the new effort for each fleet.
//      // Note that all fleets are impacted for any guild exceedance. This can and should change
//      for (area=1 ; area<=Nareas ; area++) {

//         for (int ifleet=1;ifleet<=Nfleets;ifleet++){
//          if ( mean_fishery_q(area,ifleet) < 1e-29){ // a fleet doesn't fish a particluar guild. keep effort same
//               // this will only happen at guild q not fleet q.
//               effort_updated(area,ifleet) = obs_effort(area,ifleet,yrct);
//          } else {
//               effort_updated(area,ifleet) = newExploitationLevel/mean_fishery_q(area,ifleet);
//          }
//         }
//      }
//      // now effort is used just once in initial.calcs() to obtain the Fyr terms for the whole simulation  so
//      // we need to use this new effort and create updated values for Fyr
//      for (area=1; area<=Nareas ; area++) {
//         for (spp=1; spp<=Nspecies; spp++) {
//             for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//                 for (int iy = yrct+1; iy<=Nyrs; iy++) {// set all subsequent yrs to new effort otherwise last few years will revert to original rate
//                  // this will all be updated during next assessment
// //               for (int iassess=1; iassess <= AssessmentPeriod; iassess++) {
//                   obs_effortAssess(area,ifleet,iy) = effort_updated(area,ifleet); // output only
//                   Fyr(area,spp,ifleet,iy) = fishery_q(area,spp,ifleet)*effort_updated(area,ifleet)*effortScaled(area,spp); //Andy Beet
//                 }
//             }
//         }
//      }


// //----------------------------------------------------------------------------------------
// FUNCTION calc_assessment_strategy_Step
// //----------------------------------------------------------------------------------------

// // We enter this loop every AssessmentPeriod years, during the last time period of the year.
// // We are seeing if any of the individual species or complexed have dropped below a threshold (threshold_proportion) given in data file
// // An updated level of effort is then calculated based on the threshold exceeded.
// // Each species or complex is flagged to indicate which have breached the minimum threshold (indicating big trouble).
// // Consequently landings are not allowed and all catch is considered discards. This is dealt with in catch module


//    // guild level calcs
//    // now average the guild biomass values over AssessmentPeriod yrs and then we check to see if the levels exceed some threshhold
//      for (area=1 ; area<=Nareas; area++) {
//          for (iguild=1; iguild<=Nguilds; iguild++) {
//              catchToDiscardsGuild(area,iguild) = 0;//resets flag to indicate species has not exceeded min threshold
//              maxGuildThreshold(area,iguild) = Nthresholds; // set all to maximum worst case is that no change is made to effort
//              for (iassess=1; iassess<=AssessmentPeriod;iassess++){
//                   // calculate the mean biomass and catch over the Assessment period
//                   est_survey_guild_biomass_assessment(area,iguild,yrct) += est_survey_guild_biomass(area,iguild,yrct-iassess+1)/AssessmentPeriod;
//                   for (int ifleet=1;ifleet<=Nfleets;ifleet++) {// catch by fleet over last AssessmentPeriod Years
//                       est_fleet_catch_guild_assessment(area,iguild,ifleet,yrct) += est_fleet_catch_guild_biomass(area,iguild,ifleet,yrct-iassess+1)/AssessmentPeriod;
//                   }
//              }
//              // check to see if average < threshold (threshold_proportion * biomass at equilibrium)
//              for (ithreshold=1; ithreshold<=Nthresholds; ithreshold++) {
//              // test<<yrct<<","<<iguild<<","<<ithreshold<<endl;
//                  if ((est_survey_guild_biomass_assessment(area,iguild,yrct)/B0_guilds(area,iguild)) <= threshold_proportion(ithreshold)) {

//                    maxGuildThreshold(area,iguild) = ithreshold;
//                    if (est_survey_guild_biomass_assessment(area,iguild,yrct)/B0_guilds(area,iguild) <= baseline_threshold) { // this is case where most severe threshold is breached
//                        // all catch => discards and nothing can be landed. create binary vector
//                        catchToDiscardsGuild(area,iguild) = 1;
//                     }

//                     // dont need to keep going for this guild since we've found the most severe case
//                     break;
//                   }

//               }// threshold loop
//   //                   cout<<catchToDiscardsGuild(area,iguild)<<endl;

//          } // guild loop
//        maxThreshold(area) =  min(maxGuildThreshold(area));
// //       cout<<maxThreshold(area)<<endl;
//      }  // area loop



//     // check for species falling below threshold

//    if (speciesDetection == 1) { // include species detection level in determining rate change
//      // we check for exceedances at the species level
//      // take the mean abundance over last AssessmentPeriod yrs for each species
//      for (area=1; area<=Nareas; area++){
//          for (spp=1; spp<=Nspecies; spp++){
//              catchToDiscardsSpecies(area,spp) = 0; //resets flag to indicate species has not exceeded min threshold
//              maxSpeciesThreshold(area,spp) = Nthresholds; // set all to safe level
//              for (iassess=1; iassess<=AssessmentPeriod; iassess++){
//                  // mean of last few years
//                  est_survey_biomass_assessment(area,spp,yrct) +=  est_survey_biomass(area,spp,yrct-iassess+1)/AssessmentPeriod;
//              }
//              // now check for exceedances
//              for (ithreshold=1; ithreshold<=Nthresholds; ithreshold++) {
//                  if ((est_survey_biomass_assessment(area,spp,yrct)/B0(area,spp)) <= (threshold_proportion(ithreshold)+threshold_species(spp))) {
//                     maxSpeciesThreshold(area,spp) = ithreshold;
//                     if (est_survey_biomass_assessment(area,spp,yrct)/B0(area,spp)  <= baseline_threshold ) {
//                        // all catch => discards and nothing can be landed. create binary vector
//                        catchToDiscardsSpecies(area,spp) = 1;
//                     }
//                     // dont need to keep going for this species since we've found the most severe case
//                     break;
//                   }

//              }// threshold loop
//          //    cout<<spp<<"-"<<maxSpeciesThreshold(area,spp)<<endl;
//          } // spp  loop
//          maxThreshold(area) = min(maxThreshold(area),min( maxSpeciesThreshold(area)));
//      // cout<<maxThreshold<<endl;
//      // cout<<endl;
//      } // area loop
//   } // species detection


//      // now we have checked for exceedences we need to act on them.
//      // calculate the new exploitation rate and then the new value of Effort.
//      // note that if maxThreshold = Nthresholds we revert to max exploitation.
//      // Also if the scenario is a FixedRate scenario we set exploitation rates all equal to fixed rate in data file, therefore when we
//      // encounter this phase the efort is unchanged

//      // store current exploitation level and set it for next few years until new assessment is due. used as output only
//       // if first time set exploitation_update for the first few years prior to assessment
//       for (area=1; area<=Nareas; area++) {
//          if (t==(Nstepsyr*AssessmentPeriod)) { //first assessment. assign 1st 3 yrs (not effected by assessment) to actual starting value of exploitation
//            for (int iassess=1; iassess <= AssessmentPeriod; iassess++) {
//              exploitation_update(area,iassess) =  exploitation_levels(Nthresholds);
//            }
//          }
//           // set all subsequent yrs to new exploitation otherwise last few years will revert to original rate/
//          for (int iy = yrct+1; iy<=Nyrs; iy++) {
//           exploitation_update(area,iy) = exploitation_levels(maxThreshold(area));
//          }
//      }


//      // now we calculate the new effort for each fleet.
//      // Note that all fleets are impacted for any guild exceedance. This can and should change
//      for (area=1 ; area<=Nareas ; area++) {

//         for (int ifleet=1;ifleet<=Nfleets;ifleet++){
//          if ( mean_fishery_q(area,ifleet) < 1e-29){ // a fleet doesn't fish a particluar guild. keep effort same
//               // this will only happen at guild q not fleet q.
//               effort_updated(area,ifleet) = obs_effort(area,ifleet,yrct);
//          } else {
//               effort_updated(area,ifleet) = exploitation_levels(maxThreshold(area))/mean_fishery_q(area,ifleet);
//          }
//         }
//      }
//      // now effort is used just once in initial.calcs() to obtain the Fyr terms for the whole simulation  so
//      // we need to use this new effort and create updated values for Fyr
//      for (area=1; area<=Nareas ; area++) {
//         for (spp=1; spp<=Nspecies; spp++) {
//             for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
//                 for (int iy = yrct+1; iy<=Nyrs; iy++) {// set all subsequent yrs to new effort otherwise last few years will revert to original rate
//                  // this will all be updated during next assessment
// //               for (int iassess=1; iassess <= AssessmentPeriod; iassess++) {
//                   obs_effortAssess(area,ifleet,iy) = effort_updated(area,ifleet); // output only
//                   Fyr(area,spp,ifleet,iy) = fishery_q(area,spp,ifleet)*effort_updated(area,ifleet)*effortScaled(area,spp); //Andy Beet
//                 }
//             }
//         }
//      }




// //----------------------------------------------------------------------------------------
// FUNCTION write_simout_KRAKEN
// //----------------------------------------------------------------------------------------

//   //send simulated biomass and catch data to csv for use in production model (KRAKEN)
//       ofstream simout("simKraken.csv"); // for Kraken
//       simout<<"rseed,"<<rseed<<endl;
//       simout<<"BIOMASS"<<endl;
//       for (area=1; area<=Nareas; area++){
//    	    for(spp=1; spp<=Nspecies; spp++){
//           simout<<"name_"<<spp;
//           for(yr=1; yr<=Nyrs; yr++){
//              simout<<","<<est_survey_biomass(area,spp,yr);
//           }
//         simout<<endl;
//         }
//       }
//       simout<<"CATCH"<<endl;
//       for (area=1; area<=Nareas; area++){
//    	    for(spp=1; spp<=Nspecies; spp++){
//           simout<<"name_"<<spp;
//           for(yr=1; yr<=Nyrs; yr++){
//               simout<<","<<est_catch_biomass(area,spp,yr);
//           }
//         simout<<endl;
//         }
//       }


// //----------------------------------------------------------------------------------------
// FUNCTION write_outDarwin
// //----------------------------------------------------------------------------------------
//   //send simulated biomass and catch in MSE darwinian runs
//       clock_t elapsedTime2  = clock() - startTime;
//       std::stringstream fileIndicesNames,part2Name;
//       fileIndicesNames << rseed;
//       fileIndicesNames << time(&baseTime);
//       part2Name << elapsedTime2;
//       fileIndicesNames << "_";
//       fileIndicesNames << part2Name.str();
//       fileIndicesNames << "simDarwin.text";

//       std::string fileNameIndex = fileIndicesNames.str();

//       ofstream outDarwin(fileNameIndex.c_str());

//       outDarwin<<"Nyrs\n"<<Nyrs<<endl;
//       outDarwin<<"avByr\n"<<avByr<<endl;
//       outDarwin<<"guildMembers\n"<<guildMembers<<endl;
//       outDarwin<<"est_catch_biomass\n"<<est_catch_biomass<<endl;
//       outDarwin<<"est_survey_biomass\n"<<est_survey_biomass<<endl;

//       outDarwin<<"manually exiting at end of procedure section....\n"<<endl;

// //----------------------------------------------------------------------------------------
// FUNCTION write_outIndices
// //----------------------------------------------------------------------------------------
//   //send simulated indices and metrics for use in MSE type output
//       clock_t elapsedTime2  = clock() - startTime;
//       //int rnN = (int)startTime;
//     //random number to attach to filenme
//      // random_number_generator rngInd(rnN);
//      // dvector rnFile(1,1);
//      // rnFile.fill_randu(rngInd);
//      // std::cout << rnFile << std::endl;



//       std::stringstream fileIndicesNames,part2Name;
//       fileIndicesNames << rseed;
//       fileIndicesNames << time(&baseTime);
//      // fileIndicesNames << rnFile;
//       part2Name << elapsedTime2;
//       fileIndicesNames << "_";
//       fileIndicesNames << part2Name.str();
//       fileIndicesNames << "simIndices.txt";

//       std::string fileNameIndex = fileIndicesNames.str();

//       ofstream outIndices(fileNameIndex.c_str());
      
//       // diagnose why files nort written when run in parallel
//      // std::ofstream checkFileName;
//      // checkFileName.open( "checkFile.txt",std::ios_base::app);
//      // checkFileName << fileNameIndex << endl; // test to see why not all output files names are present
      
//       outIndices<<"rseed\n"<<rseed<<endl;
//       outIndices<<"Nyrs\n"<<Nyrs<<endl;
//       outIndices<<"Nstepsyr\n"<<Nstepsyr<<endl;
//       outIndices<<"Nguilds\n"<<Nguilds<<endl;
//       outIndices<<"Nfleets\n"<<Nfleets<<endl;
//       outIndices<<"avByr\n"<<avByr<<endl;
//       outIndices<<"catch_biomass\n"<<catch_biomass<<endl;
//       outIndices<<"obs_effort\n"<<obs_effort<<endl;
//       outIndices<<"est_fleet_catch_biomass\n"<<est_fleet_catch_biomass<<endl;
//       outIndices<<"est_fleet_catch_guild_biomass\n"<<est_fleet_catch_guild_biomass<<endl;
//       outIndices<<"est_catch_guild_biomass\n"<<est_catch_guild_biomass<<endl;
//       outIndices<<"est_catch_biomass\n"<<est_catch_biomass<<endl;
//       outIndices<<"est_survey_biomass\n"<<est_survey_biomass<<endl;
//       outIndices<<"est_survey_guild_biomass\n"<<est_survey_guild_biomass<<endl;
//       outIndices<<"B0\n"<<B0<<endl;
//       outIndices<<"B0_guilds\n"<<B0_guilds<<endl;
//       outIndices<<"guildMembers\n"<<guildMembers<<endl;
//       outIndices<<"Nthresholds\n"<<Nthresholds<<endl;
//       outIndices<<"minExploitation\n"<<minExploitation<<endl;
//       outIndices<<"maxExploitation\n"<<maxExploitation<<endl;
//       //outIndices<<"threshold_proportion\n"<<threshold_proportion<<endl;
//       //outIndices<<"exploitation_levels\n"<<exploitation_levels<<endl;
//       outIndices<<"threshold_species\n"<<threshold_species<<endl;
//       outIndices<<"AssessmentPeriod\n"<<AssessmentPeriod<<endl;
//       outIndices<<"SpeciesDetection\n"<<speciesDetection<<endl;
//       outIndices<<"AssessmentOn\n"<<AssessmentOn<<endl;
//       outIndices<<"index_Simpsons_N\n"<<index_Simpsons_N<<endl;
//       outIndices<<"index_Simpsons_Nrecip\n"<<index_Simpsons_Nrecip<<endl;
//       outIndices<<"index_Simpsons_C\n"<<index_Simpsons_C<<endl;
//       outIndices<<"index_Simpsons_Crecip\n"<<index_Simpsons_Crecip<<endl;
//       outIndices<<"index_LFI_Biomass\n"<<index_LFI_Biomass<<endl;
//       outIndices<<"index_LFI_Catch\n"<<index_LFI_Catch<<endl;
//       outIndices<<"index_LFI_N\n"<<index_LFI_N<<endl;
//       outIndices<<"index_predToPreyRatio\n"<<index_predToPreyRatio<<endl;
//       outIndices<<"index_plankToPiscRatio\n"<<index_plankToPiscRatio<<endl;
//       outIndices<<"index_stdev_catch\n"<<index_stdev_catch<<endl;
//       outIndices<<"index_stdev_biomass\n"<<index_stdev_biomass<<endl;
//       outIndices<<"index_status_species\n"<<index_status_species<<endl;
//       outIndices<<"index_status_guild\n"<<index_status_guild<<endl;


//       outIndices<<"manually exiting at end of procedure section....\n"<<endl;




// //----------------------------------------------------------------------------------------
// FUNCTION write_outDiagnostics
// //----------------------------------------------------------------------------------------
//      //send all outputs to file for plotting
//       clock_t elapsedTime  = clock() - startTime;

//       std::stringstream fileNames,part1Name;
//       fileNames << rseed;
//       fileNames << time(&baseTime);
//       part1Name << elapsedTime;
//       fileNames << "_";
//       fileNames << part1Name.str();
//       fileNames << "simDiagnostics.out";

//       std::string fileName = fileNames.str();

//       ofstream outDiagnostics(fileName.c_str());

//       outDiagnostics<<"rseed\n"<<rseed<<endl;
//       outDiagnostics<<"rectype (1=gamma/'Ricker' eggprod, 2=Deriso-Schnute SSB, 3=SSB gamma, 4=SSB Ricker, 5=SSB Beverton Holt, 9=avg+dev)\n"<<rectype<<endl;
//       outDiagnostics<<"recruitment_alpha\n"<<recruitment_alpha<<endl;
//       outDiagnostics<<"recruitment_shape\n"<<recruitment_shape<<endl;
//       outDiagnostics<<"recruitment_beta\n"<<recruitment_beta<<endl;
//       outDiagnostics<<"Nyrs\n"<<Nyrs<<endl;
//       outDiagnostics<<"Nstepsyr\n"<<Nstepsyr<<endl;
//       outDiagnostics<<"stochrec\n"<<stochrec<<endl;
//       outDiagnostics<<"recsigma\n"<<recsigma<<endl;
//       outDiagnostics<<"recruitment\n"<<recruitment<<endl;
//       outDiagnostics<<"SSB\n"<<SSB<<endl;
//       outDiagnostics<<"avByr\n"<<avByr<<endl;
//       outDiagnostics<<"M2\n"<<M2<<endl;
//       outDiagnostics<<"F\n"<<F<<endl;
//       outDiagnostics<<"Z\n"<<Z<<endl;
//       outDiagnostics<<"N\n"<<N<<endl;
//       outDiagnostics<<"eaten_biomass\n"<<eaten_biomass<<endl;
//       outDiagnostics<<"discard_biomass\n"<<discard_biomass<<endl;
//       outDiagnostics<<"otherDead_biomass\n"<<otherDead_biomass<<endl;
//       outDiagnostics<<"total_biomass\n"<<total_biomass<<endl;
//       outDiagnostics<<"fleet_catch_biomass\n"<<fleet_catch_biomass<<endl;
//       outDiagnostics<<"catch_biomass\n"<<catch_biomass<<endl;
//       outDiagnostics<<"est_fleet_catch_biomass\n"<<est_fleet_catch_biomass<<endl;
//       outDiagnostics<<"est_fleet_catch_guild_biomass\n"<<est_fleet_catch_guild_biomass<<endl;
//       outDiagnostics<<"est_catch_guild_biomass\n"<<est_catch_guild_biomass<<endl;
//       outDiagnostics<<"est_catch_biomass\n"<<est_catch_biomass<<endl;
//       outDiagnostics<<"est_survey_biomass\n"<<est_survey_biomass<<endl;
//       outDiagnostics<<"est_survey_guild_biomass\n"<<est_survey_guild_biomass<<endl;
//       outDiagnostics<<"obs_survey_biomass\n"<<obs_survey_biomass<<endl;
//       outDiagnostics<<"obs_catch_biomass\n"<<obs_catch_biomass<<endl;
//       outDiagnostics<<"est_survey_guild_biomass_assessment\n"<< est_survey_guild_biomass_assessment<<endl;
//       outDiagnostics<<"predation_mortality\n"<<predation_mortality<<endl;
//       outDiagnostics<<"fishing_mortality\n"<<fishing_mortality<<endl;
//       outDiagnostics<<"predation_mortality_size\n"<<predation_mortality_size<<endl;
//       outDiagnostics<<"fishing_mortality_size\n"<<fishing_mortality_size<<endl;
//       outDiagnostics<<"B0\n"<<B0<<endl;
//       outDiagnostics<<"B0_guilds\n"<<B0_guilds<<endl;
//       outDiagnostics<<"Nguilds\n"<<Nguilds<<endl;
//       outDiagnostics<<"guildMembers\n"<<guildMembers<<endl;
//       outDiagnostics<<"Nthresholds\n"<<Nthresholds<<endl;
//       outDiagnostics<<"threshold_proportion\n"<<threshold_proportion<<endl;
//       outDiagnostics<<"exploitation_levels\n"<<exploitation_levels<<endl;
//       outDiagnostics<<"exploitation_update\n"<<exploitation_update<<endl;
//       outDiagnostics<<"threshold_species\n"<<threshold_species<<endl;
//       outDiagnostics<<"AssessmentPeriod\n"<<AssessmentPeriod<<endl;
//       outDiagnostics<<"SpeciesDetection\n"<<speciesDetection<<endl;
//       outDiagnostics<<"AssessmentOn\n"<<AssessmentOn<<endl;
//       outDiagnostics<<"index_ExploitationRate\n"<<index_ExploitationRate<<endl;
//       outDiagnostics<<"index_SystemExploitationRate\n"<<index_SystemExploitationRate<<endl;
//       outDiagnostics<<"index_Simpsons_N\n"<<index_Simpsons_N<<endl;
//       outDiagnostics<<"index_Simpsons_Nrecip\n"<<index_Simpsons_Nrecip<<endl;
//       outDiagnostics<<"index_Simpsons_C\n"<<index_Simpsons_C<<endl;
//       outDiagnostics<<"index_Simpsons_Crecip\n"<<index_Simpsons_Crecip<<endl;
//       outDiagnostics<<"index_status_species\n"<<index_status_species<<endl;
//       outDiagnostics<<"index_status_guild\n"<<index_status_guild<<endl;
//       outDiagnostics<<"LFI_threshold\n"<<LFI_threshold<<endl;
//       outDiagnostics<<"index_LFI_Biomass\n"<<index_LFI_Biomass<<endl;
//       outDiagnostics<<"index_LFI_Catch\n"<<index_LFI_Catch<<endl;
//       outDiagnostics<<"index_LFI_N\n"<<index_LFI_N<<endl;
//       outDiagnostics<<"index_predToPreyRatio\n"<<index_predToPreyRatio<<endl;
//       outDiagnostics<<"index_plankToPiscRatio\n"<<index_plankToPiscRatio<<endl;
//       outDiagnostics<<"index_stdev_catch\n"<<index_stdev_catch<<endl;
//       outDiagnostics<<"index_stdev_biomass\n"<<index_stdev_biomass<<endl;

//       outDiagnostics<<"\npin file inputs\n"<<endl;
//       outDiagnostics<<"yr1N\n"<<yr1N<<endl;
//       //cout<<"avg_F\n"<<avg_F<<endl;
//       //cout<<"F_devs\n"<<F_devs<<endl;
//       outDiagnostics<<"survey_q\n"<<survey_q<<endl;
//       outDiagnostics<<"surv_sigma\n"<<surv_sigma<<endl;
//       outDiagnostics<<"catch_sigma\n"<<catch_sigma<<endl;
//       outDiagnostics<<"\n data input time series \n"<<endl;
//       outDiagnostics<<"obs_effort\n"<<obs_effort<<endl;
//       outDiagnostics<<"obs_effortAssess\n"<<obs_effortAssess<<endl;
//       outDiagnostics<<"obs_temp\n"<<obs_temp<<endl;
//       outDiagnostics<<"manually exiting at end of procedure section....\n"<<endl;



//----------------------------------------------------------------------------------------
FUNCTION calculate_predicted_values
//----------------------------------------------------------------------------------------

// This function calculates the predicted values for the survey indices for use in the objective function
// Function by G. Fay, December 2021
// some argument for putting this in the objective function calcs, but seems less busy to do here

    // assume survey covers all areas, and uses average biomass for year, rathter than timing-specific suruvey
    pred_survey_index.initialize();

    for (int i=1;i<=Nsurvey_obs;i++) {
      int survey = obs_survey_biomass(i,1);
      int year = obs_survey_biomass(i,2);
      int spp = obs_survey_biomass(i,3);
      for (area=1; area<=Nareas; area++){
        for (int ilen=1;ilen<=Nsizebins;ilen++) {
                  //if (year == 80 && spp == 10 && survey ==2) cout << ilen << " " <<B_tot(area,spp,year,ilen) << " " << survey_sel(survey,spp,ilen) << " " << survey_q(survey,spp) << endl; 
            pred_survey_index(i) +=  B_tot(area,spp,year,ilen)*survey_sel(survey,spp,ilen)*survey_q(survey,spp)/Nstepsyr; 
        }
      }       
    }


//----------------------------------------------------------------------------------------
FUNCTION evaluate_the_objective_function
//----------------------------------------------------------------------------------------

//Original placeholder by S. Gaichas, updated by J. Boucher

// revised G. Fay April 2021
// make use of new data structures

//   resid_catch.initialize();
//     resid_bio.initialize();
//     totcatch_fit.initialize();
//     totbio_fit.initialize();
//     catchcomp_fit.initialize();
//     objfun_areaspp.initialize();

////read in survey and catch observations from .dat file
//// SKG: for now, sub in Nsurveys for the unused Nareas dimension for input survey indices and comps
//  //init_3darray obs_survey_biomass(1,Nareas,1,Nspecies,1,Nyrs)  	//spring or fall? units needed
//  init_3darray obs_survey_biomass(1,Nsurveys,1,Nspecies,1,Nyrs)  	//input spring, fall separately, in tons
//  init_3darray obs_catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)  	//total catch in tons
//  init_3darray obs_effort(1,Nareas,1,Nfleets,1,Nyrs)  	//standardized effort units needed
//  //init_4darray for survey size comp by area, species, year?
//  //init_4darray obs_survey_size(1,Nsurveys,1,Nspecies,1,Nyrs,1,Nsizebins)  //numbers, uncomment when in dat file
//  //init_5darray for catch at size by area, species, fleet, year?


  dvariable eps = 1.e-07;

  cout << "starting commercial catch nll" << endl;

//Commercial Catch
    // Ncatch_obs
    pred_catch_biomass.initialize();
    resid_catch.initialize();
    nll_catch.initialize(); 
    for (int i=1;i<=Ncatch_obs;i++) {
      //cout << "cat obs " << i << endl;
      int fleet = obs_catch_biomass(i,1);
      int area = obs_catch_biomass(i,2);
      int year = obs_catch_biomass(i,3);
      int spp = obs_catch_biomass(i,4);
      dvariable value = obs_catch_biomass(i,5)+eps;
      dvariable cv = obs_catch_biomass(i,6);
      // if (fleet==0) // add case when catch is aggregated over fleets (fleet = 0 in data file)
      pred_catch_biomass(i) = fleet_catch_biomass(area,spp,fleet,year); //predicted value for this data point
      resid_catch(i) = log(value/(pred_catch_biomass(i)+eps));
      nll_catch(i) = dlnorm(value, log(pred_catch_biomass(i)+eps), cv);
    }
   
  cout << "done commercial catch nll" << endl;

  cout << "starting commercial catch at length nll" << endl;

// Commercial catch at length 
   //Ncatch_size_obs
  int j=0;
  pred_catch_size.initialize();
  nll_catch_size.initialize();
  for (int i=1;i<=Ncatch_size_obs;i++) {
    
     //cout << "obs row " << i << endl;
  
     int fleet = obs_catch_size(i,1); //cout << "fleet" << fleet << endl;
     int area = obs_catch_size(i,2);  //cout << "area" << area << endl;
     int year = obs_catch_size(i,3);  //cout << "year" << year << endl;
     int spp = obs_catch_size(i,4);  //cout << "spp" << spp << endl;
     int type = obs_catch_size(i,5);  //cout << "type" << type << endl;//not yet used
     int effN = obs_catch_size(i,6);   //cout << "effN" << effN << endl;
     dvar_vector Lobs(1,Nsizebins);
     Lobs.initialize();
     for (int ilen=1;ilen<=Nsizebins;ilen++) Lobs(ilen) = obs_catch_size(i,6+ilen) + 0.0001;
     Lobs = Lobs/sum(Lobs);
     dvar_vector Lpred(1,Nsizebins);
     Lpred.initialize();
     for (int ilen=1;ilen<=Nsizebins;ilen++)
      Lpred(ilen) = Cfl_tot(area,spp,fleet,year,ilen); // predicted catch at length for this observation
     Lpred = (eps+Lpred)/sum(eps + Lpred);
     
     //cout << "Lpred" << Lpred << endl;
     
     for (int ilen=1;ilen<=Nsizebins;ilen++) {
      if (Lobs(ilen) > 0)
       {
        j+=1;
        // create table for data base of sizes
        // jth row of this table
        pred_catch_size(j) = Lpred(ilen);  //change for better storage table
        nll_catch_size(j) = -1.*effN*Lobs(ilen)*log(Lpred(ilen)/Lobs(ilen));
       } 
     }
  }
 
  cout << "done commercial catch at length nll" << endl;

//Survey Indices of abundance
    // Nsurvey_obs
  resid_survey.initialize();
  nll_survey.initialize();
    for (int i=1;i<=Nsurvey_obs;i++) {
      //cout << i << endl;
      //if (i==881) cout << obs_survey_biomass(881) << endl;
      int survey = obs_survey_biomass(i,1);
      int year = obs_survey_biomass(i,2);
      int spp = obs_survey_biomass(i,3);
      dvariable value = obs_survey_biomass(i,4)+eps;
      dvariable cv = obs_survey_biomass(i,5);
      //predicted value now calculated in function 'calculate_predicted_values()'
      //pred_survey_index(i) = est_survey_biomass(survey,spp,year); //survey is area here! need to change survey definitions
      resid_survey(i) = log(value/(pred_survey_index(i)+eps));
      nll_survey(i) = dlnorm(value, log(pred_survey_index(i)+eps), cv);
    }
   
  cout << "done survey abundance nll" << endl;
//  if(isinf(value(sum(nll_survey)))) {
//    cout << " INFINITE OBJ FUN" << endl;
//    gavjunk << survey_q << endl;
//    gavjunk << ln_survey_q << endl;
//    gavjunk << "survey biomass data, predicted, residual, nll" << endl;
//  for (int i=1;i<=Nsurvey_obs;i++)
//    gavjunk << obs_survey_biomass(i) << " " << pred_survey_index(i) << " " << resid_survey(i) << " " << nll_survey(i) << endl;
//  exit(-1);
//  }
  //

 //Survey Catch-at-length
    //Nsurvey_size_obs
   j=0;
   pred_survey_size.initialize();
   nll_survey_size.initialize();
   for (int i=1;i<=Nsurvey_size_obs;i++) {
      int survey = obs_survey_size(i,1);
      int year = obs_survey_size(i,2);
      int spp = obs_survey_size(i,3);
      int type = obs_survey_size(i,4);  //not yet used
      int effN = obs_survey_size(i,5);
      dvar_vector Lobs(1,Nsizebins);
      Lobs.initialize();
      for (int ilen=1;ilen<=Nsizebins;ilen++) Lobs(ilen) = obs_survey_size(i,5+ilen) + 0.0001;
      Lobs = Lobs/sum(Lobs);
      dvar_vector Lpred(1,Nsizebins);
      Lpred.initialize();
      for (int area=1;area<=Nareas;area++)
       Lpred += N_tot(area,spp,year);
      Lpred = eps + survey_q(survey,spp)*elem_prod(survey_sel(survey,spp),Lpred)/Nstepsyr; //est_survey_size(survey, year, spp, ilen);// predicted survey at length for this observation
      Lpred = Lpred/sum(Lpred);
      for (int ilen=1;ilen<=Nsizebins;ilen++) {
       if (Lobs(ilen) > 0)
        {
         j+=1;
         // create table for data base of sizes
         // jth row of this table
         pred_survey_size(j) = Lpred(ilen);  //change for better storage table
         nll_survey_size(j) = -1.*effN*Lobs(ilen)*log(Lpred(ilen)/Lobs(ilen));
         //cout << j << " " << nll_survey_size(j) << endl;
        } 
      }
   }

  cout << "done survey size comp nll" << endl;

// Survey Prey proportions 
   j=0;
   pred_dietprop.initialize();
   nll_dietprop.initialize();
   for (int i=1;i<=Ndietprop_obs;i++) {
      int survey = obs_dietprop(i,1);
      int year = obs_dietprop(i,2);
      int spp = obs_dietprop(i,3);
      int size = obs_dietprop(i,4);  
      int effN = obs_dietprop(i,5);
      dvar_vector Pobs(1,Nprey);
      Pobs.initialize();
      for (int ilen=1;ilen<=Nprey;ilen++) Pobs(ilen) = obs_dietprop(i,5+ilen) + 0.0001;
      Pobs = Pobs/sum(Pobs);
      dvar_vector Ppred(1,Nprey);
      Ppred.initialize();
      for (int ilen=1;ilen<=Nprey;ilen++)
       for (int area=1;area<=Nareas;area++)
        Ppred(ilen) += eps + biomass_prey_avail_no_size(area,spp,year,size,ilen);
      Ppred = Ppred/sum(Ppred);
     
      for (int ilen=1;ilen<=Nprey;ilen++) {
       if (Pobs(ilen) > 0)
        {
         j+=1;
         // create table for data base of sizes
         // jth row of this table
         pred_dietprop(j) = Ppred(ilen);  //change for better storage table
         nll_dietprop(j) = -1.*effN*Pobs(ilen)*log(Ppred(ilen)/Pobs(ilen));
        } 
      }
   }

  cout << "done survey prey proportions nll" << endl;


// Recruitment penalty

   j = 0;
   dvar_vector resid(1,Nareas*Nspecies*(Nyrs-1));
   dvar_vector sdrec(1,Nareas*Nspecies*(Nyrs-1));
   resid.initialize();
   sdrec.initialize();
   nll_recruit.initialize();
   for (int area=1;area <=Nareas;area++) {
    for (int spp=1;spp <=Nspecies;spp++) {
    for (int year=2;year <=Nyrs;year++) {
      j +=1;
      recdev(j) = recruitment_devs(area,spp,year);
      resid(j) = recdev(j); //+0.5*square(recsigma(area,spp));
      sdrec(j) = recsigma(area,spp);
//      nll_recruit(j) = dnorm(recdev(j)+0.5*square(recsigma(area,spp)),recsigma(area,spp),true);
    }}}
      //dvariable sigma_use = recsigma(area,spp);
   nll_recruit = dnorm(resid,sdrec);
//*/

  cout << "done recruitment nll" << endl;


// Calc objective function
   objfun = 0.;
   objfun += sum(nll_survey);
   objfun += sum(nll_survey_size);
   objfun += sum(nll_catch);
   objfun += sum(nll_catch_size);
   objfun += sum(nll_dietprop);
   objfun += nll_recruit;  

//F_devs(area,fleet,iyr)
   for (int area=1;area <=Nareas;area++) {
    for (int spp=1;spp <=Nfleets;spp++) {
     objfun += dnorm(F_devs(area,spp),1.0); }}
   
   cout << "nll_survey: " << sum(nll_survey) << endl;
   cout << "nll_survey_size: " << sum(nll_survey_size) << endl;
   cout << "nll_catch: " << sum(nll_catch) << endl;
   cout << "nll_catch_size: " << sum(nll_catch_size) << endl;
   cout << "nll_dietprop: " << sum(nll_dietprop) << endl;
   cout << "nll_recruit: " << nll_recruit << endl;
   cout << "nll_total: " << objfun << endl;


  // //est and observed survey biomass and fishery catch are 3darrays(area,spp,yr)
  // //fit matrices are area by spp

  //  resid_catch.initialize();
  //  resid_bio.initialize();
  //  totcatch_fit.initialize();
  //  totbio_fit.initialize();
  //  objfun_areaspp.initialize();

  // for (area=1; area<=Nareas; area++){
  // 	for(spp=1; spp<=Nspecies; spp++){

  //      resid_catch(area,spp) = log(obs_catch_biomass(area,spp)+o)-log(est_catch_biomass(area,spp)+o);
  //      totcatch_fit(area,spp) = norm2(resid_catch(area,spp));

  //      resid_bio(area,spp) = log(obs_survey_biomass(area,spp)+o)-log(est_survey_biomass(area,spp)+o);
  //      totbio_fit(area,spp) = norm2(resid_bio(area,spp));
  //   }
  // }
  // //cout<<"resid_catch\n"<<resid_catch<<endl;
  // //cout<<"totcatch_fit\n"<<totcatch_fit<<endl;
  // //cout<<"totbio_fit\n"<<totbio_fit<<endl;

  // objfun_areaspp = totcatch_fit + totbio_fit;
  // //cout<<"objfun_areaspp\n"<<objfun_areaspp<<endl;

  // objfun = sum(objfun_areaspp);

// //=======================================================================================
// RUNTIME_SECTION
// //=======================================================================================
//   convergence_criteria 1.e-3 ,  1.e-4
//   maximum_function_evaluations 1000

//=======================================================================================
TOP_OF_MAIN_SECTION
//=======================================================================================
 // arrmblsize = 8000000;  //Increase amount of available dvar memory
 // gradient_structure::set_CMPDIF_BUFFER_SIZE(6000000);
 // gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);

// Try to prevent *.tmp files from being created since no derivatives are needed. Purely simulation
 arrmblsize = 800000000;
// gradient_structure::set_NO_DERIVATIVES();
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1200000000);
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(600000000);

//=======================================================================================
REPORT_SECTION
//=======================================================================================

  report << "Nstepsyr model timesteps per year " << endl;  
  report << Nstepsyr << endl;
  report << "EstNsize Estimated total numbers of fish " << endl;
  report << N << endl;
  report << "EstBsize Estimated total biomass of fish " << endl;
  report << B << endl;
  report << "EstRec Estimated recruitment " << endl;
  report << recruitment << endl;
  report << "EstFsize Estimated fishing mortality " << endl;
  report << F << endl;
  report << "EstM2size Estimated predation mortality " << endl;
  report << M2 << endl;
  report << "EstM1size Estimated residual mortality " << endl;
  report << M1 << endl;  
  report << "EstOtherFood Estimated other food " << endl;
  report << otherFood << endl;    
  report << "table of fits to survey" << endl;
  report << "survey biomass data, predicted, residual, nll" << endl;
  for (int i=1;i<=Nsurvey_obs;i++)
    report << obs_survey_biomass(i) << " " << pred_survey_index(i) << " " << resid_survey(i) << " " << nll_survey(i) << endl;
  // report << "EstSurvB Estimated survey biomass of fish " << endl;
  // report << est_survey_biomass << endl;
  // report << "ObsSurvB Observed survey biomass of fish " << endl;
  // report << obs_survey_biomass << endl;
  report << "table of fits to catch" << endl;
  report << "catch data, predicted, residual, nll" << endl;
  for (int i=1;i<=Ncatch_obs;i++)
    report << obs_catch_biomass(i) << " " << pred_catch_biomass(i) << " " << resid_catch(i) << " " << nll_catch(i) << endl;
  // report << "EstCatchB Estimated catch biomass of fish " << endl;
  // report << est_catch_biomass << endl;
  // report << "ObsCatchB Observed catch biomass of fish " << endl;
  // report << obs_catch_biomass << endl;
   report << "pred_catch_size" << endl;
   report << pred_catch_size << endl;
   report << "nll_catch_size" << endl;
   report << nll_catch_size << endl;
   report << "pred_survey_size" << endl;
   report << pred_survey_size << endl;
   report << "nll_survey_size" << endl;
   report << nll_survey_size << endl;
   report << "pred_dietprop" << endl;
   report << pred_dietprop << endl;
   report << "nll_dietprop" << endl;
   report << nll_dietprop << endl;
   report << "fishsel" << endl;
   for (int i=1;i<=Nspecies;i++) {
    for (int j=1;j<=Nfleets;j++) {
     report << i << " " << j << " " << fishsel(1,i,j) << endl;
    }
   }
   report << "survey_sel" << endl;
   for (int i=1;i<=Nspecies;i++) {
    for (int j=1;j<=Nsurveys;j++) {
     report << i << " " << j << " " << survey_sel(j,i) << endl;
    }
   }
  report << "Fyr" << endl;
  for (int i=1;i<=Nspecies;i++) {
   for (int j=1;j<=Nfleets;j++) {
     report << i << " " << j << " " << Fyr(1,i,j) << endl;
  }}

  report << "Year Species SSB" << endl;

  for (int yrct = 1; yrct <= Nyrs; yrct++) {
  for (int spp = 1; spp <= Nspecies; spp++) {
    dvariable total_ssb = 0.0;
    for (int area = 1; area <= Nareas; area++) {
      total_ssb += SSB(area, spp)(yrct);
    }
    report << yrct << " " << spp << " " << total_ssb << endl;
  }
  }
//// write full time series of predicted data streams to pmse_prevals.out

//////// full table of predicted survey index
    pmse_predvals << "full time series of predicted survey index" << endl;
    pmse_predvals << "survey year spp area pred_survey" << endl;
      dvariable pred_survey_index2;
      for (int survey=1; survey<=Nsurveys;survey++) {
      for (int year=1;year<=Nyrs;year++) {
      for (int spp=1;spp<=Nspecies;spp++) {
      for (int area=1; area<=Nareas; area++){
        pred_survey_index2 = 0.;
        for (int ilen=1;ilen<=Nsizebins;ilen++) {
            //pred_survey_index2 +=  B_tot(area,spp,year,ilen)*survey_sel(survey,spp,ilen)*survey_q(survey,spp)/Nstepsyr; 
            pred_survey_index2 +=  B_tot(area,spp,year,ilen)*survey_sel(survey,spp,ilen)/Nstepsyr; 
        }
        pmse_predvals << survey << " " << year << " " << spp << " " << area << " " << pred_survey_index2 << endl;
      }}}}       

//////// full table of predicted catch
    pmse_predvals << "full time series of predicted catch" << endl;
    pmse_predvals << "fleet year spp area pred_catch" << endl;
      for (int fleet=1; fleet<=Nfleets;fleet++) {
      for (int year=1;year<=Nyrs;year++) {
      for (int spp=1;spp<=Nspecies;spp++) {
      for (int area=1; area<=Nareas; area++){
        if (indicator_fishery_q(area,fleet,spp) == 1)
         pmse_predvals << fleet << " " << year << " " << spp << " " << area << " " << fleet_catch_biomass(area,spp,fleet,year) << endl;
      }}}}       


 //full table of survey length composition
    pmse_predvals << "full time series of survey length comp" << endl;
    pmse_predvals << "survey year spp area catch-at-size" << endl;
     dvar_vector Lpred(1,Nsizebins);
      for (int survey=1; survey<=Nsurveys;survey++) {
      for (int year=1;year<=Nyrs;year++) {
      for (int spp=1;spp<=Nspecies;spp++) {
      Lpred.initialize();
      for (int area=1;area<=Nareas;area++)
       Lpred += N_tot(area,spp,year);
      Lpred = survey_q(survey,spp)*elem_prod(survey_sel(survey,spp),Lpred)/Nstepsyr; //est_survey_size(survey, year, spp, ilen);// 
        pmse_predvals << survey << " " << year << " " << spp << " " << area << " " << Lpred << endl;
      }}}

 //full table of survey length composition
    pmse_predvals << "full time series of fishery length comp" << endl;
    pmse_predvals << "fleet year spp area catch-at-size" << endl;
      for (int fleet=1; fleet<=Nfleets;fleet++) {
      for (int year=1;year<=Nyrs;year++) {
      for (int spp=1;spp<=Nspecies;spp++) {
      for (int area=1;area<=Nareas;area++) {
      Lpred.initialize();
      for (int ilen=1;ilen<=Nsizebins;ilen++)
       Lpred(ilen) = Cfl_tot(area,spp,fleet,year,ilen);
        pmse_predvals << fleet << " " << year << " " << spp << " " << area << " " << Lpred << endl;
      }}}}


//////// full table of predicted survey index
    
    report << "survey_full year spp area pred_survey" << endl;
      //dvariable pred_survey_index2;
      for (int survey=1; survey<=Nsurveys;survey++) {
      for (int year=1;year<=Nyrs;year++) {
      for (int spp=1;spp<=Nspecies;spp++) {
      for (int area=1; area<=Nareas; area++){
        pred_survey_index2 = 0.;
        for (int ilen=1;ilen<=Nsizebins;ilen++) {
            //pred_survey_index2 +=  B_tot(area,spp,year,ilen)*survey_sel(survey,spp,ilen)*survey_q(survey,spp)/Nstepsyr; 
            pred_survey_index2 +=  B_tot(area,spp,year,ilen)*survey_sel(survey,spp,ilen)/Nstepsyr; 
        }
        report << survey << " " << year << " " << spp << " " << area << " " << pred_survey_index2 << endl;
      }}}}       

//////// full table of predicted catch
    
    report << "fleet_full year spp area pred_catch" << endl;
      for (int fleet=1; fleet<=Nfleets;fleet++) {
      for (int year=1;year<=Nyrs;year++) {
      for (int spp=1;spp<=Nspecies;spp++) {
      for (int area=1; area<=Nareas; area++){
        if (indicator_fishery_q(area,fleet,spp) == 1)
         report << fleet << " " << year << " " << spp << " " << area << " " << fleet_catch_biomass(area,spp,fleet,year) << endl;
      }}}}       

