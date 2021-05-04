#############################
#CS model framework function#
### Author Linda Aulin  #####
#############################
# Including 

CS_model <- function(t_half = 5, F_Css_MIC = 1, v_Models = c("Mono A", "Mono B", "Sequential", "3 day cycling",  "1 day cycling",  "Combination"), FIT = 1, CS_A = 1, CS_B = 1, Gmin_A = -1, Gmin_B = -1, HILL_A  = 1, HILL_B = 1, u_1 = 10^-9, u_2 = 10^-9, eB0 = 6, RA0=0, RB0 =0 , RAB0 = 0,  Bmax = 9,  V = 100, n = 10, DT = 1, ST = 24){
  
  require(doParallel)
  require(doRNG)
  require(RxODE)
  require(dplyr)
  require(tidyr)
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol  # need to check for whole numbers later
  
  # argument explanation
  # t_half; half life used for both drugs [h]
  # F_Css_MIC scaling; factor to obtain Css average equal to X time MIC WT
  # FIT; Fitness cost (0-1), 1 no cost, 0 no growth 
  # CS_A ; CS factor use for change in susceptibility towards A
  # CS_B ; CS factor use for change in susceptibility towards B
  # Gmin_A; Maximal effect of drug A
  # Gmin_B; Maximal effect of drug B
  # HILL_A; Shape factor drug A
  # HILL_B; Shape factor B
  # u_1; Mutation rate of S
  # u_2; Mutation rate of RA and RB
  # eB0; Starting bacterial concentration [cfu/mL]
  # RA0; Fraction of bacterial starting in RA
  # RB0; Fraction of bacterial starting in RB
  # RAB0; Fraction of bacterial starting in RAB
  # Bmax; carrying capacity of the system [cfu/mL]
  # V; volume of infection site
  # n; number of simulations , 
  # DT; time step used [h]
  # ST; Duration of simulation [h]
  
  
  start_t <- Sys.time()  #for checking run time
  
  ##############
  ## PK Model ##
  ##############
  
  #Define PK parameters and dosing, assume same PK for both drugs
  # 1 CMT kinetics i.v adm. 
  
  
  Vd = 1                  # [L] volume of distribution. (arbitrary) 
  ke = log(2)/t_half      # [h^-1] elimination rate constant, calculated by half-life
  CL = ke*Vd              # [L/h] clearance 
  Tau = 12                # [h] dosing interval
  MIC = 1
  
  Css = MIC*F_Css_MIC     # [mg/L] Target average steady-state concentration 
  
  Dose = Css*CL*Tau       # [mg] dose giving target Css
  
  
  
  # PK model for antibiotic A and B 
  # ODE

  
  PK_mod<- RxODE({
    
    
    d/dt(A) = -ke*A   
    d/dt(B) = -ke*B
    
  });
  
  
  PK_inits <- c(A = 0, B = 0);
  
  PK_theta <- c(ke = ke)
  
  ### Dosing regimens
  # specified using argument v_Models
  # All uses same Tau and dose, except combination, which uses dose/2 for each AB
  ######
  
  
  
  
  PK_model_list <- list()
  
  if( "Mono A" %in% v_Models){  
    
    ev_PK_A <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(dosing.to = 1, dose=Dose, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST, DT))  %>%
      as.tbl()
    
    PK_A      <- as.data.frame(PK_mod$run(PK_theta, ev_PK_A,     PK_inits)) %>% 
      mutate(Model= "Mono A")
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_A  }
  
  
  if( "Mono B" %in% v_Models){
    
    ev_PK_B <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(dosing.to = 2, dose=Dose, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST, DT))  %>%
      as.tbl()
    
    PK_B      <- as.data.frame(PK_mod$run(PK_theta, ev_PK_B, PK_inits)) %>% 
      mutate(Model = "Mono B") 
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_B }
  
  if( "Sequential" %in% v_Models){
    
    ev_PK_seq<- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(start.time = 0 ,    dosing.to = 1, dose= Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau) %>%
      add.dosing(start.time = ST/2,  dosing.to = 2, dose= Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau) %>%
      add.sampling(seq(0,ST,DT))%>%
      as.tbl()
    
    PK_seq    <- as.data.frame(PK_mod$run(PK_theta, ev_PK_seq,   PK_inits))%>% 
      mutate(Model = "Sequential")
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_seq
    
  }
  
  
  if( "3 day cycling" %in% v_Models){
    
    
    ev_PK_3day <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(start.time = 0 ,  dosing.to = 1, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 72,  dosing.to = 2, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 144, dosing.to = 1, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 216, dosing.to = 2, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.dosing(start.time = 288, dosing.to = 1, dose= Dose, nbr.doses=72/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST,DT))%>%
      as.tbl()
    
    PK_3day   <- as.data.frame(PK_mod$run(PK_theta, ev_PK_3day,  PK_inits))%>% 
      mutate(Model = "3 day cycling")
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_3day
    
  }
  
  
  if( "1 day cycling" %in% v_Models){
    
    
    ev_PK_1day <- eventTable(amount.units="mg", time.units="hours") %>%   #2 dose cycling
      add.dosing(start.time = 0,     dosing.to =  1, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.dosing(start.time = Tau,   dosing.to =  1, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.dosing(start.time = Tau*2, dosing.to =  2, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.dosing(start.time = Tau*3, dosing.to =  2, dose=Dose, nbr.doses=ST/(Tau*2), dosing.interval=Tau*4) %>%
      add.sampling(seq(0,ST,DT))  %>%
      as.tbl()
    
    
    PK_1day   <- as.data.frame(PK_mod$run(PK_theta, ev_PK_1day,  PK_inits))%>% 
      mutate(Model = "1 day cycling")
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_1day
    
  }
  
  
  if( "Combination" %in% v_Models){
    
    
    ev_PK_combo <- eventTable(amount.units="mg", time.units="hours") %>%
      add.dosing(dosing.to = 1, dose=Dose/2, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.dosing(dosing.to = 2, dose=Dose/2, nbr.doses=ST/Tau, dosing.interval=Tau) %>%
      add.sampling(seq(0,ST,DT))%>%
      as.tbl()
    
    
    PK_combo  <- as.data.frame(PK_mod$run(PK_theta, ev_PK_combo, PK_inits))%>% 
      mutate(Model = "Combination")
    
    
    PK_model_list[[length( PK_model_list) + 1]] <- PK_combo
    
  }
  
  
  all_models <- c("Mono A", "Mono B",  "Sequential", "3 day cycling",  "1 day cycling", "Combination")
  dose_reg <- all_models[all_models %in% v_Models]  # selecting only models in dosing regimens in v_Models
  
  # 
  # 
  ######
  
  
  
  #------------------------------------------#
  

  ##############
  ## PD Model ##
  ##############
  
  # 4 state model including S (= WT), RA, RB, RARA (=RAB)
  # PK-PD of AB A and B
  
  #ODE
  CS_mod<- RxODE({
    
    ###PD#####
    
    #Sensitive wild type
    
    d/dt(S) = S*(1-((S+RA+RB+RARB)/10^Bmax))*KG*(1-
                                                   ((Gmax - Gmin_A/KG)*(A/MIC_S)^HILL_A/((A/MIC_S)^HILL_A - (Gmin_A/KG/Gmax)))   -
                                                   ((Gmax - Gmin_B/KG)*(B/MIC_S)^HILL_B/((B/MIC_S)^HILL_B - (Gmin_B/KG/Gmax))))  -
      GR_SRA/V -
      GR_SRB/V -
      
      ke_bac*S
    
    
    # Single resistant mutants to antibiotic A
    d/dt(RA) = RA*(1-((S+RA+RB+RARB)/10^Bmax))*KG*FIT*(1-
                                                         ((Gmax - Gmin_A/(KG*FIT))*(A/MIC_R)^HILL_A/((A/MIC_R)^HILL_A - (Gmin_A/(KG*FIT)/Gmax))) -
                                                         ((Gmax - Gmin_B/(KG*FIT))*(B/(MIC_S*CS_B))^HILL_B/((B/(MIC_S*CS_B))^HILL_B - (Gmin_B/(KG*FIT)/Gmax)))) +
      GR_SRA/V - 
      GR_RARARB/V -
      ke_bac*RA
    
    # Single resistant mutants to antibiotic B
    d/dt(RB) = RB*(1-((S+RA+RB+RARB)/10^Bmax))*KG*FIT*(1 - 
                                                         ((Gmax - Gmin_A/(KG*FIT))*(A/(MIC_S*CS_A))^HILL_A/((A/(MIC_S*CS_A))^HILL_A - (Gmin_A/(KG*FIT)/Gmax))) -
                                                         ((Gmax - Gmin_B/(KG*FIT))*(B/MIC_R)^HILL_B/((B/MIC_R)^HILL_B - (Gmin_B/(KG*FIT)/Gmax)))) +
      GR_SRB/V - 
      GR_RBRARB/V -
      ke_bac*RB
    
    # Double resistant mutants to antibiotic A and B
    d/dt(RARB) = RARB*(1-((S+RA+RB+RARB)/10^Bmax))*KG*FIT*FIT*(1 -
                                                                 ((Gmax - Gmin_A/(KG*FIT*FIT))*(A/MIC_R)^HILL_A/((A/MIC_R)^HILL_A - (Gmin_A/(KG*FIT*FIT)/Gmax)))-
                                                                 ((Gmax - Gmin_B/(KG*FIT*FIT))*(B/MIC_R)^HILL_B/((B/MIC_R)^HILL_B - (Gmin_B/(KG*FIT*FIT)/Gmax)))) +
      GR_RARARB/V + 
      GR_RBRARB/V-
      ke_bac*RARB
    
    
    
    
    ####### PK######
    
    d/dt(A) = -ke*A
    d/dt(B) = -ke*B
    
  });
  
  
  
  S0  <- (10^eB0)*(1-RA0 - RB0 - RAB0)  # starting bacterial density [cfu/ml] of sensitive bacteria (S) = WT)
  U_1 <- u_1  # Mutation rate [mutation/bacteria/hour] of first mutation step
  U_2 <- u_2  #  Mutation rate [mutation/bacteria/hour] of second mutation step (U_1 = U_2)
  
  theta <- c(
    #system specific
    Bmax = Bmax,  # maximum carrying capacity 
    KG = 0.7,   # maximal net growth, --> doubling time 1 h
    FIT = FIT,  # fitness cost per resistance
    ke_bac = 0, # Assume no non-antibiotic bacterial CL
    V     = V,    #Volume of infection (not PK related)
    
    #Drug-system hybrid parameters
    Gmax = 1, 
    Gmin_A = Gmin_A,    # Maximum drug effect , i.e. type of effect
    Gmin_B = Gmin_B,    # Maximum drug effect , i.e. type of effect
    HILL_A = HILL_A,    # shape factor, i.e. type of effect driver 
    HILL_B = HILL_B,    # shape factor, i.e. type of effect driver 
    
    MIC_S = 1,  #MIC if sensitive
    MIC_R = 10, #MIC if resistant
    
    CS_A  = CS_A, #CS factor AB A
    CS_B  = CS_B,  #CS factor AB B
    
    #Drug specific
    ke    = ke)  #PK elimination rate constant))
  
  
  n_obs  <- ST*DT+1   # number of observations ( simulation time(ST) * time step size (DT))
  n_models <- length(PK_model_list) #which models (i.e. dosing regimens) to simulate
  
  df_full_CS <-  foreach( ii = 1:n, .combine = "rbind",              #parallelization by realization
                          .packages = c("RxODE", "dplyr", "tidyr"),
                          
                          .inorder = F,
                          .options.RNG = 123) %dorng%  {            # setting seed
                            
                            
                                          #i_time <- Sys.time()
                           #Creating empty data frame to combined simulation data per realization ii
                            df_model_CS <- data.frame(time  = NA,
                                                      S     = NA,
                                                      RA    = NA,
                                                      RB    = NA,
                                                      RARB  = NA,
                                                      A     = NA,
                                                      B     = NA,
                                                      model = NA)
                            
                            for(i_mod in 1:length(PK_model_list)) {
                              #   #   
                              PK_mod = PK_model_list[[i_mod]]
                              #Creating empty data frame to  be populated with simulated data per dosing regimen i_mode for realization ii
                              
                              df_CS <- data.frame(time = rep(x = NA, times = n_obs),
                                                  S    = rep(x = NA, times = n_obs),
                                                  RA   = rep(x = NA, times = n_obs),
                                                  RB   = rep(x = NA, times = n_obs),
                                                  RARB = rep(x = NA, times = n_obs),
                                                  A    = rep(x = NA, times = n_obs),
                                                  B    = rep(x = NA, times = n_obs))
                              
                              
                              
                              
                              #population t = 0
                              df_CS[1,] <- data.frame(time = 0, S = S0, RA = RA0*S0, RB = RB0*S0, RARB = RAB0*S0,
                                                      A = PK_mod$A[1], B = PK_mod$B[1]);
                              
                              V = V
                              #########
                              # running the simulation one time step at the time (i) needed for the stochastic implementaion for mutation rate
                              
                              for (i in 1:(ST*DT)) {
                                
                               #assigning and checking previous time point for each subpopulation
                                
                                # is S NA?, set to 0
                                if(is.na(df_CS$S[i])){
                                  S_t    <- 0
                                  
                                  
                                }
                                #is the total number of S bacteria larger than 1? , if not set to 0  (protects form running into small number errors)
                                else if(df_CS$S[i]*V>=1) {
                                  S_t  <-  df_CS$S[i]
                                } else{
                                  S_t    <- 0}
                                
                                
                                # is RA NA?, set to 0
                                if(is.na(df_CS$RA[i])){
                                  RA_t    <- 0
                                  
                                }
                                #is the total number of RA bacteria larger than 1? , if not set to 0  (protects form running into small number errors)
                                else if(df_CS$RA[i]*V>=1) {
                                  RA_t    <-  df_CS$RA[i]
                                } else{
                                  RA_t    <- 0}
                                
                                
                                # is RB NA?, set to 0
                                if(is.na(df_CS$RB[i])){
                                  RB_t    <- 0
                                  
                                  #is the total number of RB bacteria larger than 1? , if not set to 0 (protects form running into small number errors)
                                } else if(df_CS$RB[i]*V>=1) {
                                  RB_t  <-  df_CS$RB[i]
                                } else{
                                  RB_t    <- 0}
                                
                                # is RAB NA?, set to 0
                                if(is.na(df_CS$RARB[i])){
                                  RARB_t    <- 0
                                  
                                }
                                
                                #is the total number of RAB bacteria larger than 1?, if not set to 0 (protects form running into small number errors)
                                else if(df_CS$RARB[i]*V>=1) {
                                  RARB_t  <-  df_CS$RARB[i]
                                } else {
                                  RARB_t    <- 0}
                                
                              #assign previous time point for A and B
                                A_t  <- PK_mod$A[i]
                                B_t  <- PK_mod$B[i]
                                
                                
            
                                #-------------------
                                # Mutations (stochastic implementation)
                                # 
                                #same approach for all subpopulations that can mutate (S, RA, RB)
                                
                                # S to RA
                                # check if the number of bacteria can be converted to integer (needed for random sampling)
                                # Workaround to avoid issues with too large numbers (> 10^9), 
  
                                if(is.na(as.integer( S_t*V ))) {
                                  
                                  n_itr <- floor(S_t*V/10^9)
                                  
                                  res_t <- floor(S_t*V - n_itr*10^9)
                                  
                                  # calculate how many S bacteria will mutant into RA over one hour with work around
                                  
                                  GR_SRA  <- sum(rbinom(n_itr, 10^9, U_1), rbinom( 1, res_t, U_1))
                                  
                                  
                                } else if (as.integer( S_t*V )>0) {
                                  # calculate how many S bacteria will mutant into RA over one hour w/o work around
                                  GR_SRA    <- rbinom(1, as.integer( S_t*V ), U_1)
                                  
                                } else {
                                  # no bacteria available for mutation
                                  GR_SRA    <- 0
                                }
                                
                                
                                
                                # S to RB
                                if(is.na(as.integer( S_t*V ))) {
                                  
                                  n_itr <- floor(S_t*V/10^9)
                                  
                                  res_t <- floor(S_t*V - n_itr*10^9)
                                  
                                  GR_SRB  <- sum(rbinom(n_itr, 10^9, U_1), rbinom( 1, res_t, U_1))
                                  
                                  
                                } else if (as.integer( S_t*V )>0) {
                                  
                                  GR_SRB    <- rbinom(1, as.integer( S_t*V ), U_1)
                                  
                                } else{
                                  GR_SRB   <- 0
                                }
                                
                                
                                # RA to RAB
                                if(is.na(as.integer( RA_t*V ))) {
                                  
                                  n_itr <- floor(RA_t*V/10^9)
                                  
                                  res_t <- floor(RA_t*V - n_itr*10^9)
                                  
                                  GR_RARARB  <- sum(rbinom(n_itr, 10^9, U_2), rbinom( 1, res_t, U_2))
                                  
                                  
                                } else if (as.integer( RA_t*V )>0) {
                                  
                                  GR_RARARB    <- rbinom(1, as.integer( RA_t*V ), U_2)
                                  
                                } else{
                                  GR_RARARB   <- 0
                                }
                                
                                
                                
                                #RB to RAB
                                if(is.na(as.integer( RB_t*V ))) {
                                  
                                  n_itr <- floor(RB_t*V/10^9)
                                  
                                  res_t <- floor(RB_t*V - n_itr*10^9)
                                  
                                  GR_RBRARB  <- sum(rbinom(n_itr, 10^9, U_2), rbinom( 1, res_t, U_2))
                                  
                                  
                                } else if (as.integer( RB_t*V )>0) {
                                  
                                  GR_RBRARB    <- rbinom(1, as.integer( RB_t*V ), U_2)
                                  
                                } else{
                                  GR_RBRARB   <- 0
                                }
                                
                                
                              # adding the mutation rates covariates to be used for the specific time step
                                
                                ev <- eventTable(amount.units="mg", time.units="hours") %>%
                                  add.sampling(seq(0,1))  %>%
                                  mutate(GR_SRA = GR_SRA,
                                         GR_SRB = GR_SRB,
                                         GR_RARARB = GR_RARARB,
                                         GR_RBRARB = GR_RBRARB) %>% 
                                  as.tbl()
                                
                                
                                #initiate model with previous time step
                                t_inits <- c(S = S_t, RA = RA_t, RB = RB_t, RARB = RARB_t,
                                             A = A_t, B = B_t);
                                
                                
                                #run the model for specific time step
                                x_CS   <- as.data.frame(CS_mod$run(theta, ev,  t_inits)) %>%
                                  mutate(time = time + i-1)
                                
                                df_CS[i:(i+1),] <- x_CS
                                
                                
                                
                              }
                              # Add name of dosing regimen to data frame
                              df_CS$model <- dose_reg[i_mod] 
                              
                              # bind all simulation for specific realization ii 
                              
                              df_model_CS <- df_model_CS %>%
                                bind_rows(df_CS)
                              
                            }
                            
                            
                            
                            # add realization identifyer
                            df_model_CS$index <- ii
                            
                              
                            #return simulated data for all dosing regimen per realization ii
                            return(df_model_CS)
                            
                            
                            
                          }
  
  
  # output from parallel model runs are bound together into df_full_CS 
  
  mean_dat <-  df_full_CS %>% 
    gather(value = "CFU", key = "Population", - time, -A, -B, -index, -model) %>% #make into long format
    group_by(time, model, index) %>%
    mutate(total = sum(CFU)) %>% #calculate  total CFU per realization, time point and dosing regimen
    ungroup() %>%
    mutate(CFU_ratio = CFU/total) %>%
    #reorder levels of bacterial subpopulations 
    mutate(Order = ifelse(Population == "S", 1,
                          ifelse(Population == "RA", 2,
                                 ifelse(Population == "RB", 3, 4)))) %>%
    mutate(Population = reorder(Population, Order)) %>% 
    ungroup() %>%
    group_by(model, index, time) %>% 
    mutate(CFU_total = sum(CFU, na.rm = T)) %>%  #calculate  total CFU per realization, time point and dosing regimen
    ungroup() %>% 
    group_by(time, model, Population) %>% 
    #calculate metrics relating to CFU per time point and dosing regimen
    mutate(CFU_MEDIAN = median(CFU, na.rm = T), 
           CFU_SD   = sd(CFU, na.rm = T), 
           CFU_95   = quantile(CFU, .95, na.rm = T),
           CFU_05   = quantile(CFU, .05, na.rm = T),
           mean_ratio = mean(CFU_ratio, na.rm = T),
           CFU_T_MEDIAN = median(CFU_total, na.rm = T),
           CFU_T_SD   = sd(CFU_total, na.rm = T), 
           CFU_T_95   = quantile(CFU_total, .95, na.rm = T),
           CFU_T_05   = quantile(CFU_total, .05, na.rm = T),
           #including specific model parameters
           CS_A = CS_A,
           CS_B = CS_B,
           HILL_A = HILL_A,
           HILL_B = HILL_B,
           GMIN_A = Gmin_A,
           GMIN_B = Gmin_B,
           U_1 = u_1,
           U_2 = u_2) %>% 
    ungroup() %>% 
    
    group_by( Population, model, index) %>% 
    #calculating time when subpopulations are over initial bacterial number [not using, remove?]
    mutate(T_10_6 = min(ifelse(CFU >= S0,
                               time, NA), na.rm = T)) %>% 
    ungroup() %>% 
    group_by(Population, model) %>% 
    mutate(End_CFU = ifelse(time == ST, CFU, NA),
           End_T_CFU = ifelse(time == ST, CFU_total, NA),
           day_id  = ifelse(is.wholenumber(time/24), index, 0),
           Day_CFU = ifelse(is.wholenumber(time/24), CFU, NA) ,
           Day_T_CFU = ifelse(is.wholenumber(time/24), CFU_total, NA)) %>% 
    #Calculation the number of realization having resistance (per subpopulation R) and total resistance subpopulations (R_T) per day *_Day_dev or end of treatment *_Dev
    mutate(R_Dev = sum(End_CFU >= 10^eB0, na.rm = T),
           R_T_Dev = sum(End_T_CFU >= 10^eB0, na.rm = T),
           R_Day_Dev = sum(Day_CFU >= 10^eB0, na.rm = T),
           R_T_Day_Dev = sum(Day_T_CFU >= 10^eB0, na.rm = T)) %>%
    ungroup() %>% 
    
    distinct(Population, model, time, day_id, .keep_all = T) %>% 
    select(time, A, B, Population, model, CFU_MEDIAN, CFU_SD, CFU_95, CFU_05, Day_CFU, Day_T_CFU, 
           CFU_T_MEDIAN, CFU_T_SD, CFU_T_95, CFU_T_05, mean_ratio , index, day_id,
           #T_10_6_MEDIAN,T_10_6_SD,T_10_6_95, T_10_6_05, 
           CS_A, CS_B, HILL_A, HILL_B, GMIN_A, GMIN_B, U_1, U_2, End_CFU, End_T_CFU, R_Dev, R_T_Dev, R_Day_Dev, R_T_Day_Dev) %>% 
    filter(!is.na(model))
  
  
  
  
  
  run_t   <-  Sys.time() - start_t 
  print(run_t)
  
  return(mean_dat)
  
  
  
}
