##################################################################################
# PULSE TRAWL PROJECY
# 8 October 2018
# IJmuiden- IMARES
##################################################################################

options(width=240)
library(RDynState5NAsigmaseason52Age5)
library(ggplot2)
library(reshape2)
library(FLCore)
library(plyr)
library(scales)

sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")

# ---------------------
# PLAICE
# ---------------------
pred <- read.csv("~/data/Plaice_Beam.csv")
pred$X <- NULL
names (pred) <- c("cat", "data", "option", "season" )
catchMean(sp1)  <- tapply(pred[,"data"], list(factor(x = pred[,"cat"], levels = unique(pred[,"cat"])), 
                                               factor(x = pred[,"season"], levels = as.character(sort(unique(pred[,"season"])))), 
                                               factor(x = pred[,"option"], levels = unique(pred[,"option"]))), sum)
load("~/data/gample_final.rdata")
catchSigma(sp1) <- array(gample$family$getTheta(), dim=dim(catchMean(sp1)),dimnames=dimnames(catchMean(sp1)))

# ---------------------
# SOLE
# ---------------------
pred <- read.csv("~/data/Sole_Beam.csv")
pred$X <- NULL
names (pred) <- c("cat", "data", "option", "season" )
catchMean(sp2)  <- tapply(pred[,"data"], list(factor(x = pred[,"cat"], levels = unique(pred[,"cat"])), 
                                              factor(x = pred[,"season"], levels = as.character(sort(unique(pred[,"season"])))), 
                                              factor(x = pred[,"option"], levels = unique(pred[,"option"]))), sum)
load("~/data/gamsol_final.rdata")
catchSigma(sp2) <- array(gample$family$getTheta(), dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))

# ---------------------
# COD
# ---------------------
pred <- read.csv("~/data/Cod_Beam.csv")
pred$X <- NULL
names (pred) <- c("cat", "data", "option", "season" )
catchMean(sp3)  <- tapply(pred[,"data"], list(factor(x = pred[,"cat"], levels = unique(pred[,"cat"])), 
                                              factor(x = pred[,"season"], levels = as.character(sort(unique(pred[,"season"])))), 
                                              factor(x = pred[,"option"], levels = unique(pred[,"option"]))), sum)
load("~/data/gamcod_final.rdata")
catchSigma(sp3) <- array(gample$family$getTheta(), dim=dim(catchMean(sp3)),dimnames=dimnames(catchMean(sp3)))

catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))
catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))

# ---------------------
# SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
# ---------------------
price <- read.csv("~/data/Visprijzen.csv")

# sp1Price <- array(c(sp1price + slope1price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
# Plaice 4 marketcat 1 discard with 0 marketvalue
sp1Price <- price[price$Spec == "Plaice" & price$Category != 0 & price$Week != 0,]
# create 5th cat with no marketvalue (discards)
sp1Price <- rbind(sp1Price,data.frame(Week = 1:52, Category = 5, Price = 0, Spec = "Plaice"))
sp1Price <- sp1Price[order(sp1Price$Week),]
sp1Price <- array(sp1Price$Price,dim=c(5,52), dimnames=dimnames(catchMean(sp1))[-3])

# SOLE = 5 cat 52 weeks
sp2Price <- price[price$Spec == "Sole" & price$Category != 0 & price$Week != 0,]
sp2Price <- sp2Price[order(sp2Price$Week),]
sp2Price <- array(sp2Price$Price,dim=c(5,52), dimnames=dimnames(catchMean(sp2))[-3])

sp3Price <- sp4Price <- sp5Price <- array(c(0), dim=c(5,52), dimnames=dimnames(catchMean(sp2))[-3])

#-------------------------------------------------------------------------------------
# EFFORT                                                       
#-------------------------------------------------------------------------------------
#URK
#effort  <- array(c(14,16,14,14,15,17,15,15,18),dim=9,2)       # kwart van een dag, so total effort divided by 4!
effort <- array(c(14,16,14,14,15,17,15,15,18), dim=c(9,52), dimnames=list(option=unique(pred[,"option"]),season=as.character(sort(unique(pred[,"season"])))))

#-------------------------------------------------------------------------------------
# CALCULATE INCS                                               
#-------------------------------------------------------------------------------------
#sp1inc <- qnbinom(0.90, mu= max(ple.input@CatchMean), size= ple.input@CatchSigma[1,1,1,1,1,1]) / (30 - 1)
sp1inc <- qnbinom(0.90, mu= max(catchMean(sp1)), size= catchSigma(sp1)[1,1,1]) / (30 - 1)
sp2inc <- qnbinom(0.90, mu= max(catchMean(sp2)), size= catchSigma(sp2)[1,1,1]) / (30 - 1)
sp1inc
sp2inc


#-------------------------------------------------------------------------------------
# Make contol and execute calculations
#-------------------------------------------------------------------------------------
# control <- DynState.control(Increments=30, PlaiceUplimit=1600000,PlaiceDiscardSteps=1, SoleDiscardSteps=0,CodDiscardSteps=0, 
#                             EffortUplimit=NA,Handling= 0.24,CrewShare = 0.33,GearCost=87,VarCost= 0.05,EffortPrice=600,
#                             FinePlaice=320,ChoiceDist=1,SimNumber=1000, NumThreads=20)
control     <- DynState.control(spp1LndQuota= 1600000,  spp2LndQuota=2000000, spp1LndQuotaFine= 320, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 600, landingCosts= 0.24,gearMaintenance= 87, addNoFishing= TRUE, increments= 30, spp1DiscardSteps= 0, spp2DiscardSteps= 0, sigma= 1, simNumber= 1000 , numThreads= 20, verbose=1)
 

z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)














source("functions.R")


ages          <- 1:6
season        <- 1:12
areas         <- c("a", "b")
stab.model    <- 10
NUMRUNS       <- 30
MPstart       <- 25

SIMNUMBER     <- 8000 #pos 8000 WORKS
SPP1DSCSTEPS  <- 0
SPP2DSCSTEPS  <- 0
endy          <- stab.model + NUMRUNS

Linf          <- 50 #50
K             <- 0.6 #0.4
alpha         <- 0.0002#0.00005
beta          <- 2.4
sages         <- array(seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)), dim=c(length(season),1,length(ages),1), dimnames=list(season=as.character(season),   year="all", cat=ages, option ="all"))
lens          <- Linf*(1-exp(-K*(sages)))
wts           <- alpha * lens ^ beta
wts           <- aperm(wts, c(3,2,1,4))

q             <- 0.000025
natmortality  <- 0.0001

migconstant   <- 0.025 # 
sp1price      <- sp2price      <- 5000
slope1price <- 0.50*sp1price
slope2price <- 0.50*sp2price#1000 # 0.50*150

FuelPrice   <- 2000

# scenario I: discarding is not allowed, YPR based in C (C=L)
# scenario II: discarding is allowed, YPR based in L, hr wanted based in catches
# scenario III: discarding ocurred but not perceived, YPR based in L, hr wanted based in landings

recs1          <- c(500000,0) 
mig1     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig1[,,,"a","a"] <- -migconstant
mig1[,,,"b","b"] <- -migconstant
mig1[,,,"a","b"] <- migconstant
mig1[,,,"b","a"] <- migconstant
aperm( mig1,c(1,3,2,4,5))

recs2          <- c(0,500000) 
mig2     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig2[,,,"a","a"] <- -migconstant
mig2[,,,"b","b"] <- -migconstant
mig2[,,,"a","b"] <- migconstant
mig2[,,,"b","a"] <- migconstant
aperm( mig2,c(1,3,2,4,5))

effort <- array(c(1), dim=c(length(areas), length(season)), dimnames=list(option =areas,season=as.character(season)))

pop1  <- pop2         <-array(0, dim=c(length(ages),endy+1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm1       <- catches.n.dsvm2       <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.n.dsvm1      <- landings.n.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm1      <- catches.wt.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.wt.dsvm1     <- landings.wt.dsvm2     <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm.tot1  <- catches.wt.dsvm.tot2  <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
landings.wt.dsvm.tot1 <- landings.wt.dsvm.tot2 <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
quota1                <- quota2                <- array(5000, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

pos_catches1 <- pos_catches2 <- pop1

hr1wanted <- hr2wanted <- array(NA, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs1, migration=mig1)
pop2 <- population_dynamics(pop=pop2, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs2, migration=mig2)


#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp1)  <- catchMean(sp2) <- array(NA, dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))

catchMean(sp3)  <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp3) <- catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

#SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
sp1Price <- array(c(sp1price + slope1price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp2Price <- array(c(sp2price + slope2price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp3Price <- sp4Price <- sp5Price <- array(c(0), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
#---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         
control     <- DynState.control(spp1LndQuota= quota1[,1,,],  spp2LndQuota=quota2[,1,,], spp1LndQuotaFine= 3e6, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = FuelPrice, landingCosts= 0,gearMaintenance= 0, 
                                addNoFishing= TRUE, increments= 25, spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER * 0.5 , numThreads= 68, verbose=1)

#this is where our loop starts, after we set up stable population
for(yy in (stab.model):(stab.model+NUMRUNS)){
  
  print("====== year yy ========")
  print(yy)

  #Calculate catch inputs for DSVM
  #for age 1, assume that catchmean over seasons and areas will be equal to the year before
  catchMean(sp1)[1,,] <- pop1[1,as.character(yy-1),,] *q*wts[1,1,,1]
  #for older ages we start from cohort i nthe previous year
  #we need to guess how much population will decline while moving into this year, and over seasons, for this wee need to look at end of two years ago and compare to last year 
  #this comes from a sweep with pop1[1:(max(ages)-1),as.character(yy-2),max(season),] and pop1[-1,as.character(yy-1),,] 
  cfracs1 <- sweep( pop1[-1,as.character(yy-1),,] ,c(1,3),  pop1[1:(max(ages)-1),as.character(yy-2),max(season),] , FUN="/")
  cfracs2 <- sweep( pop1[-1,as.character(yy-2),,] ,c(1,3),  pop1[1:(max(ages)-1),as.character(yy-3),max(season),] , FUN="/")
  cfracs <- (cfracs1 + cfracs2)/2
  #once we have the decline as a functrion of time (in cfracs) we will mulitply by cohorts in last season, and multiply by q and wts for appropriate ages. Those go in catchmean for older ages
  #it is a two step approach because we need two sweeps
  catchMean(sp1)[-1,,] <- sweep(cfracs,c(1,3),pop1[1:(max(ages)-1),as.character(yy-1),max(season),] ,FUN="*") *q
  catchMean(sp1)[-1,,] <- sweep(catchMean(sp1)[-1,,],c(1,2), wts[-1,1,,1],FUN="*")
  
  catchMean(sp2)[1,,] <- pop2[1,as.character(yy-1),,] *q*wts[1,1,,1]
  cfracs1 <- sweep( pop2[-1,as.character(yy-1),,] ,c(1,3),  pop2[1:(max(ages)-1),as.character(yy-2),max(season),] , FUN="/")
  cfracs2 <- sweep( pop2[-1,as.character(yy-2),,] ,c(1,3),  pop2[1:(max(ages)-1),as.character(yy-3),max(season),] , FUN="/")
  cfracs <- (cfracs1 + cfracs2)/2
  catchMean(sp2)[-1,,] <- sweep(cfracs,c(1,3),pop2[1:(max(ages)-1),as.character(yy-1),max(season),] ,FUN="*") *q
  catchMean(sp2)[-1,,] <- sweep(catchMean(sp2)[-1,,],c(1,2), wts[-1,1,,1],FUN="*")
  
    
  # ---No way of estimating sigma, therefore we assume that is 8% of the CPUE (note slight repetion in code for dims and dimnames of 0 catch arrays for spec 3,4,5)                                                                  
  catchSigma(sp1) <- catchMean(sp1) *0.08
  catchSigma(sp2) <- catchMean(sp2) *0.08
  
  # if we are in MP period, then set quota based on last year
  if (yy > MPstart){
    control@spp1LndQuota <-  quota1[,yy,,]
    #control@spp2LndQuota <-  quota2[,yy,,]
  }
  
  # because we wanted gradual introduction of vessels, we only put all vessels after stabmodel. First line is not strictly necessary becuase already defined when constructing control 
  if (yy == (stab.model)  ) control@simNumber <-  as.integer(SIMNUMBER*0.5)
  if (yy == (stab.model+1)) control@simNumber <-  as.integer(SIMNUMBER*0.75)
  if (yy >  (stab.model+1)) control@simNumber <-  as.integer(SIMNUMBER)
          
  
  # if we are in MPLO period, then set discardsteps =0 
  #if (yy > MPstartLO){
   # control@spp1DiscardSteps <-  0
  #  control@spp2DiscardSteps <-  0
  #}

  #run DSVM (wiht quota constraining if in MP time)   
  z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
  
  #extract DSVM results
  dsvm_res <-  extract_dsvm_res (z, control, ages, season)
  #Net revenuw from DSVM
  economics_res             <- as.data.frame(as.matrix(netRev(z)))
  names(economics_res )     <- "NetRev"
  economics_res$Grossrev    <- as.data.frame(as.matrix(grossRev(z)))$V1 
  economics_res$Fuelcosts   <- apply(effort(sim(z)) * control(z)@fuelUse * control(z)@fuelPrice, 2, sum)
  economics_res$Annualfine  <- as.data.frame(as.matrix(annualFine(z)))$V1

  if (yy == stab.model){
    dsvm_res_allyrs       <- cbind("year"= yy,dsvm_res)
    economics_res_allyrs  <- cbind("year"= yy,economics_res)
  } else {
    dsvm_res_allyrs <- rbind(dsvm_res_allyrs, (cbind("year"= yy,dsvm_res)))
    economics_res_allyrs <- rbind(economics_res_allyrs, (cbind("year"= yy,economics_res)))
  }
  
  #get catches in wts from DSVM 
  catches.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="catch.wt")
  catches.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="catch.wt") 
  
  landings.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="landings.wt")
  landings.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="landings.wt") 
  
  rm(z,dsvm_res)
  
  #calculate total catches (by summing over seasons and ages)
  catches.wt.dsvm.tot1[] <- apply(catches.wt.dsvm1,c(2),"sum")
  catches.wt.dsvm.tot2[] <- apply(catches.wt.dsvm2,c(2),"sum")
  
  landings.wt.dsvm.tot1[] <- apply(landings.wt.dsvm1,c(2),"sum")
  landings.wt.dsvm.tot2[] <- apply(landings.wt.dsvm2,c(2),"sum")
  
  #calculate numbers caught from weight caught 
  for(jj in areas){
    catches.n.dsvm1 [,yy,,as.character(jj)]  <- catches.wt.dsvm1[,yy,,as.character(jj)] / wts[,1,,1]
    catches.n.dsvm2 [,yy,,as.character(jj)]  <- catches.wt.dsvm2[,yy,,as.character(jj)] / wts[,1,,1]
    
    landings.n.dsvm1 [,yy,,as.character(jj)] <- landings.wt.dsvm1[,yy,,as.character(jj)] / wts[,1,,1]
    landings.n.dsvm2 [,yy,,as.character(jj)] <- landings.wt.dsvm2[,yy,,as.character(jj)] / wts[,1,,1]
  }

  # calculatae what happens to population based on catches
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,yy,,,drop=F], recruitment=recs1, migration=mig1)
  pop2 <- population_dynamics(pop=pop2, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm2[,yy,,,drop=F], recruitment=recs2, migration=mig2)
  
  #------------------------------------------
  #MANAGEMENT PROCEDURE
  #------------------------------------------
  
 
  #calculate the observed harvest ratios (per age) 
  hr1 <- (apply(catches.n.dsvm1,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm1,1:3,sum)+1e-20) +  apply(pop1[,1:endy,,],1:3,sum))
  hr2 <- (apply(catches.n.dsvm2,1:3,sum)+1e-20)/    (apply(catches.n.dsvm2,1:3,sum) +  apply(pop2[,1:endy,,],1:3,sum))
  
  #calculate landing ratios
  landings.ratio1 <- (apply(landings.n.dsvm1,1:3,sum)+1e-20)/ (apply(catches.n.dsvm1,1:3,sum)+1e-20)
  landings.ratio2 <- (apply(landings.n.dsvm2,1:3,sum)+1e-20)/ (apply(catches.n.dsvm2,1:3,sum)+1e-20)

  
  # if (yy == MPstartLO){
  #   landings.ratio1[,yy,]<- 1
  #   landings.ratio2[,yy,]<- 1
  # }
  
  #what sequence to use in YPR curve so that we span at least hr= 0.001 to 0.5?
  maxseq1 <- 0.5/mean(hr1[,yy,])
  maxseq2 <- 0.5/mean(hr2[,yy,])
  
  # get yield curve results given observed harvest ratios
  yc1 <- yield_curve(hr=hr1[,yy,], landings.ratio1[,yy,], wts, natmortality, R=recs1, verbose=F, sequence = seq(0.001, maxseq1, length.out=2000))
  yc2 <- yield_curve(hr=hr2[,yy,], landings.ratio2[,yy,], wts, natmortality, R=recs2, verbose=F, sequence = seq(0.001, maxseq2, length.out=2000))
  
  #store perceived HRmax    
  hr1wanted[,yy,,] <- (yc1[yc1$landings==max(yc1$landings),]$hr)[1]
  hr2wanted[,yy,,] <- (yc2[yc2$landings==max(yc2$landings),]$hr)[1]
  
  #calculate individual quota from pop, hr and hrmax 
  if (yy > (MPstart-1) & yy < endy){ #if (yy > (MPstart-1))
    prel.quota1 <- sum((hr1wanted[,yy,,]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER 
      #sum(sweep((hr1wanted[1]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum) ,1,wts[,1,,1],"*"))/SIMNUMBER
    prel.quota2 <- sum((hr2wanted[,yy,,]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER
      #sum(sweep((hr2wanted[1]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER
    
    # We constrained +- 15% TAC change, until the LO implementation, where we allow the HR wanted just for the transition
    # In the transition we assume landing ratio equal 1
    #   if (yy == MPstart){ #MPstartLO){ # landing.ratio are equal to 1 for all years before estimating them
      quota1[,yy+1,,] <- prel.quota1
      #quota2[,yy+1,,] <- prel.quota2
      
    #}else{
    #  tac.constrained1 <- c(quota1[,yy,,]*0.85, quota1[,yy,,]*1.15)
    #  tac.constrained2 <- c(quota2[,yy,,]*0.85, quota2[,yy,,]*1.15)
    #  quota1[,yy+1,,] <- max(min(prel.quota1, tac.constrained1[2]), tac.constrained1[1]) 
    #  quota2[,yy+1,,] <- max(min(prel.quota2, tac.constrained2[2]), tac.constrained2[1]) 
    #}
  }
}
ages          <- 1:6
season        <- 1:12
areas         <- c("a", "b")
stab.model    <- 10
NUMRUNS       <- 30
MPstart       <- 25

SIMNUMBER     <- 8000 #pos 8000 WORKS
SPP1DSCSTEPS  <- 0
SPP2DSCSTEPS  <- 0
endy          <- stab.model + NUMRUNS

Linf          <- 50 #50
K             <- 0.6 #0.4
alpha         <- 0.0002#0.00005
beta          <- 2.4
sages         <- array(seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)), dim=c(length(season),1,length(ages),1), dimnames=list(season=as.character(season),   year="all", cat=ages, option ="all"))
lens          <- Linf*(1-exp(-K*(sages)))
wts           <- alpha * lens ^ beta
wts           <- aperm(wts, c(3,2,1,4))

q             <- 0.000025
natmortality  <- 0.0001

migconstant   <- 0.025 # 
sp1price      <- sp2price      <- 5000
slope1price <- 0.50*sp1price
slope2price <- 0.50*sp2price#1000 # 0.50*150

FuelPrice   <- 2000

# scenario I: discarding is not allowed, YPR based in C (C=L)
# scenario II: discarding is allowed, YPR based in L, hr wanted based in catches
# scenario III: discarding ocurred but not perceived, YPR based in L, hr wanted based in landings

recs1          <- c(500000,0) 
mig1     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig1[,,,"a","a"] <- -migconstant
mig1[,,,"b","b"] <- -migconstant
mig1[,,,"a","b"] <- migconstant
mig1[,,,"b","a"] <- migconstant
aperm( mig1,c(1,3,2,4,5))

recs2          <- c(0,500000) 
mig2     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig2[,,,"a","a"] <- -migconstant
mig2[,,,"b","b"] <- -migconstant
mig2[,,,"a","b"] <- migconstant
mig2[,,,"b","a"] <- migconstant
aperm( mig2,c(1,3,2,4,5))

effort <- array(c(1), dim=c(length(areas), length(season)), dimnames=list(option =areas,season=as.character(season)))

pop1  <- pop2         <-array(0, dim=c(length(ages),endy+1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm1       <- catches.n.dsvm2       <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.n.dsvm1      <- landings.n.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm1      <- catches.wt.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.wt.dsvm1     <- landings.wt.dsvm2     <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm.tot1  <- catches.wt.dsvm.tot2  <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
landings.wt.dsvm.tot1 <- landings.wt.dsvm.tot2 <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
quota1                <- quota2                <- array(5000, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

pos_catches1 <- pos_catches2 <- pop1

hr1wanted <- hr2wanted <- array(NA, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs1, migration=mig1)
pop2 <- population_dynamics(pop=pop2, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs2, migration=mig2)


#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp1)  <- catchMean(sp2) <- array(NA, dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))

catchMean(sp3)  <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp3) <- catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

#SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
sp1Price <- array(c(sp1price + slope1price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp2Price <- array(c(sp2price + slope2price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp3Price <- sp4Price <- sp5Price <- array(c(0), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
#---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         
control     <- DynState.control(spp1LndQuota= quota1[,1,,],  spp2LndQuota=quota2[,1,,], spp1LndQuotaFine= 3e6, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = FuelPrice, landingCosts= 0,gearMaintenance= 0, 
                                addNoFishing= TRUE, increments= 25, spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER * 0.5 , numThreads= 68, verbose=1)

#this is where our loop starts, after we set up stable population
for(yy in (stab.model):(stab.model+NUMRUNS)){
  
  print("====== year yy ========")
  print(yy)

  #Calculate catch inputs for DSVM
  #for age 1, assume that catchmean over seasons and areas will be equal to the year before
  catchMean(sp1)[1,,] <- pop1[1,as.character(yy-1),,] *q*wts[1,1,,1]
  #for older ages we start from cohort i nthe previous year
  #we need to guess how much population will decline while moving into this year, and over seasons, for this wee need to look at end of two years ago and compare to last year 
  #this comes from a sweep with pop1[1:(max(ages)-1),as.character(yy-2),max(season),] and pop1[-1,as.character(yy-1),,] 
  cfracs1 <- sweep( pop1[-1,as.character(yy-1),,] ,c(1,3),  pop1[1:(max(ages)-1),as.character(yy-2),max(season),] , FUN="/")
  cfracs2 <- sweep( pop1[-1,as.character(yy-2),,] ,c(1,3),  pop1[1:(max(ages)-1),as.character(yy-3),max(season),] , FUN="/")
  cfracs <- (cfracs1 + cfracs2)/2
  #once we have the decline as a functrion of time (in cfracs) we will mulitply by cohorts in last season, and multiply by q and wts for appropriate ages. Those go in catchmean for older ages
  #it is a two step approach because we need two sweeps
  catchMean(sp1)[-1,,] <- sweep(cfracs,c(1,3),pop1[1:(max(ages)-1),as.character(yy-1),max(season),] ,FUN="*") *q
  catchMean(sp1)[-1,,] <- sweep(catchMean(sp1)[-1,,],c(1,2), wts[-1,1,,1],FUN="*")
  
  catchMean(sp2)[1,,] <- pop2[1,as.character(yy-1),,] *q*wts[1,1,,1]
  cfracs1 <- sweep( pop2[-1,as.character(yy-1),,] ,c(1,3),  pop2[1:(max(ages)-1),as.character(yy-2),max(season),] , FUN="/")
  cfracs2 <- sweep( pop2[-1,as.character(yy-2),,] ,c(1,3),  pop2[1:(max(ages)-1),as.character(yy-3),max(season),] , FUN="/")
  cfracs <- (cfracs1 + cfracs2)/2
  catchMean(sp2)[-1,,] <- sweep(cfracs,c(1,3),pop2[1:(max(ages)-1),as.character(yy-1),max(season),] ,FUN="*") *q
  catchMean(sp2)[-1,,] <- sweep(catchMean(sp2)[-1,,],c(1,2), wts[-1,1,,1],FUN="*")
  
    
  # ---No way of estimating sigma, therefore we assume that is 8% of the CPUE (note slight repetion in code for dims and dimnames of 0 catch arrays for spec 3,4,5)                                                                  
  catchSigma(sp1) <- catchMean(sp1) *0.08
  catchSigma(sp2) <- catchMean(sp2) *0.08
  
  # if we are in MP period, then set quota based on last year
  if (yy > MPstart){
    control@spp1LndQuota <-  quota1[,yy,,]
    #control@spp2LndQuota <-  quota2[,yy,,]
  }
  
  # because we wanted gradual introduction of vessels, we only put all vessels after stabmodel. First line is not strictly necessary becuase already defined when constructing control 
  if (yy == (stab.model)  ) control@simNumber <-  as.integer(SIMNUMBER*0.5)
  if (yy == (stab.model+1)) control@simNumber <-  as.integer(SIMNUMBER*0.75)
  if (yy >  (stab.model+1)) control@simNumber <-  as.integer(SIMNUMBER)
          
  
  # if we are in MPLO period, then set discardsteps =0 
  #if (yy > MPstartLO){
   # control@spp1DiscardSteps <-  0
  #  control@spp2DiscardSteps <-  0
  #}

  #run DSVM (wiht quota constraining if in MP time)   
  z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
  
  #extract DSVM results
  dsvm_res <-  extract_dsvm_res (z, control, ages, season)
  #Net revenuw from DSVM
  economics_res             <- as.data.frame(as.matrix(netRev(z)))
  names(economics_res )     <- "NetRev"
  economics_res$Grossrev    <- as.data.frame(as.matrix(grossRev(z)))$V1 
  economics_res$Fuelcosts   <- apply(effort(sim(z)) * control(z)@fuelUse * control(z)@fuelPrice, 2, sum)
  economics_res$Annualfine  <- as.data.frame(as.matrix(annualFine(z)))$V1

  if (yy == stab.model){
    dsvm_res_allyrs       <- cbind("year"= yy,dsvm_res)
    economics_res_allyrs  <- cbind("year"= yy,economics_res)
  } else {
    dsvm_res_allyrs <- rbind(dsvm_res_allyrs, (cbind("year"= yy,dsvm_res)))
    economics_res_allyrs <- rbind(economics_res_allyrs, (cbind("year"= yy,economics_res)))
  }
  
  #get catches in wts from DSVM 
  catches.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="catch.wt")
  catches.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="catch.wt") 
  
  landings.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="landings.wt")
  landings.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="landings.wt") 
  
  rm(z,dsvm_res)
  
  #calculate total catches (by summing over seasons and ages)
  catches.wt.dsvm.tot1[] <- apply(catches.wt.dsvm1,c(2),"sum")
  catches.wt.dsvm.tot2[] <- apply(catches.wt.dsvm2,c(2),"sum")
  
  landings.wt.dsvm.tot1[] <- apply(landings.wt.dsvm1,c(2),"sum")
  landings.wt.dsvm.tot2[] <- apply(landings.wt.dsvm2,c(2),"sum")
  
  #calculate numbers caught from weight caught 
  for(jj in areas){
    catches.n.dsvm1 [,yy,,as.character(jj)]  <- catches.wt.dsvm1[,yy,,as.character(jj)] / wts[,1,,1]
    catches.n.dsvm2 [,yy,,as.character(jj)]  <- catches.wt.dsvm2[,yy,,as.character(jj)] / wts[,1,,1]
    
    landings.n.dsvm1 [,yy,,as.character(jj)] <- landings.wt.dsvm1[,yy,,as.character(jj)] / wts[,1,,1]
    landings.n.dsvm2 [,yy,,as.character(jj)] <- landings.wt.dsvm2[,yy,,as.character(jj)] / wts[,1,,1]
  }

  # calculatae what happens to population based on catches
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,yy,,,drop=F], recruitment=recs1, migration=mig1)
  pop2 <- population_dynamics(pop=pop2, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm2[,yy,,,drop=F], recruitment=recs2, migration=mig2)
  
  #------------------------------------------
  #MANAGEMENT PROCEDURE
  #------------------------------------------
  
 
  #calculate the observed harvest ratios (per age) 
  hr1 <- (apply(catches.n.dsvm1,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm1,1:3,sum)+1e-20) +  apply(pop1[,1:endy,,],1:3,sum))
  hr2 <- (apply(catches.n.dsvm2,1:3,sum)+1e-20)/    (apply(catches.n.dsvm2,1:3,sum) +  apply(pop2[,1:endy,,],1:3,sum))
  
  #calculate landing ratios
  landings.ratio1 <- (apply(landings.n.dsvm1,1:3,sum)+1e-20)/ (apply(catches.n.dsvm1,1:3,sum)+1e-20)
  landings.ratio2 <- (apply(landings.n.dsvm2,1:3,sum)+1e-20)/ (apply(catches.n.dsvm2,1:3,sum)+1e-20)

  
  # if (yy == MPstartLO){
  #   landings.ratio1[,yy,]<- 1
  #   landings.ratio2[,yy,]<- 1
  # }
  
  #what sequence to use in YPR curve so that we span at least hr= 0.001 to 0.5?
  maxseq1 <- 0.5/mean(hr1[,yy,])
  maxseq2 <- 0.5/mean(hr2[,yy,])
  
  # get yield curve results given observed harvest ratios
  yc1 <- yield_curve(hr=hr1[,yy,], landings.ratio1[,yy,], wts, natmortality, R=recs1, verbose=F, sequence = seq(0.001, maxseq1, length.out=2000))
  yc2 <- yield_curve(hr=hr2[,yy,], landings.ratio2[,yy,], wts, natmortality, R=recs2, verbose=F, sequence = seq(0.001, maxseq2, length.out=2000))
  
  #store perceived HRmax    
  hr1wanted[,yy,,] <- (yc1[yc1$landings==max(yc1$landings),]$hr)[1]
  hr2wanted[,yy,,] <- (yc2[yc2$landings==max(yc2$landings),]$hr)[1]
  
  #calculate individual quota from pop, hr and hrmax 
  if (yy > (MPstart-1) & yy < endy){ #if (yy > (MPstart-1))
    prel.quota1 <- sum((hr1wanted[,yy,,]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER 
      #sum(sweep((hr1wanted[1]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum) ,1,wts[,1,,1],"*"))/SIMNUMBER
    prel.quota2 <- sum((hr2wanted[,yy,,]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER
      #sum(sweep((hr2wanted[1]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER
    
    # We constrained +- 15% TAC change, until the LO implementation, where we allow the HR wanted just for the transition
    # In the transition we assume landing ratio equal 1
    #   if (yy == MPstart){ #MPstartLO){ # landing.ratio are equal to 1 for all years before estimating them
      quota1[,yy+1,,] <- prel.quota1
      #quota2[,yy+1,,] <- prel.quota2
      
    #}else{
    #  tac.constrained1 <- c(quota1[,yy,,]*0.85, quota1[,yy,,]*1.15)
    #  tac.constrained2 <- c(quota2[,yy,,]*0.85, quota2[,yy,,]*1.15)
    #  quota1[,yy+1,,] <- max(min(prel.quota1, tac.constrained1[2]), tac.constrained1[1]) 
    #  quota2[,yy+1,,] <- max(min(prel.quota2, tac.constrained2[2]), tac.constrained2[1]) 
    #}
  }
}


##############################################################################
# POPULATION DYNAMICS
##############################################################################

population_dynamics <- function(pop, startyear, endyear, season, natmortality, catches, recruitment, migration){
  #pop[age,year, season, area]                                                                                                        
  MigToArea <- array(0, dim=dim(migration), dimnames= dimnames(migration))
  for (y in (startyear:endyear)){     # need to think this, but maybe is better change last year (y in (endy:(endy +10)))           
    for (ss in (1:length(season))){
      # move time ---------------------                                                             
      if (ss ==1){
        for(age in 2:(dim(pop)[1])){
          pop[age,as.character(y),as.character(ss),] <- pop[age-1,as.character(y-1),as.character(length(season)),];
        }
      } else {
        pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss-1),];
      }
      
      # birth/recruitment ---------------------                                                     
      if (ss ==1){
        pop[1,as.character(y),as.character(ss),] <- recruitment;
      }
      # natural mortality  ---------------------                                                             
      pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss),]*(1-natmortality)
      pop[pop < 1e-1 ] <- 1e-1
      
      # remove catches (dims of catches here is ages,season, area, just like in main )
      pop[,as.character(y),as.character(ss),] <- pop[,as.character(y),as.character(ss),] - catches[,,as.character(ss),]
      pop[pop < 1e-1 ] <- 1e-1 
      
      #migration
      for (age in (1:dim(pop)[1])){
        MigToArea[age,1,as.character(ss),,] <- 0
        for (toarea in (dimnames(pop)[4][[1]])){
          for (fromarea in (dimnames(pop)[4][[1]])){
            MigToArea[age,1,as.character(ss),toarea, fromarea] <-  MigToArea[age,1,as.character(ss),toarea,fromarea] + ( pop[age,as.character(y),as.character(ss),fromarea] * migration[age,1,as.character(ss),fromarea, toarea])
          }
        }
        for (toarea in (dimnames(pop)[4][[1]])){
          for (fromarea in (dimnames(pop)[4][[1]])){
            pop[age,as.character(y),as.character(ss),toarea]  <-  pop[age,as.character(y),as.character(ss),toarea] +  MigToArea[age,1,as.character(ss),toarea, fromarea]
          }
        }
      }
    }
  }
  return(pop)
}



##############################################################################
# YIELD CURVE
##############################################################################

yield_curve <- function(hr,lratio, wts, natmortality, R=1, sequence = seq(0.001,5,0.005), verbose=F ){
  
  # note that definition of hr is not completely correct (should be sum over seasons, and mean over ages), but as long as consistently incorect in code it should not matter
  res <- data.frame("hr"=mean(hr) *sequence,"catch"=NA, "landings"=NA) 
  iii <- 1
  sumR <- sum(R)
  if (verbose == T){ 
    print("total Recruitment")
    print(R)
  }
  for (ii in sequence){
    respop <- yld <-  matrix(0,nrow=length(ages), ncol=length(season), dimnames=list("cat"=ages,"season"=season))  
    respop[1,1] <- sumR
    for(aa in ages){
      if (aa==1){
        # respop[aa,1] <- respop[aa,1] * (1-natmortality*1)
        respop[aa,1] <- respop[aa,1] * (1-(hr[aa,1]*ii))
        yld[aa,1]    <- sumR   - respop[aa,1] 
        respop[aa,1] <- respop[aa,1] * (1-natmortality*1)
      }
      if (aa > 1){
        respop[aa,1] <- respop[aa-1,max(season)] * (1-natmortality*0)
        respop[aa,1] <- respop[aa,1]             * (1-(hr[aa,1] *ii))
        yld[aa,1]    <- respop[aa-1,max(season)] - respop[aa,1]
        respop[aa,1] <- respop[aa,1]             * (1-natmortality*1)
      }
      for (ss in 2:(max(season))){
        respop[aa,ss] <- respop[aa,ss-1] * (1-natmortality*0)
        respop[aa,ss] <- respop[aa,ss]   * (1-(hr[aa,ss]*ii))
        yld[aa,ss]    <- respop[aa,ss-1] - respop[aa,ss]
        respop[aa,ss] <- respop[aa,ss]   * (1-natmortality*1)
      }
    }
    res[iii, ]$catch    <- sum(yld*wts[,1,,1])
    res[iii, ]$landings <- sum(yld*lratio*wts[,1,,1])
    iii <- iii + 1
    if (verbose == T){ 
      print("yields (in numbers)")
      print(yld)
      print(" ")
      print("population (in numbers)")
      print(respop)
    }
  }
  return(res)
}

SIGMAS        <- seq(30000,120000,5000)   #SIGMA         <- 100 #comes from 2// chanheg from 300 to 200

for (SIGMA in SIGMAS){

ages          <- 1:6
season        <- 1:12
areas         <- c("a", "b")
stab.model    <- 10
NUMRUNS       <- 30
MPstart       <- 25

SIMNUMBER     <- 8000 #pos 8000 WORKS
SPP1DSCSTEPS  <- 0
SPP2DSCSTEPS  <- 0
endy          <- stab.model + NUMRUNS

Linf          <- 50 #50
K             <- 0.6 #0.4
alpha         <- 0.0002#0.00005
beta          <- 2.4
sages         <- array(seq(min(ages)+((1/max(season))/2), max(ages+1),1/max(season)), dim=c(length(season),1,length(ages),1), dimnames=list(season=as.character(season),   year="all", cat=ages, option ="all"))
lens          <- Linf*(1-exp(-K*(sages)))
wts           <- alpha * lens ^ beta
wts           <- aperm(wts, c(3,2,1,4))

q             <- 0.000025
natmortality  <- 0.0001

migconstant   <- 0.025 # 
sp1price      <- sp2price      <- 5000
slope1price <- 0.50*sp1price
slope2price <- 0.50*sp2price#1000 # 0.50*150

FuelPrice   <- 2000

# scenario I: discarding is not allowed, YPR based in C (C=L)
# scenario II: discarding is allowed, YPR based in L, hr wanted based in catches
# scenario III: discarding ocurred but not perceived, YPR based in L, hr wanted based in landings

recs1          <- c(500000,0) 
mig1     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig1[,,,"a","a"] <- -migconstant
mig1[,,,"b","b"] <- -migconstant
mig1[,,,"a","b"] <- migconstant
mig1[,,,"b","a"] <- migconstant
aperm( mig1,c(1,3,2,4,5))

recs2          <- c(0,500000) 
mig2     <- array(0, dim=c(length(ages),1,length(season),length(areas), length(areas)), dimnames=list(cat=ages,year="all",season=as.character(season), from =areas, to=areas)) 
mig2[,,,"a","a"] <- -migconstant
mig2[,,,"b","b"] <- -migconstant
mig2[,,,"a","b"] <- migconstant
mig2[,,,"b","a"] <- migconstant
aperm( mig2,c(1,3,2,4,5))

effort <- array(c(1), dim=c(length(areas), length(season)), dimnames=list(option =areas,season=as.character(season)))

pop1  <- pop2         <-array(0, dim=c(length(ages),endy+1,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:(endy+1)), season=as.character(season), option =areas))
catches.n.dsvm1       <- catches.n.dsvm2       <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.n.dsvm1      <- landings.n.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm1      <- catches.wt.dsvm2      <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
landings.wt.dsvm1     <- landings.wt.dsvm2     <- array(0, dim=c(length(ages),endy,length(season),length(areas)), dimnames=list(cat=ages,   year=as.character(1:endy), season=as.character(season), option =areas))
catches.wt.dsvm.tot1  <- catches.wt.dsvm.tot2  <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
landings.wt.dsvm.tot1 <- landings.wt.dsvm.tot2 <- array(0, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))
quota1                <- quota2                <- array(5000, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

pos_catches1 <- pos_catches2 <- pop1

hr1wanted <- hr2wanted <- array(NA, dim=c(1           ,endy,              1,            1), dimnames=list(cat="all", year=as.character(1:endy), season="all",                option ="all"))

#run population for 15 year
pop1 <- population_dynamics(pop=pop1, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs1, migration=mig1)
pop2 <- population_dynamics(pop=pop2, startyear=2, endyear=stab.model, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,1,,, drop=F], recruitment=recs2, migration=mig2)


#set up dsvm
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")
catchMean(sp1)  <- catchMean(sp2) <- array(NA, dim=c(length(ages), length(season),length(areas)),  dimnames=list("cat"=ages,"season"= season,"option"=areas))

catchMean(sp3)  <- catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))
catchSigma(sp3) <- catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=c(length(ages),length(season),length(areas)),dimnames=list(cat=ages,season=as.character(season),option =areas))

#SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
sp1Price <- array(c(sp1price + slope1price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp2Price <- array(c(sp2price + slope2price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
sp3Price <- sp4Price <- sp5Price <- array(c(0), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))
#---effort and prices used (note that now c is removed (but that if other runs, then make sure to fix/remove code that removes "c" option)                                                                                         
control     <- DynState.control(spp1LndQuota= quota1[,1,,],  spp2LndQuota=quota2[,1,,], spp1LndQuotaFine= 3e6, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = FuelPrice, landingCosts= 0,gearMaintenance= 0, 
                                addNoFishing= TRUE, increments= 25, spp1DiscardSteps= SPP1DSCSTEPS, spp2DiscardSteps= SPP2DSCSTEPS, sigma= SIGMA, simNumber= SIMNUMBER * 0.5 , numThreads= 68, verbose=1)

#this is where our loop starts, after we set up stable population
for(yy in (stab.model):(stab.model+NUMRUNS)){
  
  print("====== year yy ========")
  print(yy)

  #Calculate catch inputs for DSVM
  #for age 1, assume that catchmean over seasons and areas will be equal to the year before
  catchMean(sp1)[1,,] <- pop1[1,as.character(yy-1),,] *q*wts[1,1,,1]
  #for older ages we start from cohort i nthe previous year
  #we need to guess how much population will decline while moving into this year, and over seasons, for this wee need to look at end of two years ago and compare to last year 
  #this comes from a sweep with pop1[1:(max(ages)-1),as.character(yy-2),max(season),] and pop1[-1,as.character(yy-1),,] 
  cfracs1 <- sweep( pop1[-1,as.character(yy-1),,] ,c(1,3),  pop1[1:(max(ages)-1),as.character(yy-2),max(season),] , FUN="/")
  cfracs2 <- sweep( pop1[-1,as.character(yy-2),,] ,c(1,3),  pop1[1:(max(ages)-1),as.character(yy-3),max(season),] , FUN="/")
  cfracs <- (cfracs1 + cfracs2)/2
  #once we have the decline as a functrion of time (in cfracs) we will mulitply by cohorts in last season, and multiply by q and wts for appropriate ages. Those go in catchmean for older ages
  #it is a two step approach because we need two sweeps
  catchMean(sp1)[-1,,] <- sweep(cfracs,c(1,3),pop1[1:(max(ages)-1),as.character(yy-1),max(season),] ,FUN="*") *q
  catchMean(sp1)[-1,,] <- sweep(catchMean(sp1)[-1,,],c(1,2), wts[-1,1,,1],FUN="*")
  
  catchMean(sp2)[1,,] <- pop2[1,as.character(yy-1),,] *q*wts[1,1,,1]
  cfracs1 <- sweep( pop2[-1,as.character(yy-1),,] ,c(1,3),  pop2[1:(max(ages)-1),as.character(yy-2),max(season),] , FUN="/")
  cfracs2 <- sweep( pop2[-1,as.character(yy-2),,] ,c(1,3),  pop2[1:(max(ages)-1),as.character(yy-3),max(season),] , FUN="/")
  cfracs <- (cfracs1 + cfracs2)/2
  catchMean(sp2)[-1,,] <- sweep(cfracs,c(1,3),pop2[1:(max(ages)-1),as.character(yy-1),max(season),] ,FUN="*") *q
  catchMean(sp2)[-1,,] <- sweep(catchMean(sp2)[-1,,],c(1,2), wts[-1,1,,1],FUN="*")
  
    
  # ---No way of estimating sigma, therefore we assume that is 8% of the CPUE (note slight repetion in code for dims and dimnames of 0 catch arrays for spec 3,4,5)                                                                  
  catchSigma(sp1) <- catchMean(sp1) *0.08
  catchSigma(sp2) <- catchMean(sp2) *0.08
  
  # if we are in MP period, then set quota based on last year
  if (yy > MPstart){
    control@spp1LndQuota <-  quota1[,yy,,]
    #control@spp2LndQuota <-  quota2[,yy,,]
  }
  
  # because we wanted gradual introduction of vessels, we only put all vessels after stabmodel. First line is not strictly necessary becuase already defined when constructing control 
  if (yy == (stab.model)  ) control@simNumber <-  as.integer(SIMNUMBER*0.5)
  if (yy == (stab.model+1)) control@simNumber <-  as.integer(SIMNUMBER*0.75)
  if (yy >  (stab.model+1)) control@simNumber <-  as.integer(SIMNUMBER)
          
  
  # if we are in MPLO period, then set discardsteps =0 
  #if (yy > MPstartLO){
   # control@spp1DiscardSteps <-  0
  #  control@spp2DiscardSteps <-  0
  #}

  #run DSVM (wiht quota constraining if in MP time)   
  z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
  
  #extract DSVM results
  dsvm_res <-  extract_dsvm_res (z, control, ages, season)
  #Net revenuw from DSVM
  economics_res             <- as.data.frame(as.matrix(netRev(z)))
  names(economics_res )     <- "NetRev"
  economics_res$Grossrev    <- as.data.frame(as.matrix(grossRev(z)))$V1 
  economics_res$Fuelcosts   <- apply(effort(sim(z)) * control(z)@fuelUse * control(z)@fuelPrice, 2, sum)
  economics_res$Annualfine  <- as.data.frame(as.matrix(annualFine(z)))$V1

  if (yy == stab.model){
    dsvm_res_allyrs       <- cbind("year"= yy,dsvm_res)
    economics_res_allyrs  <- cbind("year"= yy,economics_res)
  } else {
    dsvm_res_allyrs <- rbind(dsvm_res_allyrs, (cbind("year"= yy,dsvm_res)))
    economics_res_allyrs <- rbind(economics_res_allyrs, (cbind("year"= yy,economics_res)))
  }
  
  #get catches in wts from DSVM 
  catches.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="catch.wt")
  catches.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="catch.wt") 
  
  landings.wt.dsvm1[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop1", catch_option="landings.wt")
  landings.wt.dsvm2[,yy,,] <- catch_dataframe_to_array(dsvm_res, ages, season, areas, "pop2", catch_option="landings.wt") 
  
  rm(z,dsvm_res)
  
  #calculate total catches (by summing over seasons and ages)
  catches.wt.dsvm.tot1[] <- apply(catches.wt.dsvm1,c(2),"sum")
  catches.wt.dsvm.tot2[] <- apply(catches.wt.dsvm2,c(2),"sum")
  
  landings.wt.dsvm.tot1[] <- apply(landings.wt.dsvm1,c(2),"sum")
  landings.wt.dsvm.tot2[] <- apply(landings.wt.dsvm2,c(2),"sum")
  
  #calculate numbers caught from weight caught 
  for(jj in areas){
    catches.n.dsvm1 [,yy,,as.character(jj)]  <- catches.wt.dsvm1[,yy,,as.character(jj)] / wts[,1,,1]
    catches.n.dsvm2 [,yy,,as.character(jj)]  <- catches.wt.dsvm2[,yy,,as.character(jj)] / wts[,1,,1]
    
    landings.n.dsvm1 [,yy,,as.character(jj)] <- landings.wt.dsvm1[,yy,,as.character(jj)] / wts[,1,,1]
    landings.n.dsvm2 [,yy,,as.character(jj)] <- landings.wt.dsvm2[,yy,,as.character(jj)] / wts[,1,,1]
  }

  # calculatae what happens to population based on catches
  pop1 <- population_dynamics(pop=pop1, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm1[,yy,,,drop=F], recruitment=recs1, migration=mig1)
  pop2 <- population_dynamics(pop=pop2, startyear=yy, endyear=yy+1, season=season, natmortality=natmortality, catches=catches.n.dsvm2[,yy,,,drop=F], recruitment=recs2, migration=mig2)
  
  #------------------------------------------
  #MANAGEMENT PROCEDURE
  #------------------------------------------
  
 
  #calculate the observed harvest ratios (per age) 
  hr1 <- (apply(catches.n.dsvm1,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm1,1:3,sum)+1e-20) +  apply(pop1[,1:endy,,],1:3,sum))
  hr2 <- (apply(catches.n.dsvm2,1:3,sum)+1e-20)/    (apply(catches.n.dsvm2,1:3,sum) +  apply(pop2[,1:endy,,],1:3,sum))
  
  #calculate landing ratios
  landings.ratio1 <- (apply(landings.n.dsvm1,1:3,sum)+1e-20)/ (apply(catches.n.dsvm1,1:3,sum)+1e-20)
  landings.ratio2 <- (apply(landings.n.dsvm2,1:3,sum)+1e-20)/ (apply(catches.n.dsvm2,1:3,sum)+1e-20)

  
  # if (yy == MPstartLO){
  #   landings.ratio1[,yy,]<- 1
  #   landings.ratio2[,yy,]<- 1
  # }
  
  #what sequence to use in YPR curve so that we span at least hr= 0.001 to 0.5?
  maxseq1 <- 0.5/mean(hr1[,yy,])
  maxseq2 <- 0.5/mean(hr2[,yy,])
  
  # get yield curve results given observed harvest ratios
  yc1 <- yield_curve(hr=hr1[,yy,], landings.ratio1[,yy,], wts, natmortality, R=recs1, verbose=F, sequence = seq(0.001, maxseq1, length.out=2000))
  yc2 <- yield_curve(hr=hr2[,yy,], landings.ratio2[,yy,], wts, natmortality, R=recs2, verbose=F, sequence = seq(0.001, maxseq2, length.out=2000))
  
  #store perceived HRmax    
  hr1wanted[,yy,,] <- (yc1[yc1$landings==max(yc1$landings),]$hr)[1]
  hr2wanted[,yy,,] <- (yc2[yc2$landings==max(yc2$landings),]$hr)[1]
  
  #calculate individual quota from pop, hr and hrmax 
  if (yy > (MPstart-1) & yy < endy){ #if (yy > (MPstart-1))
    prel.quota1 <- sum((hr1wanted[,yy,,]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER 
      #sum(sweep((hr1wanted[1]/mean(hr1[,yy,]))* hr1[,yy,]*landings.ratio1[,yy,]*apply(pop1[,yy+1,,],c(1,2), sum) ,1,wts[,1,,1],"*"))/SIMNUMBER
    prel.quota2 <- sum((hr2wanted[,yy,,]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum)*wts[,1,,1])/SIMNUMBER
      #sum(sweep((hr2wanted[1]/mean(hr2[,yy,]))* hr2[,yy,]*landings.ratio2[,yy,]*apply(pop2[,yy+1,,],c(1,2), sum) ,1,wts,"*"))/SIMNUMBER
    
    # We constrained +- 15% TAC change, until the LO implementation, where we allow the HR wanted just for the transition
    # In the transition we assume landing ratio equal 1
    #   if (yy == MPstart){ #MPstartLO){ # landing.ratio are equal to 1 for all years before estimating them
      quota1[,yy+1,,] <- prel.quota1
      #quota2[,yy+1,,] <- prel.quota2
      
    #}else{
    #  tac.constrained1 <- c(quota1[,yy,,]*0.85, quota1[,yy,,]*1.15)
    #  tac.constrained2 <- c(quota2[,yy,,]*0.85, quota2[,yy,,]*1.15)
    #  quota1[,yy+1,,] <- max(min(prel.quota1, tac.constrained1[2]), tac.constrained1[1]) 
    #  quota2[,yy+1,,] <- max(min(prel.quota2, tac.constrained2[2]), tac.constrained2[1]) 
    #}
  }
}


#what are the weights?
wts

pyrnoMP   <- MPstart- 3
pyrMP     <- endy   - 3
# pyrMP     <- MPstartLO - 4
# pyrMPLO   <- endy   - 4


# Scenarios II: full landings selectivity, when discarding is allowed
hr1 <- (apply(catches.n.dsvm1,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm1,1:3,sum)+1e-20) +    apply(pop1[,1:endy,,],1:3,sum))
hr2 <- (apply(catches.n.dsvm2,1:3,sum)+1e-20)/    ((apply(catches.n.dsvm2,1:3,sum)+1e-20) +    apply(pop2[,1:endy,,],1:3,sum))
landings.ratio1 <- (apply(landings.n.dsvm1,1:3,sum)+1e-20)/ (apply(catches.n.dsvm1,1:3,sum)+1e-20)
landings.ratio2 <- (apply(landings.n.dsvm2,1:3,sum)+1e-20)/ (apply(catches.n.dsvm2,1:3,sum)+1e-20)

yc1noMP <- yield_curve(hr=hr1[,pyrnoMP,], landings.ratio1[,pyrnoMP,], wts, natmortality, R=recs1, verbose=F, sequence=seq(0.0001,0.5/mean(hr1[,pyrnoMP,]),0.01 ))
yc2noMP <- yield_curve(hr=hr2[,pyrnoMP,], landings.ratio2[,pyrnoMP,], wts, natmortality, R=recs2, verbose=F, sequence=seq(0.0001,0.5/mean(hr2[,pyrnoMP,]),0.01 ))
#what happens in our yield curve for this hr?
#yield_curve(hr=hr1[,pyrnoMP,], landings.ratio1[,pyrnoMP,], wts, natmortality, R=recs1, sequence = 1, verbose=T)
#yield_curve(hr=hr2[,pyrnoMP,], landings.ratio2[,pyrnoMP,], wts, natmortality, R=recs2, sequence = 1, verbose=T)

yc1MP <- yield_curve(hr=hr1[,pyrMP,], landings.ratio1[,pyrMP,], wts, natmortality, R=recs1, verbose=F, sequence=seq(0.0001,0.5/mean(hr1[,pyrMP,]),0.01 ))
yc2MP <- yield_curve(hr=hr2[,pyrMP,], landings.ratio2[,pyrMP,], wts, natmortality, R=recs2, verbose=F, sequence=seq(0.0001,0.5/mean(hr2[,pyrMP,]),0.01 ))

# yc1MPLO <- yield_curve(hr=hr1[,pyrMPLO,], landings.ratio1[,pyrMPLO,], wts, natmortality, R=recs1, verbose=F)
# yc2MPLO <- yield_curve(hr=hr2[,pyrMPLO,], landings.ratio2[,pyrMPLO,], wts, natmortality, R=recs2, verbose=F)
 
Fmsy1noMP <- yc1noMP[yc1noMP$landings==max(yc1noMP$landings),]$hr
Fmsy2noMP <- yc2noMP[yc2noMP$landings==max(yc2noMP$landings),]$hr

Fmsy1MP <- yc1MP[yc1MP$landings==max(yc1MP$landings),]$hr
Fmsy2MP <- yc2MP[yc2MP$landings==max(yc2MP$landings),]$hr

# Fmsy1MPLO <- yc1MPLO[yc1MPLO$landings==max(yc1MPLO$landings),]$hr
# Fmsy2MPLO <- yc2MPLO[yc2MPLO$landings==max(yc2MPLO$landings),]$hr

#round(catches.n.dsvm1[,pyrnoMP,,],2)
#round(pop1[,pyrnoMP,,],2)

#next, what happens to theoretical pop for our harvest, what is harvest?
#hr1[,pyrnoMP,]
#mean(hr1[,pyrnoMP,])

# Effort pattern and economics
#--------------------------------------------------------------------------------------
#effort_plot_dsvm(SIMNUMBER, dsvm_res_allyrs, stab.model, economics_res_allyrs)

#effort
effort      <- subset(dsvm_res_allyrs,(spp %in% "sp1"))
effort      <- subset(effort,(cat %in% 1))
effort$year <- factor(effort$year, levels= 1:endy)
effort      <- aggregate(cbind(landings.wt, discards.wt, catch.wt, effort,trip)~ spp+cat+option+year, FUN=sum, data=effort)

trip      <- with(effort,data.frame(year, trip, option))
trip      <- dcast(trip ,option~year, value.var="trip")
trip[is.na(trip)]<- 0
rownames(trip)  <- trip$option
trip            <- trip[,-1]
trip_percentage <- matrix(apply(trip, 2, function(x){x*100/sum(x,na.rm=T)}), nrow=nrow(trip))
rownames(trip_percentage) <- rownames(trip)
colnames(trip_percentage) <- colnames(trip)

days      <- with(effort,data.frame(year, effort, option))
days<- days[!days$option=="Stay in port",]
days      <- dcast(days ,option~year, value.var="effort")
days[is.na(days)]<- 0
rownames(days)  <- days$option
days            <- days[,-1]

# Estimate mean effort days by season in the two study periods. Three year mean before the MP and after the MP
daysmen        <- subset(dsvm_res_allyrs,(spp %in% "sp1"))
daysmen        <- subset(daysmen ,(cat %in% 1))
daysmen $year   <- factor(daysmen $year, levels= 1:endy)
daysmen $season <- factor(daysmen $season, levels= 1:max(season))
daysmen         <- aggregate(effort~ spp+cat+option+year+season, FUN=sum, data=daysmen )
daysmen         <- with(daysmen ,data.frame(year, effort, option, season))
daysmen <- daysmen [!daysmen $option=="Stay in port",]
#days      <- dcast(days ,option~year, value.var="effort")
daysmen[is.na(daysmen)]<- 0

nompdays <- aggregate(effort~ option+season, FUN=mean, data=subset(daysmen,(year %in% c((MPstart-3):MPstart))))
mpdays <- aggregate(effort~ option+season, FUN=mean, data=subset(daysmen,(year %in% c((yy-3):yy))))

# First two rows are the NOMP area north and south. And 2 last rows, MP area north and south
daysnomp <- dcast(nompdays ,option~season, value.var="effort")
rownames(daysnomp) <- paste("nomp",daysnomp$option, sep="_")
                            
daysmp <- dcast(mpdays ,option~season, value.var="effort")
rownames(daysmp) <- paste("mp",daysmp$option, sep="_")
daysmen          <- rbind(daysnomp,daysmp)
daysmen[,-1]

#economics
netrev      <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("NetRev"))
grossrev    <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("Grossrev"))
annualfine  <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("Annualfine"))
fuelcosts   <-  melt(economics_res_allyrs, id.vars = "year", measure.vars = c("Fuelcosts"))

save.image(paste0("1quota12seasons_",SIGMA,".RData"))
#save.image("2quotas12seasons.RData")

rm(dsvm_res,dsvm_res_allyrs)
rm(economics_res,economics_res_allyrs)
rm(control, sp1, sp2, sp3, sp4, sp5)
}
