#################################################
#################################################
# FACE MASK USE AND OPTIMAL RESOURCE ALLOCATION
# DURING THE COVID-19 PANDEMIC
# Colin Worby & Hsiao-Han Chang
# cworby@broadinstitute.org
# hhchang@life.nthu.edu.tw
# 2020-04-26
#################################################
#################################################
library(deSolve)
###################################
# Specify basic model structure
# S - E - IP - IS/IA - R
# No mask use/interventions

basic2Imodel <- function(t, y, params) {
  
  S <- y[1]; E <- y[2]; IP <- y[3]; IS <- y[4]; IA <- y[5]; R <- y[6]
  N <- S+E+IP+IS+IA+R 
  
  with(as.list(c(params)), {

    dS <- -S*beta*(IS+asym.inf*(IA+IP))/N # susceptible
    dE <- S*beta*(IS+asym.inf*(IA+IP))/N - a1*E # exposed
    dIP <- a1*E - a2*IP # presymptomatic (but infectious)
    dIS <- a2*(1-p.asym)*IP - gS*IS # symptomatic
    dIA <- a2*p.asym*IP - gA*IA # mild/asymptomatic
    dR <- gS*IS+gA*IA # removed
    ret <- c(dS, dE, dIP, dIS, dIA, dR)
    list(ret)
  })
}

# Specify parameters
R0 <- 2.1 # basic reproduction number
asym.inf <- 0.5 # relative infectiousness of asymptomatic/mildly symptomatic cases
p.asym <- 0.3 # proportion of cases that are asymptomatic/mildly/ symptomatic
N <- 10000 # population size

# Progression rates
a1 <- 1/1 # E -> IP
a2 <- 1/5 # IP -> I
gS <- 1/14 # IS -> R
gA <- 1/14 # IA -> R

# infection rate under the parameters specified
beta <- R0*a2*gA*gS/(asym.inf*gA*gS + asym.inf*p.asym*a2*gS + (1-p.asym)*a2*gA)
# list of parameters
params <- list(beta=beta, asym.inf=asym.inf, p.asym=p.asym,
               a1=a1, a2=a2, gA=gA, gS=gS)
# Vector of starting values
start <- c(N-1,1,0,0,0,0)
# Time steps to run model
steps <- 1:500
# Run model
out <- lsoda(y=start, times=steps, func=basic2Imodel, parms=params)

# Plot epidemic curves
plot(out[,1], out[,5], type="l", las=1, bty="l", xlab="Time", ylab="Count")
lines(out[,1], out[,6], col="red")
legend("topright", legend=c("IS","IA"), col=c("black","red"), lty=1)
legend("topleft", legend=paste("Final infected: ", round(100*out[max(steps),7]/N,2), "%", sep=""), bty="n")

# Theoretical R0 based on final size
#log(out[max(steps),2]/N)/(out[max(steps),2]/N-1)

#####################################
# RESOURCE ALLOCATION MODEL
# Masks are given to a proportion of susceptibles at the start
# Masks can be provided to a proportion of infected (& detected) individuals
# Mask resources are limited to cover X% of the population

SEIRD_mask_age <- function(t,x,params) {
  
  Se=x[1:2]; Ee=x[3:4]; Ie.pre=x[5:6]; Ie.low_detect=x[7:8]; Ie.low_missed=x[9:10]; Ie.high=x[11:12]
  Sm=x[13:14]; Em=x[15:16]; Im.pre=x[17:18]; Im.low_detect=x[19:20]; Im.low_missed=x[21:22]; Im.high=x[23:24] 
  R=x[25:26]; D=x[27:28]; Nm=x[29] 

  with(as.list(c(params)), {  
    # Are masks available??
    if (Nm>mask.lim) {
      inf.mask.use <- c(0,0)
    }

    dSe <- -Se*inf.rate*(sum( c(sum(Ie.pre+Ie.low_detect+Ie.low_missed)*asym.inf,sum(Ie.high)) ) + 
                           sum( c(sum(Im.pre+Im.low_detect+Im.low_missed)*asym.inf,sum(Im.high)) )*mask.trans)/sum(N)

    dSm <- -Sm*inf.rate*mask.sus*(sum( c(sum(Ie.pre+Ie.low_detect+Ie.low_missed)*asym.inf,sum(Ie.high)) ) +
                                    sum( c(sum(Im.pre+Im.low_detect+Im.low_missed)*asym.inf,sum(Im.high)) )*mask.trans)/sum(N)

    dEe <- Se*inf.rate*(sum( c(sum(Ie.pre+Ie.low_detect+Ie.low_missed)*asym.inf,sum(Ie.high)) ) +
                          sum( c(sum(Im.pre+Im.low_detect+Im.low_missed)*asym.inf,sum(Im.high)) )*mask.trans)/sum(N) - inc.rate*Ee

    dEm <- Sm*inf.rate*mask.sus*(sum( c(sum(Ie.pre+Ie.low_detect+Ie.low_missed)*asym.inf,sum(Ie.high)) ) +
                                   sum( c(sum(Im.pre+Im.low_detect+Im.low_missed)*asym.inf,sum(Im.high)) )*mask.trans)/sum(N) - inc.rate*Em

    dIe.pre <- inc.rate*Ee - prog.rate*Ie.pre
    dIm.pre <- inc.rate*Em - prog.rate*Im.pre
    dIe.low_detect <- detection.rate*p.asym*prog.rate*(1-inf.mask.use)*Ie.pre - rem.rate[1]*Ie.low_detect
    dIm.low_detect <- detection.rate*p.asym*prog.rate*(Im.pre + inf.mask.use*Ie.pre) - rem.rate[1]*Im.low_detect
    dIe.low_missed <- (1-detection.rate)*p.asym*prog.rate*Ie.pre - rem.rate[1]*Ie.low_missed
    dIm.low_missed <- (1-detection.rate)*p.asym*prog.rate*Im.pre - rem.rate[1]*Im.low_missed
    dIe.high <- (1-p.asym)*prog.rate*(1-inf.mask.use)*Ie.pre - rem.rate[2]*Ie.high
    dIm.high <- (1-p.asym)*prog.rate*(Im.pre + inf.mask.use*Ie.pre) - rem.rate[2]*Im.high
    dR <- rem.rate[1]*(Ie.low_detect+Im.low_detect+Ie.low_missed+Im.low_missed)+rem.rate[2]*(1-death.prob)*(Ie.high+Im.high)
    dD <- rem.rate[2]*death.prob*(Ie.high+Im.high)
    dNm <- sum(prog.rate*inf.mask.use*Ie.pre*((1-p.asym) + p.asym*detection.rate))
    drecorded <- (1-p.asym)*prog.rate*(Ie.pre+Im.pre)

    list(c(dSe,dEe,dIe.pre,dIe.low_detect,dIe.low_missed,dIe.high,
           dSm,dEm,dIm.pre,dIm.low_detect,dIm.low_missed,dIm.high,dR,dD,drecorded,dNm))
  })
}

# Set parameters

R0 <- 2.5 # basic reproduction number
asym.inf <- 0.5 # relative infectiousness of mildly/asymptomatic cases
N <- 2E6*c(0.924,0.076) # population structure: young (<70), old (70+)
inc.rate <- 1/1 # incubation E -> I.pre
prog.rate <- 1/5 # progression I.pre -> I.low/high
rem.rate <- c(1/14, 1/14) # mild/symptomatic removal rates
death.prob <- c(0.013,0.097) # YOUNG/ELDERLY - proportion of recorded (severe) cases dying
detection.rate <- 0.44 # detecting low symptomatic cases
p.asym <- c(0.4,0.2) # YOUNG/ELDERLY probability of being asymptomatic
# infection rate derived from specified parameters
inf.rate <- R0*prog.rate*rem.rate/(asym.inf*rem.rate + asym.inf*p.asym*prog.rate + (1-p.asym)*prog.rate)
inf.rate <- sum(inf.rate*(N/sum(N)))
steps <- 1:500 # time steps
inf.init <- 100 # initial infected cases
sus.mask.use <- c(0.1,0.1) # proportion of susceptibles (YOUNG/ELDERLY) wearing a mask at the start of the epidemic
mask.lim <- 0.2*sum(N) # additional masks available for infected cases
inf.mask.use <- c(0.5,0.5) # proportion of detected cases (not previously wearing a mask) who adopt mask use (YOUNG/ELDERLY)
mask.sus <- 0.75 # relative susceptibility of a mask wearer (1-protection)
mask.trans <- 0.5 # relative transmissibility of a mask wearer (1-containment)

# parameter vector
params <- list(N=N, inf.rate=inf.rate, rem.rate=rem.rate, inc.rate=inc.rate, detection.rate=detection.rate,
               prog.rate=prog.rate, inf.mask.use=inf.mask.use, mask.trans=mask.trans, asym.inf=asym.inf,
               mask.lim=mask.lim, mask.sus=mask.sus, death.prob=death.prob, p.asym=p.asym)
# start vector
x=c(N*(1-sus.mask.use),0,0,inf.init,0,0,0,0,0,0,0,
    N*sus.mask.use,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0)

# Run model
out = as.data.frame(ode(x, steps, SEIRD_mask_age, params)) 
comps <- c("S","E","I.pre","I.low_detect","I.low_missed","I.high") # compartment names
# output names
nm <- c("time",paste(rep(paste(as.character(sapply(comps,rep,2)),c("y","o"), sep="_"),2), 
                     c(rep("e",12),rep("m",12)), sep=""),
        paste(as.character(sapply(c("R","D","recorded"),rep,2)),c("y","o"), sep="_"), "Nm")
names(out) <- nm # y=young, o=elderly, m=mask, e=no mask 

# Plot output
infcurve <- apply(out[,grep("I",names(out))],1,sum)/sum(N)
plot(steps, infcurve, lwd=2, type="l", bty="l", las=1, ylab="Prop. population", xlab="Time")
lines(steps, (out$D_o+out$D_y)/sum(N), lwd=2, col="red")
legend("topright", col=c("black","red"), legend=c("Cur. infected", "Cum. deaths"), lwd=2)
legend("topleft", legend=c(paste("Final infected: ", round(100*sum(out[max(steps),grep("R",nm)])/sum(N),2), "%", sep=""),
                           paste("Final dead: ", round(100*sum(out[max(steps),grep("D",nm)])/sum(N),2), "%", sep="")), bty="n")

# Compare dynamics across different mask deployment strategies
tot.masks <- 0.4 # resources to cover this fraction of the population
mask.sus <- 0.75 # relative susceptibility of a mask wearer (1-protection)
mask.trans <- 0.5 # relative transmissibility of a mask wearer (1-containment)
steps <- 1:1000
model_names <- c("S1", "S2", "S3a", "S3b", "S3c", "S4")
full_model_names <- c("100% naive", "100% elderly", "75% elderly+25% inf", 
                      "50% elderly+50% inf", "25% elderly+75% inf", "100% inf")
inf.mask.use <- rbind(c(0,0),c(0,0),c(1,1),c(1,1),c(1,1),c(1,1))
sus.mask.use <- rbind(rep(tot.masks,2),
                      c(max(c(0,sum(N)*tot.masks-N[2]))/N[1],min(c(N[2],(sum(N)*tot.masks)))/N[2]),
                      c(max(c(0,sum(N)*tot.masks*(1-0.25)-N[2]))/N[1],min(c(N[2],sum(N)*tot.masks*(1-0.25)))/N[2]),
                      c(max(c(0,sum(N)*tot.masks*(1-0.5)-N[2]))/N[1],min(c(N[2],sum(N)*tot.masks*(1-0.5)))/N[2]),
                      c(max(c(0,sum(N)*tot.masks*(1-0.75)-N[2]))/N[1],min(c(N[2],(sum(N)*tot.masks*(1-0.75))))/N[2]),
                      c(0,0))
mask.lim <- c(0,0,0.25,0.5,0.75,1)*sum(N)*tot.masks

colz <- c("pink", "lightgreen", "lightblue", "steelblue1", "royalblue", "salmon")
plot(NULL, xlim=c(0,max(time)), ylim=c(0,0.2), las=1, xlab="Time", ylab="Prop. population", bty="l")
for (j in 1:length(model_names)) {
  params <- list(N=N, inf.rate=inf.rate, rem.rate=rem.rate, inc.rate=inc.rate, detection.rate=detection.rate,
                 prog.rate=prog.rate, inf.mask.use=inf.mask.use[j,], mask.trans=mask.trans, asym.inf=asym.inf,
                 mask.lim=mask.lim[j], mask.sus=mask.sus, death.prob=death.prob, p.asym=p.asym)
  x=c(N*(1-sus.mask.use[j,]),0,0,inf.init,0,0,0,0,0,0,0,
      N*sus.mask.use[j,],0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0)
  out = as.data.frame(ode(x, steps, SEIRD_mask_age, params)) 
  names(out) <- nm
  infcurve <- apply(out[,grep("I",names(out))],1,sum)/sum(N)
  lines(steps, infcurve, col=colz[j], lwd=2)
}
legend("topright", legend=full_model_names, lwd=2, col=colz, title="Mask allocation")


##############################################
#### supply and demand model 

seirAP.ode <- function(t,x,params) {
  
  Se= x[1]        #susceptible, mask non-wearers
  Sm= x[2]        #susceptible, mask wearers
  Ee= x[3]        #exposed, mask non-wearers
  Em= x[4]        #exposed, mask wearers
  Ipe= x[5]       #pre-symptomatic, mask non-wearers
  Ipm= x[6]       #pre-symptomatic, mask wearers
  Iae= x[7]       #asymptomatic, mask non-wearers
  Iam= x[8]       #asymptomatic, mask wearers
  Ise= x[9]       #symptomatic, mask non-wearers
  Ism= x[10]      #symptomatic, mask wearers
  Re= x[11]       #recover, mask non-wearers
  Rm= x[12]       #recover, mask wearers
  M= x[13]        #mask 
  C= x[14]        #cumulative number of infections
  
  a <- params$a   #ratio of infectiouseness (beta_s/beta_a)
  N <- params$N   #population size
  alpha1 <- params$alpha1       #rate from exposed to pre-symptomatic
  alpha2 <- params$alpha2       #rate from pre-symptomatic to asymptomatic/symptomatic
  gamma_a <- params$gamma_a     #recovery rate for asymptomatic
  gamma_s <- params$gamma_s     #recovery rate for symptomatic
  R0 <- params$R0               #R0
  Pa <- params$Pa               #prop of asymptomatic infections
  Ps <- params$Ps               #relative susceptibility with mask
  Pi <- params$Pi               #relative transmissibility with mask
  omega_s_factor <- params$omega_s_factor       #ratio of rate of wearing masks (omega_s/omega_a)
  omega_a_k1 <- params$omega_a_k1               #k1, the rate of demand increase
  omega_a_k2 <- params$omega_a_k2               #k2, the timing of demand, defined as the number of reported cases when half the population seeks face masks
  mu <- params$mu                               #1/mu = average time of wearing disposable mask
  B <- params$B                                 #daily mask production
  omega_a= (1/(1+ exp(-omega_a_k1*(C-omega_a_k2))))    #the rate of healthy and asymptomatic non-wearers acquiring masks 
  omega_s= omega_a*omega_s_factor                      #the rate of symptomatic non-wearers acquiring masks 
  beta_s= R0*(2*alpha2*gamma_s/(2*alpha2 + gamma_s- Pa*alpha2)) #infection rate of symptomatic individuals 
  beta_a= beta_s*a                                              #infection rate of asymptomatic individuals 
  
  dSe <- -(beta_a*(Ipe + Iae + Ipm*Pi + Iam*Pi) + beta_s*(Ism*Pi + Ise))*Se/N + mu*Sm - omega_a*M/N*Se
  dSm <- -(beta_a*(Ipe + Iae + Ipm*Pi + Iam*Pi) + beta_s*(Ism*Pi + Ise))*Sm/N*Ps + omega_a*M/N*Se - mu*Sm
  dEe <- (beta_a*(Ipe + Iae + Ipm*Pi + Iam*Pi) + beta_s*(Ism*Pi + Ise))*Se/N - alpha1*Ee + mu*Em - omega_a*M/N*Ee
  dEm <- (beta_a*(Ipe + Iae + Ipm*Pi + Iam*Pi) + beta_s*(Ism*Pi + Ise))*Sm/N*Ps - alpha1*Em + omega_a*M/N*Ee - mu*Em
  dIpe <- alpha1*Ee - omega_a*M/N*Ipe + mu*Ipm - alpha2*Ipe
  dIpm <- alpha1*Em + omega_a*M/N*Ipe - mu*Ipm - alpha2*Ipm
  dIae <- alpha2*Ipe*Pa - Iae*gamma_a + mu*Iam - omega_a*M/N*Iae
  dIam <- alpha2*Ipm*Pa - Iam*gamma_a - mu*Iam + omega_a*M/N*Iae
  dIse <- alpha2*Ipe*(1-Pa) - Ise*gamma_s + mu*Ism - omega_s*M/N*Ise
  dIsm <- alpha2*Ipm*(1-Pa) - Ism*gamma_s - mu*Ism + omega_s*M/N*Ise
  dRe <- Iae*gamma_a + Ise*gamma_s + mu*Rm - omega_a*M/N*Re
  dRm <- Iam*gamma_a + Ism*gamma_s - mu*Rm + omega_a*M/N*Re
  dM <- B - ((Se + Ee + Ipe + Iae + Re)*omega_a+Ise*omega_s)*M/N
  dC <- alpha2*Ipe*(1-Pa) + alpha2*Ipm*(1-Pa)
  return(list(c(dSe, dSm, dEe, dEm, dIpe, dIpm, dIae, dIam, dIse, dIsm, dRe, dRm, dM, dC)))
}

#### basic demo ####
# parameters
R0= 2.5
a= 0.5                  #the ratio of asymptomatic infectiousness to symptomatic infectiousness 
Pa= 0.3                 #the probability of being asymptomatic
gamma_a= 1/14           #recovery rate for asymptomatic cases
gamma_s= 1/14           #recovery rate for symptomatic cases
alpha1= 1/1             #rate from exposed to presymptomatic
alpha2= 1/5             #rate from presymptomatic to asymptomatic/symptomatic
mu= 1                   #1/mu = average time of wearing disposable mask
k1= 0.1                 #the rate of demand increase
k2= 100                 #the timing of demand, defined as the number of reported cases when half the population seeks face masks
Pi= 0.5                 #relative transmissibility with mask
Ps= 0.75                #relative susceptibility with mask
omega_s_factor = 1000   #ratio of rate of wearing masks (omega_s/omega_a)
B = 2*10^6              #daily mask production
N = 23.6*10^6           #population size
params <- list(N=N, alpha1= alpha1, alpha2= alpha2, gamma_a= gamma_a, gamma_s= gamma_s, R0= R0, Pa= Pa, Ps= Ps, Pi= Pi, a=a,
               omega_s_factor= omega_s_factor, mu= mu, omega_a_k1= k1, omega_a_k2= k2, B= B)

x <- c(23576380, 23600, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0, B, 20) #initial conditions (Se, Sm, Ee, Em, Ipe, Ipm, Iae, Iam, Ise, Ism, Re, Rm, M, C)
out = as.data.frame(ode(x, 1:1000, seirAP.ode, params))
colnames(out)= c("time", "Se", "Sm", "Ee", "Em", "Ipe", "Ipm","Iae", "Iam","Ise" ,"Ism", "Re", "Rm", "M","C")

# make plot
plot(NULL, type = "l", xlab="Time (day)", ylab="Count",
     ylim=c(0,max(out[,2:9])), xlim=range(out$time), bty="l", las=1)
cols <- c("darkgreen", "darkgreen", "orange", "orange","purple","purple" , "pink", "pink", "red", "red")
lty <- rep(1:2,5)
for (i in 1:10) {
  lines(out[,1+i], col=cols[i], lty=lty[i], lwd=2)
}
legend("topright", legend=c(names(out)[2:11]), lty=c(lty,1,2), col=c(cols), lwd=2)
