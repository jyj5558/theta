#R function to assign a CONservation STATus based the genetic criterion
#based on het (estimated genome-wide heterozygosity), theta (estimated Watterson's theta), mu (known mutation rate), interval (known generation interval),
#which will output $Ne (effective population size calculated from 'theta'), $HetT (het after 100 years), $ConStat (Conservation Status from the genetic criterion)

# original - using genome-wide Theta
#gen_stat <- function(het,theta,mu,interval){
#    Ne = theta/(4*mu)
#    g = 100/interval
#    HetT = het*(1-(1/(2*Ne)))^g
#    if (HetT/het <= 0.9){
#        ConStat = "CR"
#    } else if (HetT/het <= 0.95){
#        ConStat = "EN"
#    } else if (HetT/het <= 0.975 || Ne < 1000 ){
#        ConStat = "VU"
#    } else if (HetT/het > 0.975 && 1000 <= Ne && Ne < 5000 ){
#        ConStat = "NT"
#    } else if (HetT/het > 0.975 && 5000 <= Ne){
#        ConStat = "LC"
#    }
#    result <- list("Ne" = Ne, "HetT" = HetT, "ConStat" = ConStat)
#    return(result)
#} 

# version 2
gen_stat <- function(het,Ne,Nc,interval){
    if (missing(Ne) && !missing(Nc)){
        Ne=0.14*Nc
    } else if (missing(Nc) && !missing(Ne)){
        Ne=Ne
    } else if (missing(Nc) && missing(Ne)){
        stop('input either one of Nc or Ne')
    } else if (!missing(Nc) && !missing(Ne)){
        stop('input either one of Nc or Ne')
    }
    g = 100/interval
    HetT = het*(1-(1/(2*Ne)))^g
    if (HetT/het <= 0.9){
        ConStat = "CR"
    } else if (HetT/het <= 0.95){
        ConStat = "EN"
    } else if (HetT/het <= 0.975 || Ne < 1000 ){
        ConStat = "VU"
    } else if (HetT/het > 0.975 && 1000 <= Ne && Ne < 5000 ){
        ConStat = "NT"
    } else if (HetT/het > 0.975 && 5000 <= Ne){
        ConStat = "LC"
    }
    result <- list("Ne" = Ne, "HetT" = HetT, "ConStat" = ConStat)
    return(result)
} 
