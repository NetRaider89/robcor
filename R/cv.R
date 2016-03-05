## Pietro Coretto, 26/02/2016: data-driven hard-thresholding  
##
## data = [patients  x  genes] 
##      = solita matrice RNA ridotta
## 
## RobCor(X) = una pseudo-funzione che denota la procedura in Python che implementa 
##             lo stimatore cormad() nel precedente scritp. 
##             RobCor(X) denota quindi la matrice di correlazione robusta ottenuta sulla 
##             matrice di dati X applicando la funzione cormad() pairwise sulle
##             colonne di X. Sarebbe opportuno salvare la matrice
##             R = RobCor(X) usando la classe "dsyMatrix" del pacchetto Matrix
## 
## Outputs che mi servono
##     R (salvata in classe "dsyMatrix" del pacchetto Matrix)
##     FLOSSES
##     CV

## Nota: se ogni istanza di RobCor prende in media 4min, la procedura 
##       dovrebbe prendere circa poco piÃ¹ di 13.5 ore.


R <- RobCor(data)

n        <- nrow(data)    ## patients
p        <- ncol(data)    ## genes
ngrid    <- 100  #
nsplits  <- 100  # 
tgrid    <- seq(1e-5, 0.99, length=ngrid)
n1       <- n - floor(n/log(n))
C1 <- C2 <- matrix(0, ncol=p, nrow=p) 
FLOSSES  <- matrix(0, nrow=nsplits, ncol=ngrid)
   
set.seed(19051977)

for (i in 1:nsplits){
   idx <- sample(1:n, size=n1, replace=FALSE)
   
   C1 <- RobCor(  data[ idx,  ]  ) 
   C2 <- RobCor(  data[-idx,  ]  )
   
   for (k in 1:ngrid){
      C1[  abs(C1) <=  tgrid[k]  ] <- 0
      FLOSSES[i,k] <- norm( C1-C2, "F")^2
   }
}















