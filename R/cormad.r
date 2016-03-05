## cormad: calcola una stima robusta della correlazione tra x ed y
cormad <- function(x,y) {

    rad2     <- sqrt(2)
    cost     <- 1.4826

    med_x    <- median(x)
    med_y    <- median(y)

    mad_x <- cost * median( abs(  x - med_x ) )
    mad_y <- cost * median( abs(  y - med_y ) )

    zx  <-  { x - med_x} / { cost * mad_x }
    zy  <-  { y - med_y} / { cost * mad_y }

    U   <-  zx + zy
    V   <-  zx - zy

    mad_U2 <- { cost * median( abs( U - median(U) )  ) }^2
    mad_V2 <- { cost * median( abs( V - median(V) )  ) }^2

    return( { mad_U2 - mad_V2 }  /  {mad_U2  + mad_V2} )

}
