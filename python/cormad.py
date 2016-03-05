import numpy as np

def cormad(x,y):
	rad2 = np.sqrt(2)
	cost = 1.4826
	
	med_x = np.median(x)
	med_y = np.median(y)
	
	mad_x =cost * np.median(np.abs(x-med_x))
	mad_y =cost * np.median(np.abs(y-med_y))
  
	zx = (x-med_x)/(cost*mad_x)
	zy = (y-med_y)/(cost*mad_y)
	
	U = zx + zy
	V = zx - zy
	mad_U2 = ( cost * np.median( np.abs( U - np.median(U) )  ) )**2
	mad_V2 = ( cost * np.median( np.abs( V - np.median(V) )  ) )**2
	
	return (mad_U2 - mad_V2)  /  (mad_U2  + mad_V2)

