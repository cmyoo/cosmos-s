##################################################################################
##   parameters for COSMOS_S code 
##   ver1.00 by Chulmoon Yoo  
##################################################################################
1000	        # maximum step of the main loop 
30.	            # maximum time to evolve
3	            # tab number of the bufer grids
20.	            # amp
400	            # maximum grid number of z
1.	            # maximum coordinate of x
1.	            # maximum coordinate of z

##################################################################################
###  parameters for Evolution
##################################################################################
0.03	        # CFL number
0.1	            # cdt number
2.0	            # etaa
5.0	            # etab(eta)
0.75	        # etabb(k)
0.	            # KO dissipation
10	            # excision grid number
0	            # excision(0:false,1:true)

##################################################################################
###  initial data parameter
##################################################################################
0	            # 0:no contnue 1:continue
out_all.dat	    # continue file
0.8	            # inr  parameter for the smoothing function
1.0	            # outr parameter for the smoothing function
0.5	            # amplitude 
10.	            # tipical wave number
0.	            # amplitude for the scalar field
20.	            # tipical wave number for the scalar field
100.0	        # Hubble

##################################################################################
###  fluid parameters
##################################################################################
0.3333333333	# fluidw
-1. 	        # Mkap
2.0	            # bminmod

##################################################################################
###  parameters for output
##################################################################################
0.	            # 1st part print interval boundary time
2.	            # 2nd part
100.	        # changing time for print interval
