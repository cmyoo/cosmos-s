##################################################################################
##   parameters for COSMOS_S code 
##   ver1.00 by Chulmoon Yoo  
##################################################################################
20              # maximum step of the main loop 
70.	            # maximum time to evolve
3	            # tab number of the bufer grids
10.	            # amp
100	            # maximum grid number of z
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
0	            # excision grid number (0:without excision)

##################################################################################
###  initial data parameter
##################################################################################
1	            # 0:no contnue 1:continue
ini_all.dat	    # continue file
0.8 	        # inr  parameter for the smoothing function
1.0	            # outr parameter for the smoothing function
0.83            # amplitude 
10.	            # typical wave number
0.	            # amplitude for the scalar field
10.	            # tipical wave number for the scalar field
50.0	        # Hubble

##################################################################################
###  fluid parameters
##################################################################################
0.3333333333	# fluidw
-1. 	        # Mkap
2.0	            # bminmod

##################################################################################
###  parameters for output
##################################################################################
10.	            # 1st part print interval boundary time
2.	            # 2nd part
100.	        # changing time for print interval
50              # horizon formation check interval
100             # constraint output interval