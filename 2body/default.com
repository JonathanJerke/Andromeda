*Body def
	*Parameters	
		#max Von Neumann entropy of mixing states
		entropy 0.1
		#power ratio of lattice change per basis-step
		attack 0.5
		# 7 = fully periodic
		# 3 = 2D periodic
		# 1 = 1D periodic
		runFlag 0
		###of cores to apply, if available	
		lanes 1
		# 2*range+1 = Basis length in 1D
		range 4
.Parameters

	*InputOutput
	#print states
	#delete to avoid
	print 1
	.InputOutput


.Body

