*Body
	*Parameters
		## quality of Monte Carlo integrals
		monteCarlo 10000
		## ## of samples
		samples 10
	.Parameters

	*InputOutput
		#change basis to suite
		read ../stage.com	
		#print radial 2-body correlation
		radial out.wave
	.InputOutput
.Body
