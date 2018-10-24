*Body  print
	*Parameters
		# set for beginning
		sectors 1
                ##length of beginning
                floorFlag 100
		# number of band-pass stages 
		cycles 4
		# number of iterations per stage, may reboot at same stage
		iterations 8
		#number of states to track
		states 12
		#expansion ratio for Symmetry Adapted Filtering per iteration
		group 2		
		# terms in first krylov subspace
		initRank 1
		##base terms in output vectors
		basisRank 8
	.Parameters
.Body
