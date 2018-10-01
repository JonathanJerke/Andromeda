*Body print
	*Parameters
	type 3
	.Parameters

	*Geometry
	3 0 0 0 
	.Geometry
	
	*InputOutput
	#read ../print.com
	#read ../run.com
	.InputOutput
.Body
