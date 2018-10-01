*Body
	*Parameters
		monteCarlo 10000
		samples 10
	.Parameters

	*InputOutput
		read ../stage.com
		vector out.wave
	.InputOutput
.Body
