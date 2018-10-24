*Body print
	*Parameters
		#basic notion of quality of Quadrature
		interval 1
		## point I switch from linear to inverse metric in Quadrature
		turn 1
		## multiplier
		scalar 1
		# 0 = no interaction
		# 1 = Pseudo-coulomb
		# 2 = Yukawa, param = mass-like
		# 3 = Coulomb
		interactionOne 3
		interactionTwo 3
	.Parameters
.Body
