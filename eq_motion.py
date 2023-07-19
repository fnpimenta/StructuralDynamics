def eq_motion(y, t , xi_l , xi_nl , w , u_ref):
	'''
	Estimates the response for a damped system under the differential equation:
	\ddot{y} = -w**2*y - 2*xi_l*w*\dot{y} - 2*xi_nl*w*abs(\dot{y}-u_ref)*(\dot{y}-u_ref).

	Parameters
	----------
	y : array_like 
		Initial conditions for the problem.
	t :array_like
		Time vector for the simulation
	xi_l : float
		Linear dimensionless damping ratio
	xi_nl : float
		Non linear dimensionless damping ratio
	w : float
		Structural natural frequency
	u_ref : float
		Environmental reference velocity
	
	Returns
	-------
	dydt : array_like
		Derivative value

	'''
	
	u,v = y

	du = v
	dv = -2*xi_l*w*v - 2*xi_nl*w*abs(v-u_ref)*(v-u_ref) - w**2*u

	dydt = [du , dv]
	
	return dydt