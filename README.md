# Free motion of damped systems
        
The free equation of motion (without any external load) is given by:
$$m\ddot{x}(t) = -kx(t) - f_d(t)$$
where $f_d(t)$ is the damping force, here taken to be a generic function of time. 
In the simplified numerical simulation above, the damping force is assumed to be such that:
$$f_d(t) = c_1x(t) + c_2|\dot{x}(t)-u_r|(\dot{x}(t)-u_r)$$
where $u_r$ is external flow velocity. 

For the linear damping model ($c_2=0$), the differential equation above has the well known solution:
$$x(t) = Ae^{-\xi\omega_0 t}\cos(w\sqrt{1-\xi^2}t+\phi) = Ae^{-\xi\omega_0 t}\cos(w_dt+\phi)$$
where $\omega_0=\sqrt{\frac{k}{m}}$ is the system undamped natural frequency and $\xi$ is the damping ratio, 
defined as the ratio between the damping coefficient and its critical value as:
$$\xi=\frac{c_1}{c_{cr}}=\frac{c_1}{2m\omega_0}$$
One may immediately see that the response is given by a periodic function modulated by a negative exponential,
implying that $\xi$ can be evaluated through the response amplitude. Here, this envelope is obtained through the peak
value in every oscillation.

Although this is no longer true if $c_2\neq0$, for low damping forces, one may still make some general considerations based on energy dissipation.
Firslty, it should be noted that the energy dissipation over a full cyle may be written as:
$$W = \int_Tf_d(t)dx$$
For the linear damping contribution, one finds:
$$W_l = c_1\int_T\dot{x}(t)dx = c_1\int_T\dot{x}^2dt \approx  c_1A^2\omega_0^2\int_T\sin(\omega_0 t)^2dt = c_1\left(A^2\omega_0\pi\right)$$
where it was assumed that to first order the motion may be approximated over a cycle as $x(t)=A\cos(\omega_0t)$.
Under the same assumption, the quadratic contribution, for $u_r=0$, may be obtained as:
$$W_q = c_2\int_T|\dot{x}(t)|\dot{x}(t)dx = 2c_2\int_{T/2}\dot{x}^3dt \approx 2c_2A^3\omega_0^3\int_{T/2}\sin(\omega_0t)^3dt = c_2 \frac{8A^3\omega_0^2}{3}$$
By comparison with the linear damping result, it follows that the linear coefficient that best approximatest the quadratic response in terms of energy dissipation is:
$$\tilde{c} = c_2 \frac{8}{3\pi} A\omega_0$$
From the expression above, it can be seen that for a purely quadratic damping force, a linear dependency on the motion amplitude is expected.