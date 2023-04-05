# reference Density

We have to proceed in reverse to the data adjustment. We not fit the $n_0$, for the model.
Instead we calculate the $n_0$ based in the radius at which the disk ($R_{disk}$) and the (Ulrick) envelope coexist. 

So, we start usign the function $n_{env} = $ res.Rho0(MRate,Rdisc,MStar).
Then, $n_{disk} = A_n * n_{env}$. 
Then we calculate 
$$n_{0} = n_{disk} \left[ \frac{R_{1AU}}{R_{disk}} \right]^p ,$$ 


Now to follow the same pipeline that in single_model.py from Roberto we need to calculate $n_0 = n_{disk}(R_{min})$, where we put $R_{min} = 1$ AU (that is small that the resolution model (2.5 AU)).
Then we calculate 
$$n_{disk} = n_0 \left[ \frac{R_{disk}}{R_{1AU}} \right]^p ,$$ 
supposing that we calculate the density in the mid plane (z = 0, in equation (1) Roberto's paper 2023) 

In this way we found that $n_0 = 1.309 \times 10^{15}$
and $n_{disk} = 3.169 \times 10^{13}$.


# about the noise 

The equation for the temperature for the noise is:
$$ T_{rms} = \frac{T_{sys}}{\sqrt{\tau \Delta\nu}}$$

where $\tau$ is the time of integration, and $\Delta\nu$ is the broadband width.

$T_{sys} \propto \theta^{-2}$, then $$T_{rms} \propto \theta^{-2} (\tau \Delta\nu)^{-1/2}.$$

In the ngVLA performance estimates (2021), they report a $T_{rms} = 4779.54$ K for the 93 GHz band, with a resolution $\theta_{1/2} = 1$ mas, $\tau = $ 1 hr, and $\Delta\nu = $ 10 km/s.

We need to calculate the $T_{rms}$ but with $\theta_{1/2} = 2$ mas, $\tau = $ 10 hr, and $\Delta\nu = $ 4 km/s.

Let $T1_{rms} = $ 4779.54 K be the noise reported by nvVLA. If we look for the $T2_{rms}$:

$$ \frac {T1_{rms}}{T2_{rms}} \propto \frac{\theta_1^{-2} (\tau_1 \Delta\nu_1)^{-1/2}}{\theta_2^{-2} (\tau_2 \Delta\nu_2)^{-1/2}}$$ 

so the $$T2_{rms} = (\frac{\theta_1}{\theta_2})^2 (\frac{\tau_1 \Delta\nu_1}{\tau_2 \Delta\nu_2})^{1/2} T1_{rms}$$   

Put the numerical values in the equation, we found that $T2_{rms} =$ 944.6395 K. 

Then, using the relation between the flux and brightness temperature  (Rayleigh-Jeans regime):

$$ S_\nu = \frac{2 k \nu^{2}}{c^2} T_B $$

The flux noise are: 5.34e-5 Jy/beam
