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
