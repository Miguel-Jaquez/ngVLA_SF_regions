# Hamburgers_piecewise disk

## Density
We have to proceed in the opposite way than in Galván-Madrid et al 2023. Here we are not going to adjust observations.

So, we need to calculate the $n_0$ based in the radius at which the disk ($R_{disk} = R_d$) and the (Ulrick) envelope coexist. 

We start usign the sf3models function $\rho_{env} = $ res.Rho0(MRate,Rdisc,MStar).
Using **MRate** = $3.33\times10^{-5}$ $M_\odot$ $yr^{-1}$ (Galvan-Madrid 2023) as the accretion rate, Rd = 250 AU is the centrifugal radius and MStar = 10 $M_\odot$ the mass of the central object.

NOTE: I don't calculated the real MRate and Rd, we hare scaling the Galván-Madrid et al. 2023 values.

The function, res.Rho0(), is the equation 11 from Keto 2007:

$\dot M = \mu n_\rho 4\pi R_d^2 V_k$.

So the Ulrich envelope density at centrifugal radius is $5.99\times10^{12} ~ m^{-3}$.

Then using the disk density scaling factor $A_\rho$, we can pass from the envelope density to the disk density at Rd using  $\rho_{disk} = A_\rho n_{\rho}$. We use the value from Galván-Madrid et al 2023 $A_\rho = 5.29$.

Roberto say that the radius of reference $R_0$ will be 1 AU. 
So, using equation 1 from Galván-Madrid et al 2023, we calculate the density $\rho_{0}$ in the plane of the disk, using the density at the centrifugal radius as reference:

$\rho_{0} = \rho_{disk} \left[ \frac{R_{1AU}}{R_{disk}} \right]^p = 3.17\times^{13} m^{-3}$. 

Now we can follow the same pipeline that in single_model.py.

To create the density, we use the function **Model.density_hamburger_piecewise(RStar, H0, R_list, P_list, rho0,Grid**, H_list = None, q_list = [0.5], rho_thres = 10.0, rho_min = 0.0, Rt = False).

Where Rstar is the stellar radius, H0 is the scaleheight normalization constant (usually H0 = shFactor * R0), R_list is a list with the polar limits, P_list list of powerlaws in R_list intervals and rho0 is the density at R_list[0].


For our model the values are:
- **$R_{Star}$** $= (M_{Star}/M_{Sun})^{0.8}$ [$R_\odot$]

- **H0** $= 1.496e11 \times sh_{factor}$ where sh_factor =0.2  Note: this value i dont know where comes from, in single_model.py I found this note: #H0 = 1.496e12 * H0_factor #m at Rmin = 10au, such that H(R=2000au) = 2000 au, valid for q = 1.0.
 
 
- **R_list** = [Rmin, Rdisc] , #Interval limits for piecewise density profile. Rmin = 1 AU and Rdisc is the centrifugal radius.
 
 
- **P_list** = [-0.674], is the: Powerlaws within intervals
 
 
- **rho0** is the density at Rmin. The "pivot" values. From create_full_model.py we found that:  rho0 = rho0	*sh_{factor}**p_list_0 	# in m^-3    #paper Izquierdo 2018.


## Temperature
We define a **constant** temperature for the model at $9\times10^3$ K. For this we use the next function:

$t_e = 9.0\times 10^3$ K


temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, backTemp=2.725480)


## Velocity

For the velocity we use the function :

vel = Model.velocity_piecewise(density, GRID,

                               R_list=R_list_rot, pR_list=pR_list, v0R=norm_vR, #polar piecewise velocity
                               
                               r_list=r_list, pr_list=pr_list, v0r=norm_vr,     #radial piecewise velocity
                               
                               )
                               
using the values :


**Rrot_pivot** = 125 * u.au # copy from the create_model_full.py

**R_list_rot** = [Rmin, Rrot_pivot, Rdisc] #Interval limits for piecewise rotation profile

**pR_list** = [-0.5, -0.5] #powerlaw exponents

**vR0** = (ct.G * MStar / R_list_rot[0])**0.5	# Circular velocity at R_list_rot[0]

**norm_vR** = [0,0,vR0] #Components of rotation velocity. Only index [2] is taken.

**rrad_int** = rrad_int * u.au

**r_list** = [rrad_int, Rdisc]

**pr_list** = [pr_list_0] #powerlaw exponents

**vr0** = vr0  = 0 # radial velocity           #m/s

**norm_vr** = [vr0,0,0]   #Components of infall. Only [0] is taken. sign determines direction, '-' for infall


# Noise 

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

# other models - density_Env_Disc function

## Only Disk 

For this model we use the function:

**density_Env_Disc**(RStar, Rd, rhoE0, Arho, GRID, 
                     exp_disc=2.25,
                     discFlag = True, envFlag = False, 
                     rdisc_max = False, renv_max = False, 
                     ang_cavity = False, average_around_Rd = None, 
                     rho_min_env = 1.0e9):

where

- **$M_{star}$** $ = 10 M_\odot$

- **MRate** = $3.33\times10^{-5}$ $M_\odot$ $yr^{-1}$ 

- **Rd** = 250 AU

- **rhoE0** = Res.Rho0(MRate, Rd, MStar) #Base density for Ulrich's model

- **Arho** = 24.1 #Disc-envelope density factor

In this way the disk density is consistently calculated. 

## Envelope plus Disk

Also, we create models with the envelope flag activated. We put Renv = 500 au, and the angle_cavity = 40 degrees (I think)

## Temperature

We use the constant temperature for both models:

$t_e = 9.0\times 10^3$ K


temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, backTemp=2.725480)

## Velocity

For this models we use the Keplerian velocity function:

**velocity**(RStar,MStar,Rd,density,GRID)

# Other model - Hamburgers function

Using the same values that the previus models, we use the hamburgers function to try to reproduce the results from the hamburgers_piecewise model.

To do this, we use the function:

**Model.density_hamburger(Rstar, H0_factor, Rdisk, rhoE0, A_rho,GRID)**.


Next, we define a **constant** temperature for the model at $9\times10^3$ K.


Finally, a keplerian velocity is defined using the function: **vel = Model.velocity(RStar, MStar, Rd, density, GRID)**.

