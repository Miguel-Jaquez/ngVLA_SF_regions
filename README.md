# Hamburgers_piecewise disk

## Density

We have to proceed in the opposite way than in Galván-Madrid et al 2023. Here we are not going to adjust observations.

Instead we calculate the $$n_0$$ based in the radius at which the disk ($$R_{disk}$$) and the (Ulrick) envelope coexist. 

The envelope base number density is defined through the acretion rate  (equation 11 from Keto (2007)) 

$$n_0=\frac{\dot{M}}{2 \mu 4 \pi R_d^2 (G M_\star/R_d)^{1/2}}$$ in units of $$m^{-3}$$.

here, $$R_d = GM/v_k^2$$ is the radius at which the gravitational force equals the centrifugal force in the midplane of the disc. 

We start usign the sf3dmodels function $$n_{env}=$$ res.Rho0(MRate,Rdisc,MStar) to calculate the density at $$R_d$$. 

We use **MRate** = $3.33\times10^{-5}$ $M_\odot$ $yr^{-1}$ (Galvan-Madrid 2023) as the accretion rate and MStar = 10 $M_\odot$ the mass of the central object. 

Then using the disk density scaling factor $A_\rho$, we can pass from the envelope density to the disk density at $$R_d$$ using  $$n_{disk}  = A_n * n_{env}$$. We use the value from Galván-Madrid et al 2023 $A_\rho = 5.29$.

The disk density is defined by (Equation 1 Galván-Madrid et al. 2023)

$$n_e(R,z)=n_0 \left[\frac{R}{R_0} \right]^p exp[-z^2/2H^2]$$.

We want to define the density at $$R=R_{min}=1 {\rm AU}$$.

Taking $$R_0 = R_{disk}$$ and $$n_0=n_{disk}$$, we solve the density $$n_{\rm 1AU}$$ at $$R = 1 {\rm AU}$$

$$n_{\rm 1AU} = n_{disk} \left[ \frac{R_{min}}{R_{disk}} \right]^p.$$

Finally, $$n_{\rm 1AU}$$ and $$R_{\rm min}$$ will be our new $$n_0$$ and $$R_0$$ for our hamburgers_piecewise models.

In this way, we found that $$n_0 = 2.79 \times 10^{15}$$ m^-3.

Now we can follow the same pipeline that in single_model.py.

To create the density, we use the function **Model.density_hamburger_piecewise(RStar, H0, R_list, P_list, rho0,Grid**, H_list = None, q_list = [0.5], rho_thres = 10.0, rho_min = 0.0, Rt = False).

Where Rstar is the stellar radius, H0 is the scaleheight normalization constant (usually H0 = shFactor * R0), R_list is a list with the polar limits, P_list list of powerlaws in R_list intervals and rho0 is the density at R_list[0].


For our model the values are:
- **$$R_{Star}$$** $$= (M_{Star}/M_{Sun})^{0.8}$$ [$$R_\odot$$]

- **H0** $$= 1.496e11 \times sh_{factor}$$ where $$sh_{factor}=0.2$$  Note: this value i dont know where comes from, in single_model.py I found this note: #H0 = 1.496e12 * H0_factor #m at Rmin = 10au, such that H(R=2000au) = 2000 au, valid for q = 1.0.
 
- **R_list** = [Rmin, Rdisc] , #Interval limits for piecewise density profile. Rmin = 1 AU and Rdisc is the centrifugal radius.
 
- **P_list** = [-0.674], is the: Powerlaws within intervals
 
- **rho0** is the density at Rmin, $$n_0$$. The "pivot" values. Note: From create_full_model.py we found that:  rho0 = rho0	*sh_{factor}**p_list_0 	# in m^-3    #paper Izquierdo 2018. but i dont found the physical motivation to do this.

## Temperature
We define a **constant** temperature for the model at $9\times10^3$ K. For this we use the next function:

$$t_e = 9.0\times 10^3$$ K


temperature = Model.temperature_Constant(density, GRID, discTemp=t_e, backTemp=2.725480)


## Velocity

For the velocity we use the function :

vel = **Model.velocity_piecewise**(density, GRID,
R_list=R_list_rot, pR_list=pR_list, v0R=norm_vR, #polar piecewise velocity
r_list=r_list, pr_list=pr_list, v0r=norm_vr,     #radial piecewise velocity)

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

# The noise 

#### From the CASA tutorials for the ngVLA:

The RMSnoise ($$\sigma_{NA}$$) in an untapered, naturally-weighted Stokes I image will be approximately

$$\sigma_{NA} \sim  \frac{\sigma_{simple}} {\sqrt{n_{ch} n_{pol} n_{baselines} n_{integrations}}}$$

where $$\sigma_{simple}$$ is the simplenoise that correspond to the noise per visibility, $$n_{ch}$$ is the total number of channels acriss akk soectrak windows, $$n_{pol}$$ is the number of polarizations used for the Stokes I (typically 2) and $$n_{integrations}$$ is the number of correlator integration times in the MS (i.e., total on-source time / integration time). For the example simulations, the track time is 4 hrs and the integration time is 60 sec, thus $$n_{integrations} = 240$$. Additionally, for these examples the total number of channels is 1 and the number of polarizations is 2. The number of baselines is N(N-1) where N is the numbers of antennas in the array. 

Sacar para la configuracion main del ngVLA.

#### For the ngVLA

We suggest the following procedure:
1)  Calculate the expected untapered, naturally weighted point source sensitivity () using one of the ngVLA performance tables. In Appendix D of ngVLA memo #55 there are key performance metrics for 6 subarrays which are tabulated as a function of frequency and resolution. For our example, we find in Table 10 of ngVLA memo #55 that the untapered, naturally weighted point source sensitivity of the Main interferometric array at 93 GHz is 0.83 uJy/beam for a 1 hour observation. 

2) Scale that number to the desired observation length, in this case $$t_{track}=4 h$$. Therefore,
$$\sigma_{NA} = 0.83 / \sqrt{t_{track}/1 hour} = 0.415 \mu Jy / beam$$.

3) Use the expected image noise ($$\sigma_{NA}$$) in the above equation (1) to solve for the scaling factor $$\sigma_{simple}$$. In this case,
$$\sigma_{simple} = 0.415 \sqrt{ 1 * 2 * 22791 * 240} = 1.4$$ mJy


# check if the description below is ok
The equation for the temperature for the noise is:
$$T_{rms} = \frac{T_{sys}}{\sqrt{\tau \Delta\nu}}$$

where $\tau$ is the time of integration(?), and $\Delta\nu$ is the broadband width.

$$T_{sys} \propto \theta^{-2}$$, then $$T_{rms} \propto \theta^{-2} (\tau \Delta\nu)^{-1/2}$$.

In the ngVLA performance estimates (2021), they report a $$T_{rms} = 4779.54$$ K for the 93 GHz band, with a resolution $$\theta_{1/2} = 1$$ mas, $$\tau =$$ 1 hr, and $$\Delta_\nu=$$ 10 km/s.

We need to calculate the $$T_{rms}$$ but with $$\theta_{1/2} = 2$$ mas, $$\tau =$$ 10 hr, and $$\Delta\nu =$$ 5 km/s.

Let $$T1_{rms} =$$ 4779.54 K be the noise reported by nvVLA. If we look for the $$T2_{rms}$$:

$$\frac {T1_{rms}}{T2_{rms}} \propto \frac{\theta_1^{-2} (\tau_1 \Delta\nu_1)^{-1/2}}{\theta_2^{-2} (\tau_2 \Delta\nu_2)^{-1/2}}$$ 

so the $$T2_{rms} = (\frac{\theta_1}{\theta_2})^2 (\frac{\tau_1 \Delta\nu_1}{\tau_2 \Delta\nu_2})^{1/2} T1_{rms}$$   

Put the values we found that $$T2_{rms} =$$ 944.6395 K. 

Then, using the relation between the flux and brightness temperature  (Rayleigh-Jeans regime):

$$S_\nu = \frac{2 k \nu^{2}}{c^2} T_B$$

The flux noise are:








# other models - density_Env_Disc function in construction

## Only Disk 

For this model we use the function:

**density_Env_Disc**(RStar, Rd, rhoE0, Arho, GRID, 
                     exp_disc=2.25,
                     discFlag = True, envFlag = False, 
                     rdisc_max = False, renv_max = False, 
                     ang_cavity = False, average_around_Rd = None, 
                     rho_min_env = 1.0e9):

where

- **$$M_{star}$$** $$ = 10 M_\odot$$

- **MRate** = $$3.33\times10^{-5}$$ $$M_\odot$$ $$yr^{-1}$$ 

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

# Other model - Hamburgers function (same that Env-disc function)

Using the same values that the previus models, we use the hamburgers function to try to reproduce the results from the hamburgers_piecewise model.

To do this, we use the function:

**Model.density_hamburger(Rstar, H0_factor, Rdisk, rhoE0, A_rho,GRID)**.


Next, we define a **constant** temperature for the model at $9\times10^3$ K.


Finally, a keplerian velocity is defined using the function: **vel = Model.velocity(RStar, MStar, Rd, density, GRID)**.


