Project 1: Armageddon - The hazard of small asteroids
=====================================================

Synopsis:
---------

Asteroids entering Earth’s atmosphere are subject to extreme drag forces
that decelerate, heat and disrupt the space rocks. The fate of an
asteroid is a complex function of its initial mass, speed, trajectory
angle and internal strength.

`Asteroids <https://en.wikipedia.org/wiki/Asteroid>`__ 10-100 m in
diameter can penetrate deep into Earth’s atmosphere and disrupt
catastrophically, generating an atmospheric disturbance
(`airburst <https://en.wikipedia.org/wiki/Air_burst>`__) that can cause
`damage on the ground <https://www.youtube.com/watch?v=tq02C_3FvFo>`__.
Such an event occurred over the city of
`Chelyabinsk <https://en.wikipedia.org/wiki/Chelyabinsk_meteor>`__ in
Russia, in 2013, releasing energy equivalent to about 520 `kilotons of
TNT <https://en.wikipedia.org/wiki/TNT_equivalent>`__ (1 kt TNT is
equivalent to :math:`4.184 \times 10^{12}` J), and injuring thousands of
people (`Popova et al.,
2013 <http://doi.org/10.1126/science.1242642>`__; `Brown et al.,
2013 <http://doi.org/10.1038/nature12741>`__). An even larger event
occurred over
`Tunguska <https://en.wikipedia.org/wiki/Tunguska_event>`__, a
relatively unpopulated area in Siberia, in 1908.

This simulator predicts the fate of asteroids entering Earth’s atmosphere,
and provides a hazard mapper for an impact over the UK.

Problem definition
------------------

Equations of motion for a rigid asteroid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dynamics of an asteroid in Earth’s atmosphere prior to break-up is
governed by a coupled set of ordinary differential equations:

.. math::
   :nowrap:

   \begin{aligned} 
   \frac{dv}{dt} & = \frac{-C_D\rho_a A v^2}{2 m} + g \sin \theta \\
   \frac{dm}{dt} & = \frac{-C_H\rho_a A v^3}{2 Q} \\
   \frac{d\theta}{dt} & = \frac{g\cos\theta}{v} - \frac{C_L\rho_a A v}{2 m} - \frac{v\cos\theta}{R_P + z} \\
   \frac{dz}{dt} & = -v\sin\theta \\
   \frac{dx}{dt} & = \frac{v\cos\theta}{1 + z/R_P}
   \end{aligned}

In these equations, :math:`v`, :math:`m`, and :math:`A` are the asteroid
speed (along trajectory), mass and cross-sectional area, respectively.
We will assume an initially **spherical asteroid** to convert from
inital radius to mass (and cross-sectional area). :math:`\theta` is the
meteoroid trajectory angle to the horizontal (in radians), :math:`x` is
the downrange distance of the meteoroid from its entry position,
:math:`z` is the altitude and :math:`t` is time; :math:`C_D` is the drag
coefficient, :math:`\rho_a` is the atmospheric density (a function of
altitude ), :math:`C_H` is an ablation efficiency coefficient, :math:`Q`
is the specific heat of ablation; :math:`C_L` is a lift coefficient; and
:math:`R_P` is the planetary radius. All terms use MKS units.

Asteroid break-up and deformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A commonly used criterion for the break-up of an asteroid in the
atmosphere is when the ram pressure of the air interacting with the
asteroid :math:`\rho_a v^2` first exceeds the strength of the asteroid
:math:`Y`.

.. math:: \rho_a v^2 = Y

Should break-up occur, the asteroid deforms and spreads laterally as it
continues its passage through the atmosphere. Several models for the
spreading rate have been proposed. In the simplest model, the fragmented
asteroid’s spreading rate is related to its along trajectory speed
`(Hills and Goda, 1993) <http://doi.org/10.1086/116499>`__:

.. math::  \frac{dr}{dt} = \left[\frac{7}{2}\alpha\frac{\rho_a}{\rho_m}\right]^{1/2} v

Where :math:`r` is the asteroid radius, :math:`\rho_m` is the asteroid
density (assumed constant) and :math:`\alpha` is a spreading
coefficient, often taken to be 0.3. It is conventional to define the
cross-sectional area of the expanding cloud of fragments as
:math:`A = \pi r^2` (i.e., assuming a circular cross-section), for use
in the above equations. Fragmentation and spreading **ceases** when the
ram pressure drops back below the strength of the meteoroid
:math:`\rho_a v^2 < Y`.

Airblast damage
~~~~~~~~~~~~~~~

The rapid deposition of energy in the atmosphere is analogous to an
explosion and so the environmental consequences of the airburst can be
estimated using empirical data from atmospheric explosion experiments
`(Glasstone and Dolan,
1977) <https://www.dtra.mil/Portals/61/Documents/NTPR/4-Rad_Exp_Rpts/36_The_Effects_of_Nuclear_Weapons.pdf>`__.

The main cause of damage close to the impact site is a strong (pressure)
blastwave in the air, known as the **airblast**. Empirical data suggest
that the pressure in this wave :math:`p` (in Pa) (above ambient, also
known as overpressure), as a function of explosion energy :math:`E_k`
(in kilotons of TNT equivalent), burst altitude :math:`z_b` (in m) and
horizontal range :math:`r` (in m), is given by:

.. math::
   :nowrap:

   \begin{equation*}
      p(r) = 3.14 \times 10^{11} \left(\frac{r^2 + z_b^2}{E_k^{2/3}}\right)^{-1.3} + 1.8 \times 10^{7} \left(\frac{r^2 + z_b^2}{E_k^{2/3}}\right)^{-0.565}
   \end{equation*}

For airbursts, we will take the total kinetic energy lost by the
asteroid at the burst altitude as the burst energy :math:`E_k`. For
low-altitude airbursts or cratering events, we will define :math:`E_k`
as the **larger** of the total kinetic energy lost by the asteroid at
the burst altitude or the residual kinetic energy of the asteroid when
it hits the ground.

The following threshold pressures can then be used to define different
degrees of damage.

+--------------+-------------------------------------+----------------+
| Damage Level | Description                         | Pressure (kPa) |
+==============+=====================================+================+
| 1            | ~10% glass windows shatter          | 1.0            |
+--------------+-------------------------------------+----------------+
| 2            | ~90% glass windows shatter          | 3.5            |
+--------------+-------------------------------------+----------------+
| 3            | Wood frame buildings collapse       | 27             |
+--------------+-------------------------------------+----------------+
| 4            | Multistory brick buildings collapse | 43             |
+--------------+-------------------------------------+----------------+

Table 1: Pressure thresholds (in kPa) for airblast damage

Additional sections
~~~~~~~~~~~~~~~~~~~

RK4 algorithm
~~~~~~~~~~~~~~~~~~~
The algorithm we use in solving the motion of asteroids in this project is the RK4 method, in which we approximate the value of
:math:`y_n+1` using the following set of equations:

.. math::
   :nowrap:

   \begin{aligned} 
   y_{n+1} & = y_n + \frac{1}{6} \Delta t (k_1 + 2k_2 +2k_3 +k_4) \\
   t_{n+1} & = t_n + \Delta t \\
   \\
   k_1 & = \Delta t f(t_n, y_n) \\
   k_2 & = \Delta t f(t_n + \frac{\Delta t}{2}, y_n + \frac{k_1}{2}) \\
   k_3 & = \Delta t f(t_n + \frac{\Delta t}{2}, y_n + \frac{k_2}{2}) \\
   k_4 & = \Delta t f(t_n + \Delta t, y_n + k_3)
   \end{aligned}

while function :math:`f` is the set of ODE functions mentioned above in the part of Equations of motion. We use RK4 as our algorithm
for the solver bacause it has a fourth order of convergence, making it a good choice when we change the time-stepping size.

Derivation of analytical solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given that the simplifying assumptions of exponential atmosphere (:math:`\rho_a = \rho_0 e^{-z/H}`), 
no gravitational acceleration (:math:`g=0`), a flat planet (:math:`R_p = \infty`), no lift (:math:`C_L = 0`), 
no mass change owing to ablation (:math:`\frac{dm}{dt}=0`) and no fragmentation (:math:`\frac{dr}{dt}=0`), we can simplify the equation 
of motion into the following set of ODEs:

.. math::
   :nowrap:

   \begin{aligned} 
   \frac{dv}{dt} & = \frac{-C_D\rho_0 A e^{-z/H} v^2}{2 m} \\
   \frac{dm}{dt} & = 0 \\
   \frac{d\theta}{dt} & = 0 \\
   \frac{dz}{dt} & = -v\sin\theta \\
   \frac{dx}{dt} & = v\cos\theta 
   \end{aligned}

What we are interested in here is the relationship between the velocity and the altitude, so by exchanging the varaible, we can obtain the following equation:

.. math::
   :nowrap:

   \begin{aligned} 
   \frac{dv}{dz} & = \frac{C_D\rho_0 A e^{-z/H} v}{2 m \sin\theta} 
   \end{aligned}

Integrate with respect to v in LHS and z in RHS:

.. math::
   :nowrap:

   \begin{aligned} 
   \ln(v) & = \frac{-C_D\rho_0 A H}{2 m \sin\theta} e^{-z/H} + Constant\\
   \\
   v & = Constant\, e^{\frac{-C_D\rho_0 A H}{2 m \sin\theta} e^{-z/H}}
   \end{aligned}

Plugging in the initial condition of :math:`z` (initial altitude) and :math:`v` (initial velocity), we can obtain a result for the constant 
in the above equation as :math:`m` and :math:`\theta` are fixed numbers in this senario.

Function API
============

.. automodule:: armageddon
  :members:
  :imported-members:
