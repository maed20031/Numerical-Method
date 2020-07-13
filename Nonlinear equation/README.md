Mathematical Theory

- Continuous function
- Intermediate Value Problem
- Convergence of a sequence



Fixed point

- Easy to implement
- absolute value of derivative should be less than 1

Bisection 

- Guarantees that the root always exists if the initial condition is satisfied

Newton method 

- Faster than bisection
- Need derivative (and second derivative)
- can fail to converge

Brent

- Faster than bisection
- Guarantees that the root always exists if the initial condition is satisfied





Stochastic ways

Robbins-Monro (from PRML 2.3.5)

Let $\theta, z$ be a pair of random variables and joint distribution $p(z,\theta)$.

Define $f(\theta)$ which denotes conditional expectation of $z$ given $\theta$,
$$
f(\theta) \equiv E[z|\theta] = \int zp(z|\theta)dz
$$
Objective : find $\theta$ which satisfies that $f(\theta) = 0$.



By iterations, 
$$
\theta_N = \theta_{N-1} + a_{N-1}z(\theta_{N-1})
$$



Assumptions

- $E[(z - f)^2] < \infty$

- $\lim_{n\to\infty} a_n = 0$

- $\sum_{n=1}^\infty a_n = \infty$

- $\sum_{n=1}^\infty a_n^2 < \infty$

  

First assumption : the sequence converges to the root with prob 1.

Second assumption : the sequence converges

Third assumption : make sure that the sequence does not converge before finding root

Last assumption : make sure that noise has only finite variance.


