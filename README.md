# A note about high-order semi-implicit differentiation: application to a numerical integration scheme with Taylor-based compensated error

This code is associated to the arXiv paper <https://arxiv.org/abs/2408.00497> and illustrates some numerical results of a third order semi-implicit differentiator, based on the implicit discretization technique [1] and last developments [2] for which a Taylor refinement algorithm allows reducing the estimation error throughout the differentiation estimation.

Below is given a very short and simplified introduction, please refer to the paper for a complete description.

## Brief presentation of the third order differentiator

Considering a real signal $x_1(t)$ to differentiate, the third order semi-implicit differentiator (SIHD3), including Taylor expansion refinement, is illustrated as:

$$\left \lbrace \begin{array}{l}
{z}_4^+= z_4+  h \left( \lambda_4 { \mu^4} \left|{e_1}\right|^{4\,\alpha-3} \mathcal{N}_4  \right)  \\
{z}_3^+= z_3+  h \left(z_4^++\lambda_3 { \mu^3} \left|{e_1}\right|^{3\,\alpha-2} \mathcal{N}_3 \right)  \\
{z_2^+}= z_2+  h \left(z_3^+ {- h \frac{1}{2} z_4^+ } + \lambda_2 { \mu^2} \left|{e_1}\right|^{2\,\alpha-1} \mathcal{N}_2  \right) \\
z_1^+= z_1+ h \left( z_2^+ {- \frac{h}{2} z_3^+ +  \frac{h^2}{3!} z_4^+ } +  \lambda_1 {\mu} |{e_1}|^{\alpha} \mathcal{N}_1 \right)  \\
\end{array} \right.$$

or, equivalently, introducing the SIHD related operator:
 
$$
 z_i = \mathfrak{D}_f ^{(i)} (z, z^+, x_1), \quad i = 1..4
$$
 
where $z_i$ is the estimation associated to the $i$ th line of the differentiator.

Denote $y = x_1$ the signal to differentiate; $\dot{y} = {x_2}$,  $\ddot{y} = x_3$ and $y^{(4)} = x_4$ are respectively the first, second, and third order differentiation;  $z_2, z_3$ and $z_4$ are respectively the first, second and third order related estimations and $e_1 = x_1 - z_1$ is the estimation error.

The projector $\mathcal{N}_q$ for $q = 1..4$ behaves either as a sign function, or as a nonlinear "slope" depending on the value of the error $e_1 = x_1 - z_1$. It is defined by:

$$
\mathcal{N}_q (\epsilon_1):= \left \lbrace
\begin{array}{l}
\epsilon_1\in SD  \rightarrow \displaystyle{  \mathcal{N}_q = {\frac{\lceil \epsilon_1\rfloor^{q(1-\alpha)}}{ \lambda_q (\mu h)^q} } } \\    
\epsilon_1 \notin SD  \rightarrow  \mathcal{N}_q = \mathrm{sign}(\epsilon_1) \\
\end{array} \right.
$$

with the convergence domain $SD=\{\epsilon_1  /  |\epsilon_1|\leq (\lambda_1 \mu h)^{\frac{1}{q(1-\alpha)}}\}$.

Figures 1 and 2 illustrate respectively the evolution of the estimated $z_1, z_2, z_3$ and $z_4$ with respect to the corresponding exact differentiation orders.


<img width="2161" alt="SIHD3_no_1" src="https://github.com/user-attachments/assets/5ddc4eee-1fe0-4530-a1a0-3be9178d5e73">
<em>Fig. 1 - Example of SIHD3-based differentiation of a sine function.</em>

<img width="2246" alt="SIHD3_no_2" src="https://github.com/user-attachments/assets/8353db69-8fba-443f-a078-9d47c5f05bab">
<em>Fig. 2 - Example of SIHD3-based differentiation of a sine function - estimation errors.</em>

$$ \, $$

The error is of order $< e_1 > = 2.96 \, 10^{-12}$  on the $z_1$ estimation, of order $< e_2 > = 4.69 \, 10^{-9}$ on the estimation of the first derivative $z_2$ and of order $< e_3 > = 4.49 \, 10^{-6}$ on the estimation of the second derivative $z_3$.

## Application to an observer-based numerical integration scheme

In the context of solving ODE, an observer-based numerical integration scheme is proposed to provide a better precision of the estimation of the ODE solution. 
Our proposed differentiation method can be associated to the nonstandard finite-difference (NSFD) methods [3] [4], described as "powerful numerical methods that preserve significant properties of exact solutions of the corresponding differential equations".

An NSFD scheme including the third order differentiator can be written:

$$
\left\{\begin{array}{l}
z_1^+  = z_1 + \psi \, \mathfrak{D}_{f}^{(1)}(z, z^+, \hat{y} ) \\
\hat{y} = z_1
\end{array}, \right. \, k \in \mathbb{N}, \, y(0) = y_0 \quad \mathrm{and} \quad \psi = h + O(h^2)
$$

where the SIHD3 differentiator is used to estimate the $f$ part of the ODE and $\psi = h + O(h^2)$ that do satisfy $h \rightarrow 0$ according to the NSFD rules.

Consider the basic first order equation:

$$
\frac{\mathrm{d}  y(t)}{\mathrm{d}  t} = -y(t) + 1, \quad y(0) = 10^{-2}
$$

whose exact solution is $y_{exact} = 1 - \exp( - t )$. Fig. 3 illustrates the comparaison of the errors between the Euler forward scheme and the Runge-Kutta scheme.


<img width="2210" alt="Result_NFSD_1" src="https://github.com/user-attachments/assets/fe4c5d2e-fcac-42c1-85e7-88bc450300a9">
<em>Fig. 3 - Illustration of the SIHD3 scheme convergence with respect to Euler forward scheme and Runge-Kutta scheme.</em>

$$ \, $$

The proposed numerical scheme provides a slightly reduced error compared with the Euler and Runge-Kutta schemes, especially when the solution is converged towards the asymptotic solution.

## References

[1] V. Acary, B. Brogliato, and Y. Orlov. Chattering-free digital sliding-mode control with state observer and disturbance rejection. IEEE Trans. on Automatic Control, 57(5):1087–1101, 2012.


[2] Loïc Michel, Malek Ghanes, Yannick Aoustin, and Jean-Pierre Barbot. An interconnected discrete time cascaded semi-implicit differentiation. 17th International Workshop on Variable Structure Systems (VSS2024), Accepted, April 2024. <https://hal.science/hal-04564290/>

[3] Ronald E Mickens. Nonstandard Finite Difference Models of Differential Equations. WORLD SCIENTIFIC, 1993.

[4] Ronald E Mickens. Advances in the Applications of Nonstandard Finite Difference Schemes. WORLD SCIENTIFIC, 2005.


## Credits and license

(c) [2024]  Nantes Université - Centrale Nantes - LS2N UMR 6004, Nantes

(c) [2024]  Quartz EA 7393, ENSEA, Cergy-Pontoise
      
Loïc MICHEL and Jean-Pierre BARBOT

All rights reserved under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.



