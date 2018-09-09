---
title: Simplifications behind thermodynamics
---

If one wants to be rigorous in the way he writes thermodynamics, the derivations become very verbose. Therefore, some implicit hypothesis on the writting style of many textbooks are done. In order to make them explicit, I give here two derivations of the activity coefficient of a pure compound. The first derivation should be fully rigorous, while the second will be done in the thermodynamic way of style.

# Rigorous derivation

The fugacity coefficient of a pure compound \(i\) is defined as :

$$
RT \ln \phi_i\left(T, P\right) = g_i^{res}\left(T,P\right) = g_i^{\phi}\left(T^{\phi} = T,P^{\phi} = P\right) - g_i^{\bullet}\left(T^{\bullet} = T,P^{\bullet} = P\right)
$$

Because pressure explicit EoS have temperature and pressure as natural variables, it is convenient to express the different terms from previous equation using such variables.

$$
RT \ln \phi_i\left(T, P\right) = a_i^{\phi}\left(T^{\phi} = T, P^{\phi} = P\right) + \underbrace{\left(Pv\right)^{\phi}}_{Pv^{\phi}\left(T,P\right)} - a_i^{\bullet}\left(T^{\bullet} = T, P^{\phi} = P\right) - \underbrace{\left(Pv\right)^{\bullet}}_{RT}
$$

Because for a real fluid, a property value is independent of the variables used to express it, one can write: \(a_i\left(T^{\phi} = T,P^{\phi} = P\right) = a_i\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right)\). Moreover, one can write:

$$
a_i^{\bullet}\left(T^{\bullet} = T, P^{\bullet} = P\right) = a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} = \frac{RT}{P}\right)
$$

Finally, one can write:

$$
RT \ln \phi_i\left(T, P\right) = a_i^{\phi}\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right) - a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} = \frac{RT}{P}\right)+ Pv^{\phi}\left(T,P\right) - RT
$$

At the infinite molar volume limit, the real fluid and the perfect gas are identical. Thus, one can write:

$$
a_i^{\phi}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right) = a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} \rightarrow +\infty\right)
$$

Then, \(RT \ln \phi\) expression rewrites:

$$
RT \ln \phi_i\left(T, P\right) = a_i^{\phi}\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right) - a_i^{\phi}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right) - \left[a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} = \frac{RT}{P}\right) - a_i^{\bullet}\left(T^{\phi} = T,v^{\bullet} \rightarrow +\infty\right)\right]+ Pv^{\bullet}\left(T,P\right) - RT
$$

The variations of state function \(a\) in variables \(T\) and \(v\) writes:

$$
da = -sdT - Pdv
$$

Applied to the equation of state and the ideal gas, this equation rewrites:

$$
\begin{cases}
da^{\phi}\left(T^{\phi}, v^{\phi}\right) = s^{\phi}\left(T^{\phi}, v^{\phi}\right)dT^{\phi} - P^{\phi}\left(T^{\phi}, v^{\phi}\right)dv^{\phi} \\
da^{\bullet}\left(T^{\bullet}, v^{\bullet}\right) = s^{\bullet}\left(T^{\bullet}, v^{\bullet}\right)dT^{\bullet} - P^{\bullet}\left(T^{\bullet}, v^{\bullet}\right)dv^{\bullet}
\end{cases}
$$

Therefore:

$$
 a_i^{\phi}\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right) - a_i^{\phi}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right) = -\int_{+\infty}^{v^{\phi}\left(T,P\right)}P^{\phi}\left(T^{\phi} = T, v^{\phi}\right) dv^{\phi}
$$

For the perfect gas:

$$
 a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet}\left(T, P\right)\right) - a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} \rightarrow +\infty\right) = -\int_{+\infty}^{v^{\bullet}\left(T,P\right)}P^{\bullet}\left(T^{\bullet} = T, v^{\bullet}\right) dv^{\bullet} = -\int_{+\infty}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet}
$$

The previous integral can be splitted into two terms:

$$
 a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet}\left(T, P\right)\right) - a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} \rightarrow +\infty\right) = -\int_{+\infty}^{v^{\phi}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} -\int_{v^{\phi}\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet}
$$

Because the integrand is a dummy variable, one can write:

$$
RT \ln \phi_i\left(T, P\right) = -\int_{+\infty}^{v^{\phi}\left(T,P\right)}\left[P^{\phi}\left(T^{\phi} = T, v\right) - \frac{RT}{v} \right] dv + \int_{v^{\phi}\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} + Pv^{\phi}\left(T, P\right) - RT
$$

Moreover:

$$
\int_{v^{\phi}\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} = RT \ln\left[\frac{v^{\bullet}\left(T,P\right)}{v^{\phi}\left(T,P\right)}\right] = RT \ln\left[\frac{RT}{P\times v^{\phi}\left(T,P\right)}\right] = - RT \ln z^{\phi}\left(T,P\right)
$$

Finally:

$$
\ln \phi_i\left(T, P\right) = -\frac{1}{RT}\int_{+\infty}^{v^{\phi}\left(T,P\right)}\left[P^{\phi}\left(T^{\phi} = T, v\right) - \frac{RT}{v} \right] dv - \ln z^{\phi}\left(T,P\right) + z^{\phi}\left(T,P\right) - 1
$$

# What hypothsesis shall we do?

From the previous derivation, the following observation can be done:

* the \(\phi\) exponent can be dropped for the real fluid state functions (\(a\), \(g\)...). From now on, a state function without exponent shall be considered as calculated using an equation of state. The \(\bullet\) exponent shall be maintained as it allows us to specify that the perfect gas state is used.

* the \(\phi\) and \(\bullet\) exponents shall be dropped for variables that are fixed by initial hypothesis to avoid writting \(T^{\phi} = T\) or \(P^{\phi} = P\). On the other hand, one should not drop the \(\phi\) and \(\bullet\) exponent for variables that are not fixed. Here, \(v^{\phi}\) and \(v^{\bullet}\) are distinct functions of temperarature and pressure and the exponents are used to recall this difference.

* the specified variables used to calculate any state function not fixed by the user shall always show the variables between parenthesis. For exemple, if the molar volume is not fixed by the user, the no one shall write \(v\) but \(v\left(T, P\right)\) because in the first case, one might possibly believe that the molar volume is fixed by the user (since we dropped the \(\phi\) exponent for the real fluid properties).

Doing so, it is possible to reduce the number of symbols in the thermodynamic derivations without adding any ambiguity.

# Second derivation

Pressure and temperature are fixed in the following derivation. The fugacity coefficient of a pure compound \(i\) is defined as :

$$
RT \ln \phi_i\left(T, P\right) = g_i^{res}\left(T,P\right) = g_i\left(T,P\right) - g_i^{\bullet}\left(T, P\right)
$$

Because pressure explicit EoS have temperature and pressure as natural variables, it is convenient to express the different terms from previous equation using such variables.

$$
RT \ln \phi_i\left(T, P\right) = a_i\left(T, P\right) + Pv\left(T,P\right) - a_i^{\bullet}\left(T, P\right) - \underbrace{\left(Pv\right)^{\bullet}}_{RT}
$$

Because for a real fluid, a property value is independent of the variables used to express it, one can write: 

$$
a_i\left(T, P\right) = a_i\left(T,v\left(T, P\right)\right)
$$

Moreover, one can write:

$$
a_i^{\bullet}\left(T, P\right) = a_i^{\bullet}\left(T,v^{\bullet} = \frac{RT}{P}\right)
$$

Finally, one can write:

$$
RT \ln \phi_i\left(T, P\right) = a_i\left(T,v\left(T, P\right)\right) - a_i^{\bullet}\left(T,v^{\bullet} = \frac{RT}{P}\right)+ Pv\left(T,P\right) - RT
$$

At the infinite molar volume limit, the real fluid and the perfect gas are identical. Thus, one can write:

$$
a_i\left(T,v \rightarrow +\infty\right) = a_i^{\bullet}\left(T,v^{\bullet} \rightarrow +\infty\right)
$$

Then, \(RT \ln \phi\) expression rewrites:

$$
RT \ln \phi_i\left(T, P\right) = a_i\left(T,v\left(T, P\right)\right) - a_i\left(T,v \rightarrow +\infty\right) - \left[a_i^{\bullet}\left(T,v^{\bullet} = \frac{RT}{P}\right) - a_i^{\bullet}\left(T,v^{\bullet} \rightarrow +\infty\right)\right]+ Pv\left(T,P\right) - RT
$$

The variations of state function \(a\) in variables \(T\) and \(v\) writes:

$$
da = -sdT - Pdv
$$

Applied to the equation of state and the ideal gas, this equation rewrites:

$$
\begin{cases}
da\left(T, v\right) = s\left(T, v\right)dT - P\left(T, v\right)dv \\
da^{\bullet}\left(T^{\bullet}, v^{\bullet}\right) = s^{\bullet}\left(T^{\bullet}, v^{\bullet}\right)dT^{\bullet} - P^{\bullet}\left(T^{\bullet}, v^{\bullet}\right)dv^{\bullet}
\end{cases}
$$

Therefore:

$$
 a_i\left(T,v\left(T, P\right)\right) - a_i\left(T,v \rightarrow +\infty\right) = -\int_{+\infty}^{v\left(T,P\right)}P\left(T, v\right) dv
$$

For the perfect gas:

$$
 a_i^{\bullet}\left(T,v^{\bullet}\left(T, P\right)\right) - a_i^{\bullet}\left(T,v^{\bullet} \rightarrow +\infty\right) = -\int_{+\infty}^{v^{\bullet}\left(T,P\right)}P^{\bullet}\left(T, v^{\bullet}\right) dv^{\bullet} = -\int_{+\infty}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet}
$$

The previous integral can be splitted into two terms:

$$
 a_i^{\bullet}\left(T,v^{\bullet}\left(T, P\right)\right) - a_i^{\bullet}\left(T,v^{\bullet} \rightarrow +\infty\right) = -\int_{+\infty}^{v\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} -\int_{v\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet}
$$

Because the integrand is a dummy variable, one can write:

$$
RT \ln \phi_i\left(T, P\right) = -\int_{+\infty}^{v\left(T,P\right)}\left[P\left(T, v\right) - \frac{RT}{v} \right] dv + \int_{v\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} + Pv\left(T, P\right) - RT
$$

Moreover:

$$
\int_{v\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} = RT \ln\left[\frac{v^{\bullet}\left(T,P\right)}{v\left(T,P\right)}\right] = RT \ln\left[\frac{RT}{P\times v\left(T,P\right)}\right] = - RT \ln z\left(T,P\right)
$$

Finally:

$$
\ln \phi_i\left(T, P\right) = -\frac{1}{RT}\int_{+\infty}^{v\left(T,P\right)}\left[P\left(T, v\right) - \frac{RT}{v} \right] dv - \ln z\left(T,P\right) + z\left(T,P\right) - 1
$$