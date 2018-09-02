---
title: Fugacity coefficient calculation
---

The fugacity coefficient of a pure compound \(i\) is defined as :

$$
RT \ln \phi_i = g_i^{res}\left(T,P\right) = g_i^{\phi}\left(T^{\phi} = T,P^{\phi} = P\right) - g_i^{\bullet}\left(T^{\bullet} = T,P^{\bullet} = P\right)
$$

Because pressure explicit EoS have temperature and pressure as natural variables, it is convenient to express the different terms from previous equation using such variables.

$$
RT \ln \phi_i = a_i^{\phi}\left(T^{\phi} = T, P^{\phi} = P\right) + \underbrace{\left(Pv\right)^{\phi}}_{Pv^{\phi}\left(T,P\right)} - a_i^{\bullet}\left(T^{\bullet} = T, P^{\phi} = P\right) - \underbrace{\left(Pv\right)^{\bullet}}_{RT}
$$

Because for a real fluid, a property value is independent of the variables used to express it, one can write: \(a_i\left(T^{\phi} = T,P^{\phi} = P\right) = a_i\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right)\). Moreover, one can write:

$$
a_i^{\bullet}\left(T^{\bullet} = T, P^{\bullet} = P\right) = a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} = \frac{RT}{P}\right)
$$

Finally, one can write:

$$
RT \ln \phi_i = a_i^{\phi}\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right) - a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} = \frac{RT}{P}\right)+ Pv^{\phi}\left(T,P\right) - RT
$$

At the infinite molar volume limit, the real fluid and the perfect gas are identical. Thus, one can write:

$$
a_i^{\phi}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right) = a_i^{\bullet}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right)
$$

Then, \(RT \ln \phi\) expression rewrites:

$$
RT \ln \phi_i = a_i^{\phi}\left(T^{\phi} = T,v^{\phi}\left(T, P\right)\right) - a_i^{\phi}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right) - \left[a_i^{\bullet}\left(T^{\bullet} = T,v^{\bullet} = \frac{RT}{P}\right) - a_i^{\bullet}\left(T^{\phi} = T,v^{\phi} \rightarrow +\infty\right)\right]+ Pv^{\phi}\left(T,P\right) - RT
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
RT \ln \phi_i = -\int_{+\infty}^{v\left(T,P\right)}\left[P^{\phi}\left(T^{\phi} = T, v\right) - \frac{RT}{v} \right] dv + \int_{v^{\phi}\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} + Pv^{\phi}\left(T, P\right) - RT
$$

Moreover:

$$
\int_{v^{\phi}\left(T,P\right)}^{v^{\bullet}\left(T,P\right)}\frac{RT}{v^{\bullet}}dv^{\bullet} = RT \ln\left[\frac{v^{\bullet}\left(T,P\right)}{v^{\phi}\left(T,P\right)}\right] = RT \ln\left[\frac{RT}{P\times v^{\phi}\left(T,P\right)}\right] = - RT \ln z^{\phi}\left(T,P\right)
$$

Finally:

$$
\ln \phi_i = -\frac{1}{RT}\int_{+\infty}^{v\left(T,P\right)}\left[P^{\phi}\left(T^{\phi} = T, v\right) - \frac{RT}{v} \right] dv - \ln z^{\phi}\left(T,P\right) + z^{\phi}\left(T,P\right) - 1
$$