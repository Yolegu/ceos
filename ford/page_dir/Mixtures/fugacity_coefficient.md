---
title: Fuagcity coefficient
---

By definition:

\[RT \ln\phi_i\left(T,P,x\right) = \bar{g}^{\phi}_i\left(T,P,x\right) - \bar{g}^{\bullet}_i\left(T,P,x\right)\]

Switching to the Helmhotz free energy:

\[RT \ln\phi_i\left(T,P,x\right) = \bar{a}^{\phi}_i\left(T,P,x\right) + P\bar{v}_i^{\phi}\left(T,P,x\right) - \bar{a}^{\bullet}_i\left(T,P,x\right) - P\bar{v}_i^{\bullet}\left(T,P,x\right)\]

The perfect gas law writes:

\[PV^{\bullet} = nRT \Rightarrow P\bar{v}_i^{\bullet} = RT\]

Moreover:

\[\begin{cases}
    \bar{a}^{\phi}_i\left(T,P,x\right) = \bar{a}^{\phi}_i\left(T,\bar{v}_i = \bar{v}_i^{\phi},x\right) \\
    \bar{a}^{\bullet}_i\left(T,P,x\right) = \bar{a}^{\bullet}_i\left(T,\bar{v}_i = \bar{v}_i^{\bullet},x\right) \\
    \lim_{ \bar{v}_i\rightarrow +\infty}\bar{a}^{\phi}_i\left(T, \bar{v}_i,x\right) = \lim_{ \bar{v}_i \rightarrow +\infty}\bar{a}^{\bullet}_i\left(T, \bar{v}_i,x\right)
\end{cases}\]

Thus:

\[RT \ln\phi_i\left(T,P,x\right) = \left[\bar{a}^{\phi}_i\left(T,\bar{v}_i = \bar{v}_i^{\phi},x\right)-\bar{a}^{\phi}_i\left(T,\bar{v}_i \rightarrow +\infty,x\right)\right] - \left[\bar{a}^{\bullet}_i\left(T,\bar{v}_i=\bar{v}_i^{\bullet},x\right)-\bar{a}^{\bullet}_i\left(T,\bar{v}_i \rightarrow +\infty,x\right)\right]  + P\bar{v}_i^{\phi}\left(T,P,x\right) - RT\]

From basic thermodynamics, one gets at constant temperature:

\[d\bar{a}_i = -Pd\bar{v}_i\]

Therefore:

\[RT \ln\phi_i\left(T,P,x\right) = -\int^{\bar{v}_i^{\phi}}_{+\infty}P^{\phi}\left(T,v,x\right)d\bar{v}_i + \int^{\bar{v}_i^{\bullet}}_{+\infty}P^{\bullet}\left(T,v,x\right)d\bar{v}_i  + P\bar{v}_i^{\phi}\left(T,P,x\right) - RT\]

Euler's theorem allows us to relate \(v\) to \(\bar{v}_i\):

\[v = \sum_i x_i \bar{v}_i\]

By splitting the second integral and injecting the perfect gas pressure definition, one gets: 

\[RT \ln\phi_i\left(T,P,x\right) = -\int^{\bar{v}_i^{\phi}}_{+\infty}\left[P^{\phi}\left(T,v,x\right)-\frac{RT}{\bar{v}_i}\right]d\bar{v}_i + \int^{\bar{v}_i^{\bullet}}_{\bar{v}_i^{\phi}}\frac{RT}{\bar{v}_i}d\bar{v}_i  + P\bar{v}_i^{\phi}\left(T,P,x\right) - RT\]

The second integral can be immediately integrated:

\[RT \ln\phi_i\left(T,P,x\right) = -\int^{\bar{v}_i^{\phi}}_{+\infty}\left[P^{\phi}\left(T,v,x\right)-\frac{RT}{\bar{v}_i}\right]d\bar{v}_i + RT \ln\left(\frac{\bar{v}_i^{\bullet}}{\bar{v}_i^{\phi}}\right) + P\bar{v}_i^{\phi}\left(T,P,x\right) - RT\]

Writting \(\bar{z}_i^{\phi} = P\bar{v}_i^{\phi} / RT\), one gets:

\[\ln\phi_i\left(T,P,x\right) = -\frac{1}{RT}\int^{\bar{v}_i^{\phi}}_{+\infty}\left[P^{\phi}\left(T,v = \sum_j{x_j \bar{v}_j},x\right)-\frac{RT}{\bar{v}_i}\right]d\bar{v}_i - \ln\left(\bar{z}_i\right) + \bar{z}_i - 1\]