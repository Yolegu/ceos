---
title: Pfaffian equations
---

From Euler's theorem:

\[G = \sum_j n_j \mu_j\]

Deriving this equation with respect to \(n_i\) at fixed temperature and pressure:

\[\bar{g}_i = \sum_{j}\left(\frac{\partial \mu_j}{\partial n_i}\right)_{T,P,n_{j \neq i}} + \mu_i\]

Remembering that \(\bar{g}_i = \mu_i\), one concludes that:

\[\sum_{j}\left(\frac{\partial \mu_j}{\partial n_i}\right)_{T,P,n_{j \neq i}}=0\]

---

\[dG = VdP - SdT + \sum_j \mu_j dn_j\]

\[d\left(\sum_j n_j \bar{g}_j\right) = \left(\sum_j n_j \bar{v}_j\right)dP - \left(\sum_j n_j \bar{s}_j\right)dT + \sum_j \mu_j dn_j\]

\[\left(\sum_j n_j d\bar{g}_j\right) + \left(\sum_j \bar{g}_j dn_j\right) = \left(\sum_j n_j \bar{v}_j\right)dP - \left(\sum_j n_j \bar{s}_j\right)dT + \sum_j \mu_j dn_j\]

\[\left(\sum_j n_j d\bar{g}_j\right) = \left(\sum_j n_j \bar{v}_j\right)dP - \left(\sum_j n_j \bar{s}_j\right)dT\]

\[\sum_j n_j \left(d\bar{g}_j - \bar{v}_j dP + \bar{s}_j dT\right) = 0\]

The molar quantites form a linearly independant set of parameters. Thus:

\[d\bar{g}_j - \bar{v}_j dP + \bar{s}_j dT = 0\]

In conlusion:

\[d\bar{g}_j = d\mu_j = \bar{v}_j dP - \bar{s}_j dT\]

---

One can write:

\[U = TS - PV + \sum_j n_j \mu_j\]

Deriving this equation with respect to \(n_i\) at fixed temperature and pressure:

\[\bar{u}_i = T\bar{s}_i - P\bar{v}_i + \underbrace{\sum_{j}\left(\frac{\partial \mu_j}{\partial n_i}\right)_{T,P,n_{j \neq i}}}_0 + \mu_i\]

Differenciating this equation:

\[d\bar{u}_i = Td\bar{s}_i + \bar{s}_idT - Pd\bar{v}_i - \bar{v}_idP + d\mu_i\]

This relation can be simplified if one writes that \(d\mu_j = \bar{v}_j dP - \bar{s}_j dT\):

\[d\bar{u}_i = Td\bar{s}_i - Pd\bar{v}_i\]