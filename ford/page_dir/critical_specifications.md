---
title: Critical Specifications
---

Apply the critical specifications to estimate the CEoS universal parameters \(\eta_c\), \(\Omega_a\), \(\Omega_b\) and \(z_c\) and the pure compound parameters \(a_c\), \(b_c\) and \(v_c\)

**CEoS universal parameters**

$$\begin{equation}
\begin{cases}
\eta_c = \frac{1}{\left[\left(1-r_1\right)\left(1-r_2\right)^2\right]^{1/3} + \left[\left(1-r_2\right)\left(1-r_1\right)^2\right]^{1/3} + 1} \\
\Omega_a =\frac{\left(1-\eta_c r_1\right)\left(1-\eta_c r_2\right)\left[2-\eta_c \left(r_1 + r_2\right)\right]}{\left(1-\eta_c\right)\left[3-\eta_c \left(1+r_1+r_2\right)\right]^2} \\
\Omega_b = \frac{\eta_c}{3-\eta_c \left(1 + r_1 + r_2\right)} \\
z_c = \frac{1}{3 - \eta_c \left(1 + r_1 + r_2\right)}
\end{cases}
 \end{equation}$$

**Pure compound universal parameter**

$$\begin{cases}
a_c = \Omega_a \times \frac{\left(RT_c\right)^2}{P_c} \\
b_c = \Omega_b \times \frac{RT_c}{P_c} \\
v_c = \frac{b_c}{\eta_c}
\end{cases}$$