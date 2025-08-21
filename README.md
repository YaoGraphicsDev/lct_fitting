A small program that fits Linearly Transformed Cosine parameters to GGX BRDF. Largely inspired by [Eric Heitz's LTC research](https://eheitzresearch.wordpress.com/415-2/).

Also accounts for implementation caveats mentioned in the lecture [Real-Time Area Lighting:
a Journey from Research to Production](https://advances.realtimerendering.com/s2016/s2016_ltc_rnd.pdf) by Stephen Hill.

The goal is to generate multiple textures that store tabulated parameters that represent $\mathbf{M}$, parameterized by viewing angle $\theta_v$ and roughness $\sqrt{\alpha}$

$$
\mathbf{M} = \mathbf{B}\mathbf{m}
$$

$\mathbf{m}$ is a transformation that scales the cosine lobe in $x$ and $y$ direction, in the meantime shearing $x$ with respect to $z$. Jacobian of this transformation will ensure $z$ component being scaled accordingly to keep the transformed distribution normalized.

$$
\mathbf{m} =
\begin{bmatrix}
m_{11} & 0 & m_{13} \\
0 & m_{22} & 0 \\
0 & 0 & 1
\end{bmatrix}
$$

$\mathbf{B}$ is a rotation matrix around the $+y$ axis:

$$
\mathbf{B} = 
\begin{bmatrix}
\cos\bar{\theta_i} & 0 & -\sin\bar{\theta_i} \\
0 & 1 & 0 \\
\sin\bar{\theta_i} & 0 & \cos\bar{\theta_i}
\end{bmatrix}
$$

$\theta_i$ is the polar angle of weighted average incident direction $\bar{\omega_i}$ given an outgoing direction $\omega_o$:

$$
\bar{\omega_i} = \int_{\Omega} f_r(\omega_i, \omega_o)\cos\theta_i\omega_id\omega_i
$$

Notice the similarity between the above weighted average and the definition of albedo:

$$
Albedo = \int_{\Omega} f_r(\omega_i, \omega_o)\cos\theta_id\omega_i
$$

Since $\bar{\omega_i}$ is typically derived from a known BRDF (in our case, GGX), only the three scaling and shearing parameters of the matrix $\mathbf{m}$ require fitting. This separation of rotation from shearing-scaling factors, as opposed to directly fitting all five parameters of $\mathbf{M}$, reduces the likelihood of converging to a local minimum. Consequently, it yields more consistent results across adjacent entries in the resulting LUTs.

Expanded form of $\mathbf{M}$:

$$
\mathbf{M}
 = \begin{bmatrix}
m_{11}\cos\bar{\theta_i} & 0 & m_{13}\cos\bar{\theta_i}-\sin\bar{\theta_i} \\
0 & m_{22} & 0 \\
m_{11}\sin\bar{\theta_i} & 0 & m_{13}\sin\bar{\theta_i}+\cos\bar{\theta_i}
\end{bmatrix}
 = \begin{bmatrix}
a & 0 & b \\
0 & c & 0 \\
d & 0 & e
\end{bmatrix}
$$

The LTC paper divides $\mathbf{M}$ by its bottom right component $e$ to save one tabulated parameter. But in practice it is seen as a "false economy" because a second texture storing albedo has to be introduced regardless. We compute and tabulate this "normalized" matrix only to compare against the results given by the author.

$$
\mathbf{\hat{M}}
 = \frac{1}{e}\mathbf{M}
 = \begin{bmatrix}
\dfrac{a}{e} & 0 & \dfrac{b}{e} \\[8pt]
0 & \dfrac{c}{e} & 0 \\[8pt]
\dfrac{d}{e} & 0 & 1 \\[8pt]
\end{bmatrix}
 = \begin{bmatrix}
\hat{a} & 0 & \hat{b} \\[8pt]
0 & \hat{c} & 0 \\[8pt]
\hat{d} & 0 & 1 \\[8pt]
\end{bmatrix}
$$

It is also convenient to write paramters of $\mathbf{M^{-1}}$ and $\mathbf{\hat{M}^{-1}}$ to an LUT for later use in real-time polygonal lighting.

$$
\mathbf{M^{-1}}
 = \begin{bmatrix}
\dfrac{e}{\Delta} & 0 & -\dfrac{b}{\Delta} \\[8pt]
0 & \dfrac{1}{c} & 0 \\[8pt]
-\dfrac{d}{\Delta} & 0 & \dfrac{a}{\Delta}\\[8pt]
\end{bmatrix}
$$

where $\Delta = \det{\mathbf{M}} = ae - bd$.

Similarly, introduce $\hat{\Delta} = \det{\mathbf{\hat{M}}} = \hat{a} - \hat{b}\hat{d}$, we can write:

$$
\mathbf{\hat{M}^{-1}}
= \begin{bmatrix}
\dfrac{1}{\hat{\Delta}} & 0 & -\dfrac{\hat{b}}{\hat{\Delta}} \\[8pt]
0 & \dfrac{1}{\hat{c}} & 0 \\[8pt]
-\dfrac{\hat{d}}{\hat{\Delta}} & 0 & \dfrac{\hat{a}}{\hat{\Delta}}\\[8pt]
\end{bmatrix}
= \frac{1}{\hat{\Delta}}\begin{bmatrix}
1 & 0 & -\hat{b} \\[8pt]
0 & \dfrac{\hat{\Delta}}{\hat{c}} & 0 \\[8pt]
-\hat{d} & 0 & \hat{a}\\[8pt]
\end{bmatrix}
$$


## $x$ Shear Versus $z$ Shear
TODO

## Preconditioning Roughness and Viewing angle
TODO

## Input Parameters

```
LCTFit.exe [OPTIONS] SUBCOMMAND


OPTIONS:
  -h,     --help              Print this help message and exit
  -d,     --output-dir TEXT:DIR [./build/out]
                              output directory

SUBCOMMANDS:
  dist-ggx                    Print values of GGX BRDF on upper hemisphere to a file
  dist-lct                    Print values of LCT BRDF on upper hemisphere to a file
  dist                        Print values of LCT and GGX on upper hemisphere to files
  gen_luts                    write M matrix and albedo to textures
```

## Compute and Write GGX BRDF data
```
LCTFit.exe dist-ggx [OPTIONS]


OPTIONS:
  -h,     --help              Print this help message and exit
  -r,     --res INT [256]     sampling resolution of theta and phi
  -t,     --theta-multiplier FLOAT REQUIRED
                              theta multiplier. Viewing angle theta_v = multiplier * 0.5 * pi
  -u,     --roughness FLOAT REQUIRED
                              roughness
  -o,     --output TEXT [ggx.dat]
                              name of output data file
```

## Compute LCT BRDF data

```
LCTFit.exe dist-lct [OPTIONS]


OPTIONS:
  -h,     --help              Print this help message and exit
  -r,     --res INT [256]     sampling resolution of theta and phi
  -m,     --m-params FLOAT x 3 REQUIRED
                              Matrix M parameters. x-scale, x-z shear, y scale
  -y,     --y-angle FLOAT [0]
                              Angle of rotation around +y axis, in degrees
  -a,     --albedo FLOAT [1]  albedo parameter
  -o,     --output TEXT [lct.dat]
                              name of output data file
```

## Generate M Matrix Parameter LUTs

subcommand `gen_luts`:

```
  -h,     --help              Print this help message and exit
  -t,     --theta-multiplier [INT,FLOAT,FLOAT] [64, 0.0, 0.99]
                              range of theta multiplier: <resolution> <min> <max>
  -r,     --roughness [INT,FLOAT,FLOAT] [64, 0.03, 1.0]
                              range of roughness: <resolution> <min> <max>
```

Data stored in each channel of each texture:

||Description|R|G|B|A|
|------|------|----|----|----|----|
|`lut_m.exr`|$\mathbf{m}$ matrix|$m_{11}$|$m_{13}$|$m_{22}$||
|`lut_b.exr`|$\mathbf{B}$ matrix|$\sin\bar{\theta_i}$|$\cos\bar{\theta_i}$|||
|`lut_bm_0.exr`|$\mathbf{Bm}$, i.e. $\mathbf{M}$|$a$|$b$|$c$|| 
|`lut_bm_1.exr`|$\mathbf{Bm}$, i.e. $\mathbf{M}$|$d$|$e$|$Albedo$||
|`lut_inv_0.exr`|$\mathbf{M^{-1}}$|$\dfrac{e}{\Delta}$|$-\dfrac{b}{\Delta}$|$\dfrac{1}{c}$|
|`lut_inv_1.exr`|$\mathbf{M^{-1}}$|$-\dfrac{d}{\Delta}$|$\dfrac{a}{\Delta}$|$Albedo$|
|`lut_norm_inv.exr`|$\mathbf{\hat{M}^{-1}}$|$\hat{a}$|$-\hat{b}$|$\dfrac{\hat{\Delta}}{\hat{c}}$|$-\hat{d}$|


