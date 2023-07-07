# Code and Erratum for *Asymptotics in finite monoidal categories*

I collected a bit of Magma and Mathematica code relevant for the paper *Asymptotics in finite monoidal categories*
<a href="https://arxiv.org/abs/2307.03044">https://arxiv.org/abs/2307.03044</a> on this page.

The code is in a **text** file or **.n** file, respectively, that can be downloaded from this site and you can run it with Mathematica.

An Erratum for the paper *Asymptotics in finite monoidal categories* can be found at the bottom of the page.

# Contact

If you find any errors in the paper *Asymptotics in finite monoidal categories* **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.

# Background

Let $\mathbb{K}=\mathbb{R}_{\geq 0}$, our ground semiring. The slogan is: **Nonnegativity is crucial**.

A finite based $\mathbb{K}$-algebra $R=(R,C)$ is an algebra with a fixed basis $C$ such that 
the structure constants with respect to $C$ are in $\mathbb{K}$.

Fix $c\in R$.
We write $m_{n}^{\ast}(c)$ for these coefficients as 
they appear in $c^{n}$ where $\ast\in\{0,\dots,r-1\}$.
Define
$$b_{n}^{R,c}:=\text{total sum of coefficients $m_{n}^{\ast}(c)$}.$$
Moreover, we define the function
$$b^{R,c}(n)\colon\mathbb{N}\to\mathbb{R},n\mapsto b_{n}^{R,c}.$$
We are interested in the **asymptotic behavior** of the function 
$b^{R,c}(n)$.

For $a_{i}\in\mathbb{K}$, the (transposed) **action matrix** of 
$c=a_{0}\cdot c_{0}+\dots+a_{r-1}\cdot c_{r-1}\in R$ is the matrix $(\sum_{i}a_{i}m_{i,j}^{k})_{k,j}$. Abusing language, we 
will call the 
submatrix of it corresponding to the connected component of $1$ 
also the action matrix and use this below.

Assume (only for the sake of this page - the note does not need that assumption) that the **Perron--Frobenius theorem** holds, that is
the action matrix of $c\in R$ has a leading eigenvalue $\lambda_{0}=\mathrm{PFdim}c$ of multiplicity one that we call 
the **Perron--Frobenius dimension** of $c$. 
Moreover, the action matrix 
has some period $h\in\mathbb{N}$ such that there are precisely $h-1$ other eigenvalues 
$\lambda_{i}=\zeta^{i}\mathrm{PFdim}c$, and all of these are of multiplicity one, where 
$\zeta=\exp(2\pi i/h)$.

Take the left and right eigenvectors $v_{i}$ and $w_{i}$, normalized such that the transpose of $w_i$ times $v_i$ is one.

Let $v_{i}w_{i}^{T}[1]$ denote taking the sum of the first column of the matrix 
$v_{i}w_{i}^{T}$. Define

$$a(n)=\big(v_{0}w_{0}^{T}[1]\cdot 1+v_{1}w_{1}^{T}[1]\cdot\zeta^{n}+v_{2}w_{2}^{T}[1]\cdot(\zeta^{2})^{n}+\dots+v_{h-1}w_{h-1}^{T}[1]\cdot(\zeta^{h-1})^{n}\big)
\cdot(\mathrm{PFdim}c)^{n}.$$

Let $\lambda^{sec}$ be the second largest eigenvalue of the action matrix of $c$.

**Theorem**

>We have
$$b^{R,c}(n)\sim a(n)$$
>and the convergence is geometric with ratio $|\lambda^{sec}/\mathrm{PFdim} c|$. In particular,
$$\beta^{R,c}:=\lim_{n\to\infty}\sqrt[n]{b_{n}^{R,c}}=\mathrm{PFdim} c.$$

This theorem is interesting from a **categorical point of view** since Grothendieck rings of certain nice enough monoidal categories are 
finite based $\mathbb{K}$-algebras. In particular, we also get asymptotics for the
$$\text{number of indecomposable summands in $\mathtt{X}^{\otimes n}$ counted with multiplicities},$$
where $\mathtt{X}$ is an object in a fixed nice enough monoidal category.

The code attached to this page computes the function $a(n)$ for various examples.

# The Magma code

The code can, for example, be run in the online calculator of Magma:

<a href="https://magma.maths.usyd.edu.au/calc/">Magma online</a>

The code we collected is about the dihedral group of order $2m$ and can be adjusted to 
other groups easily.

Let us go through the code. 

The first part is about comupting $b(n)$ -- **the numbers we what to approximate**.  

```
m:=4;
G:=DihedralGroup(m);
X:=CharacterTable(G);
M:=X[#X];
n:=20;
out:=0;
outlist:=[];

for k in [1..n] do
for j in [1..#X] do
out+:=InnerProduct(M^k,X[j]);
end for;
outlist:=Append(outlist,out);
out:=0;
end for;
outlist;
```

which sets up the dihedral group, its character table and takes $M$ to be the last representation in the character table and what will be eventually the output.
It then runs a loop to check the numbers $b(n)$ and outputs them up the fixed n. 

The output in this example is

```
[ 1, 4, 4, 16, 16, 64, 64, 256, 256, 1024, 1024, 4096, 4096, 16384, 16384,
65536, 65536, 262144, 262144, 1048576 ]
```

# The Mathematica code

The code is hopefully pretty straightforward. All that is happening is that 
the data created with the Magma code is plotted.

For example, here is the dihedral group of order ten:

![The dimension 3 case](https://github.com/dtubbenhauer/growth-pfdim/blob/main/growth-pfdim-mathematica.png)

Let me know if there are nay questions!

# Erratum

Empty so far.
