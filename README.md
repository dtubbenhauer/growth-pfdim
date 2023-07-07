# Code and Erratum for *Asymptotics in finite monoidal categories*

I collected a bit of Magma and Mathematica code relevant for the paper *Asymptotics in finite monoidal categories*
<a href="https://arxiv.org/abs/2307.03044">https://arxiv.org/abs/2307.03044</a> on this page.

The code is in a **.odt** file or **.n** file, respectively, that can be downloaded from this site and you can run it with Mathematica.

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
m:=5; //The order
G:=DihedralGroup(m); //The dihedral group
X:=CharacterTable(G);
M:=X[#X]; //Take the last representation
n:=20; //Mamimal tensor power
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
[ 1, 3, 4, 11, 17, 42, 71, 163, 292, 639, 1189, 2522, 4811, 9999, 19388, 39763,
77929, 158442, 312703, 632171 ]
```

Second, the **asymptotic formula** for $a(n)$:

```
m:=5;
G:=DihedralGroup(m); //set group: dihedral group of order 2m
X:=CharacterTable(G);
V:=X[#X]; //set module

fus:=[[0 : i in [1..#X]] : j in [1..#X]];

for i in [1..#X] do
for j in [1..#X] do
fus[j][i]:=InnerProduct(X[i]*V,X[j]); //computes multiplicity of X[j] in X[i]\otimes V
end for;
end for;

M:=Matrix(Rationals(),fus); //use a cyclotomic field in order to have every complex eigenvalues

w0:=Eigenspace(M,2).1; //computes the vector w_0 (warning: magma computes the eigenspace of the linear map given by right multiplication by M)
vv0:=Eigenspace(Transpose(M),2).1;
v0:=vv0/ScalarProduct(w0,vv0); //computes v_0 with the normalization so that w_0^{T}v_0=1

&+[(Matrix(Rationals(),#X,1,ElementToSequence(v0))*Matrix(Rationals(),1,#X,ElementToSequence(w0)))[i][1] : i in [1..#X]]; //computes v_0w_0^{T}[1]


if m mod 2 eq 0 then //if m is even, we have a coefficient in front of (-1)^n
w1:=Eigenspace(M,-2).1; //computes the vector w_1 (warning: magma computes the eigenspace of the linear map given by right multiplication by M)
vv1:=Eigenspace(Transpose(M),-2).1;
v1:=vv1/ScalarProduct(w1,vv1); //computes v_1 with the normalization so that w_1^{T}v_1=1

&+[(Matrix(Rationals(),#X,1,ElementToSequence(v1))*Matrix(Rationals(),1,#X,ElementToSequence(w1)))[i][1] : i in [1..#X]]; //computes v_1w_1^{T}[1]
end if; 
```

The code computes the action matrix and its eigenspaces and pieces them together as explained above. The output in this case is

```
3/5
```

which are the (at most) two constants in the formula

$$a(n)=
\begin{cases}
\frac{m+1}{2m}\cdot 2^{n} & \text{if } m \text{ is odd},
\\
\frac{m+2}{2m}\cdot 2^{n} & \text{if } m \text{ is even and }m^{\prime}\text{ is odd},
\\
\left(\frac{(m+2)}{2m}\cdot 1+\frac{1}{m}\cdot(-1)^{n}\right)\cdot 2^{n} & \text{if }m\text{ is even and } m^{\prime}\text{ is even}.
\end{cases}$$

where $m^{\prime}=m/2$, if $m$ is even, and $m^{\prime}=(m-1)/2$, if $m$ is odd.

# The Mathematica code

The code is hopefully pretty straightforward. All that is happening is that 
the data created with the Magma code is plotted.

For example, here is the dihedral group of order ten:

![The dimension 3 case](https://github.com/dtubbenhauer/growth-pfdim/blob/main/growth-pfdim-mathematica.png)

Let me know if there are any questions!

# Erratum

Empty so far.
