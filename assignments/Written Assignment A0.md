# Discrete Differential Geometry: Written Assignment A0

- Exercise 2.1

My proof sketch is to triangulate the original shape, then show the changes of V, E and F

- Exercise 2.2

https://www.mathsisfun.com/geometry/platonic-solids-why-five.html

Suppose we have $s$ sides for each faces, $m=\text{val}(v_i)$ for each vertex, then we decompose the sphere, got:
$$
sF=2E=mV
$$
Then we put the equation into Euler formula:
$$
\frac{2E}{m}-E+\frac{2E}{S}=2
$$
Derive it:
$$
\frac{1}{m}+\frac{1}{s}=\frac{1}{E}+\frac{1}{2}
$$
And we have $m \ge 3$, $s \ge 3, E \ge 6$, and valid solution for m, s pairs are $(3,3), (3,4), (3,5), (5,3), (4,3)$, Q.E.D.

- Exercise 2.3

Proof Sketch is similar to the topology solution to ex 2.2

First we know $V-E+F=2-2g$, then we decompose the surfaces, then we got two extra constraints:

1. For vertices, we have $\text{val}(v)=6$, so $6V=2E$
2. For new edges, we got $e'F=2E$, since the surface is connected and orientable and finite, we could decompose it to a triangular mesh, then $e'=3$, so $3F=2E$

Putting together, we got $\frac{1}{3}E-E+\frac{2}{3}E=2-2g$, then $g=1$, so the only solution is a torus.

- Exercise 2.4

http://yblin.is-programmer.com/posts/43217

Decompose vertices into $n$ irregular vertices and $V-n$ regular vertices, then we got:
$$
2E=3F=6(V-n)+\sum_{i=1}^{n}val(v_i)
$$
Put it into Euler formula:
$$
V-\frac{6(V-n)+\sum_{i=1}^{n}val(v_i)}{6}=2-2g
$$

$$
n-\frac{1}{6}\sum_{i=1}^{n} val(v_i) = 2-2g
$$

$$
n = \frac{1}{6}\sum_{i=1}^{n} val(v_i) + 2-2g
$$

Since $val(v_i)\ge3$:
$$
n\geq2-2g+\frac{n}{2}
$$
So
$$
n \geq 4-4g
$$
Case 1: $g=0$, then we got $n \geq 4$

Case 2: $g = 1$, then we got $n \geq 0$ 

Case 3: $g \geq 2$, we cannot use equation (9) anymore, when $n=0$ does not satisfy equation (7), and $n=1$ is a valid solution.

- Exercise 2.5

We already have the equation
$$
2E=3F=\sum val(v_i) = \frac{\sum val(v_i)}{V}V
$$
And now that $\frac{\sum val(vi)}{V} = 6$, Then $V:E:F=1:3:2$. 

- Exercise 2.6

Consider grid-like quad mesh, then $2E=4F$, and each vertex is connected to 4 faces, so $V:E:F=1:2:1$

- Exercise 2.7

Consider the combination of icosahedrons, each face corresponds to 2 tets, so $2F=4T$. And In each icosahedrons, one vertex corresponds to 20 tets, and each tet correspond to 4 vertices, so $V:T=1:5$, Now we got $V:E:F:T=1:6:10:5$ due to the Euler formula.

- Exercise 2.8

Star of edges are edges themselves and triangles that contains that edge (in 2D.

Star of triangles are themselves (in 2D).

$St(S)$:

Vertices: e, k, o

Edges: ea, eb, ef, ek, ej, ed, of, ol, oq, op, ok, kp, kj

Faces: aed, aeb, ebf, ekf, ejk, ejd, eda, ofl, ols, opq, oks, oks, jkp, njp



$Cl(S)$:

Vertices: e, f, o, k, j, p, n

Edges: jk, kp, jp, jn, np, of

Faces: jnp, jkp



$Lk(S)=Cl(St(S)) \backslash St(Cl(S)) $

Edges: ab, lq, ad

Vertices: a, b, d, l, q

(This should be verified with code... Too complicated)

- Exercise 2.9

Bd(K)=abflmsqpjea

Int(K)=everything inside cl(K) besides bd(K), this is easy.

- Exercise 2.10

Twin of 0-9: 4215037698

Next of 0-9: 1204563978

- Exercise 2.11

Draw it on your paper :)

- Exercise 2.12

For $A_0$, mark corresponding vertices to be 1 for each face, and do the same job for edges in $A_1$

- Exercise 2.13

Explain Simplicity 1-manifold looks like a closed loop or path. (This is not a proof)

An important to mention is that a simplicial 0-sphere is a set of two points, since we have the definition: 
$$
\mathbb{S}^n:=\{x\in \mathbb R^{n+1} : |x|=1\}
$$
When n = 0, it reduced to ${-1,1}$

Branches are not allowed in the link of simplicial 1-manifold since it would lead to at least 3 points, which are not homeomorphic to 0-sphere.

So when we expand the link, it could not generate the higher dimensional simplifies like faces and tets. So we just got path and closed loops

- Exercise 2.14

Proof sketch is the same. The link of a vertex in simplicial 2-manifold is homeomorphic to a circle, which looks like a closed loop. So we do not have the "suddenly ending" points(like y-forks, whose link is a 0-sphere)

Personally I cannot write a strict proof, but those two questions give me a feeling of a well defined simplicial k-manifold with links :)

- Exercise 2.15

The explanation is based on the definition of boundary on simplicial k-manifold:

all (k−1)-simplices in K that are contained in exactly one k-simplex

And we cannot find a (k-2)-simplices in K that belongs to exactly one (k-1)-simplices.

(Take k=2 as an example, boundary of simplicial 2-manifold is the outside edges, but each vertex on the boundary belongs to 2 edges, so there is no boundary.)

There is a proof with a sequence definition of boundary maps: https://algebrology.github.io/simplicial-complexes-and-boundary-maps/

