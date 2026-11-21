---
title: "Perelman's Proof, Read Through Addition: Heat Carries the Answer, Integers Decide It"
published: true
description: "Why a shape problem needed heat and entropy, and a ledger that climbs one whole step at a time — Poincaré, Ricci flow, and Fibonacci on a 13-hour clock."
tags: math, geometry, python, showdev
cover_image:
# canonical_url: (Zenn URL of the Japanese original, if you want it)
---

<!-- Cover / Figure 1: upload perelman_dumbbell.gif and paste its URL in cover_image and below -->

![perelman_dumbbell.gif](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/08i1cu85mgzwrbjbgtqm.gif)


## The triple disappointment

A story circulates about the moment the Poincaré conjecture fell. Mathematicians listening to Perelman were disappointed three times: first that the conjecture had been solved at all, second that it had been solved by differential geometry instead of topology, and third that they could not follow the explanation. (It is an anecdote passed from mouth to mouth; I do not vouch for the details.)

The proof is full of words that seem to come from a physics textbook: heat, entropy, free energy. Why would a question about shape need heat?

This article follows one reading of that "why", and it fits in one line:

**Geometry is continuous. The decision is an integer. What carries the shape from one to the other is addition — averaging.**

Part 1 retells the proof in that light. Part 2 builds a tiny device that uses nothing but addition and asks what "entropy" means inside it. All of the mathematics in Part 1 is known; what I add is only the arrangement. Part 2 is computed, and every number in it can be rerun with the short script at the end.

---

## Part 1 — The heat that carried the answer

## 1. The question is about integers

The Poincaré conjecture (1904):

> Every simply connected, closed 3-manifold is homeomorphic to the 3-sphere {% katex inline %}S^3{% endkatex %}.

**Closed** means finite in extent and without an edge. The surface of the Earth is a closed 2-dimensional manifold: you can walk forever and never fall off a cliff.

**Simply connected** means that any loop you draw can be shrunk to a point without snagging. On a sphere every rubber band slides off. On the surface of a doughnut, a band that winds around the hole cannot be removed without cutting it.

How badly a band is stuck can be counted: how many times it winds around each hole. On a torus the pair {% katex inline %}(m, n) \in \mathbb{Z}^2{% endkatex %} is a complete ledger of loops. This ledger is the fundamental group {% katex inline %}\pi_1{% endkatex %}. "Simply connected" means **the ledger is blank**: {% katex inline %}\pi_1 = 0{% endkatex %}.

So the conjecture says: **from integer data alone — a blank winding ledger — the shape of a closed 3-dimensional space is decided. It is the sphere and nothing else.** Both the question and the answer are discrete.

The order of the solutions was strange. Smale settled dimensions 5 and up in the 1960s, and Freedman settled dimension 4 in 1982. Higher dimensions have more room to untangle loops and discs. Dimension 3 has too little room for cut-and-paste tricks, and it held out for a century.

## 2. Curvature is already a sum

Hamilton and then Perelman changed the toolbox. Give the space a **metric** {% katex inline %}g{% endkatex %}, a ruler at every point and in every direction. Then you can measure **curvature**.

The most concrete definition: draw a small ball of radius {% katex inline %}r{% endkatex %} around a point {% katex inline %}p{% endkatex %} and measure its volume. In flat space it is {% katex inline %}\tfrac{4}{3}\pi r^3{% endkatex %}. If the real volume is smaller, the curvature is positive (sphere-like); if it is larger, the curvature is negative (saddle-like):

{% katex %}
\mathrm{Vol}(B_r(p)) = \tfrac{4}{3}\pi r^3\left(1 - \frac{R(p)}{30}\,r^2 + \cdots\right)
{% endkatex %}

{% katex inline %}R{% endkatex %} is the scalar curvature. Ricci curvature averages the bending over the planes that contain a given direction; scalar curvature adds that up over all directions. The very material of the proof is **made by adding up and averaging what surrounds a point**. Addition is in the ingredients before it is in the method.

## 3. Ricci flow: let heat run through the shape

Hamilton's equation (1982):

{% katex %}
\frac{\partial g}{\partial t} = -2\,\mathrm{Ric}(g)
{% endkatex %}

Read it as: **where the bending is strongest, change the ruler fastest.** Positively curved, bulging parts shrink; negatively curved, saddle-shaped parts spread. In local coordinates the leading part is

{% katex %}
\frac{\partial g_{ij}{\partial t} = \Delta g_{ij} + (\text{lower-order terms})
{% endkatex %}

and {% katex inline %}\Delta{% endkatex %} moves each value toward the average of its neighbours. That is the **heat equation**. Heat evens out temperature; Ricci flow evens out curvature. An average is a sum, so the engine of the proof is addition made continuous.

**Example 1: the round sphere.** Put a round {% katex inline %}S^3{% endkatex %} of radius {% katex inline %}r_0{% endkatex %} into the flow. By symmetry it stays round, and only the radius changes:

{% katex %}
r(t)^2 = r_0^2 - 4t
{% endkatex %}

The squared radius falls in a straight line, and **the sphere vanishes to a point at the finite time {% katex inline %}t = r_0^2/4{% endkatex %}**, still round. This simplest case turns out to be the whole conclusion in miniature.

**Example 2: the dumbbell.** Two big balls joined by a thin neck (Figure 1 computes this case). The neck is strongly curved around its circumference, so **the neck pinches first**. Its radius hits zero in finite time and the equation breaks there. That is a **singularity**.


<!-- Figure 1: perelman_dumbbell.gif -->
![A wireframe dumbbell whose neck pinches under Ricci flow and is then cut at its thinnest point, leaving two bulbs with thin horns and tiny caps](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/3g47an4kasyk412j4vl7.gif)


*The neck pinch is computed from the flow equations for a 3-manifold whose cross-sections are round 2-spheres (each drawn as a circle; the two outer ends are closed with round caps). Frames are spaced evenly in the neck radius {% katex inline %}h{% endkatex %} on a log scale, from 0.26 down to 0.01. The zoom at bottom left always shows {% katex inline %}\pm 15h{% endkatex %} by {% katex inline %}\pm 4h{% endkatex %}, so it shows the neck's shape: it gets straighter as it thins. At {% katex inline %}h = 0.01{% endkatex %} the neck is cut at its thinnest cross-section and each side gets a cap of radius {% katex inline %}h{% endkatex %} (drawn as a hemisphere; Perelman's standard cap is not reproduced). What is left is two large bulbs, each carrying a thin horn. The flow after the cut is not computed.*


*Interactive version (the page computes the flow in your browser): [CodePen — perelman_dumbbell.html](https://codepen.io/editor/morcb13-bit/pen/01a1fb70-ab2e-7a60-83ab-6a6e0587c477)


At the moment of the cut the numbers are lopsided. The neck radius is {% katex inline %}h = 0.01{% endkatex %} while the bulbs have radius about 0.94; the curvature of the neck's 2-spheres is about 10,000, against about 1 in the bulbs. The part thinner than {% katex inline %}2h{% endkatex %} is only 12 neck radii long, and within {% katex inline %}\pm 4h{% endkatex %} of the thinnest point the slope {% katex inline %}|\psi_s|{% endkatex %} still reaches 0.22. The neck straightens as it thins, but slowly (from 0.75 at {% katex inline %}h = 0.26{% endkatex %} to 0.22 at {% katex inline %}h = 0.01{% endkatex %}). A real surgery waits for a neck much thinner and straighter than this one. What it removes is tiny, and what it leaves is a horn with a small cap, not two clean spheres.

## 4. Where Hamilton got stuck

For twenty years Hamilton pushed this programme. He proved that positive Ricci curvature rounds the space into a sphere. In general, though, the flow tears partway through, like the dumbbell. To handle the tear you zoom in on the singularity, rescale, and classify the limiting shapes.

One nightmare lived there: the **cigar soliton**, a self-similar cigar-shaped solution. If the zoomed-in limit could be a cigar, the tearing would be uncontrollable and there would be no way to operate. Hamilton could not rule it out.

## 5. Perelman's three tools

**Tool 1: an entropy for geometry.** Perelman defined

{% katex %}
\mathcal{W}(g, f, \tau) = \int_M \Big[\tau\big(R + |\nabla f|^2\big) + f - n\Big]\,(4\pi\tau)^{-n/2} e^{-f}\, dV
{% endkatex %}

{% katex inline %}e^{-f}\,dV{% endkatex %} is a probability distribution, and the skeleton is that of free energy in statistical mechanics (Perelman devotes a section of his paper to the analogy). He proved that **{% katex inline %}\mathcal{W}{% endkatex %} never decreases along the flow** — a second law of thermodynamics for geometry. A monotone entropy **forbids** behaviour that looks like time running backwards, such as a region collapsing flat and losing its volume. This is the *no local collapsing* theorem, and it rules out the cigar. Hamilton's twenty-year block came loose with a single entropy.

**Tool 2: the catalogue of singularities is finite.** Once collapse is forbidden, the zoomed-in limits (κ-solutions) can be listed. In three dimensions the answer is almost anticlimactic: **a pinching neck always looks like a round cylinder {% katex inline %}S^2 \times \mathbb{R}{% endkatex %}**, plus round spheres and capped cylinders. The dumbbell picture was the general answer. Infinitely many conceivable ways of tearing fall into a finite catalogue.

**Tool 3: surgery, with a guarantee that it ends.** When curvature passes a threshold, **cut the neck along a round {% katex inline %}S^2{% endkatex %}, glue a round cap on each side, and keep flowing.** Each surgery removes a definite amount of volume, so surgeries cannot pile up infinitely in finite time. Cut, cap, flow — a whole number of times.

## 6. The shortcut: the space evaporates completely

For a general 3-manifold you run the flow with surgery for all time and get to Thurston's geometrization. The Poincaré conjecture itself has a shortcut. Perelman's third paper (with a separate proof by Colding and Minicozzi) shows:

**If the space is simply connected, the flow with surgery makes it vanish completely in finite time.**

Every piece rounds off and evaporates. Run the surgeries backwards and the pieces just before extinction are round spheres, cylinders and caps sewn together — that is, {% katex inline %}S^3{% endkatex %}. Rewind all the way, and the original space was {% katex inline %}S^3{% endkatex %} too.

The round sphere vanishing at {% katex inline %}t = r_0^2/4{% endkatex %} in Section 3 was the conclusion in miniature: **blank winding ledger = a shape heat can evaporate completely = the sphere.** The ledger (integers) and the flow (continuous) lie on top of each other here.

## 7. What is left when the flow runs out: eight geometries

Going all the way to geometrization, every closed 3-manifold, after standard cuts, splits into pieces that each wear exactly one of **eight** homogeneous geometries:

{% katex %}
S^3,\ \mathbb{E}^3,\ \mathbb{H}^3,\ S^2\times\mathbb{R},\ \mathbb{H}^2\times\mathbb{R},\ \widetilde{SL_2},\ \mathrm{Nil},\ \mathrm{Sol}
{% endkatex %}

Not nine, not twelve: eight. The input is infinite-dimensional (the space of metrics), the process is continuous (heat), and the output is a list you can count on your fingers.

| Stage | Content | Character |
|---|---|---|
| Question | is the winding ledger {% katex inline %}\pi_1{% endkatex %} blank? | **integer** |
| Material | metric, curvature (sums of volume) | continuous, built from addition |
| Process | Ricci flow = heat = averaging | continuous **addition** |
| Safety device | entropy is monotone | a real number whose job is to **forbid** |
| Singularities | necks are round cylinders, **finite** catalogue | discrete |
| Surgery | cut and cap, **finitely many** times | integer (a count) |
| Output | eight geometries / extinction if simply connected | **discrete list** |

**The continuous is allowed. The decision is an integer.** Geometry ripples freely, but the name of the space is called by the countable side: eight, a winding number, a number of surgeries. What carried the space to the place where its name is read was averaging — addition.

## 8. The disappointments, reread

The second disappointment — geometry instead of topology — looks inevitable in this light. Perelman exposed the continuous process running **underneath** the ledger, and showed that this process always lands on a discrete shelf. Heat and entropy are not decoration in the proof. **The decision lives on the integer side, but the carrier is addition**, and its official names are heat and entropy.

The third disappointment has a structure too. The proof speaks discrete (topology) and continuous (analysis and geometry) at once, and its safety device is imported from statistical mechanics. A native speaker of any one of the three hears the other two as foreign languages. What was hard was crossing three language borders — and crossing them was what it took.

---

## Part 2 — Entropy as a ledger that only climbs

Part 1 ended with a pair: a continuous flow that is cut a **whole number of times**, and an entropy that **never goes down**. Here is a device built only from addition that has both.

## 9. The device: firing by addition

Advance the Fibonacci numbers by addition:

{% katex %}
r_{n+1} = r_n + r_{n-1}
{% endkatex %}

and read them on an {% katex inline %}m{% endkatex %}-hour clock. If the sum reaches {% katex inline %}m{% endkatex %}, subtract {% katex inline %}m{% endkatex %} **once**.

Both {% katex inline %}r_n{% endkatex %} and {% katex inline %}r_{n-1}{% endkatex %} are below {% katex inline %}m{% endkatex %}, so their sum is below {% katex inline %}2m{% endkatex %}. One subtraction is always enough. Whether it happened is a clean 0 or 1:

{% katex %}
f_n \in \{0, 1\}
{% endkatex %}

Call this a **firing**. It is not a metaphor; it is the fact that one subtraction took place. On a 13-hour clock, one full turn looks like this:

```plaintext
residue  1  1  2  3  5  8  0  8  8  3 11  1 12  0
        12 12 11 10  8  5  0  5  5 10  2 12  1  0

firing   0  0  0  0  1  0  0  1  0  1  0  1  0  0
         1  1  1  1  1  0  0  0  1  0  1  1  0  0
```

The period is 28, and there are 12 firings in it. After 14 steps the residues come back with the sign flipped ({% katex inline %}F_{n+14} \equiv -F_n \bmod 13{% endkatex %}); after 28 they are back where they started.

## 10. The ledger

Add up the firings:

{% katex %}
L_n = \sum_{k \le n} f_k
{% endkatex %}

{% katex inline %}L_n{% endkatex %} is an integer. It never decreases. It climbs one whole step at a time. **Inside this device, this is what I call entropy.** Nothing smooth was discretized here; there were only steps from the start.

The ledger has a closed form. Write the step with the subtraction showing — {% katex inline %}r_{n+2} = r_{n+1} + r_n - m f_n{% endkatex %} — and add it up over one period. The sum of residues {% katex inline %}S{% endkatex %} appears once on the left and twice on the right, so:

{% katex %}
S = S + S - mF \quad\Longrightarrow\quad m \times F = S
{% endkatex %}

where {% katex inline %}F{% endkatex %} is the number of firings per period. It is an integer identity, checked without a single division: it holds for every {% katex inline %}m{% endkatex %} from 2 to 399. For {% katex inline %}m = 13{% endkatex %}: {% katex inline %}S = 156{% endkatex %}, {% katex inline %}F = 12{% endkatex %}, and {% katex inline %}13 \times 12 = 156{% endkatex %}.


<!-- Figure 2: perelman_firing_ledger.gif -->
![Fibonacci residues stepping around a 13-hour clock; each subtraction lights up and the ledger climbs to 12](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/jdh6m4e17rwg4rlp0poi.gif)


*Each step adds the two previous residues. When the sum reaches 13, one subtraction happens (FIRE) and the ledger on the right climbs one step. After one full turn of 28 steps it reads 12, and 13 × 12 = 156 is the sum of the residues.*


*Interactive version (integer arithmetic only): [CodePen — perelman_firing_ledger.html](https://codepen.io/editor/morcb13-bit/pen/01a1fb7a-6ab5-7508-bafe-43b24e04637d)*


## 11. Change the clock

Now read **the same sequence** on different clocks.

| clock {% katex inline %}m{% endkatex %} | period | firings {% katex inline %}F{% endkatex %} | residue sum {% katex inline %}S{% endkatex %} | {% katex inline %}m \times F{% endkatex %} |
|---|---|---|---|---|
| 7 | 16 | 7 | 49 | 49 |
| 11 | 10 | 3 | 33 | 33 |
| **13** | **28** | **12** | **156** | **156** |
| 19 | 18 | 7 | 133 | 133 |
| 23 | 48 | 23 | 529 | 529 |
| 29 | 14 | 4 | 116 | 116 |
| 41 | 40 | 19 | 779 | 779 |
| 61 | 60 | 28 | 1708 | 1708 |

Over the clocks {% katex inline %}m = 2{% endkatex %} to {% katex inline %}31{% endkatex %}, the firing count takes 18 different values.

Meanwhile **the integers themselves have not changed at all.** {% katex inline %}1, 1, 2, 3, 5, 8, 13, 21, 34, 55, \dots{% endkatex %} is the same row whichever clock you read it on.


<!-- Figure 3: perelman_four_clocks.gif -->
![Four ledger staircases for clocks 11, 13, 23 and 29 climbing at different rates from the same Fibonacci row](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/fwpangkg9jk16w6idt48.gif)


*The row of integers at the top never changes. Read on four clocks, it produces four different staircases. The tick marks on the bottom axis show where each clock completes a period.*


*Interactive version (integer arithmetic only): [CodePen — perelman_four_clocks.html]
(https://codepen.io/editor/morcb13-bit/pen/01a1fb87-239e-726f-a97a-cfa396bbee27)*


So the number of firings is **not a property of the sequence. It belongs to the clock** — to the address system you read it in.

## 12. Crossing zero by turning the sign over

Go back to the round sphere of Section 3. Counted in the radius {% katex inline %}r{% endkatex %}, the story ends at {% katex inline %}t = r_0^2/4{% endkatex %}: a radius cannot be negative, so zero is a wall. Counted in {% katex inline %}u = r^2{% endkatex %}, the same motion is the straight line {% katex inline %}u(t) = r_0^2 - 4t{% endkatex %}, and nothing in that line stops at zero. It simply carries on into negative values. Whether zero is an end or a point where the sign turns over depends on which variable you count in.

The neck does the same. For a rotationally symmetric neck of radius {% katex inline %}\psi(s,t){% endkatex %}, the flow reads

{% katex %}
\psi_t = \psi_{ss} - (n-1)\,\frac{1-\psi_s^2}{\psi}
{% endkatex %}

The {% katex inline %}1/\psi{% endkatex %} is what breaks when the neck pinches. Switch to {% katex inline %}u = \psi^2{% endkatex %}. For {% katex inline %}n = 2{% endkatex %} everything collapses to

{% katex %}
u_t = u_{ss} - 2
{% endkatex %}

Here {% katex inline %}n{% endkatex %} is the dimension of the neck's cross-section. Perelman's space is three-dimensional and its necks look like {% katex inline %}S^2 \times \mathbb{R}{% endkatex %}, so {% katex inline %}n = 2{% endkatex %}: this is his case, not a toy. No division is left, and {% katex inline %}u = 0{% endkatex %} is not forbidden. I checked the substitution symbolically, and also in the computed dumbbell of Figure 1: at the neck, {% katex inline %}du/dt{% endkatex %} and {% katex inline %}u_{ss} - 2{% endkatex %} agree at every stage I checked ({% katex inline %}-1.6248{% endkatex %} and {% katex inline %}-1.6248{% endkatex %} at {% katex inline %}h = 0.10{% endkatex %}; {% katex inline %}-1.8366{% endkatex %} and {% katex inline %}-1.8347{% endkatex %} at {% katex inline %}h = 0.013{% endkatex %}). That is why the neck's {% katex inline %}u{% endkatex %} in the inset falls almost in a straight line, with a slope approaching {% katex inline %}-2{% endkatex %}.

The device of Part 2 does this on every turn. On the 13-hour clock the residue reaches 0 at the seventh term, and the sequence does not stop there: 8, 0, 8, 8, … After 14 steps every residue has come back as its negative (12 is −1 on this clock), and after 28 it is positive again. Zero is a place the sequence passes through, and the sign keeps the information that would otherwise be lost there. On the real line, the Fibonacci ratios never arrive at {% katex inline %}\varphi{% endkatex %}; on the clock the sequence closes, carrying a sign.

In an earlier piece I drew this as a dial {% katex inline %}\mu{% endkatex %} stepping through the integers {% katex inline %}-4, \dots, +4{% endkatex %}: a torus whose hole shrinks, closes into a disc at {% katex inline %}\mu = 0{% endkatex %}, and comes out the other side turned over — nine integer steps, and no cut. Read in {% katex inline %}r{% endkatex %}, the same passage needs one surgery at zero. Read in {% katex inline %}\mu = \pm r^2{% endkatex %}, it needs none.


<!-- Figure 4: perelman_mu_dial.gif -->
![A wireframe torus stepping through the integers mu = -4 to +4: the hole shrinks, closes into a disc at zero, and the surface comes out the other side turned over](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/6fy1yjd1eewgr20r6r2z.gif)


*The dial {% katex inline %}\mu{% endkatex %} moves only through whole numbers, from −4 to +4 and back. Below zero the hole has radius {% katex inline %}\sqrt{-\mu}{% endkatex %}; at {% katex inline %}\mu = 0{% endkatex %} it closes (red); above zero the centre opens to a thickness {% katex inline %}2\sqrt{\mu}{% endkatex %}. Nothing is cut anywhere along the way.*


*Interactive version (drag to turn it): [CodePen — perelman_mu_dial.html](https://codepen.io/editor/morcb13-bit/pen/01a1fb89-e4c1-7695-9d13-842462e12b3f)*


So with discrete values, the limit is crossed by turning the sign over instead of by cutting. This is the sharpest form of Section 11: **the number of cuts belongs to the reading, and in a suitable reading it can be zero.**

How far this reaches is a separate question. The division has moved rather than vanished: the length along the neck evolves by {% katex inline %}\varphi_t = 2(\psi_{ss}/\psi)\,\varphi{% endkatex %}, and in {% katex inline %}u{% endkatex %} that factor reads {% katex inline %}\frac{u_{ss}}{2u} - \frac{u_s^2}{4u^2}{% endkatex %} (checked symbolically). And crossing zero in a variable is not yet the same as the geometry flowing on: what lives at {% katex inline %}u < 0{% endkatex %} is not decided by the equation.

## 13. Back to Perelman (my labelling, not a result)

From here on I am labelling borrowed knowledge, not deriving anything.

- **Surgery**: cut once, at the instant a threshold is crossed. A discrete event — like a firing.
- **Monotone {% katex inline %}\mathcal{W}{% endkatex %}**: the guarantee that the cuts stay finite. A ledger that only climbs.

In Perelman's proof, entropy did not earn its place by growing smoothly. It earned it by **guaranteeing that the number of cuts is a finite whole number.** The device in Part 2 has the same two pieces: a 0/1 event and a monotone integer ledger. What it shows by computation is narrower: in the device, the count belongs to the reading, not to the thing being read.

## 14. One step further, stated softly

If the count of cuts belongs to the clock, then "how many times was it cut?" is partly a question about how you read it. The old question of whether the universe is continuous or discrete seems to have the same shape: is continuity out there, or is it on the side of the reading? "Nature knows only addition" narrows that question to a single point — can it be written with addition alone?

Several famous open problems also have the shape "a continuous process, landing on a discrete shelf". The Poincaré conjecture was actually solved in that shape: eight geometries, finitely many surgeries. That is an observer's arrangement of known facts; the device above does not reach that far.

## 15. What I have not checked

- **Whether firings and surgeries are the same quantity.** Both are one-off events at a threshold, both counted by a monotone ledger. That is a match in form. Turning it into a claim would require recounting surgeries by addition alone, and I have not done that.
- **Whether the sign flip removes surgery in Perelman's 3-dimensional flow.** It removes the division from the radius equation in Perelman's own case ({% katex inline %}n = 2{% endkatex %}), but the length equation still divides by {% katex inline %}u{% endkatex %}, and I have not run the flow through {% katex inline %}u = 0{% endkatex %}.
- **The flow after surgery.** Figure 1 stops at the cut; how the two horned pieces flow on is not computed.
- **Whether the universe is discrete.** The device does not touch that. Section 14 is about the shape of the question, not its answer.

---

## Tests

Each test was written down before the run, together with what a failure would mean. Two of them are negative controls, there to show the checks can say no.

| test | content | result |
|---|---|---|
| E1 | one subtraction is always enough, so firing is 0/1 ({% katex inline %}m = 2..399{% endkatex %}) | OK |
| E2 | {% katex inline %}m \times F = S{% endkatex %} ({% katex inline %}m = 2..399{% endkatex %}) | OK |
| E3 | changing the clock changes {% katex inline %}F{% endkatex %} (18 values for {% katex inline %}m = 2..31{% endkatex %}) | OK |
| E5 | control: break the rule to {% katex inline %}r_{n+1} = r_n + 2r_{n-1}{% endkatex %} → one subtraction is no longer enough | fails as it should (31 of 43 clocks) |
| E6 | clock 13: period 28, sign flips after 14 | OK |
| N1 | neck pinch: grid 1601 vs 3201 at {% katex inline %}h = 0.01{% endkatex %} (slope within {% katex inline %}\pm 4h{% endkatex %}; length of {% katex inline %}\psi < 2h{% endkatex %}) | 0.2226 / 0.2236; 12.3h / 11.8h |
| N2 | at the neck, {% katex inline %}du/dt = u_{ss} - 2{% endkatex %} in the computed flow | agrees to 3–4 digits at three stages |
| E7 | is there a ceiling of 1/2 on firings per step (integer test {% katex inline %}2F > \text{period}{% endkatex %})? | **NG** — exceeded 73 times below {% katex inline %}m = 2000{% endkatex %} (first at {% katex inline %}m = 148{% endkatex %}) |

## What is known, and what this adds

**Known:** Ricci flow (Hamilton 1982); the {% katex inline %}\mathcal{W}{% endkatex %}-entropy, no local collapsing, κ-solutions and surgery (Perelman 2002–03); finite extinction for simply connected spaces (Perelman's third paper; Colding–Minicozzi); geometrization (Thurston's conjecture, a theorem through Perelman); the higher dimensions (Smale, Freedman); Pisano periods of Fibonacci numbers.

**What this adds:** the arrangement (addition carries, integers decide), and the device in Part 2 — firings, a monotone integer ledger, the identity {% katex inline %}m \times F = S{% endkatex %}, the measurement that the count changes with the clock while the sequence does not, and the reading of zero as a sign flip ({% katex inline %}u = r^2{% endkatex %}, the half-turn of the 13-hour clock) instead of an end that has to be cut.

## The code

Each figure is a single HTML page (linked under the figure); the dumbbell page computes the flow itself, and the GIFs here were recorded frame by frame from those same pages. The neck computation is also in Python, `neck2.py` ([GitHub](PASTE_GITHUB_URL)). The ledger needs only this:

Standard library only, integers only.

```python
def ledger(m):
    a, b, r, f = 1, 1, [1, 1], []
    while True:
        s = a + b
        fired = s >= m
        if fired:
            s -= m              # one subtraction is always enough
        f.append(int(fired))
        a, b = b, s
        r.append(s)
        if (a, b) == (1, 1):
            break
    r = r[:-2]                  # one full period of residues
    return len(r), sum(f), sum(r)

for m in (7, 11, 13, 19, 23, 29, 41, 61):
    P, F, S = ledger(m)
    print(m, P, F, S, m * F == S)

print(all(m * ledger(m)[1] == ledger(m)[2] for m in range(2, 400)))
```

---

## Summary

Perelman's proof answered an integer question — is the winding ledger blank? — by letting heat run through the shape. The heat is averaging, which is addition; the entropy forbids collapse and keeps the number of surgeries finite; and the answer arrives as something you can count: eight geometries, or complete extinction for the sphere.

A device made only of addition shows the same two pieces in miniature. A firing is one subtraction, a clean 0 or 1. Its running total is a ledger that only climbs, one whole step at a time, and {% katex inline %}m \times F = S{% endkatex %} holds exactly. Read the same Fibonacci integers on a different clock and the ledger changes, while the integers do not. Count in {% katex inline %}u = r^2{% endkatex %} instead of {% katex inline %}r{% endkatex %}, and zero stops being a wall: the sign turns over and the line goes on.

**Heat carried the shape. The integers decided its name. And the count of cuts belongs to the clock you read it on — with discrete values, the limit can be crossed by turning the sign over.**
