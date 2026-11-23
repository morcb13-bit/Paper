---
title: "BepiColombo: Two Eyes on Mercury — MPO and Mio, After Nine Flybys"
published: false
description: "BepiColombo entered orbit around Mercury on 21 November 2026. What the two orbiters will do, why there are two, and how nine flybys got them there, with animated figures and a small integer-only orbit check."
tags: space, astronomy, simulation, python
---

On 21 November 2026, BepiColombo, the joint mission of the European Space Agency (ESA) and the Japan Aerospace Exploration Agency (JAXA), entered orbit around Mercury. It was launched on 20 October 2018, so the trip took a little over eight years.

What is circling Mercury right now is a single stack of two spacecraft: ESA's Mercury Planetary Orbiter (MPO) and JAXA's Mercury Magnetospheric Orbiter, nicknamed **Mio**. They are due to separate on 9–10 December and move into orbits of their own. It will be the first time two spacecraft orbit Mercury at the same time.

This post looks at what each orbiter is going to study, why the mission needs two of them, and what it took to get there, with two animated figures along the way.

## Mercury is close, and very hard to reach

Seen from Earth, Mercury sits right next to the Sun. By distance alone it is often closer to us than Mars. Even so, it is one of the hardest planets to send a spacecraft to.

A spacecraft leaving Earth speeds up as it falls toward the Sun. If it simply arrived at Mercury like that, it would be going far too fast for Mercury's weak gravity to capture it, and it would fly straight past. To get into orbit, the spacecraft has to match its speed around the Sun to Mercury's almost exactly.

To do that, BepiColombo made nine planetary flybys: one of Earth, two of Venus and six of Mercury, more than any planetary mission before it. Between flybys, its electric (ion) propulsion trimmed the orbit little by little.

The figure below shows the eight-year path, looking down on the Sun from ecliptic north. Grey lines are the orbits of Mercury, Venus and Earth, the coloured line is the spacecraft, and red rings mark the flybys.



*fig1:BepiColombo's eight-year path around the Sun*　(bepicolombo_path_en.gif)

![BepiColombo's eight-year path around the Sun]
(https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/w4scce8uehwd5p1dmss8.gif)

<!-- CodePen: bepicolombo_path_en.html (interactive version of the figure above) -->
{% codepen https://codepen.io/editor/morcb13-bit/pen/01a205b8-11bf-72a7-a76f-96da42be7625 %}


The positions are JPL Horizons values (the ephemeris service of NASA's Jet Propulsion Laboratory), taken every 20 days, interpolated to daily steps, and drawn together with the planets.

Here is how the spacecraft's orbit around the Sun changed across the flybys, computed from Horizons osculating elements:

| After | Perihelion (AU) | Aphelion (AU) | Period (days) |
|---|---|---|---|
| Launch | 0.887 | 1.165 | 379.5 |
| Earth | 0.624 | 1.013 | 270.3 |
| Venus 1 | 0.520 | 0.833 | 203.1 |
| Venus 2 | 0.335 | 0.731 | 142.1 |
| Mercury 1 | 0.321 | 0.695 | 132.2 |
| Mercury 2 | 0.304 | 0.632 | 116.8 |
| Mercury 3 | 0.306 | 0.595 | 110.4 |
| Mercury 4 | 0.306 | 0.470 | 88.3 |
| Mercury 5 | 0.307 | 0.467 | 87.9 |
| Mercury 6 | 0.308 | 0.456 | 86.2 |
| (Mercury itself) | 0.307 | 0.467 | 88.0 |

The Earth and Venus flybys pulled the perihelion down to Mercury's orbit. The six Mercury flybys then pulled the aphelion down, until the spacecraft's ellipse matched Mercury's.

People call this "braking", but the spacecraft did not get slower around the Sun. Its speed at perihelion actually rose from 33.7 km/s after launch to 60.2 km/s after the second Venus flyby. What it shed was the *difference* between its speed and Mercury's.

## Nine hyperbolas

Each flyby, seen from the planet, is a hyperbola. The spacecraft comes in from far away, gets its direction bent by the planet's gravity, and leaves again.


*fig2:The nine flyby hyperbolas, seen from each planet*　(bepicolombo_flybys_en.gif)
![The nine flyby hyperbolas, seen from each planet](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/cxyq8gc4g3gex1op0lq6.gif)


<!-- CodePen: bepicolombo_flybys_en.html (interactive version of the figure above) -->
{% codepen https://codepen.io/editor/morcb13-bit/pen/01a205c5-45cc-7e6f-a779-dfffa25427d4 %}


The bars at the bottom show the speed at which the spacecraft approached each planet (the relative speed "at infinity", v∞). Over the six Mercury flybys it drops: 6.62 → 6.33 → 3.64 → 2.79 → 2.66 → 1.86 km/s.

The slower and closer the pass, the more the path bends. The last one, Mercury 6 on 8 January 2025, passed 296 km above the surface and turned the spacecraft through 88.6°. Mercury 5, by contrast, passed far out at 37,626 km and turned it by only 8.3°. v∞ shrinking toward zero in this figure and the ellipse settling onto Mercury's in the first figure are the same event, seen from two places.

## Becoming a satellite of Mercury

Even after the last flyby (Mercury 6), the spacecraft did not immediately become a satellite of Mercury. For about two years, from January 2025 to November 2026, it rode almost the same ellipse as Mercury and travelled around the Sun side by side with it.

The first part of the next figure shows the spacecraft during those two years, seen from Mercury, turning with the planet so that the Sun always stays below (computed from Horizons values).


*fig3:The spacecraft seen from Mercury, then captured*　(bepicolombo_capture_en.gif)
![The spacecraft seen from Mercury, then captured](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/jtvystlkqqm2t2u3b647.gif)


<!-- CodePen: bepicolombo_capture_en.html (interactive version of the figure above) -->
{% codepen https://codepen.io/editor/morcb13-bit/pen/01a205cb-76f4-765e-9827-cd1e0b996651 %}


The two ellipses differ slightly in shape and orientation. So once every Mercury year (88 days) the spacecraft traces a loop around Mercury, coming closer and drifting away. Right after Mercury 6, the spacecraft's period was 86.2 days against Mercury's 88, so it gained about two days per lap, and the whole loop slowly shifted.

The distance from Mercury, half a year at a time:

| Period | Closest | Farthest |
|---|---|---|
| Early 2025 | 0.015 AU | 0.163 AU |
| Late 2025 | 0.096 AU | 0.220 AU |
| Early 2026 | 0.051 AU | 0.198 AU |
| Jul–Oct 2026 | 0.001 AU | 0.064 AU |

At its farthest it was more than 0.2 AU, over 30 million km, away. During this time the spacecraft was not orbiting Mercury; it was simply running alongside it around the Sun.

In 2026 the ion engine shrank the loops, and from July on the distance clearly closes. On 21 November, with its speed relative to Mercury close to zero, the spacecraft was captured by Mercury's gravity and entered orbit. (The distances come from 20-day values that were interpolated, so they carry an error of about 0.005 AU.)

The second part of the figure is a projection of what follows. Separation and the final orbits are still ahead, so they are drawn from the planned values. The orbit right after capture is not yet the one used for separation: a series of burns takes the stacked pair into Mio's orbit (590–11,640 km altitude, 9.3 hours), where they are due to separate on 9–10 December. The MPO will then use its own engine to lower itself to its orbit (480–1,500 km, 2.3 hours), planned for March 2027. Mercury and both orbits are drawn to scale using the planned values.

From now on, the two spacecraft orbit Mercury while travelling around the Sun together with it, the same way the Moon orbits Earth while going around the Sun with Earth. "Running alongside" before capture and "orbiting as a satellite" after it may look similar, but they are different states.

## MPO: eyes on the surface and the interior

The MPO, built by ESA, is three-axis stabilised and keeps the face carrying its instruments pointed at Mercury. It will fly a polar orbit 480 × 1,500 km above the surface, once every 2.3 hours, mapping the whole planet from close in.

Sunlight at Mercury is more than ten times stronger than at Earth, and the hot dayside surface reflects heat back up. The MPO dumps its heat through a 3.7 m radiator that always faces away from the Sun.

It carries 11 instrument suites. Grouped by role:

- **Topography and surface composition:** the BELA laser altimeter, the SIMBIO-SYS cameras and spectrometer, the MERTIS thermal-infrared spectrometer, the MIXS X-ray spectrometer and the MGNS gamma-ray and neutron spectrometer. Together they measure heights, minerals and elements.
- **Interior, gravity and relativity:** the MORE radio-science experiment, which tracks the orbit very precisely, and the ISA accelerometer. The gravity field reveals the interior, and the same measurements test Einstein's general relativity right next to the Sun.
- **Thin atmosphere and magnetic field:** the PHEBUS ultraviolet spectrometer, the SERENA particle package and the MPO-MAG magnetometer.
- **Watching the Sun:** SIXS, which monitors X-rays and particles coming from the Sun.

## Mio: eyes on the magnetic field and plasma

Mio was built by JAXA, and Japan's Institute of Space and Astronautical Science runs it from design through operations.

It is an octagonal prism that fits inside a 1.8 m circle and weighs about 255 kg. It spins once every 4 seconds (15 rpm), with its spin axis held at right angles to the Sun. Spinning spreads the Sun's heat evenly over the body, and lets the instruments sweep every direction in turn. At Mercury it will unfold four 15 m wire antennas and two 5 m masts.

Its orbit is a long, thin polar ellipse, 590 × 11,640 km, with a period of 9.3 hours. On every orbit it passes from close to the planet out to the boundary where Mercury's magnetic field meets the solar wind.

It has five instruments:

- **MGF** (magnetometer): Mercury's magnetic field and how it pushes against the solar wind.
- **MPPE** (plasma particles): ions and electrons in the magnetosphere, with seven sensors.
- **PWI** (plasma waves and radio): electric fields, plasma waves and radio waves.
- **MSASI** (sodium imager): the distribution and motion of sodium in Mercury's thin atmosphere.
- **MDM** (dust monitor): dust grains along Mercury's orbit.

Mio could not spin during the cruise, so it travelled wrapped in MOSIF, a sunshield built by ESA. MOSIF will be released on 16 December, after the two orbiters separate.

## Why two spacecraft?

Mercury's magnetic field is only about one percent as strong as Earth's, and because the planet is so close to the Sun, the solar wind presses hard on it. Earth's magnetosphere extends to about ten Earth radii on the Sun-facing side; Mercury's reaches only about one and a half Mercury radii. It is the most compressed magnetosphere of any planet in the Solar System.

It also changes fast: when the solar wind shifts, its shape can change within minutes. With one spacecraft, you cannot tell whether a change you measured is a difference between two *places* or a change over *time*.

So the MPO, low and close to the planet, and Mio, high and reaching the outer edge, measure at the same moment in different places. The MPO takes the surface and interior, Mio takes the magnetic field and plasma, and together their data let scientists read Mercury as a single system.

## Background: it starts with Bepi Colombo's calculation

The mission is named after Giuseppe "Bepi" Colombo (1920–1984), a mathematician and engineer at the University of Padua in Italy.

In 1974, NASA's Mariner 10 used a Venus flyby to reach Mercury. Colombo worked out that if the spacecraft were put on an orbit around the Sun with a period of 176 days, exactly twice Mercury's year, it would meet Mercury again and again. Thanks to that idea, Mariner 10 passed Mercury three times instead of once.

Colombo also explained why Mercury turns exactly three times on its axis for every two trips around the Sun (a 3:2 spin–orbit resonance). BepiColombo's nine flybys are a direct descendant of his calculation half a century ago.

Mariner 10 discovered that Mercury has a global magnetic field, about one percent of Earth's. NASA's MESSENGER, which orbited Mercury from 2011 to 2015, left even more puzzles behind:

- The centre of the magnetic field is shifted well to the north of the planet's centre.
- Despite sitting next to the Sun, Mercury has water ice in the floors of polar craters that never see sunlight.
- The surface holds more sulfur, potassium and other easily evaporated elements than expected.
- Shallow, rimless depressions called "hollows" are scattered across the surface, apparently formed as volatile material escaped.

BepiColombo is going to look at all of these with two pairs of eyes.

On the Japanese side, the story began with a Mercury exploration working group formed in 1997. ESA proposed a joint study in 1999, and in 2000 ESA selected BepiColombo as a Cornerstone (large) mission. Japan approved its part in 2003.

The name **Mio** was chosen in 2018 from 6,494 public entries. In Japanese, *mio* (澪) is the trail of water a ship leaves behind as it moves. The channel markers that guide ships along it are called *miotsukushi*, a word that classical Japanese poems use as a pun on *mi wo tsukushi*, "giving oneself completely". The name joins the image of a craft sailing through the stream of solar-wind plasma with the persistence of the team that spent decades building it.

## A quick orbit check on an ordinary computer

As a small side trip, I propagated BepiColombo's orbit for 20 days on my own computer. I used the B13 library from my series *Nature Knows Only Addition*: a small Python library in which multiplication, division and square roots are built from integer addition alone, using balanced base-13 digits.

The interval was 8 → 28 December 2024, between Mercury flybys 5 and 6. The orbital elements barely change over these 20 days, so the engine can be taken as off.

Before running anything, I fixed what would count as a pass:

- Start from the Horizons position and velocity on 8 December, include only the Sun's gravity, and take 320 leapfrog steps of 1/16 day each (half a velocity kick, a position drift, half a kick).
- **Pass** if the result lands within 0.0005 AU (about 75,000 km) of the Horizons position on 28 December.
- **Negative control:** the same run with gravity switched off must miss by more than 0.05 AU. If it does not, the comparison says nothing.

Positions are integers with 1 AU = 10¹², and every operation runs on balanced base-13 digit lists.

```python
def accel(P, GMS, Rg):
    r2 = A(A(M(P[0], P[0]), M(P[1], P[1])), M(P[2], P[2]))   # r^2  (scale S^2)
    R = isqrt(r2, Rg)                                         # r    (Newton iteration using divide only)
    r3 = M(R, r2)                                             # r^3
    den = M(r3, SIXTEEN)                                      # kick for 1/16 day
    return [neg(Q(M(M(GMS, c), S2), den)) for c in P], R
# M = b13_arith.mul, Q = b13_arith.divide, A = b13_balanced.add (all built from addition)
```

The output:

```plaintext
Position on 12-28, Horizons : [-0.395217, -0.057955, 0.020093]
Position on 12-28, B13 int  : [-0.395223, -0.057954, 0.020050]   run time 4.0 s
Difference (B13 - Horizons) : 0.000043 AU = 6501 km
Negative control, GM = 0    : 0.500059 AU
Integer vs floating-point   : 1.24e-10 AU
Pass (< 0.0005): True / control misses (> 0.05): True
```

Over those 20 days the spacecraft travels about 100° around the Sun, and the integer run lands within 6,500 km of where Horizons puts it. Most of what is left over probably comes from forces this simple run leaves out, such as the pull of Mercury and Venus and the pressure of sunlight. The integer arithmetic itself contributes about 10⁻¹⁰ AU, which is negligible. The whole run took about four seconds on an ordinary PC.

Real mission operations use far more detailed models, with the gravity of every planet and the spacecraft's attitude included. This is only the doorway. Still, for checking the skeleton of an orbit, you don't need a big machine: an addition-only toolkit on a laptop gets you there in seconds.

## The journey on a ladder of powers of 13

To finish, a look at the journey in B13 terms. Take the Planck length as address 0, and step up one address for every factor of 13 (13^k). Reading the size of the spacecraft's ellipse (its semi-major axis) on this ladder:

- Just after launch: address 41.275
- After Mercury 1: address 41.000
- After Mercury 6: address 40.889
- Mercury itself: address 40.895

Eight years and nine flybys amount to a move of about 0.4 on the ladder, crossing the boundary at 41 just once. On this cosmic ladder, Earth's orbit and Mercury's sit within neighbouring steps, and stepping down that one rung took eight years and nine encounters.

## Closing

The MPO will watch Mercury's surface and interior from a low orbit; Mio will watch the magnetic field and plasma from a high one. Measuring at the same time from different places will let scientists separate differences in space from changes in time inside the small, tightly squeezed magnetosphere right next to the Sun.

Getting there took eight years of a strategy Bepi Colombo showed half a century ago: meet the planet again and again, and match its speed a little more each time. The two orbiters separate on 9–10 December, and full science operations begin in April 2027.

---

*A Japanese version of this article is on Zenn.*

**Data and sources**

- Spacecraft positions and orbital elements: JPL Horizons (spacecraft −121; Sun-centred and planet-centred, ecliptic J2000)
- ESA, "Latest updates: BepiColombo's arrival at Mercury": https://www.esa.int/articles/about-us/directorates/science/bepicolombo/latest-updates-bepicolombo-s-arrival-at-mercury
- ESA, "Mercury Planetary Orbiter": https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/Mercury_Planetary_Orbiter
- ESA, "Mercury Magnetospheric Orbiter": https://m.esa.int/Our_Activities/Space_Science/BepiColombo/Mercury_Magnetospheric_Orbiter
- JAXA, Mio (MMO) fact sheet: https://mio.isas.jaxa.jp/topics/uploads/2026/MMO_FACT_SHEET_260603.pdf
- JAXA, "Mio" nickname announcement: https://mio.isas.jaxa.jp/en/topics/news/000022.html
- NASA, "45 Years Ago: Mariner 10 First to Explore Mercury": https://www.nasa.gov/solar-system/45-years-ago-mariner-10-first-to-explore-mercury/
- ESA, "Giuseppe 'Bepi' Colombo (1920–1984)": https://www.esa.int/spaceinimages/Images/2002/01/Giuseppe_Bepi_Colombo_1920-1984
- Benkhoff et al., "BepiColombo – Mission Overview and Science Goals", *Space Science Reviews* 217:90 (2021)
