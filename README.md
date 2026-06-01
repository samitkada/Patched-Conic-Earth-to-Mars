## Physics & Engineering Models

This project simulates an Earth-to-Mars interplanetary mission using the patched conic approximation, a classical astrodynamics technique commonly used during preliminary spacecraft mission design. The trajectory is divided into three gravitational environments: Earth departure, heliocentric transfer around the Sun, and Mars arrival. By treating each region independently and connecting them at planetary spheres of influence, the simulator produces realistic first-order estimates of mission trajectory, transfer time, and propulsion requirements.

### Patched Conic Approximation

The patched conic method assumes that only one celestial body dominates the spacecraft's motion at a given time. Near Earth, the spacecraft is influenced primarily by Earth's gravity. During interplanetary cruise, the Sun becomes the dominant gravitational body. Upon arrival, Mars becomes the primary influence. This approach significantly simplifies trajectory design while maintaining good accuracy for early mission planning.

### Heliocentric Transfer Orbit

The transfer between Earth and Mars is modeled as a minimum-energy Hohmann transfer orbit. This transfer ellipse connects the nearly circular orbits of Earth and Mars around the Sun and provides an efficient trajectory requiring relatively low propulsion energy. The simulation calculates transfer geometry, spacecraft velocity changes throughout the transfer, and total time of flight.

### Hyperbolic Earth Departure

The mission begins in a low Earth parking orbit. A Trans-Mars Injection (TMI) maneuver increases spacecraft velocity enough to escape Earth's gravitational influence and enter the heliocentric transfer orbit. The simulator models this departure as a hyperbolic escape trajectory and computes the required departure delta-V and characteristic launch energy.

### Hyperbolic Mars Arrival

As the spacecraft approaches Mars, it enters the planet's sphere of influence with excess velocity relative to Mars. A Mars Orbit Insertion (MOI) burn is then performed to reduce spacecraft energy and achieve orbital capture. The arrival phase is modeled as a hyperbolic approach trajectory, allowing calculation of capture requirements and insertion delta-V.

### Orbital Mechanics Modeling

The heliocentric transfer trajectory is propagated using classical two-body orbital mechanics. Spacecraft motion is governed by the Sun's gravitational field during the cruise phase, while planetary gravitational fields govern departure and arrival phases. Numerical integration is performed using MATLAB's adaptive ODE45 solver to generate a complete Earth-to-Mars trajectory.

### Mission Performance Analysis

The simulation evaluates several key mission design parameters, including:

- Earth departure hyperbolic excess velocity
- Mars arrival hyperbolic excess velocity
- Characteristic launch energy (C3)
- Trans-Mars Injection (TMI) delta-V
- Mars Orbit Insertion (MOI) delta-V
- Total mission delta-V
- Transfer orbit geometry
- Earth-to-Mars time of flight

These metrics are fundamental to launch vehicle sizing, spacecraft propulsion design, and interplanetary mission architecture.

### Numerical Methods

The simulator combines analytical astrodynamics with numerical trajectory propagation, including:

- Patched conic trajectory design
- Hohmann transfer analysis
- Hyperbolic escape and capture modeling
- Two-body orbital propagation
- Adaptive Runge-Kutta numerical integration
- Three-dimensional mission visualization

This methodology mirrors the preliminary trajectory design process used in many historical and modern robotic Mars mission studies and provides a strong foundation for more advanced mission analysis techniques such as Lambert targeting, gravity assists, and low-thrust trajectory optimization.
