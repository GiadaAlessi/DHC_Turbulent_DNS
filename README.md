# High-Rayleigh Simulation of a Rectangular Differentially Heated Cavity

This project presents a 2D Direct Numerical Simulation (DNS) of a rectangular cavity heated differentially, developed as part of a Master's thesis in Thermal Engineering. The cavity features a high aspect ratio (1:4) and is subject to high Rayleigh numbers up to \( \mathrm{Ra} \sim 10^9 \), in order to investigate the emergence of quasi-turbulent behavior in buoyancy-driven natural convection.

## Physical Problem

- **Geometry**: Rectangular cavity with \( H = 4L \)
- **Boundary Conditions**:
  - Left wall: \( T_\text{hot} = 1 \)
  - Right wall: \( T_\text{cold} = 0 \)
  - Top & Bottom: adiabatic
- **Fluid**: Air (Pr = 0.71)
- **Rayleigh numbers tested**: \( 10^8 \), \( 6.4 \times 10^8 \), \( 1.07 \times 10^9 \)

## Numerical Method

- **Discretization**: Finite Volume Method (FVM) on staggered mesh  
- **Time Integration**: Adams-Bashforth (2nd order explicit)  
- **Pressure-Velocity Coupling**: Fractional Step Method (FSM)  
- **Solver**: Gauss–Seidel for Poisson equation  
- **Convective Schemes**: UDS/CDS switching based on flow development  
- **Mesh**: Uniform, with Nx × (4·Nx) cells for square CVs  

# Post-Processing

All post-processing was done in **Python**, using the posted script.

It generates:
- Time-averaged colormaps for \( T(x, y) \), \( u(x, y) \), \( v(x, y) \)
- Streamfunction \( \psi(x, y) \) from time-averaged \( u, v \)
- Nusselt number and velocity profiles at mid-sections
- Kinetic energy budget breakdown
- TKE and Reynolds stresses fields

## Quantities of Interest

- Time-averaged temperature and velocity fields  
- Streamfunction contours and streamlines  
- Average Nusselt number \( \overline{Nu} \)  
- Maximum \( u \) and \( v \) values and locations  
- Profiles at mid-height and mid-width  
- Global and local kinetic energy balance  
- Turbulent kinetic energy (TKE)  
- Reynolds stress components  

## References

- Trias et al., “Direct numerical simulations of the differentially heated cavity at high Ra”, *Int. J. Heat Mass Transfer*, 2007.  
- G. de Vahl Davis, "Natural Convection of Air in a Square Cavity", *Int. J. Numer. Methods Fluids*, 1983.  
- Bejan, “Convection Heat Transfer”, 4th Edition.

---

*Developed by Giada Alessi for the Master's Degree in Thermal Engineering (UPC Barcelona, 2025).*

