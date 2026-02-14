/**
 * ClassicalDiracElectron
 *
 * This project simulates a classical Dirac electron in time, following the
 * Barut–Zanghi model (Phys. Rev. Lett. 52, 2009 (1984)).
 *
 * ---
 *
 * CDE_4thRK.c — Main simulator
 *
 *   Evolves the electron's position x^μ, momentum p^μ, and spinor z = z_r + i z_i
 *   (4 complex components). Time integration uses 4th-order Runge–Kutta.
 *
 *   The electromagnetic field is set via ConstantEB() (constant E and B) or
 *   OneOverr2E() (Coulomb-like 1/r^2 field). Parameters at the top of the file
 *   include λ, charge q, E_x, E_z, B_z, time step dt0, final time T, and initial
 *   angles theta[]. Output is written to a file such as CE_trajectory_RK_2_.dat.
 *
 * ---
 *
 * Plotter.c — Post-processing and plotting
 *
 *   Writes a gnuplot script to plot 3D trajectories from the .dat files
 *   (columns 4, 5, 6 are the spatial coordinates). It also reads an "out1" file,
 *   bins radii r = sqrt(x^2 + y^2) into a histogram, and prints bin centres,
 *   counts, and a cumulative sum for radial distribution analysis.
 *
 * ---
 *
 * Theory reference:
 *
 *   A. O. Barut, Nino Zanghi
 *   Phys. Rev. Lett. 52, 2009–2012 (1984)
 *   Classical Model of the Dirac Electron
 *   http://prl.aps.org/abstract/PRL/v52/i23/p2009_1
 *
 **/
