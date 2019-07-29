# StochasticBurst
This is the code associated with my article "Timing and Shape of Stochastic Autocatalytic Burst Formation"
A link to the article will follow shortly, once uploaded.

As a very brief and light introduction to the code:
* The Various "back-forward" solvers are designed to solve the PDE problems described in section four of the paper.
* DerpyNonLinSimulation and DerpyNonLinSimulationIterate are used to produce the simulations of the non-linear equations which I compare my analytic results to in section 5.2. In particular, DerpyNonLinSimulation simulates once, and was used to make figure 10 (right), and DerpyNonLinSimulationIterate simulates many times, and was used to create figure 9, right. The simulations are derpy simply because they use a simplistic finite element method with normal valued noise.
* LinearUnstable1dFourierRepresentation and LinearUnstable1dFourierRepresentationRepeat were used to make the left hand side of figures 9 and 10. They use fourier representations of the system, and are thus much faster, but can not handle non-linearity.
*MoriModelInfDiffuseForPaper and CalculateMoriParamsGivenStartingValue were used to create the simulations in Section 2, figure 3. Please note that because Noise does not work in 2d, these simulations are sensitive to mesh size.


All figure numbers currently refer to the ArXiv article, and may need to be updated depending how ArXiv numbers things.
