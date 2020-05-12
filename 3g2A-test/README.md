# 3g2A-test
In this folder the results of evaluating 3g2A helicity amplitudes are compared against the numerical results of NJet.

* Read_3g2A: the analytical expressions computed with FiniteFlow are re-written in a conveninet basis of functions. Bubble, triangle and box integrals are replaced with their known expressions and the resulting logs are combined into Lhat functions described in [this paper.](https://arxiv.org/abs/1605.02172)

* 3g2A_Comparison: the analytical expressions are compared with the numerical squared amplitudes computed by NJet. The Polylog functions coming from bubble and boxes integrals need to be defined and properly analytically continued in the complex plane. The helicity factor of the analytical amplitude is properly reconstructed.
