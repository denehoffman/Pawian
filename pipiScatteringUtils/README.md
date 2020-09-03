# pipi scattering utilities

This directory contains the specific tools for the analysis for pi+ pi-
scattering data. It should be possible to consider all two particle final
states induces from pi pi reactions like pi+ pi-, K+ K-, eta eta, eta etaprime,
etc. In addition one has to specify well defined intermediate states ( in
general defined by IG(JPC), like the (pipi)S wave with 0^+(0^++) or the f_2
wave with 0^+(2^++) or the rho-wave with 1^+(1^--) ). The supported dynamics is
so far only the T-Matrix representation. These scattering reaction tools are
also embadded in the general PAWIAN PWA framework. The main differences
compared to all other reactions are treated by the implementation of specific
virtual methods and are as follows:

1. no production amplitude is required here
2. exactly one decay channel with the specific intermediate state and the two
   final state particles
3. the T-matrix dynamics is required and part of the decay amplitude
4. no decay angular distributions required
5. measurements are a scan of sqrt(s) or invariant mass of the final state
   paricles, the phase motion and the inelasticity including errors
6. the measuremet at one defined mass is technically assosiated with one event
7. log likelihood fit is replaced by a chi^2 fit
8. each individual channel and intermediate state is treated as an individal
   application. For combining it, one has to perform "coupled channel fits" in
   the server/client mode
