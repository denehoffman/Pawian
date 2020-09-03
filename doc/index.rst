Welcome to PAWIAN
=================

.. note::

  These instructions were extracted `from the P̅ANDA wiki
  <https://panda-wiki.gsi.de/foswiki/bin/view/PWA/PawianPwaSoftware>`_.

PAWIAN (PArtial Wave Interactive ANalysis Software) is a partial wave analysis
software package with the aim to support analyses dedicated to all physics
cases at P̅ANDA and partly to specific items at other hadron spectroscopy
experiments. Real data originated from the pbar p annihilation processes and
from the e+ e- reactions are currently under investigation with PAWIAN.
Essential parts of the powerful, user-friendly and highly modular software are
already available. It is possible to define complicated decay trees, to choose
different spin formalisms and dynamics as well as to set up numerous other
parameters via a simple configuration file without modifying and compiling any
code. The key features of PAWIAN which can be setup via such a configuration
file with a user-friendly syntax are as follows:

* The amplitudes can be described with different formalisms. The widely-used
  helicity, canonical and the Lorentz-invariant Rarita-Schwinger formalisms are
  supported so far.

* A couple of different descriptions of the dynamics can be chosen which are
  among other things the Breit-Wigner parametrization with or without barrier
  factors, the Flatté-formalism and the K-matrix formalism.

* The minimization is realized with an event-based maximum likelihood fit. It
  makes use of the MINUIT2 minimization package. In addition one can start with
  a pre-fitter based on generic algorithms.

* To improve the performance for the very time-consuming fit procedure the
  applications can be run in parallel on multi-core CPUs and/or on computer
  clusters.

* Coupled channel analyses are supported.

* The analyst can generate events based on the user-defined decay model or on
  the fit result obtained with real data.

* Tools for quality assurance (QA), histogramming and for extracting values of
  different goodness-of-fit criteria are also available.

.. toctree::
  :maxdepth: 2
  :caption: Table of Contents

  install
  usage
  contribute
