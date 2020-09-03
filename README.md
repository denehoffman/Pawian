# PAWIAN (**PA**rtial **W**ave **I**nteractive **AN**alysis)

**PAWIAN** is a powerful, user-friendly and highly modular partial wave
analysis software package with the aim to support analyses for a multitude of
different physics cases at hadron physics experiments. Real data originated
from the pbar p annihilation process and from e+ e- reactions are currently
under investigation with **PAWIAN**.

It is possible to define complicated decay trees, to choose different spin
formalisms and resonance dynamics as well as to set up numerous other
parameters via a simple plain text configuration file without modifying and
compiling any code. Parameters that can be set in the configuration file
include, among many more:

- The amplitudes can be described with different formalisms. The widely-used
  helicity, canonical and the Lorentz-invariant Rarita-Schwinger formalisms are
  supported so far
- A couple of different descriptions for the dynamics of resonances can be
  chosen, e.g. the Breit-Wigner parametrization with or without barrier
  factors, the Flatté-formalism and the K-matrix formalism
- The minimization is realized with an event-based maximum likelihood fit. It
  makes use of the MINUIT2 minimization package. In addition one can start with
  a pre-fitter based on a genetic algorithm
- To improve the performance for the very time-consuming fit procedure the
  applications can be run in parallel on multi-core CPUs (multi-threading)
  and/or on computer clusters (networking)
- Coupled channel analyses are supported
- The analyst can generate events based on the user-defined decay model or on
  the fit result obtained with real data
- Tools for quality assurance, histogramming and for extracting values of
  different goodness-of-fit criteria are also available

**PAWIAN** is written in C++ and follows an object-oriented approach with a
wide range of flexibility. The code therefore allows to be easily extended for
further decay models, complete amplitudes and also another descriptions for the
dynamics.

**PAWIAN** and all its associated files are licensed under the GNU General
Public License v3 (GPLv3), which you can find
[here](https://www.gnu.org/licenses/gpl.txt) or in the `COPYING` file in the
**PAWIAN** root directory.

**For information about the dependencies, the installation procedure, and
usage, see the
[Wiki pages](https://panda-wiki.gsi.de/foswiki/bin/view/PWA/PawianPwaSoftware).
Issues can be discussed through the
[PawianTalk user forum](https://pawiantalk.ep1.rub.de/).**
