Set-up of the dynamics
======================

.. note::
  See original Wiki page `here
  <https://panda-wiki.gsi.de/foswiki/bin/view/PWA/PawianDynamics>`__.

Several different descriptions for the dynamics of the decay of a resonance are
supported. It has to be set-up in the main configuration-file with the option
"addDynamics" followed by in general one argument. For specific descriptions
more than one argument is needed (e.g. for the Flatte formalism) or even a
separate configuration file (for the generic :math:`K`-matrix description) is
mandatory. The user can specify between:


* :option:`addDynamics particleName BreitWigner <addDynamics>` |br|
  a simple Breit-Wigner parametrization

* :option:`addDynamics particleName BreitWignerRel <addDynamics>` |br|
  a relativistic Breit-Wigner parametrization

* :option:`addDynamics particleName BreitWignerBlattWRel
  {qR-factor} <addDynamics>` |br|
  a relativistic Breit-Wigner parametrization with Blatt-Weisskopf
  centrifugal-barrier factors and optional an argument for setting th​e
  qR-barrier constant (default: 0.197)

* :option:`addDynamics particleName Voigt <addDynamics>` |br|
  a convolution of a Lorentz profile and a Gaussian distribution which is
  suitable for very narrow resonances

* :option:`addDynamics particleName Flatte {particle1}
  {particle2} <addDynamics>` |br|
  a generic Flatte parametrization. The second decay channel has to be defined
  by the specification of :command:`particle1` and :command:`particle2`.

* :option:`addDynamics particleName K0star1430Lass <addDynamics>` |br|
  the LASS parametrization for the :math:`K^{0*}(1430)`

* :option:`addDynamics KMatrixParticleName PiPiSWaveAS <addDynamics>` |br|
  (pi pi)S-wave parametrization by Anisovich and Sarantsev ("Eur. Phys. J. A16,
  229(2013)") for the mass range up to 1900 MeV which contains 5 different
  poles and 5 decay channels

* :option:`addDynamics KMatrixParticleName KpiSWaveIso12 <addDynamics>` |br|
  :math:`(K pi)S`-wave parametrization for the isospin 1/2 component used by
  the FOCUS experiment ("Phys. Lett. B653 (2007) 1-11").

* :option:`addDynamics particleName KpiSWaveIso32 <addDynamics>` |br|
  :math:`(K pi)S`-wave parametrization for the isospin 3/2 component used by
  the FOCUS experiment ("Phys. Lett. B653 (2007) 1-11").

* :option:`addDynamics particleName BlattWBarrier {qR-factor} <addDynamics>`
  |br|
  simple Blatt-Weisskopf factors without any description of the mass shape; the
  argument for setting th​e qR-barrier constant is optional (default: 0.197).

* :option:`addDynamics KMatrixParticleName KMatrix {name of the K-matrix
  configuration file} <addDynamics>` |br|
  generic :math:`K`-matrix defined by a separate configuration file. A summary
  about the options and different arguments are listed under
  :doc:`/usage/arguments`.
