Parser arguments
================

.. note::
  See original Wiki page `here
  <https://panda-wiki.gsi.de/foswiki/bin/view/PWA/ParserArguments>`__.

Exactly one PWA application is provided for each individual reaction (e.g.
:code:`pbarpReactionApp` and :code:`epemReactionApp`). The choice of the
application mode, the definition of the decay trees, and a lot of other parser
options can be set by a configuration file. An example of such a configuration
file can be found in the tutorial. All supported parser options are listed
below.

Arguments required in the Command Line
--------------------------------------

.. option:: -c <path>, --config <path>

  Specify path to the configuration (:file:`.cfg`) file


.. option:: -C <path>, --coupledChannelConfigFile <path>

  Specify path to another configuration file for a coupled channel


Arguments to define the PWA model
---------------------------------

Commons
~~~~~~~

.. cspell:ignore noiso

.. option:: addDynamics {particle} {dynamic} [argument1 ...]

  Adds a dynamic amplitude to the decay of a particle. Arguments have to be
  added for specific dynamics. When the k-matrix formalism is used a
  configuration file has to be given.


.. option:: cloneParticle {particle1} {particle2}

  Create a new particle :command:`{particle2}` in the particle table that has
  the same properties as particle :command:`{particle1}`


.. option:: decay {formalism} [noiso] {particle1} To {particle2} {particle3}

  Adds the decay of particle1 to the amplitude using a certain spin formalism.
  The :command:`[noiso]` option is optional and disables isospin couplings


.. option:: finalStateParticle {name}

  Name of one final state particle. The order of final state particles in the
  configuration file has to be the same as in the data/monte carlo files


.. option:: production {particle1} {particle2}

  Adds the production of two particles to the amplitude


.. option:: replaceParamSuffix {string1} {string2}

  Replace a part of a parameter name :command:`{string1}` with
  :command:`{string2}`. This option can be used to merge parameters.


.. option:: replaceMassKey {particle name} {new name}

  Replace the name of a particle inside the Mass/Width parameter name. This
  option can be used to merge parameters.


.. option:: useProductionBarrier {bool} {double)

  Enable or disable barrier factors for the production system and optional an
  argument for setting th​e :math:`qR`-barrier constant (default: 0.197) can be
  added.


.. option:: massRange {min} {max} {particle1} {particle2} ... {particleN}

  Set mass interval for the invariant mass of :command:`particle1 particle2 ...
  particleN`.


.. option:: productionFormalism {string}

  Choose formalism for the production process: :code:`"Heli"`, :code:`"Cano"`
  or :code:`"Tensor"` are supported.


Antiproton-Proton-Annihilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: max {value}

  Maximal contributing orbital antiproton-proton momentum


.. option:: dropPbarpLForParticle {particle name} {value}

  Drop a pbarp L for certain particle


.. option:: pbarmom {momentum in GeV}

  Antiproton momentum in lab frame. Required for histogramming


.. option:: spinDensity {particle name}

  Calculate spin density matrix elements of this particle


Other Arguments
---------------

General
~~~~~~~

.. option:: datFile {path}

  Path to file containing data events


.. option:: errLogMode {string}

  Choose mode for error logger: "debug", "trace", "routine", "warning", "error"
  and "alert" are supported.


.. option:: mode {mode name}

  Set execution mode. See list of modes.


.. option:: mcFile {path}

  Path to file containing monte carlo events


.. option:: orderInFile {string1} {string2} {string3} {string4}

  Order of the 4-vector elements in input files: :code:`"Px Py Pz E"` and
  :code:`"E Px Py Pz"` are supported.


.. option:: serializationFile {path}

  Path to serialization file (for qa, sdm modes).


.. option:: name {string}

  Set a name that is attached to all output file names.


.. option:: noOfDataEvents {value}

  Set maximum number of data events for PWA- or QA-mode.


.. option:: paramFile {path}

  Path to file containing start parameters.


.. option:: pdgTableFile {path}

  Path of the PDG-table file, containing particle definitions and properties.


.. option:: preFactor {value}

  Set prefactor for an individual amplitude.


.. option:: ratioMcToData {value}

  Limit number of read monte carlo events by setting a ratio.


.. option:: unitInFile {string: GEV or MEV}

  Chosen unit in the ASCII input files.


.. option:: useEventWeight {bool}

  Read a (leading) weight for each set of four-vectors in the :option:`datFile`
  (**not** :option:`mcFile` currently).


Minimization
------------

.. option:: cacheAmps {bool}

  Enable/disable caching of amplitudes.


.. option:: mnParFix {parameter name}

  Fix parameter in Minuit.


.. option:: doScaling {bool}

  Enable/disable pwa fit just for scaling (scaling parameter free; remaining
  parameter fixed).


.. option:: contributionError {none}

  Calculate error of wave contribution in QA mode.


.. option:: evoIterations {int}

  Max number of iterations for evolutionary minimizer.


.. option:: evoPopulation {int}

  Population size for evolutionary minimizer.


.. option:: minimumTolerance {float}

  Sets the migrad tolerance argument. The tolerance that can optionally be
  passed to Migrad defines the convergence criterion for the minimization by
  comparing it to the current estimated distance to minimum (:code:`edm`). The
  minimization is stopped when :code:`edm < 0.0002*\mathrm{tolerance}` (see
  :file:`VariableMetricBuilder.cxx` of Minuit). The default value (in Migrad as
  well as in Pawian) is :code:`tolerance = 0.1`.


Histogramming
-------------

.. option:: histAngles {particle name 1} from {particle name 1} {particle name 2} ...

  Create histogram of angles of particle 1 in the cms frame of (particle
  1)+(particle 2)+....

.. option:: histMass {particle name 1} {particle name 2} ...

  Create histogram of the invariant mass of the specified particles.


Options for the Client/Server-Mode
----------------------------------

.. option:: serverAddress {value}

  IP or hostname of server.

.. option:: serverPort {value}

  Port for server. Default is 22222.

.. option:: noOfClients {value}

  Number of clients the server expects.


Options for Coupled Channel Mode
--------------------------------

.. option:: clientNumberWeights {weight channel 1} {weight channel 2} ...

  Weights to assign the number of clients to the respective channels. If the
  option is not given the clients are auto-assigned. The weights become
  normalized so :option:`clientNumberWeights = 0.2 0.8 <clientNumberWeights>`
  is the same as :option:`clientNumberWeights = 5 20 <clientNumberWeights>`.


Options for Event Generator
---------------------------

.. option:: generateWithModel {bool}

  Use model or phase space only.

.. option:: randomSeed {int}

  Start seed for random generator.

.. option:: noOfGenEvents {int}

  Number of events to be generated.


Options for :code:`qaMode`
--------------------------

.. option:: serializationFile {path}

  The serialization file contains the covariance matrix of the pwa result. This
  file is needed for later error calculations.

.. option:: contributionError {bool}

  Calculates the error of the total contribution and the individual
  contributions defined by :option:`calcContribution`. The --serializationFile
  argument needs to be set.

.. option:: calcContribution {name} withZeroAmp {string} {string} ...

  Calculates the number of events for a part of the full amplitude (i.e. only
  for a certain resonance). This is done by setting all parameters that include
  :command:`{string}` as a substring in their name to zero. The calculated
  contribution ist tagged with :command:`{name}`.

.. option:: saveContributionHistos {bool}

  Generates an additional root-File with histograms for each contribution
  defined with the option :option:`calcContribution` in the configuration file.
  The name of the ROOT file is composed by the general option :command:`{name}`
  and the :command:`{name}` of the contribution as set in the calcContribution
  option.
