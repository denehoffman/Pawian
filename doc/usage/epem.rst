Example: :math:`e^+e^-` data
============================

.. note::
  See original Wiki page `here
  <https://panda-wiki.gsi.de/foswiki/bin/view/PWA/Benchmarkepem>`__.

For testing and benchmarking purposes, an example PWA of Monte Carlo toy data
is located in the directory :file:`Examples/benchmark/epem`. It is an analysis
of the :math:`e^+e^-` annihilation to a :math:`J/\psi` charmonium with a
radiative decay to the final state particles :math:`\omega\omega\gamma` via two
intermediate resonances, the :math:`f_2(1910)` and the :math:`\eta(2225)`.

This is a step by step tutorial on running this example.

Relevant Files
--------------

The following files should be located in the :file:`benchmark` folder:

====================================  =========================================
File                                  Description
====================================  =========================================
:file:`JpsiGamOmegaOmega.cfg`         The configuration file. The fitted model
                                      is defined and all options set here
:file:`toyDataSignal.dat`             The "experimental" toy data as an ASCII
                                      file containing four-vectors
:file:`toyMC.dat`                     The phase space distributed toy Monte
                                      Carlo data as an ASCII file containing
                                      four-vectors
:file:`startparams.dat`               A list of default start parameters for
                                      the fit
:file:`benchmarkGam2Omega_runscript`  A shell script to run the PWA in
                                      client/server-mode
====================================  =========================================

Setting the Environment
-----------------------

If not yet done, the relevant environment variables have to be set. An example
for the Bochum site is given in the file :file:`SetEnv_rub` in the top
directory of the PAWIAN installation.

Starting the Server and the Clients
-----------------------------------

This example PWA is performed using the client/server mode with two clients,
each one using two threads. The relevant configuration file options are

.. code-block:: ini

  serverPort = 50003
  serverAddress = localhost
  noOfClients = 2
  noOfThreads = 2

To start the server run:

.. code-block:: bash

  $TOP_DIR/bin/singleChannelApp --epemFiles $TOP_DIR/Examples/benchmark/epem/JpsiGamOmegaOmega.cfg --mode server !>& server.log &

Wait some seconds and then start the two clients on the same machine:

.. code-block:: bash

  $TOP_DIR/bin/singleChannelApp --epemFiles $TOP_DIR/Examples/benchmark/epem/JpsiGamOmegaOmega.cfg --mode client !>& client1.log &
  $TOP_DIR/bin/singleChannelApp --epemFiles $TOP_DIR/Examples/benchmark/epem/JpsiGamOmegaOmega.cfg --mode client !>& client2.log &

In the server.log, lines similar to the following should appear:

.. code-block:: none

  INFO: ********* Server mode ****************
  INFO: Listening on port 50003
  INFO: Number of clients for channel 0 : 2
  INFO: Waiting for 2 clients ...
  INFO: Client pc14 logged in (ID 0).
  INFO: Client pc14 logged in (ID 1).
  INFO: All clients ready.
  INFO: start migrad

During the minimization, the current values of the likelihood are printed:

.. code-block:: none

  INFO: current LH = -3298.143159127945

while the current parameter values are frequently dumped into a file named
currentResult.dat. After the minimization has finished, the log file should
show:

.. code-block:: none

  **************** Minuit FunctionMinimum information ******************
  INFO: Function minimum is valid.
  INFO: Covariance matrix is valid.
  INFO: Final LH: -6184.421575
  INFO: # of function calls: 1329
  INFO: minimum edm: 1.288794568e-06

The file :file:`finalResult.dat` now holds the list of final fit parameters
with proper errors. Another file :file:`serializedOutput.dat` is created in
which (currently only) the covariance matrix is stored for later use.

Plotting the Result and Calculating the Contributions of the Intermediate Resonances
------------------------------------------------------------------------------------

For this purpose it is needed to start the QA mode by running:

.. code-block:: bash

  $TOP_DIR/bin/singleChannelApp -c $TOP_DIR/Examples/benchmark/epem/JpsiGamOmegaOmega.cfg --mode qaMode --paramFile finalResult_epem.dat --serializationFile serializedOutput_epem.dat !>& logQA.log &
