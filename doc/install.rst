Installation
============

.. note::

  These instructions were extracted from the `original Wiki page
  <https://panda-wiki.gsi.de/foswiki/bin/view/PWA/FirstStepInstructions2>`__.

This documentation briefly summarizes the first steps on how to install and
compile the PAWIAN PWA software and how to run a test application.

Software Repository and Access to the Source Code
-------------------------------------------------

Navigate to `gitlab.ep1.rub.de <http://gitlab.ep1.rub.de>`_. You should be
greeted by the “Sign In” page of our GitLab server. As an external user, please
select the “Register” tab right above the User/Password fields. By supplying
your name, user name and an e-mail address you can now create a user account
for this GitLab. After having done so, you should be able to browse the PAWIAN
sources by navigating to :menuselection:`Explore public projects --> All -->
EP1/Pawian` ― or, even easier, by clicking `this link
<http://gitlab.ep1.rub.de/EP1/Pawian>`__.

You can view the code in the browser, or use the link at the top of the PAWIAN
project page to clone the repository to your local machine. Try a

.. code-block:: bash

  git clone http://jollyj.ep1.rub.de/EP1/Pawian.git

on your machine; it should ask you for your credentials. Supply the
username/password for the newly created GitLab account and it should start
cloning.

After having done this, make sure you have all the external packages installed
(Minuit2, boost-build, root). You may want to have a look at one of the
:file:`SetEnv_*` files, modify one of them to fit your needs and source it.

Some useful information and instructions about Git can be found `here
<https://panda-wiki.gsi.de/foswiki/bin/view/PWA/Git>`__.

External Packages
-----------------

Before you can compile and link the code the following external packages must
be installed:

* `Boost.Build <https://boostorg.github.io/build/>`_ (:code:`bjam`/:code:`b2`)

* `Boost <https://www.boost.org/users/download>`_, at least version 1.48.0

* `Minuit2 <https://www.cern.ch/minuit>`_, at least
  `v5.24.00
  <http://www.cern.ch/mathlibs/sw/5_24_00/Minuit2/Minuit2-5.24.00.tar.gz>`_.
  This has to be a standalone version with parallelization flags disabled
  (:command:`./configure --disable-openmp`).

* `ROOT <https://root.cern.ch>`_ *with Minuit2 disabled*


Compilation and Linking
-----------------------

To compile and run PAWIAN some environment variables have to be set. An example
is given in the file :code:`SetEnv_rub` in the top directory. The
:file:`Jamroot` file most likely has to be modified to let the paths
(libraries, includes) match your local environment. Used path constants like
:code:`${ROOTSYS}` or :code:`${extern}` are best defined in the
:file:`site-config.jam` configuration file in the boost-build directory, e.g
:file:`/usr/share/boost-build/site-config.jam` (example directory; create file
if not existing).

.. code-block:: tcsh

  path-constant ROOTSYS : <path to rootsys>
  path-constant extern : <path to folder with external libraries>
  ...

Then you are ready for compiling and linking the code. By typing:

* :code:`bjam` executables with debugging options will be created. These
  executables and the relevant libraries are then located at sub-directory
  :file:`/bin/<compile-version>/debug/link-static/`.

* :code:`bjam release` executables with optimized performance will be created
  (recommended for the actual PWA). Here you can find the relevant executables
  in sub-directory :file:`/bin/<compile-version>/release/link-static/`.

By typing :code:`bjam clean` or :code:`bjam release clean` you clean up your
created executables and libraries and you can start to compile the code from
scratch.

Starting an Example PWA
-----------------------

To test your PAWIAN installation, please have a look at the :doc:`example PWA
tutorials </usage>`.
