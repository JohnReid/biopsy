..  This is a reStructuredText document
    (http://docutils.sourceforge.net/rst.html) Although you can read it as text, it
    will be easier to read when processed with a tool like rst2html.


Build instructions
==================

Requirements
------------

To build the biopsy C++ code you will need to install:

- the Boost_ C++ source libraries as described in `their documentation
  <http://www.boost.org/doc/libs/1_54_0/more/getting_started/index.html>`_.
- version 2.7.7 (or similar) of the ANTLR_ parser generator and make
  sure that the ``antlr`` tool is in your PATH
- version 2 of the `Xerces C++ XML parser`_ headers and library files
- the `GNU Scientific Library`_ headers and library files
- a working `SSL library`_ (will be installed by default on most major OSes)

Note that easiest way to install these packages on most Linux OSes will be through their
package manager. For example, on Ubuntu/Debian you can install all of them with one command::

    sudo apt-get install libboost1.53-all-dev antlr libxerces-c2-dev libgsl0-dev libssl-dev

.. _Boost: http://www.boost.org/
.. _ANTLR: http://www.antlr.org/
.. _`Xerces C++ XML parser`: http://xerces.apache.org/xerces-c/
.. _`GNU Scientific Library`: http://www.gnu.org/software/gsl/
.. _`SSL library`: http://www.openssl.org/



Configuration
-------------

Create entries in your ``site-config.jam`` named ``antlr``, ``gsl``, ``ssl``
and ``xerces-c`` that describe where to find the header files and library files
of the requirements you have installed. Instructions on where to put your
``site-config.jam`` can be found `here
<http://www.boost.org/boost-build2/doc/html/bbv2/overview/configuration.html>`_
and instructions on the format of the `lib` rule you will need to use are `here
<http://www.boost.org/boost-build2/doc/html/bbv2/overview/configuration.html>`_.
For reference, the entries in my ``site-config.jam`` look like:

::

    # ANTLR
    lib antlr : : <name>antlr-pic ;

    # Xerces-C++
    lib xerces-c : : <name>xerces-c <link>shared ;

    # GSL
    lib gslcblas : : <file>/usr/lib/libgslcblas.a ;
    lib gsl : : <file>/usr/lib/libgsl.a : : <library>gslcblas ;

    # SSL
    lib ssl : : <name>crypto ;



Build
-----

To build the biopsy C++ code, first run the ``bash`` script

    ./grammars-run-antlr

to parse the ANTLR_ grammar files and create C++ lexers and parsers.
Now you are ready to build the biopsy C++ libraries

    b2 -q install

This command will install the libraries into the ``lib`` and ``../Python``
directories.

