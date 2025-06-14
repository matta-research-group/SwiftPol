Getting Started with Swiftpol
=============================

Welcome to the Swiftpol documentation! This guide will help you get started with using Swiftpol for polymer molecular dynamics simulations.





Installation
------------

To install Swiftpol dependencies, you can use `conda`. Run the following commands in your terminal:

.. code-block:: bash

    git clone https://github.com/matta-research-group/SwiftPol
    cd swiftpol
    conda env create -f Dev_tools/swiftpol.yml -n swiftpol -y 
    conda activate swiftpol

Swiftpol is in development. To install an editable version, clone the repository and install it locally:

.. code-block:: bash

    git clone https://github.com/matta-research-group/SwiftPol
    cd swiftpol
    pip install -e .

Installing With a Docker Container
----------------------------------

To install SwiftPol as a docker container, install `docker desktop <https://docs.docker.com/desktop/>`_ and run the following command in the terminal

.. code-block:: bash
    docker build --no-cache --platform=linux/amd64 -f Dockerfile -t my-image -t setup/swiftpol:latest .


to open a jupyter notebook run

.. code-block:: bash
    docker run --rm -it -p 8888:8888 setup/swiftpol:latest


Basic Usage
-----------

Example usage can be found in our `example notebooks <https://github.com/matta-research-group/SwiftPol/tree/main/Example_Notebooks>`_.

Contributing
------------

We welcome contributions from the community. All contributors to SwiftPol should follow our `Code of Conduct <https://github.com/matta-research-group/SwiftPol/tree/main/CODE_OF_CONDUCT.md>`_.
To submit your feature to be incorporated to the main branch, you should submit a Pull Request. The repository maintainers will review your pull request before accepting your changes.

Support
-------

If you encounter any issues or have any questions regarding the use of SwiftPol, please open an issue on our `GitHub repository <https://github.com/matta-research-group/SwiftPol/issues>`_.

License
-------

Swiftpol is licensed under the `BSD-3-Clause License <https://github.com/matta-research-group/SwiftPol/tree/main/LICENSE>`_.