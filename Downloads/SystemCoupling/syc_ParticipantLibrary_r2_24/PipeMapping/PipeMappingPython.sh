#!/bin/sh

# Set version-specific environment variables
export ANSYS_INSTALL_DIR="$AWP_ROOT242"
export FLUENT_MULTIPORT_VERSION=24.2.0

# Setup run-time environment
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/bin/compiler:"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/bin:"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/cnlauncher/fluent/fluent${FLUENT_MULTIPORT_VERSION}/multiport/mpi_wrapper/lnamd64/stub:"$LD_LIBRARY_PATH"
export PYTHONPATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/bin:"$PYTHONPATH"
export SYSC_PYTHON="${ANSYS_INSTALL_DIR}"/commonfiles/CPython/3_10/linx64/Release/python/bin

# Run the application
"${SYSC_PYTHON}"/python PipeMapping.py "$@"

# Exit.
exit $?

