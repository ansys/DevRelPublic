#!/bin/sh

# Set version-specific environment variables
export ANSYS_INSTALL_DIR="$AWP_ROOT242"
export FLUENT_MULTIPORT_VERSION=24.2.0

NPROCS=1
while [ $# -gt 0 ]; do
    case "$1" in
    -np)
        NPROCS=$2
        shift
        shift
        ;;
    -hostlist)
        MACHINES=$2
        shift
        shift
        ;;
    *)
        if [ -z "${newargs}" ]; then
          newargs=$1
        else
          newargs="${newargs} ${1}"
        fi
        shift
        ;;
    esac
done
set -- ${newargs}
unset newargs

# Set up Intel MPI environment
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/bin/compiler:"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/bin:"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/cnlauncher/fluent/fluent${FLUENT_MULTIPORT_VERSION}/multiport/mpi_wrapper/lnamd64/intel:"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/cnlauncher/fluent/fluent${FLUENT_MULTIPORT_VERSION}/multiport/mpi/lnamd64/intel/lib:"$LD_LIBRARY_PATH"
export I_MPI_ROOT="${ANSYS_INSTALL_DIR}"/SystemCoupling/runTime/linx64/cnlauncher/fluent/fluent${FLUENT_MULTIPORT_VERSION}/multiport/mpi/lnamd64
export TMI_CONFIG="${I_MPI_ROOT}"/intel/etc/tmi.conf
if [ -z "$MACHINES" ]; then
    # Local parallel
    "$I_MPI_ROOT"/intel/bin/mpiexec.hydra -n "$NPROCS" ./ChannelFlowMockSolverMPI "$@"
else
    # Write hosts to Intel MPI-compatible machines file
    INTEL_MPI_MACHINESFILE="intel_mpi_machines.txt"
    if [ -f $INTEL_MPI_MACHINESFILE ]; then
        echo "${INTEL_MPI_MACHINESFILE} exists - removing it"
        rm $INTEL_MPI_MACHINESFILE
    fi
   (IFS=','; for HOST in $MACHINES; do echo "$HOST" >> $INTEL_MPI_MACHINESFILE; done)
    # Distributed parallel
    "$I_MPI_ROOT"/intel/bin/mpiexec.hydra -machinefile $INTEL_MPI_MACHINESFILE ./ChannelFlowMockSolverMPI "$@"
fi
exit $?
