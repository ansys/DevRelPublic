@ECHO OFF
SETLOCAL

REM Set version-specific environment variables
set ANSYS_INSTALL_DIR=%AWP_ROOT242%
set FLUENT_MULTIPORT_VERSION=24.2.0

REM parse parallel execution command line arguments (if any)
set NPROCS=1
set MACHINES=
:ARGLOOP
  if [%1]==[] goto :ARGDONE
  if [%1]==[-np] (
    set NPROCS=%2
    shift
    shift
    goto :ARGLOOP
  )
  if [%1]==[-hostlist] (
    set MACHINES=%~2
    shift
    shift
    goto :ARGLOOP
  )
  set ARGS=%ARGS% %1
  shift
  goto :ARGLOOP
:ARGDONE

REM Set up Intel MPI environment.
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin\compiler;%PATH%
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin;%PATH%;
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\cnlauncher\fluent\fluent%FLUENT_MULTIPORT_VERSION%\multiport\mpi_wrapper\win64\intel;%PATH%;
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\cnlauncher\fluent\fluent%FLUENT_MULTIPORT_VERSION%\multiport\mpi\win64\intel\bin;%PATH%;
set I_MPI_ROOT=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\cnlauncher\fluent\fluent%FLUENT_MULTIPORT_VERSION%\multiport\mpi\win64
set INTEL_MPI_MACHINESFILE=intel_mpi_machines.txt
if "%MACHINES%" == "" (
  REM Local parallel
  set SYSC_MPI_CMD="%I_MPI_ROOT%\intel\bin\mpiexec.exe" -n %NPROCS% -localonly -noprompt -hosts localhost
) else (
  REM Distributed parallel
  if exist "%INTEL_MPI_MACHINESFILE%" (
    REM file exists
    echo "%INTEL_MPI_MACHINESFILE% exists - removing it"
    del %INTEL_MPI_MACHINESFILE%
  )
  :MACHINELOOP
  for /F "tokens=1* delims=," %%a in ("%MACHINES%") do (
    echo %%a >> %INTEL_MPI_MACHINESFILE%
    set MACHINES=%%b
  )
  if not "%MACHINES%" == "" goto :MACHINELOOP
  set SYSC_MPI_CMD="%I_MPI_ROOT%\intel\bin\mpiexec.exe" -machinefile %INTEL_MPI_MACHINESFILE% -noprompt
)

REM Run the application.
%SYSC_MPI_CMD% ChannelFlowMockSolverMPI.exe %ARGS%
exit /B %errorlevel%
