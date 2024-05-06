@ECHO OFF
SETLOCAL

REM Set version-specific environment variables
set ANSYS_INSTALL_DIR=%AWP_ROOT242%
set FLUENT_MULTIPORT_VERSION=24.2.0

REM Setup run-time environment
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin\compiler;%PATH%
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin;%PATH%
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\cnlauncher\fluent\fluent%FLUENT_MULTIPORT_VERSION%\multiport\mpi_wrapper\win64\stub;%PATH%

REM Run the application
OscillatingPlateDamping.exe %*
exit /B %errorlevel%
