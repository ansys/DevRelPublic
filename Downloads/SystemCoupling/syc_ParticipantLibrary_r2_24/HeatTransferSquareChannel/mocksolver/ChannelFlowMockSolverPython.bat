@ECHO OFF
SETLOCAL

REM Set version-specific environment variables
set ANSYS_INSTALL_DIR=%AWP_ROOT242%
set FLUENT_MULTIPORT_VERSION=24.2.0

REM Setup run-time environment
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin\compiler;%PATH%
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin;%PATH%
set PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\cnlauncher\fluent\fluent%FLUENT_MULTIPORT_VERSION%\multiport\mpi_wrapper\win64\stub;%PATH%
set PYTHONPATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin;%PYTHONPATH%;
set PYTHON_DLL_PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin;%PYTHON_DLL_PATH%
set PYTHON_DLL_PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\cnlauncher\fluent\fluent%FLUENT_MULTIPORT_VERSION%\multiport\mpi_wrapper\win64\stub;%PYTHON_DLL_PATH%
set PYTHON_DLL_PATH=%ANSYS_INSTALL_DIR%\SystemCoupling\runTime\winx64\bin\compiler;%PYTHON_DLL_PATH%

REM Run the application
"%ANSYS_INSTALL_DIR%\commonfiles\CPython\3_10\winx64\Release\python\python.exe" ChannelFlowMockSolver.py %*
exit /B %errorlevel%
