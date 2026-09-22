
 ****************************************************************
 *                                                              *
 *                     ANSYS MOTION 2025 R1                     * 
 *                                                              *
 ****************************************************************
 * RELEASE                  = 2025 R1
 * BUILD                    = 25.10.41029 (build.88dd8299)
 * VERSION                  = WINDOWS x64
 ================================================================
 < START OF ANSYS WINDOWS SOLVER ( ANSYS WINDOWS SOLVER VERSION = 25.10.41029 (build.88dd8299) ) >
 * START TIME = Tue Mar 18 11:12:57 2025 AM
-----------------------------------------------------------------
 * DFS FILE    = C:\Users\dyun\AppData\Local\Temp\DFS1D0F.dfs
 * SSC FILE    = C:\Users\dyun\AppData\Local\Temp\DFS1D0F.ssc
 * RESULT FILE = D:\ansysdev\2025R1_ToturalModel\Modal Flex\Modal Flex Result\Tutorial_ModalFlex
-----------------------------------------------------------------
 * OS  = Windows 10 Enterprise
 * CPU = 13th Gen Intel(R) Core(TM) i7-13850HX
 * MEM = 65228 MB
 * MACHINE NAME = AAPJOI4NEPEKZLP
-----------------------------------------------------------------

 < DF_SOLVER:: START READING A SSC FILE > 
 < DF_SOLVER:: END READING A SSC FILE > 


 < DF_SOLVER:: START READING A DFS FILE > 
========================================================
 * GRAVITY = 0.0000e+00 -9.8066e+03 0.0000e+00
========================================================
 < DF_SOLVER:: END READING A DFS FILE > 


 < DF_SOLVER:: START READING A GINF FILE > 
 < DF_SOLVER:: END READING A GINF FILE > 


 < DF_SOLVER:: START CHECKING LICENSE > 
 < DF_SOLVER:: END CHECKING LICENSE > 


 < DF_SOLVER:: START PRE-ADDRESSING > 
========================================================
 * NO. OF PROCESSORS = 1
 * MY PROCESSOR      = 0
 * NO. OF THREADS    = 8
========================================================

 < DF_SOLVER:: PRE-ADDRESSING - Moving Reference Frame > 

 < DF_SOLVER:: PRE-ADDRESSING - Frontal Solver > 

 < DF_SOLVER:: PRE-ADDRESSING - EasyFlex > 

 < DF_SOLVER:: PRE-ADDRESSING - Moving Reference Frame > 

 < DF_SOLVER:: PRE-ADDRESSING - Super Solver > 

 < DF_SOLVER:: PRE-ADDRESSING - Constraint > 

 < DF_SOLVER:: PRE-ADDRESSING - Frontal Solver > 

 < DF_SOLVER:: PRE-ADDRESSING - Thermal > 

 < DF_SOLVER:: PRE-ADDRESSING - Modal Data > 

 < DF_SOLVER:: PRE-ADDRESSING - Super Solver > 

 < DF_SOLVER:: PRE-ADDRESSING - Body > 

 < DF_SOLVER:: END PRE-ADDRESSING > 


 < DF_SOLVER:: START ADDRESSING > 

 < DF_SOLVER:: ADDRESSING - Simulation Scenario > 

 < DF_SOLVER:: ADDRESSING - EasyFlex > 

 < DF_SOLVER:: ADDRESSING - Threading > 

 < DF_SOLVER:: ADDRESSING - Function > 

 < DF_SOLVER:: ADDRESSING - Contact > 

 < DF_SOLVER:: ADDRESSING - Geometry > 

 < DF_SOLVER:: ADDRESSING - Finite Element > 

 < DF_SOLVER:: ADDRESSING - EasyFlex Element > 

 < DF_SOLVER:: ADDRESSING - Material > 

 < DF_SOLVER:: ADDRESSING - Constraint > 

 < DF_SOLVER:: ADDRESSING - Force > 

 < DF_SOLVER:: ADDRESSING - Equation > 

 < DF_SOLVER:: ADDRESSING - Sub-entity > 

 < DF_SOLVER:: ADDRESSING - Body > 

-----------------------------------------------------------------
 < MASTER NODE ID AND POSITION >
-----------------------------------------------------------------
 FBODY   MASTER
    ID   NODE ID          MASTER NODE POSITION         FBODY NAME
-----------------------------------------------------------------
     0         0  -1.5661e+01  1.0011e+02  1.2817e+01  /Tutorial_ModalFlex/FEBody_01
                   0.0000e+00 -0.0000e+00 -2.2105e-01
-----------------------------------------------------------------

 < DF_SOLVER:: ADDRESSING - Thermal > 

 < DF_SOLVER:: END ADDRESSING > 


 < DF_SOLVER:: START POST-ADDRESSING > 

 < DF_SOLVER:: POST-ADDRESSING - Super Solver > 

 < DF_SOLVER:: POST-ADDRESSING - Constraint > 

 < DF_SOLVER:: POST-ADDRESSING - EasyFlex > 

 < DF_SOLVER:: POST-ADDRESSING - Force > 

 < DF_SOLVER:: POST-ADDRESSING - Body > 
-----------------------------------------------------------------
 < DF_SOLVER:: MASS OF DUMMY BODIES >
-----------------------------------------------------------------
 TYPE      ID   COMP         VALUE   BODY NAME     
-----------------------------------------------------------------
 RBODY     3   Mass    3.8651e-07   /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02
 RBODY     3   JXX     1.0000e-04   /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02
 RBODY     3   JYY     1.0000e-04   /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02
 RBODY     3   JZZ     1.0000e-04   /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02
 RBODY     4   Mass    3.8651e-07   /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
 RBODY     4   JXX     1.0000e-04   /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
 RBODY     4   JYY     1.0000e-04   /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
 RBODY     4   JZZ     1.0000e-04   /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
------------------------------------------------------------------

 < DF_SOLVER:: END POST-ADDRESSING > 


-----------------------------------------------------------------
 * DF_SOLVER:: WORKING DIRECTORY = D:\ansysdev\2025R1_ToturalModel\Modal Flex\Modal Flex Result
-----------------------------------------------------------------

=================================================================
 < START ANALYSIS >
 * ANALYSIS START TIME = Tue Mar 18 11:12:58 2025 AM
-----------------------------------------------------------------

 ---------------------------------------------------------------
 DFSOLVER:: CO-SIMULATION MODE = 0
 ---------------------------------------------------------------

---------------------------------------------------------------

 SIMULATION PROGRESS = 0.00 % 

---------------------------------------------------------------

 < DF_SOLVER:: START POSITION ANALYSIS > 
 < THE NO. OF REDUNDANT CONSTRAINT = 3 > 
-----------------------------------------------------------------
 REVOLUTE JOINT ID = 1 | ROT.Y | /Tutorial_ModalFlex/RJ_02
 REVOLUTE JOINT ID = 2 | TRA.Z | /Tutorial_ModalFlex/RJ_03
 TRANSLATIONAL JOINT ID = 0 | ROT.X | /Tutorial_ModalFlex/TJ_01
-----------------------------------------------------------------

-----------------------------------------------------------------------------
 < POSITION ANALYSIS > GLOBAL ITERATION(FOR LARGE MOTION) = 1
-----------------------------------------------------------------------------

--------------------------------------------------
 < POSITION ANALYSIS > ITERATION = 0
 -------------------------------------------------- 
 * < POSITION ANALYSIS > MAXIMUM DELTA OF FBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = -2.1643e-08
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 1.2538e-10
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = EULER. Z
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < POSITION ANALYSIS > MAXIMUM DELTA OF RBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = -3.2191e-07
 * BODY ID    = 2
 * BODY NAME  = /Tutorial_ModalFlex/Crank_002
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 1.2538e-10
 * BODY ID    = 3
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02
 * COORDINATE = EULER. Z
 -------------------------------------------------- 

 -------------------------------------------------- 
 * < POSITION ANALYSIS > CONSTRAINT VIOLATION 
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 3.0224e-07
 * JOINT TYPE      = REVOLUTE JOINT
 * JOINT ID        = 2
 * JOINT NAME      = /Tutorial_ModalFlex/RJ_03
 * CONSTRAINT TYPE = TRA.Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = -9.3913e-11
 * JOINT TYPE      = MOTION
 * JOINT ID        = 0
 * JOINT NAME      = /Tutorial_ModalFlex/RJ_01/DriverFunction
 * CONSTRAINT TYPE = ROT.Z
 -------------------------------------------------- 

--------------------------------------------------
 < POSITION ANALYSIS > NORM = 3.2307e-06
--------------------------------------------------

--------------------------------------------------
 < POSITION ANALYSIS > ITERATION = 1
 -------------------------------------------------- 
 * < POSITION ANALYSIS > MAXIMUM DELTA OF FBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = -1.9740e-15
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 2.7071e-17
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = EULER. Z
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < POSITION ANALYSIS > MAXIMUM DELTA OF RBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 5.5295e-14
 * BODY ID    = 4
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = -1.4043e-16
 * BODY ID    = 1
 * BODY NAME  = /Tutorial_ModalFlex/Crank
 * COORDINATE = EULER. Z
 -------------------------------------------------- 

 -------------------------------------------------- 
 * < POSITION ANALYSIS > CONSTRAINT VIOLATION 
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 1.6983e-09
 * JOINT TYPE      = REVOLUTE JOINT
 * JOINT ID        = 2
 * JOINT NAME      = /Tutorial_ModalFlex/RJ_03
 * CONSTRAINT TYPE = TRA.Z
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = -1.3259e-14
 * JOINT TYPE      = TRANSLATIONAL JOINT
 * JOINT ID        = 0
 * JOINT NAME      = /Tutorial_ModalFlex/TJ_01
 * CONSTRAINT TYPE = ROT.X
 -------------------------------------------------- 

--------------------------------------------------
 < POSITION ANALYSIS > NORM = 3.8877e-13
--------------------------------------------------
--------------------------------------------------
 < DF_SOLVER:: POSITION ANALYSIS SUCCESSFULLY FINISHED  > 

-----------------------------------------------------------------
 < DF_SOLVER:: SIMULATION TIME = 0.0000e+00 >
 < DF_SOLVER:: SOLVED FBODY POSITION >
-----------------------------------------------------------------
 FBODY   MASTER
    ID   NODE ID          MASTER NODE POSITION         FBODY NAME
-----------------------------------------------------------------
     0         0  -1.5661e+01  1.0011e+02  1.2817e+01  /Tutorial_ModalFlex/FEBody_01
                  -1.9018e-19  8.4626e-19 -2.2105e-01
-----------------------------------------------------------------
 < DF_SOLVER:: END POSITION ANALYSIS > 

---------------------------------------------------------------

 SIMULATION PROGRESS = 20.00 % 

---------------------------------------------------------------

 < DF_SOLVER:: START VELOCITY ANALYSIS > 
 < THE NO. OF VELOCITY CONSTRAINTS = 33 > 
 < THE NO. OF USER DEFINED INITIAL VELOCITIES = 0 > 
 < THE NO. OF TOTAL REDUNDANT CONSTRAINT = 5 > 
 < THE NO. OF REDUNDANT CONSTRAINT = 5 > 
 < THE NO. OF REDUNDANT INITIAL VELOCITY = 0 > 
-----------------------------------------------------------------
 REVOLUTE JOINT ID = 1 | TRA.Y | /Tutorial_ModalFlex/RJ_02
 REVOLUTE JOINT ID = 1 | TRA.Z | /Tutorial_ModalFlex/RJ_02
 REVOLUTE JOINT ID = 2 | TRA.Z | /Tutorial_ModalFlex/RJ_03
 TRANSLATIONAL JOINT ID = 0 | ROT.X | /Tutorial_ModalFlex/TJ_01
 TRANSLATIONAL JOINT ID = 0 | ROT.Z | /Tutorial_ModalFlex/TJ_01
-----------------------------------------------------------------

 <WARNING> A constraint close to a redundant constraint was detected during velocity analysis.
 Check the constraints below carefully.
 -----------------------------------------------------------------
  REVOLUTE JOINT ID = 1 | TRA.X | /Tutorial_ModalFlex/RJ_02
  REVOLUTE JOINT ID = 2 | TRA.X | /Tutorial_ModalFlex/RJ_03
 -----------------------------------------------------------------

 < THE NO. OF INDEPENDENT COORDINATE = 2 > 
 --------------------------------------------------------- 
 ..... 1ST FLEXIBLE BODY TRANSLATION VELOCITY IN Y
 ..... 1ST FLEXIBLE BODY TRANSLATION VELOCITY IN Z
 --------------------------------------------------------- 

--------------------------------------------------
 < VEOLOCITY ANALYSIS > ITERATION = 0
 -------------------------------------------------- 
 * < VELOCITY ANALYSIS > MAXIMUM DELTA OF FBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = TRA. X
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = EULER. X
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < VELOCITY ANALYSIS > MAXIMUM DELTA OF RBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = TRA. X
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = EULER. X
 -------------------------------------------------- 

 -------------------------------------------------- 
 * < VELOCITY ANALYSIS > CONSTRAINT VIOLATION 
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * JOINT TYPE      = NULL
 * JOINT ID        = -1
 * JOINT NAME      = NULL
 * CONSTRAINT TYPE = NULL
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * JOINT TYPE      = NULL
 * JOINT ID        = -1
 * JOINT NAME      = NULL
 * CONSTRAINT TYPE = NULL
 -------------------------------------------------- 

--------------------------------------------------
 < VELOCITY ANALYSIS > NORM = 0.0000e+00
--------------------------------------------------

--------------------------------------------------
 < VEOLOCITY ANALYSIS > ITERATION = 1
 -------------------------------------------------- 
 * < VELOCITY ANALYSIS > MAXIMUM DELTA OF FBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = TRA. X
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = EULER. X
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < VELOCITY ANALYSIS > MAXIMUM DELTA OF RBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = TRA. X
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = EULER. X
 -------------------------------------------------- 

 -------------------------------------------------- 
 * < VELOCITY ANALYSIS > CONSTRAINT VIOLATION 
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * JOINT TYPE      = NULL
 * JOINT ID        = -1
 * JOINT NAME      = NULL
 * CONSTRAINT TYPE = NULL
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * JOINT TYPE      = NULL
 * JOINT ID        = -1
 * JOINT NAME      = NULL
 * CONSTRAINT TYPE = NULL
 -------------------------------------------------- 

--------------------------------------------------
 < VELOCITY ANALYSIS > NORM = 0.0000e+00
--------------------------------------------------
--------------------------------------------------
 < DF_SOLVER:: VELOCITY ANALYSIS SUCCESSFULLY FINISHED  > 

-----------------------------------------------------------------
 < DF_SOLVER:: SIMULATION TIME = 0.0000e+00 >
 < DF_SOLVER:: SOLVED FBODY VELOCITY >
-----------------------------------------------------------------
 FBODY   MASTER
    ID   NODE ID          MASTER NODE VELOCITY         FBODY NAME
-----------------------------------------------------------------
     0         0   0.0000e+00  0.0000e+00  0.0000e+00  /Tutorial_ModalFlex/FEBody_01
                   0.0000e+00  0.0000e+00  0.0000e+00
-----------------------------------------------------------------
 < DF_SOLVER:: END VELOCITY ANALYSIS > 

---------------------------------------------------------------

 SIMULATION PROGRESS = 40.00 % 

---------------------------------------------------------------

 < DF_SOLVER:: START ACCELERATION ANALYSIS > 
 < THE NO. OF REDUNDANT CONSTRAINT = 0 > 

 ========================================================= 
                    SYSTEM INFORMATION                     
 ========================================================= 
 < LINEAR SOLVER TYPE        = SPARSE SOLVER > 
 < JACOBIAN CALCULATE TYPE >
   - FINITE ELEMENT          = ANALYTIC
   - CONTACT (FLEX TO RIGID) = SIMPLIFY
   - CONTACT (FLEX TO FLEX)  = SIMPLIFY
 < STATIC ANALYSIS TYPE      = LINEAR, WITHOUT RIGID > 
 < HIGH-SPEED ROTATION       = PROGRAM CONTROLLED (DEPENDING ON SPEED) > 
 --------------------------------------------------------- 
 < SEARCHING CONTACT POINT = FREQUENTLY >
 < EVALUATION OF CONTACT ENTITY IN JACOBIAN = ALL >
 --------------------------------------------------------- 
 < THE DEGREE OF FREEDOM             = 23 > 
 < THE NO. OF VARIABLE EQUATION      = 0 > 
 < THE NO. OF 1ST DIFF. EQUATION     = 0 > 
   - IN TIRE MODULE                  = 0 > 
 < THE NO. OF 2ND DIFF. EQUATION     = 0 > 
 --------------------------------------------------------- 
 < THE NO. OF GENERALIZED COORDINATE = 62 > 
 < THE NO. OF CONSTRAINED EQUATION   = 39 > 
 --------------------------------------------------------- 
 < THE NO. OF RIGID BODY             = 5 > 
 < THE NO. OF FLEXIBLE BODY          = 1 > 
   - NODAL FE BODY                   = 0 > 
   - MODAL FE BODY                   = 1 > 
   - NODAL EasyFlex                  = 0 > 
   - MODAL EasyFlex                  = 0 > 
 < THE NO. OF NODE                   = 0 > 
   - FE NODE                         = 0 > 
   - OUTPUT NODE                     = 0 > 
 < THE NO. OF MODE                   = 26 > 
   - FE MODE                         = 26 > 
   - EasyFlex MODE                   = 0 > 
 < THE NO. OF EasyFlex NODE          = 0 > 
 --------------------------------------------------------- 
 < THE NO. OF FBEAM2                 = 0 > 
 < THE NO. OF FSHELL3                = 0 > 
   - FE FSHELL3                      = 0 > 
   - EasyFlex FSHELL3                = 0 > 
 < THE NO. OF FSHELL4                = 0 > 
 < THE NO. OF FSHELL6                = 0 > 
 < THE NO. OF FSHELL8                = 0 > 
 < THE NO. OF FSOLID4                = 0 > 
   - FE FSOLID4                      = 0 > 
   - EasyFlex FSOLID4                = 0 > 
 < THE NO. OF FSOLID5                = 0 > 
 < THE NO. OF FSOLID6                = 0 > 
 < THE NO. OF FSOLID8                = 0 > 
 < THE NO. OF FSOLID10               = 0 > 
 < THE NO. OF FSOLID13               = 0 > 
 < THE NO. OF FSOLID15               = 0 > 
 < THE NO. OF FSOLID20               = 0 > 
 --------------------------------------------------------- 
 < THE NO. OF MATERIAL               = 2 > 
 < THE NO. OF HYPER-ELASTIC MAT.     = 0 > 
 < THE NO. OF PLASTICITY MAT.        = 0 > 
 < THE NO. OF FOAM MAT.              = 0 > 
 --------------------------------------------------------- 
 < THE NO. OF R-T-R CONTACT          = 0 > 
 < THE NO. OF F-T-R CONTACT          = 0 > 
 < THE NO. OF F-T-F CONTACT          = 0 > 
 < THE NO. OF GENERAL CONTACT        = 0 > 
   - RIGID TO RIGID CASE             = 0 > 
   - FLEX TO RIGID CASE              = 0 > 
   - FLEX TO FLEX CASE               = 0 > 
 < THE NO. OF 2.5D CONTACT           = 0 > 
 < THE NO. OF TIE CONTACT            = 0 > 
 < THE NO. OF NTSF CONTACT           = 0 > 
 < THE NO. OF PTCV CONTACT           = 0 > 
 < THE NO. OF CVCV CONTACT           = 0 > 
 < THE NO. OF PTSF CONTACT           = 0 > 
 < THE NO. OF BC                     = 0 > 
 --------------------------------------------------------- 
 < THE NO. OF NURBS SURFACE          = 0 > 
 < THE NO. OF PIECEWISE SURFACE      = 0 > 
 < THE NO. OF PIECEWISE PATCH        = 0 > 
 < THE NO. OF PIECEWISE VERTEX       = 0 > 
 --------------------------------------------------------- 
 < THE NO. OF FE PLOAD               = 0 > 
 < THE NO. OF EasyFlex PLOAD         = 0 > 
 < THE NO. OF FE CLOAD               = 0 > 
 < THE NO. OF EasyFlex CLOAD         = 0 > 
 ========================================================= 

--------------------------------------------------
 < ACCELERATION ANALYSIS > ITERATION = 0
 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > MAXIMUM DELTA OF MBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = -1.9829e+03
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = MODE
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = EULER. X
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > MAXIMUM DELTA OF FBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 6.0936e+03
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 1.4374e+01
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = EULER. Z
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > MAXIMUM DELTA OF RBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 9.8223e+03
 * BODY ID    = 4
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 9.8820e+00
 * BODY ID    = 3
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02
 * COORDINATE = EULER. Z
 -------------------------------------------------- 

 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > CONSTRAINT VIOLATION 
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 0.0000e+00
 * JOINT TYPE      = NULL
 * JOINT ID        = -1
 * JOINT NAME      = NULL
 * CONSTRAINT TYPE = NULL
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = -3.7694e-01
 * JOINT TYPE      = MOTION
 * JOINT ID        = 0
 * JOINT NAME      = /Tutorial_ModalFlex/RJ_01/DriverFunction
 * CONSTRAINT TYPE = ROT.Z
 -------------------------------------------------- 

--------------------------------------------------
 < ACCELERATION ANALYSIS > NORM = 1.6923e+03
 < ACCELERATION ANALYSIS > CST NORM = 1.6234e+03
--------------------------------------------------

--------------------------------------------------
 < ACCELERATION ANALYSIS > ITERATION = 1
 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > MAXIMUM DELTA OF MBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 1.7924e-13
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = MODE
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 0.0000e+00
 * BODY ID    = -1
 * BODY NAME  = NULL
 * COORDINATE = EULER. X
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > MAXIMUM DELTA OF FBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 3.0115e-13
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = -2.4946e-15
 * BODY ID    = 0
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01
 * COORDINATE = EULER. Z
 -------------------------------------------------- 
 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > MAXIMUM DELTA OF RBODY
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = 9.4886e-13
 * BODY ID    = 2
 * BODY NAME  = /Tutorial_ModalFlex/Crank_002
 * COORDINATE = TRA. Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = -3.2385e-14
 * BODY ID    = 4
 * BODY NAME  = /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03
 * COORDINATE = EULER. Z
 -------------------------------------------------- 

 -------------------------------------------------- 
 * < ACCELERATION ANALYSIS > CONSTRAINT VIOLATION 
 -------------------------------------------------- 
 * MAX. VALUE IN TRANSLATIONAL TERM = -5.8283e+03
 * JOINT TYPE      = MTR JOINT
 * JOINT ID        = 0
 * JOINT NAME      = /Tutorial_ModalFlex/FEBody_01/FRBE_001^/Tutorial_ModalFlex/RJ_02/Tutorial_ModalFlex/FEBody_01
 * CONSTRAINT TYPE = TRA.Y
 -------------------------------------------------- 
 * MAX. VALUE IN ROTATIONAL TERM = 3.5527e-14
 * JOINT TYPE      = MTR JOINT
 * JOINT ID        = 1
 * JOINT NAME      = /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03/Tutorial_ModalFlex/FEBody_01
 * CONSTRAINT TYPE = ROT.Z
 -------------------------------------------------- 

--------------------------------------------------
 < ACCELERATION ANALYSIS > NORM = 3.6206e-13
 < ACCELERATION ANALYSIS > CST NORM = 1.6222e-11
--------------------------------------------------
----------------------------------------------------------------
 < DF_SOLVER:: ACCELERATION ANALYSIS SUCCESSFULLY FINISHED  > 

-----------------------------------------------------------------
 < DF_SOLVER:: SIMULATION TIME = 0.0000e+00 >
 < DF_SOLVER:: SOLVED FBODY ACCELERATION >
-----------------------------------------------------------------
 FBODY   MASTER
    ID   NODE ID          MASTER NODE ACCELERATION     FBODY NAME
-----------------------------------------------------------------
     0         0   2.4539e+02 -6.0936e+03 -1.9555e+00  /Tutorial_ModalFlex/FEBody_01
                   6.9041e-03  3.3693e-01 -1.4374e+01
-----------------------------------------------------------------
----------------------------------------------------------------
 < DF_SOLVER:: MAXIMUM RIGID BODY ACCELERATION > 
 < DF_SOLVER:: RIGID BODY ID   = 4 >
 < DF_SOLVER:: RIGID BODY NAME = /Tutorial_ModalFlex/FEBody_01/FRBE_002^/Tutorial_ModalFlex/RJ_03 >
 < DF_SOLVER:: ACCELERATION = 9.8223e+03 >
 < DF_SOLVER:: RIGID BODY POSITION  > 
  5.5704e-02  1.6903e+02  1.2500e+01
----------------------------------------------------------------
 < DF_SOLVER:: END ACCELERATION ANALYSIS > 

---------------------------------------------------------------

 SIMULATION PROGRESS = 60.00 % 

---------------------------------------------------------------

 < DF_SOLVER:: START EIGENVALUE ANALYSIS > 
 < THE NO. OF REDUNDANT CONSTRAINT = 1 > 
-----------------------------------------------------------------
 REVOLUTE JOINT ID = 0 | MOTION.ROT | /Tutorial_ModalFlex/RJ_01
-----------------------------------------------------------------
 < DF_SOLVER:: THE INITIAL NUMBER OF MODE REQUESTED IS 20.> 
 < DF_SOLVER:: THE NUMBER OF DEGREE OF FREEDOM IS 23.> 
 < DF_SOLVER:: THE NUMBER OF NORMAL MODES MUST BE LESS THAN OR EQUAL TO DEGREE OF FREEDOM.> 
-----------------------------------------------------------------
 < THE NUMBER OF NORMAL MODE REQUESTED IS CHANGED TO 23 > 
-----------------------------------------------------------------


 < DF_SOLVER:: START GET NORMAL MODE >

** LANCZOS :: Active Normal Mode = 23
** LANCZOS :: EIGENVALUE (NATURAL FREQUENCY) AFTER SORTING 0
-------------------------------------------------
       1    3.5527e-15 (9.4864e-09)
       2    1.2125e+08 (1.7525e+03)
       3    3.2765e+08 (2.8809e+03)
       4    1.2092e+09 (5.5343e+03)
       5    1.6134e+09 (6.3928e+03)
       6    6.7145e+09 (1.3042e+04)
       7    6.8865e+09 (1.3207e+04)
       8    8.1335e+09 (1.4354e+04)
       9    1.8363e+10 (2.1567e+04)
      10    1.9259e+10 (2.2087e+04)
      11    2.4865e+10 (2.5097e+04)
      12    3.2318e+10 (2.8612e+04)
      13    4.2794e+10 (3.2924e+04)
      14    5.5520e+10 (3.7501e+04)
      15    7.5279e+10 (4.3667e+04)
      16    7.9972e+10 (4.5008e+04)
      17    8.2883e+10 (4.5820e+04)
      18    9.6581e+10 (4.9461e+04)
      19    9.7419e+10 (4.9675e+04)
      20    1.0231e+11 (5.0907e+04)
      21    1.2050e+11 (5.5247e+04)
      22    1.2498e+11 (5.6266e+04)
      23    1.4676e+11 (6.0972e+04)
-------------------------------------------------
-----------------------------------------------------------------
< THE NUMBER OF NORMAL MODE REQUESTED IS 23 > 
< THE NUMBER OF CONVERGED NORMAL MODE IS 20 > 
-----------------------------------------------------------------

 < DF_SOLVER:: START RECOVER EIGENVALUE & EIGENVECTOR >

 < DF_SOLVER:: START MAKE DFEIG FILE >

----------------------------------------------------------
 < DF_SOLVER:: EIGENVALUE ANALYSIS SUCCESSFULLY FINISHED  > 

 < DF_SOLVER:: END EIGENVALUE ANALYSIS > 

---------------------------------------------------------------

 SIMULATION PROGRESS = 80.00 % 

---------------------------------------------------------------

 < DF_SOLVER:: START DYNAMICS ANALYSIS > 
 < NUMERICAL DAMPING = 1.0000e+00 > 
------------------------------------------------------------------------------------------
       TIME       STEPSIZE       ANORM     NRES   NJAC   STEP  NFAIL  IFAIL  RFAIL   IPCG
------------------------------------------------------------------------------------------
     0.0000e+00  1.0000e-04  0.0000e+00      11     11      0      0      0      0      0
     3.0000e-04  1.0000e-04  3.4928e-05      20     14      3      0      0      0      0
     * 3.4915e-05   * 1.2317e-08  
     5.0000e-04  2.0000e-04  3.9302e-05      23     15      4      0      0      0      0
     * 3.9253e-05   * 4.9254e-08  
     4.9000e-03  4.0000e-04  1.1584e-04      56     26     15      0      0      0      0
     * 1.1565e-04   * 1.8949e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 2.1000 % 
------------------------------------------------------------------------------------------
     1.0500e-02  8.0000e-04  9.0306e-04      84     33     22      0      0      0      0
     * 9.0236e-04   * 7.0688e-07   * 4.0983e-13  
     1.3700e-02  8.0000e-04  8.8928e-04     100     37     26      0      0      0      0
     * 8.8857e-04   * 7.0996e-07   * 5.0846e-13  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 4.1400 % 
------------------------------------------------------------------------------------------
     2.0700e-02  1.0000e-03  1.6733e-03     128     44     33      0      0      0      0
     * 1.6722e-03   * 1.1341e-06   * 1.4909e-12  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 6.1400 % 
------------------------------------------------------------------------------------------
     3.0700e-02  1.0000e-03  1.5824e-03     168     54     43      0      0      0      0
     * 1.5810e-03   * 1.3943e-06   * 2.6736e-12  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 8.1400 % 
------------------------------------------------------------------------------------------
     4.0700e-02  1.0000e-03  1.5885e-03     208     64     53      0      0      0      0
     * 1.5866e-03   * 1.9358e-06   * 4.5762e-12  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 10.1400 % 
------------------------------------------------------------------------------------------
     5.0700e-02  1.0000e-03  1.9413e-03     248     74     63      0      0      0      0
     * 1.9386e-03   * 2.7204e-06   * 6.4913e-12  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 12.1400 % 
------------------------------------------------------------------------------------------
     6.0700e-02  1.0000e-03  2.7218e-03     288     84     73      0      0      0      0
     * 2.7182e-03   * 3.5673e-06   * 5.6239e-12  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 14.1400 % 
------------------------------------------------------------------------------------------
     7.0700e-02  1.0000e-03  3.6995e-03     328     94     83      0      0      0      0
     * 3.6952e-03   * 4.2686e-06   * 4.4044e-12  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 16.1400 % 
------------------------------------------------------------------------------------------
     8.0700e-02  1.0000e-03  4.9326e-03     368    104     93      0      0      0      0
     * 4.9277e-03   * 4.9019e-06   * 2.5343e-11  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 18.1400 % 
------------------------------------------------------------------------------------------
     9.0700e-02  1.0000e-03  7.7060e-03     408    114    103      0      0      0      0
     * 7.7000e-03   * 5.9924e-06   * 3.8038e-11  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 20.1400 % 
------------------------------------------------------------------------------------------
     1.0070e-01  1.0000e-03  1.0095e-02     448    124    113      0      0      0      0
     * 1.0086e-02   * 8.6838e-06   * 3.4362e-11  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 22.1400 % 
------------------------------------------------------------------------------------------
     1.1070e-01  1.0000e-03  1.2206e-02     488    134    123      0      0      0      0
     * 1.2191e-02   * 1.4357e-05   * 2.9195e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 24.1400 % 
------------------------------------------------------------------------------------------
     1.2070e-01  1.0000e-03  2.4763e-02     528    144    133      0      0      0      0
     * 2.4734e-02   * 2.9335e-05   * 1.0073e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 26.1400 % 
------------------------------------------------------------------------------------------
     1.3070e-01  1.0000e-03  5.0738e-02     568    154    143      0      0      0      0
     * 5.0691e-02   * 4.7335e-05   * 1.7948e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 28.1400 % 
------------------------------------------------------------------------------------------
     1.4070e-01  1.0000e-03  4.3835e-02     608    164    153      0      0      0      0
     * 4.3782e-02   * 5.2654e-05   * 4.5930e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 30.1400 % 
------------------------------------------------------------------------------------------
     1.5070e-01  1.0000e-03  7.2838e-02     648    174    163      0      0      0      0
     * 7.2747e-02   * 9.1237e-05   * 2.8486e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 32.1400 % 
------------------------------------------------------------------------------------------
     1.6070e-01  1.0000e-03  6.6376e-02     688    184    173      0      0      0      0
     * 6.6249e-02   * 1.2729e-04   * 8.1121e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 34.1400 % 
------------------------------------------------------------------------------------------
     1.7070e-01  1.0000e-03  1.9092e-01     728    194    183      0      0      0      0
     * 1.9073e-01   * 1.8808e-04   * 2.1046e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 36.1400 % 
------------------------------------------------------------------------------------------
     1.8070e-01  1.0000e-03  1.5008e-01     768    204    193      0      0      0      0
     * 1.4993e-01   * 1.4174e-04   * 3.9747e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 38.1400 % 
------------------------------------------------------------------------------------------
     1.9070e-01  1.0000e-03  2.4755e-01     808    214    203      0      0      0      0
     * 2.4732e-01   * 2.2954e-04   * 3.2527e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 40.1400 % 
------------------------------------------------------------------------------------------
     2.0070e-01  1.0000e-03  2.9184e-01     848    224    213      0      0      0      0
     * 2.9150e-01   * 3.4352e-04   * 3.7490e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 42.1400 % 
------------------------------------------------------------------------------------------
     2.1070e-01  1.0000e-03  2.9129e-01     888    234    223      0      0      0      0
     * 2.9092e-01   * 3.7737e-04   * 3.9688e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 44.1400 % 
------------------------------------------------------------------------------------------
     2.2070e-01  1.0000e-03  3.4797e-01     928    244    233      0      0      0      0
     * 3.4741e-01   * 5.6323e-04   * 1.3779e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 46.1400 % 
------------------------------------------------------------------------------------------
     2.3070e-01  1.0000e-03  8.8405e-01     968    254    243      0      0      0      0
     * 8.8327e-01   * 7.7641e-04   * 1.7767e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 48.1400 % 
------------------------------------------------------------------------------------------
     2.4070e-01  1.0000e-03  8.8512e-01    1008    264    253      0      0      0      0
     * 8.8439e-01   * 7.3085e-04   * 1.4932e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 50.1400 % 
------------------------------------------------------------------------------------------
     2.5070e-01  1.0000e-03  1.7993e+00    1048    274    263      0      0      0      0
     * 1.7979e+00   * 1.3859e-03   * 1.2176e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 52.1400 % 
------------------------------------------------------------------------------------------
     2.6070e-01  1.0000e-03  1.2934e+00    1088    284    273      0      0      0      0
     * 1.2922e+00   * 1.2563e-03   * 4.9991e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 54.1400 % 
------------------------------------------------------------------------------------------
     2.7070e-01  1.0000e-03  1.8367e+00    1135    294    283      0      0      0      0
     * 1.8352e+00   * 1.5595e-03   * 9.3538e-07   * 3.3874e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 56.1400 % 
------------------------------------------------------------------------------------------
     2.8070e-01  1.0000e-03  9.8481e-01    1181    304    293      0      0      0      0
     * 9.8358e-01   * 1.2327e-03   * 6.3951e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 58.1400 % 
------------------------------------------------------------------------------------------
     2.9070e-01  1.0000e-03  1.1985e+00    1221    314    303      0      0      0      0
     * 1.1969e+00   * 1.5903e-03   * 9.6560e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 60.1400 % 
------------------------------------------------------------------------------------------
     3.0070e-01  1.0000e-03  1.8760e+00    1270    324    313      0      0      0      0
     * 1.8739e+00   * 2.0871e-03   * 2.6739e-06   * 4.1784e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 62.1400 % 
------------------------------------------------------------------------------------------
     3.1070e-01  1.0000e-03  2.3445e+00    1318    334    323      0      0      0      0
     * 2.3424e+00   * 2.1038e-03   * 4.8288e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 64.1400 % 
------------------------------------------------------------------------------------------
     3.2070e-01  1.0000e-03  1.8253e+00    1363    344    333      0      0      0      0
     * 1.8230e+00   * 2.3678e-03   * 2.0837e-06   * 2.8309e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 66.1400 % 
------------------------------------------------------------------------------------------
     3.3070e-01  1.0000e-03  3.7656e+00    1412    354    343      0      0      0      0
     * 3.7619e+00   * 3.7231e-03   * 1.9559e-06   * 1.0122e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 68.1400 % 
------------------------------------------------------------------------------------------
     3.4070e-01  1.0000e-03  1.9578e+00    1462    364    353      0      0      0      0
     * 1.9552e+00   * 2.5913e-03   * 2.9852e-06   * 5.5029e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 70.1400 % 
------------------------------------------------------------------------------------------
     3.5070e-01  1.0000e-03  3.0680e+00    1509    374    363      0      0      0      0
     * 3.0641e+00   * 3.8523e-03   * 2.0554e-06   * 6.1097e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 72.1400 % 
------------------------------------------------------------------------------------------
     3.6070e-01  1.0000e-03  2.6336e+00    1559    384    373      0      0      0      0
     * 2.6300e+00   * 3.5598e-03   * 1.4702e-05   * 2.7115e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 74.1400 % 
------------------------------------------------------------------------------------------
     3.7070e-01  1.0000e-03  5.7263e+00    1606    394    383      0      0      0      0
     * 5.7210e+00   * 5.3251e-03   * 2.4379e-06  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 76.1400 % 
------------------------------------------------------------------------------------------
     3.8070e-01  1.0000e-03  5.2613e+00    1655    404    393      0      0      0      0
     * 5.2557e+00   * 5.5685e-03   * 6.6697e-06   * 3.6134e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 78.1400 % 
------------------------------------------------------------------------------------------
     3.9070e-01  1.0000e-03  3.3728e+00    1705    414    403      0      0      0      0
     * 3.3684e+00   * 4.4723e-03   * 5.4590e-06   * 2.9035e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 80.1400 % 
------------------------------------------------------------------------------------------
     4.0070e-01  1.0000e-03  3.5274e+00    1753    424    413      0      0      0      0
     * 3.5231e+00   * 4.2657e-03   * 1.9942e-05   * 4.1268e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 82.1400 % 
------------------------------------------------------------------------------------------
     4.1070e-01  1.0000e-03  4.6131e+00    1803    434    423      0      0      0      0
     * 4.6083e+00   * 4.7779e-03   * 6.5231e-06   * 8.4745e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 84.1400 % 
------------------------------------------------------------------------------------------
     4.2070e-01  1.0000e-03  6.5085e+00    1851    444    433      0      0      0      0
     * 6.5020e+00   * 6.5516e-03   * 3.8095e-07  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 86.1400 % 
------------------------------------------------------------------------------------------
     4.3070e-01  1.0000e-03  4.8408e+00    1901    454    443      0      0      0      0
     * 4.8342e+00   * 6.5689e-03   * 8.9499e-06   * 8.7697e-10  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 88.1400 % 
------------------------------------------------------------------------------------------
     4.4070e-01  1.0000e-03  4.3226e+00    1950    464    453      0      0      0      0
     * 4.3160e+00   * 6.5498e-03   * 1.3791e-05   * 2.7480e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 90.1400 % 
------------------------------------------------------------------------------------------
     4.5070e-01  1.0000e-03  4.3518e+00    2000    474    463      0      0      0      0
     * 4.3443e+00   * 7.5017e-03   * 2.6946e-05   * 6.1022e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 92.1400 % 
------------------------------------------------------------------------------------------
     4.6070e-01  1.0000e-03  4.2414e+00    2050    484    473      0      0      0      0
     * 4.2348e+00   * 6.4963e-03   * 1.9041e-05   * 5.6352e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 94.1400 % 
------------------------------------------------------------------------------------------
     4.7070e-01  1.0000e-03  4.5029e+00    2100    494    483      0      0      0      0
     * 4.4951e+00   * 7.8186e-03   * 3.7777e-05   * 1.0546e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 96.1400 % 
------------------------------------------------------------------------------------------
     4.8070e-01  1.0000e-03  4.4017e+00    2150    504    493      0      0      0      0
     * 4.3950e+00   * 6.6204e-03   * 2.1078e-05   * 7.1127e-09  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 98.1400 % 
------------------------------------------------------------------------------------------
     4.9070e-01  1.0000e-03  4.6318e+00    2200    514    503      0      0      0      0
     * 4.6238e+00   * 7.9892e-03   * 4.1816e-05   * 1.2325e-08  
     4.9870e-01  1.0000e-03  8.2952e+00    2240    522    511      0      0      0      0
     * 8.2867e+00   * 8.4015e-03   * 1.2258e-05   * 1.2552e-08  
------------------------------------------------------------------------------------------
 ANALYSIS PROGRESS = 100.0000 % 
------------------------------------------------------------------------------------------
     5.0000e-01  1.3000e-03  1.1739e+01    2245    523    512      0      0      0      0
     * 1.1727e+01   * 1.1967e-02   * 7.3022e-05   * 8.1272e-08  

 < DF_SOLVER:: END DYNAMICS ANALYSIS > 

-----------------------------------------------------------------

 SIMULATION PROGRESS = 100.00 % 

-----------------------------------------------------------------

-----------------------------------------------------------------
 < DF_SOLVER:: FUNCTION EVALUATION RESULT >
-----------------------------------------------------------------
                             MEMORY(MB)      TIME(SEC)      COUNT
 TOTAL                           129.32           3.07          1

         * PRE PROCESSING *
 READ_DFSFILE                     11.80           0.14          1
 READ_GINFFILE                     0.00           0.00          1
 PRE_ADDRESSING                    6.98           0.01          1
 ADDRESSING                       62.27           0.00          1
 POST_ADDRESSING                   0.00           0.00          1

         * ANALYSIS *
 POSITION_ANALYSIS                 5.61           0.13          1
 VELOCITY_ANALYSIS                 0.00           0.00          1
 ACCELERATION_ANALYSIS             0.00           0.01          1
 STATIC_ANALYSIS                   0.00           0.00          0
 DYNAMICS_ANALYSIS                 0.45           2.20          1
 EIGENVALUE_ANALYSIS               5.08           0.01          1
 THERMAL_ANALYSIS                  0.00           0.00          0

         * SPARSE LINEAR SOLVER *
 INITIALIZATION                    0.08           0.00          6
 SYMBOLIC_FACTORIZATION            0.67           0.11          5
 NUMERICAL_FACTORIZATION           4.91           0.09        519
 SOLVE                             0.06           0.17       1728
 CALC_RES_ONLY                     0.00           0.05       1210
 CALC_RES_JAC                      0.06           0.06        519

         * ITERATION LINEAR SOLVER *
 INITIALIZATION                    0.00           0.00          0
 CALC_PRECONDITIONER               0.00           0.00          0
 SYMBOLIC_FACTORIZATION            0.00           0.00          0
 NUMERICAL_FACTORIZATION           0.00           0.00          0
 SOLVE                             0.00           0.00          0

         * FRONTAL & SUPER SOLVER *
 SYMBOLIC_FACTORIZATION            0.00           0.00          0
 - ALL                             0.00           0.00          0
 - REDUCED                         0.00           0.00          0
 NUMERICAL_FACTORIZATION           0.00           0.00          0
 - STEP_1                          0.00           0.00          0
 - STEP_2                          0.00           0.00          0
 CALC_RES_ONLY                     0.00           0.00          0
 SOLVE                             0.00           0.00          0

         * GENERAL PROCESS *
 PRE-PROCESSING(PREGET)            0.00           0.13       1724
 - RECOVER                         0.00           0.12       1723
 - PRE_CALC_RES_JAC                0.00           0.00       1723
 POST-PROCESSING(RECOVER_RES)      0.19           0.19        513
 - RECOVER                         0.00           0.04        514
 - PRE_CALC_RESIDUAL               0.00           0.00        514
 - CALC_RESIDUAL                   0.06           0.02        513
 - POST_RECOVER                    0.12           0.13        513
 CALC_VEL_ACC_CONST(CSTCOR)        1.06           1.15       1722
 - RECOVER                         0.12           0.40       5677
 - PRE-PROCESSING                  0.00           0.00          0
 - CALC_RESIDUAL                   0.06           0.22       7153
 - SOLVE                           0.19           0.36       3955
 CHECK_REDUNDANT                   0.01           0.00          4

         * EIGEN PROCESS *
 INITIALIZATION                    0.14           0.00          1
 CALC_NORMAL_MODES                 4.10           0.01          1
 CALC_CORRECTION_MODES             0.00           0.00          0
 CALC_FINAL_MODES                  0.00           0.00          1
 - STEP_1(MODGRAM)                 0.00           0.00          1
 - STEP_2(REDUCEDK)                0.00           0.00          0
 - STEP_3(SOLVE)                   0.00           0.00          0
 - STEP_4(RECOVER)                 0.00           0.00          0
 RECOVER                           0.00           0.00          1
 MAKE_DFMF_FILE                    0.00           0.00          0
 - CALC_INVARIANT                  0.00           0.00          0
 - WRITE                           0.00           0.00          0
 - REPORT                          0.00           0.00          0
 MAKE_DFEIG_FILE                   0.00           0.00          1

         * CONTACT *
 CONTACT                           0.00           0.00          1
 - UPDATE_INFORMATION              0.00           0.00          0
 - FIND_POINT_CONTACT              0.00           0.00          0
   - GENERAL_CONTACT               0.00           0.00          0
   - 3D_CONTACT                    0.00           0.00          0
   - 25D_CONTACT                   0.00           0.00          0
 - FIND_POINT_GAP                  0.00           0.00          0
 - PREGET_RESJAC                   0.00           0.00          1

         * REPORT *
 REPORT                            0.04           0.07        513
-----------------------------------------------------------------
 < USED TIME = 0 H/ 0 M/ 2 S ( 2.35 SECONDS ) >
-----------------------------------------------------------------

--------------------------------------------------------------------------------
 < DF_SOLVER:: RECOMMENDED SETTINGS FOR ANALYSIS >
--------------------------------------------------------------------------------
 ===============================================================================
 =                                                                             =
 = 1. Improving Analysis Performance by Managing Velocity and Acceleration     =
 =    Constraints                                                              =
 =                                                                             =
 ===============================================================================
  Calculations to satisfy velocity and acceleration constraints have taken more 
  than a quarter of the dynamic analysis time. When the number of constraint 
  equations is excessive relative to the system's degrees of freedom, the 
  computation time related to constraints can be lengthy. To improve analysis 
  performance, consider the following:
  1. Model with Force Entities: Use force entities like bushings instead of 
     constraints like joints. This approach allows for some minor errors while 
     reducing the size of the system matrix and decreasing the computation time 
     related to constraints. However, be aware that penalties may occur 
     depending on stiffness, or if the stiffness is too high relative to the 
     system's rigidity, it can cause high-frequency noise or even NR failure.
  2. Disable Velocity and Acceleration Constraints: Enter "cstcortp, 0" in 
     Solver Commands(Beta) to bypass the calculation of velocity and 
     acceleration constraints. For more detailed instructions on applying Solver
     Commands, refer to the Help documentation. This will enhance analysis 
     performance but may introduce errors in the velocity and acceleration 
     solutions. Static problems may be less affected, but high-speed rotational 
     problems might experience increased errors, so review the results 
     carefully.
 ===============================================================================
--------------------------------------------------------------------------------

-----------------------------------------------------------------
 * END TIME = Tue Mar 18 11:13:00 2025 AM
 < END OF ANSYS WINDOWS SOLVER >
=================================================================
