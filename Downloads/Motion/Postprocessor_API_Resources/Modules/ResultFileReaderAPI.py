import clr, sys
from pathlib import Path

clr.AddReference("mscorlib")
clr.AddReference("System")
clr.AddReference("System.Collections")
clr.AddReference("VM.Enums.Post")
clr.AddReference("VM.Models")
clr.AddReference("VM.Models.OutputReader")
clr.AddReference("VM.Models.Post")
clr.AddReference("VM.Models.Post.ChartMathLib")
clr.AddReference("VM.Post.API.OutputReader")

from System import *
from System import Array, Action, Double, Int32
from System.IO import *
from System.Collections.Generic import IList, List

from VM import *
from VM.Enums.Post import *
from VM.Models import *
from VM.Models.OutputReader import *
from VM.Models.Post import *
from VM.Post.API.OutputReader import *
from VM.Models.Post.ChartMathLib import *
