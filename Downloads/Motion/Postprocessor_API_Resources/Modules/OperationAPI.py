import clr

clr.AddReference("PresentationCore")
clr.AddReference("System")
clr.AddReference("System.Collections")
clr.AddReference("System.Runtime")
clr.AddReference("SciChart.Core")

from System import *
from System.Collections.Generic import List
from System.Windows.Media import *
from SciChart.Core import *
from System.IO import Path
from System.IO import DirectoryInfo
from VM import *
from VM.API.Post.Operations import *
from VM.Models import *
from VM.Models.OutputReader import *
from VM.Models.Post import *
from VM.Models.Post.ManagedMathLib import *
from VM.Operations.Post import *
from VM.Operations.Post.Interfaces import *
from VM.Operations.Post.Models import *
from VM.Post.API.OutputReader import *
from VM.ViewModels.Post import *
from VM.ViewModels.Post.Entities.Charts import *
from VM.Windows.Post.Controls.Model import *
