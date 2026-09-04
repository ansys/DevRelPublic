'''Hello, Mechanical!'''

#1- Rename an object in the tree
Model.Name="Hello, Mechanical!"

#2- Change a material 
Model.Geometry.Children[1].Material = "Aluminum Alloy"

if Model.Geometry.Children[1].Material!="Aluminum Alloy": # Optional: to check if the material was found in the tree
    print("The material was not found")

#3- Mesh with default settings
Model.Mesh.GenerateMesh() #Calling a function always ends with "()"

#4- Add a Fixed Support
myNS=DataModel.GetObjectsByName("Fixed_support_face")[0] #Select the firt item in the tree wth the name "Fixed_support_face"
fixed_support=Model.Analyses[0].AddFixedSupport() #Add fixed support
fixed_support.Location=myNS #Scope named selection

#5- Insert Deformation
Model.Analyses[0].Solution.AddTotalDeformation()

#6- Solve Solution
Model.Analyses[0].Solve()

#7- Save results in a file
result=Model.Analyses[0].Solution.Children[1] #Take the Total Deformation Result
result.ExportToTextFile(Model.Analyses[0].WorkingDir+"result.txt") #Save the data in the chosen path with chosen name

