#Model associated: Model1.wbpz
#In this Script, the goal is to get the Ids of the faces scoped to the existing 2 bonded contacts (defined inside Connections)

contact_group=DataModel.GetObjectsByName("Contacts_to_change")[0] #Look for the first ([0]) object in the tree with the name "Contacts_to_change"

my_contact_faces=[] #Initiate a list
my_target_faces=[] #Initiate a list

for contact in contact_group.Children:#Loop for all the Children of "Contacts_to_Change" (2 here)
    my_contact_faces.append(contact.SourceLocation)
    my_target_faces.append(contact.TargetLocation)


IDs_contact_faces=[face.Ids[0] for face in my_contact_faces]
IDs_target_faces=[face.Ids[0] for face in my_target_faces]

print(IDs_contact_faces)
print(IDs_target_faces)
