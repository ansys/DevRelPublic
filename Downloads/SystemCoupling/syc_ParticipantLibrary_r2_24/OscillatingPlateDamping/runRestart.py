#
# Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
#

# Participant Library #
# Transient Oscillating Plate Damping Tutorial - Restart #

# Open the analysis at the end of coupling step 50.
Open(CouplingStep = 50)

# (Re-)Solve the coupled analysis until the end.
Solve()
