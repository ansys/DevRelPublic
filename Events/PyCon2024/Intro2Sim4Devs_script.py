from ansys.mapdl.core import launch_mapdl

# start mapdl and clear it
mapdl = launch_mapdl()
mapdl.clear()  # optional as MAPDL just started
mapdl.units("SI")  # SI - International system (m, kg, s, K).
mapdl.prep7()

mapdl.antype("STATIC")

mapdl.et(1, "BEAM188")
mapdl.sectype(1, "BEAM", "RECT")
mapdl.secdata(0.01, 0.01)

# ice
# Young's Modulus = 10 GPa
# Poisson's ratio = 0.33
# Yield strength = 6 MPa
mapdl.mp("EX", 1, 10.e9)
mapdl.mp("PRXY", 0.33)
mapdl.mp("DENS", 917.)
yield_strength = 6.0e6

coordinates = [[0, 0],
               [0.1, 0],
               [0.2, 0],
               [0.3, 0],
               [0.4, 0],
               [0.5, 0],
               [0.6, 0],
               [0.7, 0],
               [0.8, 0],
               [0.9, 0],
               [1.0, 0]]

nodes = [mapdl.n(i+1, *coord, 0) for i, coord in enumerate(coordinates)]
# Uncomment to plot/list nodes
# mapdl.nlist()
# mapdl.nplot(cpos="xy")

node_pairings = [(1, 2), (2, 3), (3, 4), (4, 5),
                 (5, 6), (6, 7), (7, 8), (8, 9),
                 (9, 10), (10, 11)]
    
elements = [mapdl.e(*pair) for pair in node_pairings]
# Uncomment to plot/list elements
# mapdl.elist()
# mapdl.eplot(show_node_numbering=True, cpos="xy")


# constrain nodes at fixed end
mapdl.nsel("ALL")
for n in [1, 11]:
    mapdl.d(n, "ALL")
mapdl.f(6, "FY", -60)
mapdl.finish()

mapdl.run("/SOLU")
mapdl.solve()
mapdl.finish()
mapdl.post1()

yield_strength = 6.0e6
mapdl.set('LAST')
for beam in elements:
    max_equivalent_stress = mapdl.get_value('secr', beam, 's', 'eqv', 'max')
    percentage_diff = 100. * max_equivalent_stress / yield_strength
    print(f"Bridge Element {beam} experienced {percentage_diff}% the max stress it can handle!")

mapdl.exit()