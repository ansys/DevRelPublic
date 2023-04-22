from ansys.mapdl.core import launch_mapdl
from helper import Beam, Node, assess_for_breaks, plot


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
# Yield strength = 6 MPa (similar UTS also)
mapdl.mp("EX", 1, 10.e9)
mapdl.mp("PRXY", 0.33)
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
               [1.0, 0],
               [0, -0.1],
               [1.0, -0.1]]

nodes = [Node(mapdl.n(i+1, *coord, 0), *coord)
         for i, coord in enumerate(coordinates)]

mapdl.nplot(cpos="xy")

node_pairings = [(1, 2), (2, 3), (3, 4), (4, 5),
                 (5, 6), (6, 7), (7, 8), (8, 9),
                 (9, 10), (10, 11), (12, 2), (13, 10)]

elements = [Beam(*pair, mapdl.e(*pair)) for pair in node_pairings]

mapdl.eplot(show_node_numbering=True, cpos="xy")

# constrain nodes at fixed end
mapdl.nsel("ALL")
for n in [1, 12, 11, 13]:
    mapdl.d(n, "ALL")
mapdl.f(6, "FY", -60)
mapdl.finish()

mapdl.run("/SOLU")
mapdl.solve()
mapdl.finish()
mapdl.post1()

does_it_break, output = assess_for_breaks(mapdl, elements, yield_strength)
plot(elements, nodes, yield_strength)
print(output)
mapdl.exit()
