from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.reader.rst import Result


def assess_for_breaks(result, node_list, yield_str):
    nodes_array, stresses = result.principal_nodal_stress(0, node_list)
    equiv_stresses = [s[-1] for s in stresses]
    failed_nodes = []
    for n_, equiv in zip(nodes_array, equiv_stresses):
        percentage_diff = 100. * equiv / yield_str
        if abs(equiv) > yield_str:
            failed_nodes.append((n_, True, equiv, percentage_diff))
        else:
            failed_nodes.append((n_, False, equiv, percentage_diff))
    num_failures = len([i for i in failed_nodes if i[1]])
    print(f'Node Failures: {num_failures} out of {len(failed_nodes)}')
    print('Details:')
    for node in failed_nodes:
        if node[1]:
            addendum = '- BREAK'
        else:
            addendum = ''
        print(f'{node[0]}: {node[-1]:g} % of yield {addendum}')
    if num_failures > 0:
        print(f'{"BRIDGE BREAKS":*^30}')


# start mapdl and clear it
mapdl = launch_mapdl()
mapdl.clear()  # optional as MAPDL just started
mapdl.units("SI")  # SI - International system (m, kg, s, K).
mapdl.prep7()

mapdl.antype("STATIC")
mapdl.et(1, "BEAM188")
mapdl.sectype(1, "BEAM", "CSOLID")
mapdl.secdata(.01/2.)
# ice
# Young's Modulus = 10 GPa
# Poisson's ratio = 0.33
# Yield strength = 6 MPa (similar UTS also)
mapdl.mp("EX", 1, 10.e9)
mapdl.mp("PRXY", 0.33)
yield_strength = 6.0e6
print(mapdl.slist())

nodes = [mapdl.n(1, 0, 0, 0), mapdl.n(2, 0.1, 0, 0), mapdl.n(3, 0.2, 0, 0),
         mapdl.n(4, 0.3, 0, 0), mapdl.n(5, 0.4, 0, 0), mapdl.n(6, 0.5, 0, 0),
         mapdl.n(7, 0.6, 0, 0), mapdl.n(8, 0.7, 0, 0), mapdl.n(9, 0.8, 0, 0),
         mapdl.n(10, 0.9, 0, 0), mapdl.n(11, 1.0, 0, 0),
         mapdl.n(12, 0, -0.1, 0), mapdl.n(13, 1.0, -0.1, 0)]

mapdl.nplot(cpos="xy")

elements = [mapdl.e(1, 2), mapdl.e(2, 3), mapdl.e(3, 4), mapdl.e(4, 5),
            mapdl.e(5, 6), mapdl.e(6, 7), mapdl.e(7, 8), mapdl.e(8, 9),
            mapdl.e(9, 10), mapdl.e(10, 11), mapdl.e(12, 2), mapdl.e(13, 10)]

mapdl.eplot(show_node_numbering=True, cpos="xy")

# constrain nodes at fixed end
mapdl.nsel("ALL")
for n in [1, 12, 11, 13]:
    mapdl.d(n, "ALL")
mapdl.f(5, "FY", -20)
mapdl.f(6, "FY", -20)
mapdl.f(7, "FY", -20)
mapdl.finish()

mapdl.run("/SOLU")
# mapdl.outres("ALL", "ALL")
mapdl.solve()
mapdl.finish()
mapdl.post1()

simulation_result = mapdl.result  # type: Result

assess_for_breaks(simulation_result, nodes, yield_strength)
simulation_result.plot_principal_nodal_stress(
    0,
    "SEQV",
    show_edges=True,
    cmap="viridis",
    cpos="xy",
    render_lines_as_tubes=True,
    line_width=20.
)

mapdl.finish()
mapdl.exit()
