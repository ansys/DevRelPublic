from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.reader.rst import Result


# start mapdl and clear it
mapdl = launch_mapdl()
mapdl.clear()  # optional as MAPDL just started
mapdl.units("SI")  # SI - International system (m, kg, s, K).
mapdl.prep7()

mapdl.antype("STATIC")
mapdl.et(1, "BEAM188")
mapdl.sectype(1, "BEAM", "RECT")
mapdl.secdata(0.01, 0.01)


mapdl.mp("EX", 1, 10.e9)
mapdl.mp("PRXY", 0.33)


nodes = [mapdl.n(1, 0, 0, 0),
         mapdl.n(2, 0.1, 0, 0),
         mapdl.n(3, 0.2, 0, 0),
         mapdl.n(4, 0.3, 0, 0),
         mapdl.n(5, 0.4, 0, 0),
         mapdl.n(6, 0.5, 0, 0),
         mapdl.n(7, 0.6, 0, 0),
         mapdl.n(8, 0.7, 0, 0),
         mapdl.n(9, 0.8, 0, 0),
         mapdl.n(10, 0.9, 0, 0),
         mapdl.n(11, 1.0, 0, 0),
         mapdl.n(12, 0, -0.1, 0),
         mapdl.n(13, 1.0, -0.1, 0)]


elements = [mapdl.e(1, 2), mapdl.e(2, 3), mapdl.e(3, 4), mapdl.e(4, 5),
            mapdl.e(5, 6), mapdl.e(6, 7), mapdl.e(7, 8), mapdl.e(8, 9),
            mapdl.e(9, 10), mapdl.e(10, 11), mapdl.e(12, 2), mapdl.e(13, 10)]

mapdl.nplot(True,cpos='xy')
# constrain nodes at fixed end
mapdl.nsel("ALL")
for n in [1, 12, 11, 13]:
    mapdl.d(n, "ALL")
    mapdl.f(5, "FY", -20)
    mapdl.f(6, "FY", -20)
    mapdl.f(7, "FY", -20)

mapdl.finish()

mapdl.run("/SOLU")
mapdl.solve()
mapdl.finish()
mapdl.post1()
mapdl.eshape(1)
simulation_result = mapdl.result  # type: Result

# mapdl.set('LAST')
# mapdl.esel("ALL")
# mapdl.etable("MYSTRS", 'S', 'EQV')
# for e in elements:
#    try:
#        stress = mapdl.get_value("ELEM", e, "S", "EQV")
#    except Exception:
#        pass

simulation_result.plot_principal_nodal_stress(
    0,
    "SEQV",
    show_edges=True,
    cmap="viridis",
    cpos="xy",
    render_lines_as_tubes=True,
    line_width=20.
)
nodes_, stresses = simulation_result.principal_nodal_stress(0)
stresses = [s[-1] for s in stresses]

for n, s in zip(nodes_, stresses):
    if n not in [12, 13]:
        print(n, s)

mapdl.set('LAST')
mapdl.view(1, 1, 1, 1)
mapdl.graphics("power")
mapdl.eshape(1, 1)
mapdl.show("PNG")
mapdl.plnsol('S', 'EQV')
mapdl.replot()

mapdl.finish()
mapdl.exit()
