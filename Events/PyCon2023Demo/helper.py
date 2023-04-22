from dataclasses import dataclass
from collections import namedtuple
import matplotlib.pyplot as plt
import matplotlib as mpl
from ansys.mapdl.core.mapdl import _MapdlCore
from mpl_toolkits.axes_grid1 import make_axes_locatable


Node = namedtuple('Node', ['number', 'x', 'y'])


@dataclass
class Beam:
    start: int
    end: int
    number: int = None
    stress: float = None

    def __str__(self):
        return f'Beam {self.number: <3}: ' \
               f'{self.start: ^3} - {self.end: ^3}'


def assess_for_breaks(mapdl: _MapdlCore,
                      element_list: list[Beam],
                      yield_str: float):
    mapdl.set('LAST')
    failed_elements = []
    breakage = False
    for beam in element_list:
        equiv = mapdl.get_value('secr', beam.number, 's', 'eqv', 'max')
        beam.stress = equiv
        percentage_diff = 100. * equiv / yield_str
        if abs(equiv) > yield_str:
            failed_elements.append((beam, True, equiv, percentage_diff))
        else:
            failed_elements.append((beam, False, equiv, percentage_diff))
    num_failures = len([i for i in failed_elements if i[1]])
    string_result = [
        f'Beam Failures: {num_failures} out of {len(failed_elements)}'
    ]
    for element in failed_elements:
        if element[1]:
            addendum = '- BREAK'
        else:
            addendum = ''
        string_result.append(f'{element[0]}: {element[-1]:g} % '
                             f'of yield {addendum}')
    if num_failures > 0:
        string_result.append(f'{"BRIDGE BREAKS":*^30}')
        breakage = True
    return breakage, '\n'.join(string_result)


def plot(beams: list[Beam], nodes: list[Node], yield_stress: float):
    nodes_dict = {n.number: n for n in nodes}
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect(1)
    cmap = plt.cm.get_cmap('viridis')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="25%", pad=0.5)
    max_stress = max([b.stress for b in beams])
    for beam in beams:
        start_node = nodes_dict[beam.start]
        end_node = nodes_dict[beam.end]
        xx = [start_node.x, end_node.x]
        yy = [start_node.y, end_node.y]
        ax.plot(xx, yy,
                marker='o',
                linestyle='-',
                color=cmap(beam.stress/max_stress),
                lw=5., mfc='w')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    # ax.legend(loc='best')

    cb = fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap,
                              norm=plt.Normalize(vmin=0.,
                                                 vmax=max_stress)),
        cax=cax, orientation='horizontal',
        label='Stress [MPa]'
    )
    cb.ax.plot([yield_stress, yield_stress], [0, 1], color='w')
    ax.set_ylim(-0.15, 0.15)
    plt.show()
