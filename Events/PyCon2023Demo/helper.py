from dataclasses import dataclass
from collections import namedtuple


Node = namedtuple('Node', ['number', 'x', 'y'])


@dataclass
class Beam:
    start: Node
    end: Node
    number: int = None

    def __str__(self):
        return f'Beam {self.number: <3}: ' \
               f'{self.start.number: ^3} - {self.end.number: ^3}'


def assess_for_breaks(mapdl, element_list, yield_str):
    mapdl.set('LAST')
    failed_elements = []
    breakage = False
    for beam in element_list:
        equiv = mapdl.get_value('secr', beam, 's', 'eqv', 'max')
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
