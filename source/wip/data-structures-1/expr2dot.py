from sympy import *

def top_sort(expression):
    seen = set()
    sorted_expressions = []

    def _recurse(e):
        if e in seen:
            return
        else:
            for e2 in e.args:
                _recurse(e2)
            seen.add(e)
            sorted_expressions.append(e)

    _recurse(expression)
    return sorted_expressions

def make_dot(expression):
    exprs = top_sort(expression)
    indices = {e: n for n, e in enumerate(exprs)}
    nodename = lambda e: f'node{indices[e]}'

    atoms = []

    lines = []
    lines.append('digraph g {')
    #lines.append('concentrate=True')
    for e in exprs:
        lines.append(f'{nodename(e)} [')
        if not e.args:
            lines.append(f'label = "{e}"')
            atoms.append(e)
        else:
            boxes = [f'<head> {e.func.__name__}']
            boxes += [f'<f{n}>' for n in range(len(e.args))]
            line = ' | '.join(boxes)
            lines.append(f'label = "{line}"')
            lines.append('shape = "record"')
        lines.append('];')
        for n, a in enumerate(e.args):
            lines.append(f'"{nodename(e)}":f{n} -> "{nodename(a)}";')

    rankstr = ' '.join(f'{nodename(e)};' for e in atoms)
    lines.append(f'{{rank = same; {rankstr} }}')

    lines.append('}')
    return '\n'.join(lines)


e = diff(x*cos(sin(exp(tan(x)))), x)
print(make_dot(e))

with open('diff_graph_2.dot', 'w') as fout:
    fout.write(make_dot(e))
