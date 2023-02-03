class UserExpression:

    def __init__(self, rep):
        self._rep = rep

    def __repr__(self):
        return to_str(self._rep)

    def __neg__(self):
        return UserExpression(Mul(negone, self._rep))

    def __add__(self, other):
        other_expr = expressify(other)
        if other_expr is None:
            return NotImplemented
        return UserExpression(Add(self._rep, other_expr._rep))

    def __radd__(self, other):
        other_expr = expressify(other)
        if other_expr is None:
            return NotImplemented
        return UserExpression(Add(other_expr._rep, self._rep))

    def __sub__(self, other):
        other_expr = expressify(other)
        if other_expr is None:
            return NotImplemented
        return UserExpression(Add(self._rep, Mul(negone, other_expr._rep)))

    def __mul__(self, other):
        other_expr = expressify(other)
        if other_expr is None:
            return NotImplemented
        return Mul(self._rep, other_expr._rep)

    def __rmul__(self, other):
        other_expr = expressify(other)
        if other_expr is None:
            return NotImplemented
        return UserExpression(Mul(other_expr._rep, self._rep))

    def __pow__(self, other):
        other_expr = expressify(other)
        if other_expr is None:
            return NotImplemented
        return UserExpression(Pow(self._rep, other_expr._rep))

    def eval_f64(self, values=None):
        if values is None:
            values = {}
        return eval_f64(self._rep, values)

    def diff(self, symbol, ntimes=1):
        derivative = self._rep
        for n in range(ntimes):
            derivative = diff(derivative, symbol._rep)
        return UserExpression(derivative)


expressify_dict = {
    int: Integer,
}

def expressify(obj) -> UserExpression:
    """Convert a general Python object to an expression"""
    if isinstance(obj, UserExpression):
        return obj
    converter = expressify_dict.get(type(obj), None)
    if converter is not None:
        return UserExpression(converter(obj))
    raise ValueError("Cannot convert %s to expression" % type(obj))

def symbols(names: str) -> list[UserExpression]:
    return [UserExpression(Symbol(name)) for name in names.split(', ')]
