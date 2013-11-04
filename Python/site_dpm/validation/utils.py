
def matrices_by_name(name):
    "@return: Those matrices/consensus sequences for a factor of the given name"
    import biopsy.transfac as T
    name = name.upper()
    for m in T.Matrix.all():
        if -1 != m.name.find(name):
            yield m
    for s in T.Site.all():
        if -1 != s.name.find(name):
            yield s
