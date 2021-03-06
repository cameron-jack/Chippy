"""simple test to catch syntax errors from incorrect editing of conflicts"""
import sys
sys.path.extend(['../src', '..'])

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.
    
    __import__ only imports the top-level module.
    
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod


_modules_to_import = [
    'chippy.core.read_count',
    'chippy.parse.sam',
    'chippy.ref.genome'
    ]

for mod in _modules_to_import:
    my_import(mod)
