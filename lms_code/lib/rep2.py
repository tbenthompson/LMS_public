import os
import dill as pickle

cfg = dict()
cfg['data_dir'] = 'data'
cfg['repeat_enabled'] = True

def full(filename):
    return cfg['data_dir'] + '/' + filename

# Note that pickling a numpy array is much faster than pickling the associated
# list because the numpy array is known to be of a homogeneous type.
def save(filename, data):
    full_path = full(filename)
    dirname = os.path.dirname(full_path)
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else: raise
    with open(full_path, 'w') as f:
        pickle.dump(data, f, protocol = -1)

def load(filename, data = None):
    full_path = full(filename)
    new_data = None
    with open(full_path, 'r') as f:
        new_data = pickle.load(f)
    if data is None:
        return new_data
    else:
        data.update(new_data)

def run_or_load(filename, function, *args):
    full_path = full(filename)
    if os.path.exists(full_path):
        return load(filename)
    else:
        result = function(*args)
        save(filename, result)
        return result

def repeat(module, func, *args):
    """
    Repeatedly reload a module and re-run the specified function until the
    user uses Ctrl-C to exit. Useful for plotting data that takes a long time
    to load.
    """
    if cfg['repeat_enabled'] is False:
        return module.__dict__[func](*args)
    while True:
        try:
            module.__dict__[func](*args)
        except Exception, e:
            print("Exception: " + str(e))
        raw_input("Press enter to run again. Ctrl-C to quit.")
        module = reload(module)
