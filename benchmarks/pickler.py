# Helper functions for easy saving and loading of named results
#
# Amit Moscovich, Tel Aviv University, 2022.

import os
import pickle
import gzip
import errno
import inspect
import time


CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PICKLES_PATH = os.path.join(CURRENT_DIR, 'pickles')


def _pickled_name(name):
    return os.path.join(PICKLES_PATH, name) + '.pickle.gz'


def dump(filename, **kwargs):
    os.makedirs(PICKLES_PATH, exist_ok=True)
    filename = _pickled_name(filename)

    metadata = {'date': time.ctime()}

    print('Saving to', filename)
    print("Saved fields: ", ', '.join(sorted(kwargs.keys())))

    with gzip.GzipFile(filename, 'wb') as f:
        pickle.dump({'metadata': metadata, 'data': kwargs}, f, 2)


class StructFromDict(object):
    def __init__(self, d): 
        self.__dict__.update(d)
    def __repr__(self):
        return repr(self.__dict__)


def load(name):
    filename = _pickled_name(name)
    print('Loading', filename)
    with gzip.GzipFile(filename, 'rb') as f:
        d = pickle.load(f)
    print('Creation time:', d['metadata']['date'])
    return StructFromDict(d['data'])


