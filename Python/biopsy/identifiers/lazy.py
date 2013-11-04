#
# Copyright John Reid 2007
#

"""
Code for lazy cached initialisation.
"""
import cPickle

class LazyInitialiser(object):
    "Calls an initialiser function on first use and thereafter returns its value."

    def __init__(self, initialiser):
        self.initialiser = initialiser
        self.instance = None

    def __call__(self):
        if None == self.instance:
            self.instance = self.initialiser()
        return self.instance

class PersistableLazyInitialiser(LazyInitialiser):
    """
    A PersistableLazyInitialiser is a proxy for an object that is unpickled or
    initialised on first use. To get the value of the object do: self().
    """
    def __init__(self, initialiser, pickle_file):
        LazyInitialiser.__init__(self, initialiser)
        self.pickle_file = pickle_file

    def __call__(self):
        if None == self.instance:
            try:
                print 'Unpickling: %s' % self.pickle_file
                self.instance = cPickle.load(open(self.pickle_file))
            except:
                print 'Unpickling failed: %s' % self.pickle_file
                self.instance =  LazyInitialiser.__call__(self)
        return self.instance

    def persist(self):
        "Pickle object to disk"
        print 'Pickling: %s' % self.pickle_file
        cPickle.dump(self.instance, open(self.pickle_file, 'w'))
