class Orf:    
    def __init__(self, start, stop, value):
        self.start = start
        self.stop = stop
        self.value = value
    
    def __str__(self):
        return 'ORF <Start: {}, Stop: {}>'.format(self.start, self.stop)