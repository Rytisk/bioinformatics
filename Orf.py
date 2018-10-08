class Orf:    
    def __init__(self, start, stop, frame, value, length):
        self.start = start
        self.stop = stop
        self.value = value
        self.frame = frame
        self.length = length
    
    def __str__(self):
        return 'ORF <Start: {}, Stop: {}, Length: {}, Frame: {}>'.format(self.start, self.stop, self.length, self.frame)