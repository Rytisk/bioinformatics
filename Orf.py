class Orf:    
    def __init__(self, start, stop, frame, value, length, name):
        self.start = start
        self.stop = stop
        self.value = value
        self.frame = frame
        self.length = length
        self.name = name
    
    def __str__(self):
        return '{} <Start: {}, Stop: {}, Length: {}, Frame: {}>'.format(self.name, self.start, self.stop, self.length, self.frame)