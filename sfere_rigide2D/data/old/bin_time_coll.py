from matplotlib import pyplot

filename= 'time_collision.dat'
lines = open(filename).readlines()
x = [int(line.strip()) for line in lines]
bins = [i*1000 for i in range(10)]

pyplot.hist(x,bins=bins, facecolor='green', alpha=0.75)
pyplot.grid(True)

pyplot.save( 'tempi_collisione.png')

