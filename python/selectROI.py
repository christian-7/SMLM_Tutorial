import matplotlib.pyplot as plt

	def selectROI(data):

		def line_select_callback(eclick, erelease):
    			global x1, y1, x2, y2
    			x1, y1 = eclick.xdata, eclick.ydata
    			x2, y2 = erelease.xdata, erelease.ydata
    			print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))

		fig, current_ax = plt.subplots() 

		plt.hexbin(data.x_nm,data.y_nm, gridsize=25,bins='log', cmap='inferno');

		# drawtype is 'box' or 'line' or 'none'
		toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True)
