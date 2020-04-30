import matplotlib.pyplot as plt
from matplotlib.widgets  import RectangleSelector

def selectROI(data):

	click = [None,None]
	release = [None,None]

	def line_select_callback(eclick, erelease):
		click[:] = eclick.xdata, eclick.ydata
		release[:] = erelease.xdata, erelease.ydata
		x1, y1 = eclick.xdata, eclick.ydata
		x2, y2 = erelease.xdata, erelease.ydata
		print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))

	def toggle_selector(event):
		print(' Key pressed.')
		if event.key in ['Q', 'q'] and toggle_selector.RS.active:
			print(' RectangleSelector deactivated.')
			toggle_selector.RS.set_active(False)
		if event.key in ['A', 'a'] and not toggle_selector.RS.active:
			print(' RectangleSelector activated.')
			toggle_selector.RS.set_active(True)

	fig, current_ax = plt.subplots() 

	plt.hexbin(data.x_nm,data.y_nm, gridsize=25,bins='log', cmap='inferno');

# drawtype is 'box' or 'line' or 'none'
	toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True)
	return click, release