import numpy
from scipy import interpolate

from matplotlib import pyplot

#from sklearn.metrics import roc_curve, auc

'''
Created on Fri Oct 1, 2010

@author: Bastiaan van den Berg
'''


class RocCollection(object):

    def __init__(self):
        self.roc_list = []
        self.roc_names = []

    def add(self, r, name=''):
        self.roc_list.append(r)
        self.roc_names.append(name)

    def is_empty(self):
        return len(self.roc_list) == 0

    def save_roc_plot(self, f, colors, linestyles, transparent=False, dpi=300):

        # TODO turn this into create empty ROC axes function
        # create figure axes
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)

        # plot the random classification line
        ax.plot([0, 1], [0, 1], color='#babdb6', linestyle='--')

        # plot all rocs in the collection
        for i, r in enumerate(self.roc_list):
            r.add_to_roc_axes(ax, label=self.roc_names[i], color=colors[i],
                              linestyle=linestyles[i])

        # general plot settings
        ax.grid()
        ax.set_xlabel('false positive rate')
        ax.set_ylabel('true positive rate')
        ax.legend(loc="lower right", prop={'size': 8})

        fig.savefig(f, transparent=transparent, bbox_inches='tight')
        pyplot.close(fig)

    def save_avg_roc_plot(self, f, color='#3465a4', linestyle='-',
                          subcolor='#d3d7cf', sublinestyle='-',
                          transparent=False):
        fig = self.get_avg_roc_plot(color, linestyle, subcolor, sublinestyle)
        fig.savefig(f, transparent=transparent, bbox_inches='tight')
        pyplot.close(fig)

    def get_avg_roc_plot(self, color='#3465a4', linestyle='-',
                         subcolor='#d3d7cf', sublinestyle='-', avg_only=False):

        # TODO turn this into create empty ROC axes function
        # create figure axes
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)

        # plot the random classification line
        ax.plot([0, 1], [0, 1], color='#babdb6', linestyle='--')

        if not(avg_only):
            # plot all rocs in the collection
            for index, r in enumerate(self.roc_list):
                r.add_to_roc_axes(ax, color=subcolor, linestyle=sublinestyle)

        # plot the average roc
        x, y = self.avg_roc()
        avg, std = self.avg_auc()
        ax.plot(x, y, color=color, label='avg-auc = %0.3f (std = %0.3f))' %
                (avg, std))

        # general plot settings
        ax.grid()
        ax.set_xlabel('false positive rate')
        ax.set_ylabel('true positive rate')
        ax.legend(loc="lower right", prop={'size': 8})
        return fig

    def avg_roc(self):

        step = 0.01
        x_grid = numpy.arange(0.0, 1.0 + 0.5 * step, step)

        y_values = []
        for x in x_grid:
            y_values.append(numpy.mean([r.get_y(x) for r in self.roc_list]))

        return (x_grid, y_values)

    def avg_auc(self):
        '''
        Returns mean and std of the area under the curves of all ROC-curves in
        the collection.

        NOTE: the area under the average roc-curve might be slightly different.
        '''
        aucs = [r.auc() for r in self.roc_list]
        return (numpy.mean(aucs), numpy.std(aucs))


class ROC(object):
    '''
    '''

    def __init__(self, class_labels, predictions, class0=-1, class1=1):
        '''
        ROC-curve for 2-class classification
        The predictions and the class labels must be a list with numbers,
        negatives for one class, positives for the other class.
        predictions: list with real numbers, negative values for one class
                     positive values for the other
        class_labels: list with class labels, -1.0 for one class, 1.0 for
                      the other class. both values must be present in the
                      list at least once. A ValueError will be raised
                      otherwise.
        A value error will be raised when the lists have different sizes
        '''
        self._check_args(class_labels, predictions, class0, class1)
        self.predictions = predictions
        self.class_labels = class_labels
        self.class0 = class0
        self.class1 = class1
        # TODO make roc monotonically increasing, so that
        # scipy.interpolate.interp1d can be used to for interpolation
        # UPDATE not sure if possible, vertical lines impossible...
        self.roc = self._roc()

    def _check_args(self, class_labels, predictions, class0, class1):
        if not(len(predictions) == len(class_labels)):
            raise ValueError('Unequal number of predictions and class labels.')
        if not(type(class0) == int and type(class1) == int):
            raise ValueError('Class labels must be int values.')
        if(class1 == class0):
            raise ValueError('Classes must be different.')
        if(class1 < class0):
            raise ValueError('Class 0 must be smaller than class 1.')
        if not(set(numpy.unique(class_labels)) == set([class0, class1])):
            raise ValueError('Labels contain other than class0 and class1')
        # TODO check classes
        #class_set = set(class_labels)
        #if not(val_set == set(classes)):
        #    raise ValueError('Class labels may only be -1.0 and 1.0,' +
        #                     'and both values must be at least once present' +
        #                     'in the class labels list.')

    def _roc(self):
        '''
        >>> from bio import roc
        >>> pred = [-1.0, 1.0]
        >>> clab = [-1.0, 1.0]
        >>> r = roc.ROC(pred, clab)
        >>> r.roc
        ([0.0, 0.0, 1.0], [0.0, 1.0, 1.0])

        '''

        # list with x and y values roc-curve
        x = []
        y = []

        bins = [self.class0, self.class1, self.class1 + 1]
        hist = numpy.histogram(self.class_labels, bins=bins)
        class_counts = dict(zip(hist[1], hist[0]))

        # determine N and P (number of negative and positive class labels)
        n = class_counts[self.class0]
        p = class_counts[self.class1]

        assert(n + p == len(self.class_labels))

        pred_sorted = sorted(self.predictions[:])

        first = pred_sorted[0] - 0.1
        last = pred_sorted[-1] + 0.1
        middle = []
        for i in range(len(pred_sorted) - 1):
            if not(pred_sorted[i] == pred_sorted[i + 1]):
                midpoint = pred_sorted[i] + 0.5 *\
                    (pred_sorted[i + 1] - pred_sorted[i])
                middle.append(midpoint)

        thresholds = []
        thresholds.append(first)
        thresholds.extend(middle)
        thresholds.append(last)
        thresholds.reverse()

        for threshold in thresholds:

            # determine number of false and true positives
            fp = 0
            tp = 0
            for i in range(len(self.predictions)):
                if(self.predictions[i] > threshold):
                    if(self.class_labels[i] == self.class0):
                        fp += 1
                    elif(self.class_labels[i] == self.class1):
                        tp += 1
                    else:
                        raise ValueError('Labels can only be %i or %i.' %
                                         (self.class0, self.class1))

            # calculate false and true positive rate
            fpr = float(fp) / n
            tpr = float(tp) / p
            # add to list
            x.append(fpr)
            y.append(tpr)

        return (x, y)

    def get_y(self, x):
        '''
        pre: 0.0 <= x <= 1.0
        '''

        y_indices = self._get_y_indices(x)

        # there is exactly one y-value for given x, return that
        if(len(y_indices) == 1):
            return self.roc[1][y_indices[0]]

        # there are multiple y-values for giver x, return max - min / 2
        # TODO: check if this is the correct approach
        elif(len(y_indices) > 1):
            miny = self.roc[1][y_indices[0]]
            maxy = self.roc[1][y_indices[-1]]
            return (maxy + miny) / 2.0

        # else, no y_value available, use interpolation surrounding y-values
        else:
            x_pre_i = self._get_before_index(x)
            x_pre = self.roc[0][x_pre_i]
            x_post = self.roc[0][x_pre_i + 1]
            y_pre = self.roc[1][x_pre_i]
            y_post = self.roc[1][x_pre_i + 1]
            f = interpolate.interp1d([x_pre, x_post], [y_pre, y_post])
            return f(x)

    def _get_y_indices(self, x):
        '''
        Return roc list indices where roc x-value is x.
        '''
        return [i for i, roc_x in enumerate(self.roc[0]) if roc_x == x]

    def _get_before_index(self, x):
        for index in xrange(len(self.roc[0]) - 1):
            x_before = self.roc[0][index]
            x_after = self.roc[0][index + 1]
            if(x_before < x and x_after > x):
                return (index)

    def save_roc_plot(self, f, label='', color='#3465a4', linestyle='-',
                      transparent=False):
        fig = self.get_roc_plot(f, label, color, linestyle)
        fig.savefig(f, transparent=transparent, bbox_inches='tight')
        pyplot.close(fig)

    def get_roc_plot(self, f, label='', color='#3465a4', linestyle='-'):
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot([0, 1], [0, 1], color='#babdb6', linestyle='--')
        self.add_to_roc_axes(ax, label, color, linestyle)
        ax.grid()
        ax.set_xlabel('false positive rate')
        ax.set_ylabel('true positive rate')
        ax.legend(loc="lower right", prop={'size': 8})
        return fig

    def add_to_roc_axes(self, ax, label='', color='#d3d7cf', linestyle='-'):
        x, y = self.roc
        if(len(label) > 0):
            label = '%s (auc = %0.3f)' % (label, self.auc())
            ax.plot(x, y, color=color, label=label)
        else:
            ax.plot(x, y, color=color)

    #def auc_roc(self, limit=1.0):
    def auc(self):

        #if(limit < 0.0 or limit > 1.0):
        #    raise ValueError('Limit must be in range [0.0, 1.0].')

        (xvalues, yvalues) = self.roc

        '''
        area = 0.0
        roci = 0

        while roci + 1 < len(xvalues) and xvalues[roci + 1] <= limit:
            x0 = xvalues[roci]
            x1 = xvalues[roci + 1]
            y0 = yvalues[roci]
            y1 = yvalues[roci + 1]
            miny = min(y0, y1)
            binarea = (x1 - x0) * miny + 0.5 * (x1 - x0) * abs(y1 - y0)
            area += binarea
            roci += 1

        # add left part of the last bin (x0 to limit)
        if(roci + 1 < len(xvalues) and xvalues[roci + 1] > limit):
            x0 = xvalues[roci]
            x1 = xvalues[roci + 1]
            y0 = yvalues[roci]
            y1 = yvalues[roci + 1]
            miny = min(y0, y1)
            xdist = limit - x0
            ydist = (xdist / (x1 - x0)) * abs(y1 - y0)
            binarea = xdist * miny + 0.5 * xdist * ydist
            area = area + binarea

        return area
        '''
        return numpy.trapz(yvalues, xvalues)
